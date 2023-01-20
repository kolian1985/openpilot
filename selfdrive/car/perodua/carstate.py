from cereal import car
from collections import deque
from math import ceil
from opendbc.can.parser import CANParser
from opendbc.can.can_define import CANDefine
from common.numpy_fast import mean
from selfdrive.config import Conversions as CV
from selfdrive.car.interfaces import CarStateBase
from selfdrive.car.perodua.values import DBC, CAR, ACC_CAR, HUD_MULTIPLIER
from time import time

# todo: clean this part up
pedal_counter = 0
pedal_press_state = 0
PEDAL_COUNTER_THRES = 35
PEDAL_UPPER_TRIG_THRES = 0.125
PEDAL_NON_ZERO_THRES = 0.01

SEC_HOLD_TO_STEP_SPEED = 0.6

class CarState(CarStateBase):
  def __init__(self, CP):
    super().__init__(CP)
    can_define = CANDefine(DBC[CP.carFingerprint]['pt'])
    self.shifter_values = can_define.dv["TRANSMISSION"]['GEAR']
    if CP.carFingerprint in ACC_CAR:
      self.set_distance_values = can_define.dv['ACC_CMD_HUD']['FOLLOW_DISTANCE']
    self.is_cruise_latch = False
    self.cruise_speed = 30 * CV.KPH_TO_MS
    self.cruise_speed_counter = 0
    self.acttrGas = 0

    self.is_plus_btn_latch = False
    self.is_minus_btn_latch = False
    # shared by both + and - button, since release of another button will reset this
    self.rising_edge_since = 0
    self.last_frame = time() # todo: existing infra to reuse?
    self.dt = 0

  def update(self, cp):
    ret = car.CarState.new_message()

    # there is a backwheel speed, but it will overflow to 0 when reach 60kmh
    # perodua vehicles doesn't have a good standard for their wheelspeed scaling
    ret.wheelSpeeds = self.get_wheel_speeds(
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_F'],
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_F'],
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_F'],
      cp.vl["WHEEL_SPEED"]['WHEELSPEED_F'],
    )
    ret.vEgoRaw = mean([ret.wheelSpeeds.rr, ret.wheelSpeeds.rl, ret.wheelSpeeds.fr, ret.wheelSpeeds.fl])

    # unfiltered speed from CAN sensors
    ret.vEgo, ret.aEgo = self.update_speed_kf(ret.vEgoRaw)
    ret.standstill = ret.vEgoRaw < 0.01

    # safety checks to engage
    can_gear = int(cp.vl["TRANSMISSION"]['GEAR'])

    ret.doorOpen = any([cp.vl["METER_CLUSTER"]['MAIN_DOOR'],
                     cp.vl["METER_CLUSTER"]['LEFT_FRONT_DOOR'],
                     cp.vl["METER_CLUSTER"]['RIGHT_BACK_DOOR'],
                     cp.vl["METER_CLUSTER"]['LEFT_BACK_DOOR']])

    ret.seatbeltUnlatched = cp.vl["METER_CLUSTER"]['SEAT_BELT_WARNING'] == 1
    if self.CP.carFingerprint in ACC_CAR:
      ret.seatbeltUnlatched |= cp.vl["METER_CLUSTER"]['SEAT_BELT_WARNING2'] == 1
    ret.gearShifter = self.parse_gear_shifter(self.shifter_values.get(can_gear, None))
    disengage = ret.doorOpen or ret.seatbeltUnlatched
    if disengage:
      self.is_cruise_latch = False

    # gas pedal
    ret.gas = cp.vl["GAS_PEDAL"]['APPS_1']
    # todo: let gas pressed legit
    if self.CP.carFingerprint in ACC_CAR:
      ret.gasPressed = not bool(cp.vl["GAS_PEDAL_2"]['GAS_PEDAL_STEP'])
    else:
      ret.gasPressed = False

    self.acttrGas = (cp.vl["GAS_SENSOR"]['INTERCEPTOR_GAS']) # KommuActuator gas, read when stock pedal is being intercepted
    if self.acttrGas < 0:
      self.acttrGas = 0

    # brake pedal
    ret.brake = cp.vl["BRAKE"]['BRAKE_PRESSURE']

    # perodua bezza has a lower resolution brake pressure sensor
    if self.CP.carFingerprint == CAR.BEZZA:
      ret.brakePressed = ret.brake > 1.2
    elif self.CP.carFingerprint in ACC_CAR:
      ret.brakePressed = bool(cp.vl["BRAKE"]['BRAKE_ENGAGED'])
    else:
      ret.brakePressed = ret.brake > 1e-5

    # steer
    if self.CP.carFingerprint in ACC_CAR:
      ret.steeringAngleDeg = cp.vl["STEERING_MODULE"]['STEER_ANGLE']
      ret.steeringTorque = cp.vl["STEERING_MODULE"]['MAIN_TORQUE']
      ret.steeringTorqueEps = cp.vl["EPS_SHAFT_TORQUE"]['STEERING_TORQUE']
    else:
      ret.steeringAngleDeg = cp.vl["STEERING_ANGLE_SENSOR"]['STEER_ANGLE']
      steer_dir = 1 if (ret.steeringAngleDeg >= 0) else -1
      ret.steeringTorque = cp.vl["STEERING_TORQUE"]['MAIN_TORQUE'] * steer_dir
      ret.steeringTorqueEps = ret.steeringTorque

    if self.CP.carFingerprint == CAR.AXIA:
      ret.steeringPressed = bool(abs(ret.steeringTorque) > 18)
    elif self.CP.carFingerprint in ACC_CAR:
      ret.steeringPressed = bool(abs(ret.steeringTorque) > 20)
    else:
      ret.steeringPressed = bool(abs(ret.steeringTorque) > 70)

    ret.steerWarning = False    # since Perodua has no LKAS, make it always no warning
    ret.steerError = False        # since Perodua has no LKAS, make it always no warning

    if self.CP.carFingerprint not in ACC_CAR:

      ret.vEgoCluster = ret.vEgoRaw * HUD_MULTIPLIER
      ret.stockAeb = cp.vl["ADAS_AEB"]['BRAKE_REQ'] != 0
      ret.stockFcw = cp.vl["ADAS_HUD"]['AEB_ALARM'] != 0
      ret.cruiseState.available = True

      if self.is_cruise_latch:
        # pedal disengage
        if self.check_pedal_engage(self.acttrGas, pedal_press_state):
          self.is_cruise_latch = False

        # increase cruise_speed using pedal when engage
        self.cruise_speed_counter += 1
        if self.cruise_speed_counter % 100 == 0 and self.acttrGas > 0.2:
          self.cruise_speed += (5 * CV.KPH_TO_MS)
          self.cruise_speed_counter = 0

      # latching cruiseState logic
      if not self.is_cruise_latch:
        if self.check_pedal_engage(ret.gas, pedal_press_state):
          self.cruise_speed = max(30 * CV.KPH_TO_MS, ret.vEgoCluster)
          self.is_cruise_latch = True

      # set distance as SetDistance.aggresive
      ret.cruiseState.setDistance = 3
    else:

      ret.vEgoCluster = cp.vl["BUTTONS"]["UI_SPEED"] * CV.KPH_TO_MS
      if self.CP.carFingerprint == CAR.MYVI_PSD:
          ret.vEgoCluster *= 1.04
      elif self.CP.carFingerprint == CAR.ALZA:
          ret.vEgoCluster *= 1.04
      elif self.CP.carFingerprint == CAR.ATIVA:
          ret.vEgoCluster *= 1.04
      ret.stockAdas.frontDepartureHUD = bool(cp.vl["LKAS_HUD"]["FRONT_DEPART"])
      ret.stockAdas.laneDepartureHUD = bool(cp.vl["LKAS_HUD"]["LDA_ALERT"])
      ret.stockAdas.ldpSteerV = cp.vl["STEERING_LKAS"]['STEER_CMD']

      ret.stockAeb = bool(cp.vl["LKAS_HUD"]['AEB_BRAKE'])
      ret.stockFcw = bool(cp.vl["LKAS_HUD"]['AEB_ALARM'])

      ret.cruiseState.available = cp.vl["PCM_BUTTONS"]["ACC_RDY"] != 0
      distance_val = int(cp.vl["ACC_CMD_HUD"]['FOLLOW_DISTANCE'])
      ret.cruiseState.setDistance = self.parse_set_distance(self.set_distance_values.get(distance_val, None))

      # set speed logic
      # todo: check if the logic needs to be this complicated

      minus_button = bool(cp.vl["PCM_BUTTONS"]["SET_MINUS"])
      plus_button = bool(cp.vl["PCM_BUTTONS"]["RES_PLUS"])

      if self.is_cruise_latch:
        cur_time = time()
        self.dt += cur_time - self.last_frame
        self.last_frame = cur_time

        if self.is_plus_btn_latch != plus_button: # rising or falling
          if not plus_button: # released, falling
            if cur_time - self.rising_edge_since < 1:
              self.cruise_speed += CV.KPH_TO_MS
          else: # pressed, rising, init
            self.rising_edge_since = cur_time
            self.dt = 0
        elif plus_button: # is holding
          while self.dt >= SEC_HOLD_TO_STEP_SPEED:
            kph = self.cruise_speed * CV.MS_TO_KPH
            kph += 5 - (kph % 5)  # step up to next nearest 5
            self.cruise_speed = kph * CV.KPH_TO_MS
            self.dt -= SEC_HOLD_TO_STEP_SPEED

        if self.is_minus_btn_latch != minus_button: # rising or falling
          if not minus_button: # released, falling
            if cur_time - self.rising_edge_since < 1:
              self.cruise_speed -= CV.KPH_TO_MS
          else: # pressed, rising
            self.rising_edge_since = cur_time
            self.dt = 0
        elif minus_button: # is holding
          while self.dt >= SEC_HOLD_TO_STEP_SPEED:
            kph = self.cruise_speed * CV.MS_TO_KPH
            kph = ((kph / 5) - 1) * 5  # step down to next nearest 5
            kph = max(30, kph)
            self.cruise_speed = kph * CV.KPH_TO_MS
            self.dt -= SEC_HOLD_TO_STEP_SPEED

      if not self.is_cruise_latch:
        # activate cruise onReleased
        if self.is_plus_btn_latch and not plus_button:
          self.is_cruise_latch = True

        elif self.is_minus_btn_latch and not minus_button:
          self.cruise_speed = max(30 * CV.KPH_TO_MS, ret.vEgoCluster)
          self.is_cruise_latch = True


      self.is_plus_btn_latch = plus_button
      self.is_minus_btn_latch = minus_button

      if bool(cp.vl["PCM_BUTTONS"]["CANCEL"]):
        self.is_cruise_latch = False

    if ret.brakePressed:
      self.is_cruise_latch = False

    # set speed in range of 30 - 130kmh only
    self.cruise_speed = max(min(self.cruise_speed, 130 * CV.KPH_TO_MS), 30 * CV.KPH_TO_MS)
    ret.cruiseState.speedCluster = self.cruise_speed
    ret.cruiseState.speed = ret.cruiseState.speedCluster / HUD_MULTIPLIER
    if self.CP.carFingerprint == CAR.MYVI_PSD:
      ret.cruiseState.speed *= 1.04
    ret.cruiseState.standstill = False
    ret.cruiseState.nonAdaptive = False
    ret.cruiseState.enabled = self.is_cruise_latch
    if not ret.cruiseState.available:
      self.is_cruise_latch = False

    # button presses
    ret.leftBlinker = bool(cp.vl["METER_CLUSTER"]["LEFT_SIGNAL"])
    ret.rightBlinker = bool(cp.vl["METER_CLUSTER"]["RIGHT_SIGNAL"])
    ret.genericToggle = bool(cp.vl["RIGHT_STALK"]["GENERIC_TOGGLE"])

    # blindspot sensors
    if self.CP.enableBsm:
      # used for lane change so its okay for the chime to work on both side.
      ret.leftBlindspot = bool(cp.vl["BSM"]["BSM_CHIME"])
      ret.rightBlindspot = bool(cp.vl["BSM"]["BSM_CHIME"])
    else:
      ret.leftBlindspot = False
      ret.rightBlindspot = False

    return ret

  @staticmethod
  def check_pedal_engage(gas,state):
    ''' Pedal engage logic '''
    global pedal_counter
    global pedal_press_state
    if (state == 0):
      if (gas > PEDAL_UPPER_TRIG_THRES):
        pedal_counter += 1
        if (pedal_counter == PEDAL_COUNTER_THRES):
          pedal_counter = 0
          return False
      if (pedal_counter > 2 and gas <= PEDAL_NON_ZERO_THRES):
        pedal_press_state = 1
        pedal_counter = 0
      return False
    if (state == 1):
      pedal_counter += 1
      if (pedal_counter == PEDAL_COUNTER_THRES):
        pedal_counter = 0
        pedal_press_state = 0
        return False
      if (gas > PEDAL_UPPER_TRIG_THRES):
        pedal_press_state = 2
        pedal_counter = 0
      return False
    if (state == 2):
      pedal_counter += 1
      if (pedal_counter == PEDAL_COUNTER_THRES):
        pedal_counter = 0
        pedal_press_state = 0
        return False
      if (gas <= PEDAL_NON_ZERO_THRES):
        pedal_counter = 0
        pedal_press_state = 0
        return True
    return False


  @staticmethod
  def get_can_parser(CP):
    signals = [
      # sig_name, sig_address, default
      ("WHEELSPEED_F", "WHEEL_SPEED", 0.),
      ("GEAR", "TRANSMISSION", 0),
      ("APPS_1", "GAS_PEDAL", 0.),
      ("BRAKE_PRESSURE", "BRAKE", 0.),
      ("BRAKE_ENGAGED", "BRAKE", 0),
      ("INTERCEPTOR_GAS", "GAS_SENSOR", 0),
      ("GENERIC_TOGGLE", "RIGHT_STALK", 0),
      ("FOG_LIGHT", "RIGHT_STALK", 0),
      ("LEFT_SIGNAL", "METER_CLUSTER", 0),
      ("RIGHT_SIGNAL", "METER_CLUSTER", 0),
      ("SEAT_BELT_WARNING", "METER_CLUSTER", 0),
      ("MAIN_DOOR", "METER_CLUSTER", 1),
      ("LEFT_FRONT_DOOR", "METER_CLUSTER", 1),
      ("RIGHT_BACK_DOOR", "METER_CLUSTER", 1),
      ("LEFT_BACK_DOOR", "METER_CLUSTER", 1)
    ]
    checks = []

    if CP.carFingerprint in ACC_CAR:
      signals.append(("BSM_CHIME","BSM", 0))
      signals.append(("SEAT_BELT_WARNING2","METER_CLUSTER", 0))
      signals.append(("STEER_ANGLE", "STEERING_MODULE", 0.))
      signals.append(("MAIN_TORQUE", "STEERING_MODULE", 0.))
      signals.append(("STEERING_TORQUE", "EPS_SHAFT_TORQUE", 0.))
      signals.append(("ACC_RDY", "PCM_BUTTONS", 0))
      signals.append(("SET_MINUS", "PCM_BUTTONS", 0))
      signals.append(("RES_PLUS","PCM_BUTTONS", 0))
      signals.append(("CANCEL","PCM_BUTTONS", 0))
      signals.append(("PEDAL_DEPRESSED","PCM_BUTTONS", 0))
      signals.append(("LKAS_ENGAGED", "LKAS_HUD", 0))
      signals.append(("LDA_ALERT", "LKAS_HUD", 0))
      signals.append(("ACC_CMD", "ACC_CMD_HUD", 0))
      signals.append(("STEER_CMD", "STEERING_LKAS", 0))
      signals.append(("STEER_REQ", "STEERING_LKAS", 0))
      signals.append(("FRONT_DEPART", "LKAS_HUD", 0))
      signals.append(("AEB_BRAKE", "LKAS_HUD", 0))
      signals.append(("AEB_ALARM", "LKAS_HUD", 0))
      signals.append(("SET_SPEED", "ACC_CMD_HUD", 0))
      signals.append(("FOLLOW_DISTANCE", "ACC_CMD_HUD", 0))
      signals.append(("LDA_ALERT", "LKAS_HUD", 0))
      signals.append(("GAS_PEDAL_STEP", "GAS_PEDAL_2", 0))
      signals.append(("UI_SPEED", "BUTTONS", 0))
    else:
      signals.append(("MAIN_TORQUE", "STEERING_TORQUE", 0))
      signals.append(("STEER_ANGLE", "STEERING_ANGLE_SENSOR", 0.))
      signals.append(("AEB_ALARM", "ADAS_HUD", 0))
      signals.append(("BRAKE_REQ", "ADAS_AEB", 0))
      signals.append(("WHEELSPEED_B", "WHEEL_SPEED", 0.))

    # todo: make it such that enforce_checks=True
    return CANParser(DBC[CP.carFingerprint]['pt'], signals, checks, 0, enforce_checks=False)
