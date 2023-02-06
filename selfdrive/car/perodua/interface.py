#!/usr/bin/env python3
from cereal import car
from selfdrive.swaglog import cloudlog
from selfdrive.config import Conversions as CV
from selfdrive.car import STD_CARGO_KG, scale_rot_inertia, scale_tire_stiffness, gen_empty_fingerprint, get_safety_config
from selfdrive.car.interfaces import CarInterfaceBase
from selfdrive.car.perodua.values import CAR, ACC_CAR
from selfdrive.controls.lib.desire_helper import LANE_CHANGE_SPEED_MIN

from common.features import Features
from common.params import Params

EventName = car.CarEvent.EventName

class CarInterface(CarInterfaceBase):

  @staticmethod
  def get_params(candidate, fingerprint=gen_empty_fingerprint(), has_relay=False, car_fw=None):
    ret = CarInterfaceBase.get_std_params(candidate, fingerprint)
    ret.carName = "perodua"
    ret.safetyConfigs = [get_safety_config(car.CarParams.SafetyModel.perodua)]
    ret.safetyConfigs[0].safetyParam = 1
    ret.transmissionType = car.CarParams.TransmissionType.automatic
    ret.radarOffCan = True
    ret.enableApgs = False                 # advanced parking guidance system
    ret.enableDsu = False                  # driving support unit

    ret.steerRateCost = 0.7                # Lateral MPC cost on steering rate, higher value = sharper turn
    ret.steerLimitTimer = 0.1              # time before steerLimitAlert is issued
    ret.steerControlType = car.CarParams.SteerControlType.torque
    ret.steerActuatorDelay = 0.48          # Steering wheel actuator delay in seconds

    ret.lateralTuning.init('pid')
    ret.lateralTuning.pid.kiBP, ret.lateralTuning.pid.kpBP = [[0.], [0.]]
    ret.longitudinalTuning.kpV = [0.9, 0.8, 0.8]

    ret.enableGasInterceptor = 0x201 in fingerprint[0] or 0x401 in fingerprint[0]
    ret.openpilotLongitudinalControl = True

    f = Features()
    if f.has("StockAcc"):
      ret.safetyConfigs[0].safetyParam = 2
      ret.openpilotLongitudinalControl = False

    if candidate == CAR.AXIA:
      ret.wheelbase = 2.455                # meter
      ret.steerRatio = 16.54
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.8371
      ret.mass = 850. + STD_CARGO_KG
      ret.wheelSpeedFactor = 0.927

      ret.lateralTuning.pid.kf = 0.0000715
      ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.08], [0.32]]
      ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[11, 26], [480, 830]]

      ret.longitudinalTuning.kpV = [0.8, 0.9, 1.0]

    elif candidate == CAR.MYVI:
      ret.wheelbase = 2.5
      ret.steerRatio = 12.14
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.8371
      ret.mass = 1015. + STD_CARGO_KG
      ret.wheelSpeedFactor = 0.94

      ret.lateralTuning.pid.kf = 0.000091
      ret.longitudinalTuning.kpV = [1.0, 1.2, 1.4]

      f = Features()
      if f.has("MyviAzri"):
        ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.16], [0.41]]
        ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[11, 22], [380, 670]]
      elif f.has("MyviKevin"):
        ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.10], [0.32]]
        ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[11, 28], [390, 650]]
      else:
        ret.lateralTuning.pid.kiBP, ret.lateralTuning.pid.kpBP = [[0.], [0., 18]]
        ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.10], [0.31, 0.32]]
        ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[11, 28], [390, 580]]

    elif candidate == CAR.BEZZA:
      ret.wheelbase = 2.455
      ret.steerRatio = 12.14
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.8371
      ret.mass = 940. + STD_CARGO_KG
      ret.wheelSpeedFactor = 1.0

      ret.lateralTuning.pid.kf = 0.0000918
      ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.05], [0.45]]
      ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[11, 22], [400, 630]]

    elif candidate == CAR.ARUZ:
      ret.wheelbase = 2.685
      ret.steerRatio = 16.14
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.68371
      ret.mass = 1310. + STD_CARGO_KG
      ret.wheelSpeedFactor = 1.0

      ret.lateralTuning.pid.kf = 0.0000917
      ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.098], [0.135]]
      ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[11, 22], [380, 600]]

      ret.longitudinalTuning.kpV = [1.6, 1.1, 1.1]

    elif candidate == CAR.MYVI_PSD:
      ret.wheelbase = 2.5
      ret.steerRatio = 17.44
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.9871
      ret.mass = 1025. + STD_CARGO_KG
      ret.wheelSpeedFactor = 1.32

      ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.12], [0.20]]
      ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[0.], [255]]
      ret.lateralTuning.pid.kf = 0.0000007

      ret.longitudinalTuning.kpBP = [0., 5., 20.]
      ret.longitudinalTuning.kpV = [0.6, 0.55, 0.4]
      ret.longitudinalActuatorDelayLowerBound = 0.42
      ret.longitudinalActuatorDelayUpperBound = 0.60

    elif candidate == CAR.ATIVA:
      ret.wheelbase = 2.525
      ret.steerRatio = 16.74
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.9871
      ret.mass = 1035. + STD_CARGO_KG
      ret.wheelSpeedFactor = 1.54

      ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.12], [0.22]]
      ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[0.], [255]]
      ret.lateralTuning.pid.kf = 0.0000007

      ret.longitudinalTuning.kpBP = [0., 5., 20.]
      ret.longitudinalTuning.kpV = [0.6, 0.55, 0.3]
      ret.longitudinalActuatorDelayLowerBound = 0.42
      ret.longitudinalActuatorDelayUpperBound = 0.60

    elif candidate == CAR.ALZA:
      ret.wheelbase = 2.750
      ret.steerRatio = 17.00
      ret.centerToFront = ret.wheelbase * 0.44
      tire_stiffness_factor = 0.9871
      ret.mass = 1170. + STD_CARGO_KG
      ret.wheelSpeedFactor = 1.42

      ret.lateralTuning.pid.kiV, ret.lateralTuning.pid.kpV = [[0.13], [0.23]]
      ret.lateralParams.torqueBP, ret.lateralParams.torqueV = [[0.], [255]]
      ret.lateralTuning.pid.kf = 0.0000007

      ret.longitudinalTuning.kpBP = [0., 5., 20.]
      ret.longitudinalTuning.kpV = [0.6, 0.55, 0.45]
      ret.longitudinalActuatorDelayLowerBound = 0.42
      ret.longitudinalActuatorDelayUpperBound = 0.60

    else:
      ret.dashcamOnly = True
      ret.safetyModel = car.CarParams.SafetyModel.noOutput

    if candidate in ACC_CAR:
      ret.longitudinalTuning.kiBP = [5, 7, 28]
      ret.longitudinalTuning.kiV = [0.15, 0.12, 0.06]

      if candidate == CAR.ALZA:
        ret.longitudinalTuning.kiBP = [0.]
        ret.longitudinalTuning.kiV = [0.1]

      ret.minEnableSpeed = -1
      ret.steerActuatorDelay = 0.30           # Steering wheel actuator delay in seconds
      ret.enableBsm = True
      ret.stoppingDecelRate = 0.05 # reach stopping target smoothly
    else:
      ret.longitudinalTuning.kiBP = [0.]
      ret.longitudinalTuning.kiV = [0.]

    ret.rotationalInertia = scale_rot_inertia(ret.mass, ret.wheelbase)
    ret.tireStiffnessFront, ret.tireStiffnessRear = scale_tire_stiffness(ret.mass, ret.wheelbase, ret.centerToFront, tire_stiffness_factor=tire_stiffness_factor)

    return ret

  # returns a car.CarState
  def update(self, c, can_strings):
    # to receive CAN Messages
    self.cp.update_strings(can_strings)

    ret = self.CS.update(self.cp)
    ret.canValid = self.cp.can_valid
    ret.steeringRateLimited = self.CC.steer_rate_limited if self.CC is not None else False

    # events
    events = self.create_common_events(ret)

    # create events for auto lane change below allowable speed
    if ret.vEgo < LANE_CHANGE_SPEED_MIN and (ret.leftBlinker or ret.rightBlinker):
      events.add(EventName.belowLaneChangeSpeed)

    # events for non ACC cars
    if self.CP.carFingerprint not in ACC_CAR:
      # create events to warn user that their vehicle doesn't have brakes
      if self.CC.brake_pressed and (ret.cruiseState.speed >= ret.vEgo):
        events.add(EventName.promptDriverBrake)

    ret.events = events.to_msg()

    self.CS.out = ret.as_reader()
    return self.CS.out

  # pass in a car.CarControl to be called at 100hz
  def apply(self, c):

    isLdw = c.hudControl.leftLaneDepart or c.hudControl.rightLaneDepart

    can_sends = self.CC.update(c.enabled, self.CS, self.frame, c.actuators, c.hudControl.leadVisible, c.hudControl.rightLaneVisible, c.hudControl.leftLaneVisible, c.cruiseControl.cancel, isLdw)

    self.frame += 1
    return can_sends
