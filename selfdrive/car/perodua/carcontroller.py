from cereal import car
from selfdrive.car import make_can_msg, apply_std_steer_torque_limits, create_gas_command
from selfdrive.car.perodua.peroduacan import create_steer_command, perodua_create_gas_command, perodua_aeb_brake
from opendbc.can.packer import CANPacker
from selfdrive.car.perodua.values import DBC
import cereal.messaging as messaging
from common.numpy_fast import clip, interp

# livetuner import
from selfdrive.livetune_conf import livetune_conf

class CarControllerParams():
  def __init__(self):

    self.STEER_MAX = round(float((livetune.conf['maxSteer'])))                          # dac write value
    self.STEER_STEP = 1                                                      # how often we update the steer cmd
    self.STEER_DELTA_UP = int(livetune.conf['steerDeltaUp'])                 # torque increase per refresh, 0.8s to max
    self.STEER_DELTA_DOWN = int(livetune.conf['steerDeltaDown'])             # torque decrease per refresh
    self.STEER_DRIVER_ALLOWANCE = int(livetune.conf['steerDriverAllowance']) # allowed driver torque before start limiting
    self.STEER_DRIVER_MULTIPLIER = int(livetune.conf['steerDriverMult'])     # weight driver torque heavily
    self.STEER_DRIVER_FACTOR = 1                                             # from dbc


class CarController():
  def __init__(self, dbc_name, CP, VM):
    self.last_steer = 0
    self.steer_rate_limited = False
    self.steering_direction = False
    #self.params = CarControllerParams()
    self.packer = CANPacker(DBC[CP.carFingerprint]['pt'])
    self.brake_pressed = True

  def update(self, enabled, CS, frame, actuators, visual_alert, pcm_cancel):
   
    can_sends = []

    livetune = livetune_conf()
    self.STEER_MAX = round(float((livetune.conf['maxSteer'])))               # dac write value
    self.STEER_STEP = 1                                                      # how often we update the steer cmd
    self.STEER_DELTA_UP = round(float(livetune.conf['steerDeltaUp']))        # torque increase per refresh, 0.8s to max
    self.STEER_DELTA_DOWN = round(float(livetune.conf['steerDeltaDown']))    # torque decrease per refresh
    self.STEER_DRIVER_ALLOWANCE = round(float(livetune.conf['steerDriverAllowance'])) # allowed driver torque before start limiting
    self.STEER_DRIVER_MULTIPLIER = round(float(livetune.conf['steerDriverMult']))     # weight driver torque heavily
    self.STEER_DRIVER_FACTOR = 1                                             # from dbc
    
    # generate steering command
    if (frame % self.STEER_STEP) == 0 and enabled:
      new_steer = int(round(actuators.steer * self.STEER_MAX))                # actuator.steer is from range from -1.0 - 1.0
      apply_steer = apply_std_steer_torque_limits(new_steer, self.last_steer, CS.out.steeringTorqueEps, self)
      self.steer_rate_limited = ( new_steer != apply_steer ) and (apply_steer != 0)
    else:
      apply_steer = 0
    # limit steering based on angle
    if (abs(CS.out.steeringAngle) > 220):
      apply_steer = 0

    # lower the fighting torque during manual steer
    if (CS.out.steeringPressed):
      apply_steer = apply_steer / 100
    
    self.last_steer = apply_steer
    if apply_steer >= 0:
      self.steering_direction = True
    else:
      self.steering_direction = False
    
    if (frame % self.STEER_STEP) == 0:
      can_sends.append(create_steer_command(self.packer, apply_steer, self.steering_direction, enabled, frame))

    FRAME_DIVIDER = 2

    # Gas Command Interceptor
    if (frame % FRAME_DIVIDER) == 0:
      idx = frame // FRAME_DIVIDER
      apply_gas = clip(actuators.gas, 0., 1.)
      apply_brake = clip(actuators.brake, 0., 1.)
      apply_gas = abs(apply_gas * round(float(livetune.conf['maxGas'])))

      # gas
      if CS.CP.enableGasInterceptor:
        # create_gas_command inherited from car
        can_sends.append(perodua_create_gas_command(self.packer, apply_gas, enabled, idx))

      # brakes
      if apply_brake > 0.2:
        if not self.brake_pressed:
          can_sends.append(perodua_aeb_brake(self.packer, apply_brake))
          self.brake_pressed = True
      else:
        self.brake_pressed = False

    return can_sends
