"""
Luis Mercado
lurmerca 1658336
ECE 163
LAB #4
"""


import math
import sys
import pickle
import enum

import ece163.Containers.Inputs as Inputs
import ece163.Containers.States as States
import ece163.Containers.Controls as Controls
import ece163.Constants.VehiclePhysicalConstants as VPC
import ece163.Modeling.VehicleAerodynamicsModel as VAM

class PDControl():
    def __init__(self, kp=0.0, kd=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return

    def setPDGains(self, kp = 0.0, kd = 0.0, trim = 0.0, lowLimit = 0.0, highLimit = 0.0):
        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        error = command - current
        u = (self.kp * error) - (self.kd * derivative) + self.trim
        if u > self.highLimit:
            u = self.highLimit
        elif u < self.lowLimit:
            u = self.lowLimit
        return u


class PIControl():
    def __init__(self, dT=VPC.dT, kp=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        self.prevError = 0.0
        self.accumulator = 0.0
        return

    def Update(self, command, current):
        error = command - current
        self.accumulator += self.dT * (error + self.prevError)
        u = (self.kp * error) + (self.ki * self.accumulator) + self.trim
        if u > self.highLimit:
            u = self.highLimit
        elif u < self.lowLimit:
            u = self.lowLimit
        return u

    def resetIntegrator(self):
        self.prevError = 0.0
        self.accumulator = 0.0
        return

    def setPIGains(self, dT=VPC.dT, kp = 0.0, ki = 0.0, trim = 0.0, lowLimit = 0.0, highLimit = 0.0):
        self.dT = dT
        self.kp = kp
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return


class PIDControl():
    def __init__(self, dT=VPC.dT, kp=0.0, kd=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.ki = kd
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit

        self.prevError = 0.0
        self.accumulator = 0.0
        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        error = command - current
        self.accumulator += 0.5 * self.dT*(error + self.prevError)
        u = (self.kp * error) - (self.kd * derivative) + (self.ki * self.accumulator) + self.trim

        if u > self.highLimit:
            u = self.highLimit
            self.accumulator -= 0.5 * self.dT * (error + self.prevError)
        elif u < self.lowLimit:
            u = self.lowLimit
            self.accumulator -= 0.5 * self.dT*(error + self.prevError)

        self.prevError = error
        return u

    def resetIntegrator(self):
        self.prevError = 0.0
        self.accumulator = 0.0
        return

    def setPIDGains(self, dT=VPC.dT, kp=0.0, kd=0.0, ki=0.0, trim=0.0, lowLimit=0.0, highLimit=0.0):
        self.dT = dT
        self.kp = kp
        self.kd = kd
        self.ki = ki
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return



class VehicleClosedLoopControl():
    def __init__(self, dT = VPC.dT):
        self.controlGains = Controls.controlGains()
        self.trimInputs = Inputs.controlInputs()
        self.controlSurface = Inputs.controlInputs()
        self.VAM = VAM.VehicleAerodynamicsModel()
        self.VAM.vehicleDynamics.dT = dT
        self.dT = self.VAM.vehicleDynamics.dT
        self.aileronFromRoll = PIDControl()
        self.rollFromCourse = PIControl()
        self.rudderFromSideSlip = PIControl()
        self.elevatorFromPitch = PDControl()
        self.throttleFromAirpseed = PIControl()
        self.pitchFromAltitude = PIControl()
        self.pitchFromAirspeed = PIControl()
        self.climbMode = Controls.AltitudeStates.HOLDING
        return

    def setControlGains(self, controlGains = Controls.controlGains()):
        self.controlGains = controlGains

        self.aileronFromRoll.setPIDGains(self.dT, controlGains.kp_roll, controlGains.kd_roll, controlGains.ki_roll,
                                         self.trimInputs.Aileron, VPC.minControls.Aileron, VPC.maxControls.Aileron)

        self.rollFromCourse.setPIGains(self.dT, controlGains.kp_course, controlGains.ki_course, 0.0,
                            -math.radians(VPC.bankAngleLimit), math.radians(VPC.bankAngleLimit))

        self.rudderFromSideSlip.setPIGains(self.dT, controlGains.kp_sideslip, controlGains.ki_sideslip,
                                           self.trimInputs.Rudder, VPC.minControls.Rudder, VPC.maxControls.Rudder)

        self.elevatorFromPitch.setPDGains(controlGains.kp_pitch, controlGains.kd_pitch, 0.0,
                               VPC.minControls.Elevator, VPC.maxControls.Elevator)

        self.throttleFromAirpseed.setPIGains(self.dT, controlGains.kp_SpeedfromThrottle, controlGains.kp_SpeedfromThrottle,
                                             self.trimInputs.Throttle, VPC.minControls.Throttle, VPC.maxControls.Throttle)

        self.pitchFromAltitude.setPIGains(self.dT, controlGains.kp_altitude, controlGains.ki_altitude, 0.0,
                               -math.radians(VPC.pitchAngleLimit), math.radians(VPC.pitchAngleLimit))

        self.pitchFromAirspeed.setPIGains(self.dT, controlGains.kp_SpeedfromElevator, controlGains.ki_SpeedfromElevator, 0.0,
                               -math.radians(VPC.pitchAngleLimit), math.radians(VPC.pitchAngleLimit))
        return


    def setTrimInputs(self, trimInputs = Inputs.controlInputs(Throttle=0.5, Aileron=0.0, Elevator=0.0, Rudder=0.0)):
        self.trimInputs.Throttle = trimInputs.Throttle
        self.trimInputs.Aileron = trimInputs.Aileron
        self.trimInputs.Elevator = trimInputs.Elevator
        self.trimInputs.Rudder = trimInputs.Rudder
        return


    def setVehicleState(self, state):
        self.VAM.vehicleDynamics.state = state
        return


    def getControlGains(self):
        return self.controlGains


    def getVehicleState(self):
        return self.VAM.vehicleDynamics.state


    def getVehicleAerodynamicsModel(self):
        return self.VAM


    def getVehicleControlSurfaces(self):
        return self.controlSurface

    def getTrimInputs(self):
        return self.trimInputs


    def reset(self):
        self.VAM.reset()
        self.aileronFromRoll.resetIntegrator()
        self.rollFromCourse.resetIntegrator()
        self.rudderFromSideSlip.resetIntegrator()
        self.throttleFromAirpseed.resetIntegrator()
        self.pitchFromAltitude.resetIntegrator()
        self.pitchFromAirspeed.resetIntegrator()
        return

    def Update(self, referenceCommands = Controls.referenceCommands()):

        courseErr = referenceCommands.commandedCourse - self.VAM.vehicleDynamics.state.chi

        if courseErr >= math.pi:
            self.VAM.vehicleDynamics.state.chi += (2 * math.pi)
        if courseErr <= -math.pi:
            self.VAM.vehicleDynamics.state.chi -= (2 * math.pi)

        roll_cmd = self.rollFromCourse.Update(referenceCommands.commandedCourse, self.VAM.vehicleDynamics.state.chi)
        aileron_cmd = self.aileronFromRoll.Update(roll_cmd, self.VAM.vehicleDynamics.state.roll, self.VAM.vehicleDynamics.state.p)
        rudder_cmd = self.rudderFromSideSlip.Update(0.0, self.VAM.vehicleDynamics.state.beta)

        altitude = -self.VAM.vehicleDynamics.state.pd
        if altitude > (referenceCommands.commandedAltitude + VPC.altitudeHoldZone):
            if self.climbMode is not Controls.AltitudeStates.DESCENDING:
                self.climbMode = Controls.AltitudeStates.DESCENDING
                self.pitchFromAirspeed.resetIntegrator()
            throttle_cmd = VPC.minControls.Throttle
            pitch_cmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, self.VAM.vehicleDynamics.state.Va)
        elif altitude < (referenceCommands.commandedAltitude - VPC.altitudeHoldZone):
            if self.climbMode is not Controls.AltitudeStates.CLIMBING:
                self.climbMode = Controls.AltitudeStates.CLIMBING
                self.pitchFromAirspeed.resetIntegrator()
            throttle_cmd = VPC.maxControls.Throttle
            pitch_cmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, self.VAM.vehicleDynamics.state.Va)
        else:
            if self.climbMode is not Controls.AltitudeStates.HOLDING:
                self.climbMode = Controls.AltitudeStates.HOLDING

            throttle_cmd = self.throttleFromAirpseed.Update(referenceCommands.commandedAirspeed, self.VAM.vehicleDynamics.state.Va)
            pitch_cmd = self.pitchFromAltitude.Update(referenceCommands.commandedAltitude, altitude)
        elevator_cmd = self.elevatorFromPitch.Update(pitch_cmd, self.VAM.vehicleDynamics.state.pitch, self.VAM.vehicleDynamics.state.q)

        self.controlSurface.Throttle = throttle_cmd
        self.controlSurface.Elevator = elevator_cmd
        self.controlSurface.Aileron = aileron_cmd
        self.controlSurface.Rudder = rudder_cmd
        referenceCommands.commandedPitch = pitch_cmd
        referenceCommands.commandedRoll = roll_cmd

        self.VAM.Update(self.controlSurface)
        return