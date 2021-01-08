"""
Luis Mercado, Pauline Mae F. Cuaresma, and Sabrina Fong
ECE 163 Final Project
Due 12/17/20

Added to this module for the final project are lines to ascend and descend the UAV when a collision procedure is
initiated by a crash signal. It takes precedence over the holding zone.
Lines 255 - 280 contain the changes to this module for the final project.
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
        #feel like Im missing stuff here
        return

    def setPDGains(self, kp = 0.0, kd = 0.0, trim = 0.0, lowLimit = 0.0, highLimit = 0.0):
        self.kp = kp
        self.kd = kd
        self.trim = trim
        self.lowLimit = lowLimit
        self.highLimit = highLimit
        return

    def Update(self, command=0.0, current=0.0, derivative=0.0):
        err = command - current # calculate path error
        u = (self.kp * err) - (self.kd * derivative) + self.trim
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
        # feel like Im missing stuff here
        return

    def Update(self, command, current):
        err = command - current
        self.accumulator += self.dT * (err + self.prevError) #accumulator -- ?? previous error is it actually 0?
        u = (self.kp * err) + (self.ki * self.accumulator) + self.trim
        if u > self.highLimit:
            u = self.highLimit
        elif u < self.lowLimit:
            u = self.lowLimit
        return u

    def resetIntegrator(self):
        self.prevError = 0.0  # previous error
        self.accumulator = 0.0  # accumulator
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
        err = command - current
        self.accumulator += 0.5 * self.dT*(err + self.prevError)  # accumulator -- ?? previous error is it actually 0?
        u = (self.kp * err) - (self.kd * derivative) + (self.ki * self.accumulator) + self.trim

        if u > self.highLimit:
            u = self.highLimit
            self.accumulator -= 0.5 * self.dT * (err + self.prevError)
        elif u < self.lowLimit:
            u = self.lowLimit
            self.accumulator -= 0.5 * self.dT*(err + self.prevError)

        self.prevError = err
        return u

    def resetIntegrator(self):
        self.prevError = 0.0  # previous error
        self.accumulator = 0.0  # accumulator
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
        self.controlSurfaceOutputs = Inputs.controlInputs()
        self.VAM = VAM.VehicleAerodynamicsModel()
        self.VAM.vehicleDynamics.dT = dT
        self.dT = self.VAM.vehicleDynamics.dT                 #self.dT must be separate from VAM dT
        self.aileronFromRoll = PIDControl()       #PID
        self.rollFromCourse = PIControl()         #PI
        self.rudderFromSideSlip = PIControl()     #PI
        self.elevatorFromPitch = PDControl()      #PD
        self.throttleFromAirpseed = PIControl()   #PI
        self.pitchFromAltitude = PIControl()      #PI
        self.pitchFromAirspeed = PIControl()      #PI
        self.climbMode = Controls.AltitudeStates.HOLDING
        return


    """ Setting gains for all controllers listed in Init """
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
                               -math.radians(VPC.pitchAngleLimit), math.radians(VPC.pitchAngleLimit)) #?? unsure about this
        return


    def setTrimInputs(self, trimInputs = Inputs.controlInputs(Throttle=0.5, Aileron=0.0, Elevator=0.0, Rudder=0.0)):
        # self.trimInputs.Throttle = trimInputs.Throttle
        # self.trimInputs.Aileron = trimInputs.Aileron
        # self.trimInputs.Elevator = trimInputs.Elevator
        # self.trimInputs.Rudder = trimInputs.Rudder
        self.trimInputs = trimInputs
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
        return self.controlSurfaceOutputs

    def getTrimInputs(self):
        return self.trimInputs


    def reset(self):
        self.VAM.reset() #reset aero dynamic states
        """ reset all (non-PD) controller related attributes """
        self.aileronFromRoll.resetIntegrator()       # PID
        self.rollFromCourse.resetIntegrator()        # PI
        self.rudderFromSideSlip.resetIntegrator()    # PI
        self.throttleFromAirpseed.resetIntegrator()  # PI
        self.pitchFromAltitude.resetIntegrator()     # PI
        self.pitchFromAirspeed.resetIntegrator()     # PI
        return


    """ Taking course error, altitude, climb mode and calculating throttle, pitch, elevator command """
    def Update(self, referenceCommands = Controls.referenceCommands()):
        """ renaming states to make it shorter reference """
        chi = self.VAM.vehicleDynamics.state.chi # check this later, might cause problems
        roll = self.VAM.vehicleDynamics.state.roll
        p = self.VAM.vehicleDynamics.state.p
        beta = self.VAM.vehicleDynamics.state.beta
        Va = self.VAM.vehicleDynamics.state.Va
        pitch = self.VAM.vehicleDynamics.state.pitch
        q = self.VAM.vehicleDynamics.state.q

        courseErr = referenceCommands.commandedCourse - chi

        if courseErr >= math.pi:
            chi += (2 * math.pi)
        if courseErr <= -math.pi:
            chi += -(2 * math.pi)

        """ recursive """
        roll_cmd = self.rollFromCourse.Update(referenceCommands.commandedCourse, chi)
        aileron_cmd = self.aileronFromRoll.Update(roll_cmd, roll, p)
        rudder_cmd = self.rudderFromSideSlip.Update(0.0, beta)

        """----working on longitudinal air pilot---
        For final project: accounted longitudinal air pilot control if crash is detected (lines 255 - 280)
        """

        altitude = -self.VAM.vehicleDynamics.state.pd #recall neg of down positive of state
        if self.VAM.crash == True:
            dir = self.VAM.avoidance()
            # set the first direction to the previous direction to be remembered
            if self.VAM.preDir == None:
                self.VAM.preDir = dir
            # if a previous direction is determined, use that direction
            if dir != self.VAM.preDir:
                dir = self.VAM.preDir
            if dir == True: #up
                if self.climbMode is not Controls.AltitudeStates.CLIMBING:
                    self.climbMode = Controls.AltitudeStates.CLIMBING
                    # reset integrator/ pitch from airspeed (which cntrls altitude) if too much air accumulated vertically
                    self.pitchFromAirspeed.resetIntegrator()
                """get throttle & pitch cmds to deal which change made"""
                throttle_cmd = VPC.maxControls.Throttle
                pitch_cmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, Va)
            if dir == False: #down
                if self.climbMode is not Controls.AltitudeStates.DESCENDING:
                    self.climbMode = Controls.AltitudeStates.DESCENDING
                    # reset integrator/ pitch from airspeed (which cntrls altitude) if too much air accumulated vertically
                    self.pitchFromAirspeed.resetIntegrator()
                """get throttle & pitch cmds to deal which change made"""
                throttle_cmd = VPC.minControls.Throttle
                pitch_cmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, Va)
        else:
            self.VAM.preDir = None
            if altitude > (referenceCommands.commandedAltitude + VPC.altitudeHoldZone):
                # if not already descending, descend
                if self.climbMode is not Controls.AltitudeStates.DESCENDING:
                    self.climbMode = Controls.AltitudeStates.DESCENDING
                    # reset integrator/ pitch from airspeed (which cntrls altitude) if too much air accumulated vertically
                    self.pitchFromAirspeed.resetIntegrator()
                """get throttle & pitch cmds to deal which change made"""
                throttle_cmd = VPC.minControls.Throttle
                pitch_cmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, Va)
            elif altitude < (referenceCommands.commandedAltitude - VPC.altitudeHoldZone):
                if self.climbMode is not Controls.AltitudeStates.CLIMBING:
                    self.climbMode = Controls.AltitudeStates.CLIMBING
                    # reset integrator/ pitch from airspeed (which cntrls altitude) if too much air accumulated vertically
                    self.pitchFromAirspeed.resetIntegrator()
                """get throttle & pitch cmds to deal which change made"""
                throttle_cmd = VPC.maxControls.Throttle
                pitch_cmd = self.pitchFromAirspeed.Update(referenceCommands.commandedAirspeed, Va)
            else:
                if self.climbMode is not Controls.AltitudeStates.HOLDING:
                    self.climbMode = Controls.AltitudeStates.HOLDING   #make it hold
                    self.pitchFromAirspeed.resetIntegrator()

                throttle_cmd = self.throttleFromAirpseed.Update(referenceCommands.commandedAirspeed, Va)
                pitch_cmd = self.pitchFromAltitude.Update(referenceCommands.commandedAltitude, altitude)
            # ------------------------------------------
        elevator_cmd = self.elevatorFromPitch.Update(pitch_cmd, pitch, q)

        """ setting controlSurfaceOutputs to corresponding commands """
        self.controlSurfaceOutputs.Throttle = throttle_cmd
        self.controlSurfaceOutputs.Elevator = elevator_cmd
        self.controlSurfaceOutputs.Aileron = aileron_cmd
        self.controlSurfaceOutputs.Rudder = rudder_cmd
        """setting reference commands to correlating commands """
        referenceCommands.commandedPitch = pitch_cmd
        referenceCommands.commandedRoll = roll_cmd

        """ updating """
        self.VAM.Update(self.controlSurfaceOutputs)
        return









