import math
import sys
import pickle
import enum
from ..Containers import Inputs as Inputs
from ..Containers import States as States
from ..Containers import Controls as Controls
from ..Containers import Linearized as Linearized
from ..Constants import VehiclePhysicalConstants as VPC
from ..Modeling import VehicleAerodynamicsModel as VehicleAerodynamicsModule

def computeGains(tuningParameters=Controls.controlTuning(), linearizedModel=Linearized.transferFunctions()):

    # lateral Gains
    controlGains = Controls.controlGains()

    controlGains.kp_roll = ((tuningParameters.Wn_roll**2) / linearizedModel.a_phi2)
    controlGains.kd_roll = (((2 * tuningParameters.Zeta_roll*tuningParameters.Wn_roll) - linearizedModel.a_phi1)
                            / linearizedModel.a_phi2)
    controlGains.ki_roll = 0.001

    controlGains.kp_sideslip = (((2 * tuningParameters.Zeta_sideslip * tuningParameters.Wn_sideslip)
                                 - linearizedModel.a_beta1) / linearizedModel.a_beta2)
    controlGains.ki_sideslip = (tuningParameters.Wn_sideslip**2) / linearizedModel.a_beta2

    controlGains.kp_course = ((2 * tuningParameters.Zeta_course * tuningParameters.Wn_course * linearizedModel.Va_trim)
                              / VPC.g0)
    controlGains.ki_course = ((tuningParameters.Wn_course**2) * linearizedModel.Va_trim) / VPC.g0

    # longitudinal Gains
    controlGains.kp_pitch = ((tuningParameters.Wn_pitch**2) - linearizedModel.a_theta2) / linearizedModel.a_theta3
    controlGains.kd_pitch = (((2 * tuningParameters.Zeta_pitch * tuningParameters.Wn_pitch) - linearizedModel.a_theta1)
                             / linearizedModel.a_theta3)

    K_thetaDC = ((controlGains.kp_pitch * linearizedModel.a_theta3)
                 / (linearizedModel.a_theta2 + (controlGains.kp_pitch * linearizedModel.a_theta3)))

    controlGains.kp_altitude = ((2 * tuningParameters.Zeta_altitude * tuningParameters.Wn_altitude) /
                                (K_thetaDC * linearizedModel.Va_trim))
    controlGains.ki_altitude = (tuningParameters.Wn_altitude**2) / (K_thetaDC * linearizedModel.Va_trim)

    controlGains.kp_SpeedfromThrottle = (((2 * tuningParameters.Zeta_SpeedfromThrottle
                                           * tuningParameters.Wn_SpeedfromThrottle) - linearizedModel.a_V1)
                                         / linearizedModel.a_V2)
    controlGains.ki_SpeedfromThrottle = (tuningParameters.Wn_SpeedfromThrottle ** 2) / linearizedModel.a_V2

    controlGains.kp_SpeedfromElevator = ((linearizedModel.a_V1 - (2 * tuningParameters.Zeta_SpeedfromElevator
                                                                  * tuningParameters.Wn_SpeedfromElevator))
                                         / (K_thetaDC * VPC.g0))
    controlGains.ki_SpeedfromElevator = (-tuningParameters.Wn_SpeedfromElevator**2) / (K_thetaDC * VPC.g0)
    return controlGains

def computeTuningParameters(controlGains=Controls.controlGains(), linearizedModel=Linearized.transferFunctions()):
    tuningParameters = Controls.controlTuning()

    # tuning for lateral control
    tuningParameters.Wn_roll = math.sqrt(controlGains.kp_roll * linearizedModel.a_phi2)
    tuningParameters.Zeta_roll = ((linearizedModel.a_phi1 + (controlGains.kd_roll * linearizedModel.a_phi2))
                                  / (2 * tuningParameters.Wn_roll))

    tuningParameters.Wn_course = math.sqrt((controlGains.ki_course * VPC.g0) / linearizedModel.Va_trim)
    tuningParameters.Zeta_course = ((controlGains.kp_course * VPC.g0)
                                    / (2 * tuningParameters.Wn_course * linearizedModel.Va_trim))

    tuningParameters.Wn_sideslip = math.sqrt(linearizedModel.a_beta2 * controlGains.ki_sideslip)
    tuningParameters.Zeta_sideslip = ((linearizedModel.a_beta1 + (linearizedModel.a_beta2 * controlGains.kp_sideslip))
                                      / (2 * tuningParameters.Wn_sideslip))

    # tuning knows for longitudinal control
    tuningParameters.Wn_pitch = math.sqrt((controlGains.kp_pitch * linearizedModel.a_theta3) + linearizedModel.a_theta2)
    tuningParameters.Zeta_pitch = (((controlGains.kd_pitch * linearizedModel.a_theta3) + linearizedModel.a_theta1)
                                   / (2 * tuningParameters.Wn_pitch))

    K_thetaDC = ((controlGains.kp_pitch * linearizedModel.a_theta3)
                 / (linearizedModel.a_theta2 + (controlGains.kp_pitch * linearizedModel.a_theta3)))

    tuningParameters.Wn_altitude = math.sqrt(controlGains.ki_altitude * K_thetaDC * linearizedModel.Va_trim)
    tuningParameters.Zeta_altitude = ((controlGains.kp_altitude * K_thetaDC * linearizedModel.Va_trim)
                                      / (2 * tuningParameters.Wn_altitude))

    tuningParameters.Wn_SpeedfromThrottle = math.sqrt(controlGains.ki_SpeedfromThrottle * linearizedModel.a_V2)
    tuningParameters.Zeta_SpeedfromThrottle = ((linearizedModel.a_V1 + (controlGains.kp_SpeedfromThrottle
                                                                        * linearizedModel.a_V2))
                                               / (2 * tuningParameters.Wn_SpeedfromThrottle))


    tuningParameters.Wn_SpeedfromElevator = math.sqrt(-(controlGains.ki_SpeedfromElevator * K_thetaDC * VPC.g0))
    tuningParameters.Zeta_SpeedfromElevator = ((linearizedModel.a_V1 - (controlGains.kp_SpeedfromElevator * K_thetaDC
                                                                        * VPC.g0))
                                               / (2 * tuningParameters.Wn_SpeedfromElevator))
    return tuningParameters