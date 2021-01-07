"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create VehicleControlGains and VehicleClosedLoopControl.
"""

import math
import pickle
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Controls
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations

def computeGains (tuningParameters = Controls.controlTuning(), linearizedModel = Linearized.transferFunctions()):
    gainz = Controls.controlGains()

    gainz.kp_roll = (tuningParameters.Wn_roll**2)/linearizedModel.a_phi2
    gainz.ki_roll = 0.001
    gainz.kd_roll = ((2*tuningParameters.Zeta_roll*tuningParameters.Wn_roll)-linearizedModel.a_phi1)/linearizedModel.a_phi2

    gainz.kp_course = (2*tuningParameters.Zeta_course*tuningParameters.Wn_course*linearizedModel.Va_trim)/VPC.g0
    gainz.ki_course = ((tuningParameters.Wn_course**2)*linearizedModel.Va_trim)/VPC.g0

    gainz.ki_sideslip = (tuningParameters.Wn_sideslip**2)/linearizedModel.a_beta2
    gainz.kp_sideslip = ((2*tuningParameters.Zeta_sideslip*tuningParameters.Wn_sideslip)-linearizedModel.a_beta1)/linearizedModel.a_beta2

    gainz.kp_pitch = ((tuningParameters.Wn_pitch**2)-linearizedModel.a_theta2)/linearizedModel.a_theta3
    gainz.kd_pitch = ((2*tuningParameters.Zeta_pitch*tuningParameters.Wn_pitch)-linearizedModel.a_theta1)/linearizedModel.a_theta3

    K_DC = (gainz.kp_pitch*linearizedModel.a_theta3)/(linearizedModel.a_theta2 + (gainz.kp_pitch*linearizedModel.a_theta3))

    gainz.ki_altitude = (tuningParameters.Wn_altitude**2)/(K_DC*linearizedModel.Va_trim)
    gainz.kp_altitude = (2*tuningParameters.Zeta_altitude*tuningParameters.Wn_altitude)/(K_DC*linearizedModel.Va_trim)

    gainz.ki_SpeedfromThrottle = (tuningParameters.Wn_SpeedfromThrottle**2)/linearizedModel.a_V2
    gainz.kp_SpeedfromThrottle = ((2*tuningParameters.Zeta_SpeedfromThrottle*tuningParameters.Wn_SpeedfromThrottle)-linearizedModel.a_V1)/linearizedModel.a_V2

    gainz.ki_SpeedfromElevator = -(tuningParameters.Wn_SpeedfromElevator**2)/(K_DC*VPC.g0)
    gainz.kp_SpeedfromElevator = (linearizedModel.a_V1 - (2*tuningParameters.Zeta_SpeedfromElevator*tuningParameters.Wn_SpeedfromElevator))/(K_DC*VPC.g0)

    return gainz

def computeTuningParameters (controlGains = Controls.controlGains(), lineraizedModel = Linearized.transferFunctions()):
    tunez = Controls.controlTuning()

    try:
        tunez.Wn_roll = math.sqrt(controlGains.kp_roll*lineraizedModel.a_phi2)
        tunez.Zeta_roll = (lineraizedModel.a_phi1+lineraizedModel.a_phi2*controlGains.kd_roll)/(2*tunez.Wn_roll)

        tunez.Wn_course = math.sqrt((VPC.g0*controlGains.ki_course)/lineraizedModel.Va_trim)
        tunez.Zeta_course = (controlGains.kp_course*VPC.g0)/(lineraizedModel.Va_trim*tunez.Wn_course*2)

        tunez.Wn_sideslip = math.sqrt(lineraizedModel.a_beta2*controlGains.ki_sideslip)
        tunez.Zeta_sideslip = (lineraizedModel.a_beta1 + (lineraizedModel.a_beta2*controlGains.kp_sideslip))/(2*tunez.Wn_sideslip)

        tunez.Wn_pitch = math.sqrt(lineraizedModel.a_theta2+(controlGains.kp_pitch*lineraizedModel.a_theta3))
        tunez.Zeta_pitch = (lineraizedModel.a_theta1+(controlGains.kd_pitch*lineraizedModel.a_theta3))/(2*tunez.Wn_pitch)

        K_DC = (controlGains.kp_pitch*lineraizedModel.a_theta3)/(lineraizedModel.a_theta2 + (controlGains.kp_pitch*lineraizedModel.a_theta3))

        tunez.Wn_altitude = math.sqrt(K_DC*lineraizedModel.Va_trim*controlGains.ki_altitude)
        tunez.Zeta_altitude = (K_DC*lineraizedModel.Va_trim*controlGains.kp_altitude)/(2*tunez.Wn_altitude)

        tunez.Wn_SpeedfromThrottle = math.sqrt(lineraizedModel.a_V2*controlGains.ki_SpeedfromThrottle)
        tunez.Zeta_SpeedfromThrottle = (lineraizedModel.a_V1+(lineraizedModel.a_V2*controlGains.kp_SpeedfromThrottle))/(2*tunez.Wn_SpeedfromThrottle)

        tunez.Wn_SpeedfromElevator = math.sqrt(-K_DC*VPC.g0*controlGains.ki_SpeedfromElevator)
        tunez.Zeta_SpeedfromElevator = (lineraizedModel.a_V1 - (K_DC*VPC.g0*controlGains.kp_SpeedfromElevator))/(2*tunez.Wn_SpeedfromElevator)
    except ValueError:
        return tunez

    return tunez