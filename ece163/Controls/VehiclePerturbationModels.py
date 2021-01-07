"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create VehiclePerturbationModels.py and WindModels.
"""

import math
import pickle
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Constants import VehiclePhysicalConstants as VPC
from ece163.Containers import States
from ece163.Containers import Inputs
from ece163.Containers import Linearized
from ece163.Utilities import MatrixMath
from ece163.Utilities import Rotations

def dThrust_dVa(Va, Throttle, epsilon = 0.5):
    VAM = VehicleAerodynamicsModel.VehicleAerodynamicsModel()

    prop_force = VAM.CalculatePropForces(Va, Throttle)
    EP_force = VAM.CalculatePropForces((Va+epsilon), Throttle)

    delta = (EP_force[0] - prop_force[0])/epsilon

    return delta

def dThrust_dThrottle(Va, Throttle, epsilon = 0.01):
    VAM = VehicleAerodynamicsModel.VehicleAerodynamicsModel()

    prop_force = VAM.CalculatePropForces(Va, Throttle)
    EP_force = VAM.CalculatePropForces(Va, (Throttle+epsilon))

    delta = (EP_force[0] - prop_force[0])/epsilon

    return delta

def CreateTransferFunction(trimState, trimInputs):
    tf = Linearized.transferFunctions()

    if(trimState.Va != 0):
        tf.beta_trim = math.sin(trimState.v/trimState.Va)
    else:
        tf.beta_trim = math.copysign((math.pi/2), trimState.u)

    tf.Va_trim = trimState.Va
    tf.alpha_trim = trimState.alpha
    tf.beta_trim = trimState.beta
    tf.theta_trim = trimState.pitch
    tf.gamma_trim = trimState.pitch - trimState.alpha
    tf.phi_trim = trimState.roll

    phi = 0.5*VPC.rho*(trimState.Va**2)*VPC.S*VPC.b
    beta = (VPC.rho*trimState.Va*VPC.S)/(2*VPC.mass)
    theta = (VPC.rho*(trimState.Va**2)*VPC.c*VPC.S)/(2*VPC.Jyy)
    V = (VPC.rho*trimState.Va*VPC.S)/(VPC.mass)
    V_coeff = VPC.CD0 + (VPC.CDalpha*trimState.alpha) + (VPC.CDdeltaE*trimInputs.Elevator)

    tf.a_phi1 = -phi*VPC.Cpp*(VPC.b/(2*trimState.Va))
    tf.a_phi2 = phi*VPC.CpdeltaA

    tf.a_beta1 = -beta*VPC.CYbeta
    tf.a_beta2 = beta*VPC.CYdeltaR

    tf.a_theta1 = -theta*VPC.CMq*(VPC.c/(2*trimState.Va))
    tf.a_theta2 = -theta*VPC.CMalpha
    tf.a_theta3 = theta*VPC.CMdeltaE

    dThrottle = dThrust_dThrottle(trimState.Va, trimInputs.Throttle)
    dVa = dThrust_dVa(trimState.Va, trimInputs.Throttle)
    tf.a_V1 = (V*V_coeff)-((1/VPC.mass)*dVa)
    tf.a_V2 = (1/VPC.mass)*dThrottle
    tf.a_V3 = VPC.g0*math.cos(trimState.pitch-trimState.alpha)

    return tf