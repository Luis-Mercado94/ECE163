import math
import pickle
from ..Modeling import VehicleAerodynamicsModel
from ..Constants import VehiclePhysicalConstants as VPC
from ..Controls import VehicleTrim
from ..Containers import States
from ..Containers import Inputs
from ..Containers import Linearized
from ..Utilities import MatrixMath
from ..Utilities import Rotations


def CreateTransferFunction(trimState, trimInputs):
    # create instance of transfer functions
    tf = Linearized.transferFunctions()

    # determine transfer function parameters
    trimState.Va = math.hypot(trimState.u, trimState.w, trimState.v)
    trimState.alpha = math.atan2(trimState.w, trimState.u)
    if trimState.Va != 0:
        trimState.beta = math.asin(trimState.v/trimState.Va)
    else:
        trimState.beta = math.copysign(math.pi / 2, trimState.u)

    dTdVa = dThrust_dVa(trimState.Va, trimInputs.Throttle)
    dTdThrottle = dThrust_dThrottle(trimState.Va, trimInputs.Throttle)

    tf.Va_trim = trimState.Va
    tf.alpha_trim = trimState.alpha
    tf.beta_trim = trimState.beta

    tf.phi_trim = trimState.roll
    tf.gamma_trim = trimState.pitch - trimState.alpha
    tf.theta_trim = trimState.pitch

    tf.a_beta1 = -((VPC.rho * trimState.Va * VPC.S) / (2 * VPC.mass)) * VPC.CYbeta
    tf.a_beta2 = ((VPC.rho * trimState.Va * VPC.S) / (2 * VPC.mass)) * VPC.CYdeltaR

    tf.a_phi1 = -(0.5 * VPC.rho * (trimState.Va**2) * VPC.S * VPC.b * VPC.Cpp * (VPC.b / (2 * trimState.Va)))
    tf.a_phi2 = (0.5 * VPC.rho * (trimState.Va**2) * VPC.S * VPC.b * VPC.CpdeltaA)

    tf.a_theta1 = ((-((VPC.rho * (trimState.Va**2) * VPC.c * VPC.S) / (2 * VPC.Jyy)))
                   * VPC.CMq * (VPC.c / (2 * trimState.Va)))
    tf.a_theta2 = (-((VPC.rho * (trimState.Va**2) * VPC.c * VPC.S) / (2 * VPC.Jyy))) * VPC.CMalpha
    tf.a_theta3 = ((VPC.rho * (trimState.Va**2) * VPC.c * VPC.S) / (2 * VPC.Jyy)) * VPC.CMdeltaE

    tf.a_V1 = ((((VPC.rho * trimState.Va * VPC.S) / VPC.mass) * (VPC.CD0 + (VPC.CDalpha * trimState.alpha)
                                                                 + (VPC.CDdeltaE * trimInputs.Elevator)))
               - ((1 / VPC.mass) * dTdVa))
    tf.a_V2 = ((1 / VPC.mass) * dTdThrottle)
    tf.a_V3 = VPC.g0 * math.cos(tf.gamma_trim)

    return tf


def dThrust_dThrottle(Va, Throttle, epsilon=0.01):
    # find dThrottle
    V = VehicleAerodynamicsModel.VehicleAerodynamicsModel()

    Fx, Mx = V.CalculatePropForces(Va, Throttle)
    Fx_hat, Mx_hat = V.CalculatePropForces(Va, Throttle + epsilon)

    dTdDeltaT = (Fx_hat - Fx) / epsilon
    return dTdDeltaT


def dThrust_dVa(Va, Throttle, epsilon=0.5):
    # find dVa
    V = VehicleAerodynamicsModel.VehicleAerodynamicsModel()

    Fx, Mx = V.CalculatePropForces(Va, Throttle)
    Fx_hat, Mx_hat = V.CalculatePropForces(Va + epsilon, Throttle)
    
    dTdVa = (Fx_hat - Fx) / epsilon
    return dTdVa