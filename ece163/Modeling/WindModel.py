"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create VehiclePerturbationModels.py and WindModels.
"""

import math
import random
from ..Containers import States
from ..Utilities import MatrixMath
from ..Constants import VehiclePhysicalConstants as VPC

class WindModel:
    def __init__(self, dT = VPC.dT, Va = VPC.InitialSpeed, drydenParameters=VPC.DrydenNoWind):
        self.Wind = States.windState()

        self.x_u = [[0.0]]
        self.x_v = [[0.0], [0.0]]
        self.x_w = [[0.0], [0.0]]

        self.CreateDrydenTransferFns(dT, Va, drydenParameters)

        # self.Phi_u = [[1.0]]
        # self.Gamma_u = [[0.0]]
        # self.H_u = [[1.0]]
        #
        # self.Phi_v = [[1.0, 1.0], [1.0, 1.0]]
        # self.Gamma_v = [[0.0], [0.0]]
        # self.H_v = [[1.0, 1.0]]
        #
        # self.Phi_w = [[1.0, 1.0], [1.0, 1.0]]
        # self.Gamma_w = [[0.0], [0.0]]
        # self.H_w = [[1.0, 1.0]]

    def CreateDrydenTransferFns(self, dT, Va, drydenParameters):
        # print("\n",drydenParameters)

        if (drydenParameters == VPC.DrydenNoWind):
            # print("\n", drydenParameters)
            self.Phi_u = [[1.0]]
            self.Gamma_u = [[0.0]]
            self.H_u = [[1.0]]

            self.Phi_v = [[1.0, 0.0], [0.0, 1.0]]
            self.Gamma_v = [[0.0], [0.0]]
            self.H_v = [[1.0, 1.0]]

            self.Phi_w = [[1.0, 0.0], [0.0, 1.0]]
            self.Gamma_w = [[0.0], [0.0]]
            self.H_w = [[1.0, 1.0]]

            return

        if (math.isclose(Va, 0.0)):
            raise ArithmeticError("No speed error")
        # print("\nhere\n")
        exp_u = math.exp((-Va/drydenParameters.Lu)*dT)
        exp_v = (Va/drydenParameters.Lv)*dT
        exp_w = (Va/drydenParameters.Lw)*dT

        PV_matrix = [[1-(Va/drydenParameters.Lv)*dT, -(Va/drydenParameters.Lv)**2*dT],
                    [dT, 1+(Va/drydenParameters.Lv)*dT]]
        GV_matrix = [[dT],
                    [((drydenParameters.Lv/Va)**2*(math.exp(exp_v) - 1))-((drydenParameters.Lv/Va)*dT)]]
        HV_matrix = [[1, Va/(drydenParameters.Lv*math.sqrt(3))]]

        Pw_matrix = [[1 - (Va / drydenParameters.Lw) * dT, -(Va / drydenParameters.Lw) ** 2 * dT],
                     [dT, 1 + (Va / drydenParameters.Lw) * dT]]
        Gw_matrix = [[dT],
                     [(drydenParameters.Lw / Va) ** 2 * (math.exp(exp_w) - 1) - (drydenParameters.Lw / Va) * dT]]
        Hw_matrix = [[1, Va / (drydenParameters.Lw * math.sqrt(3))]]

        u_sqrt = (2*Va)/(drydenParameters.Lu*math.pi)
        v_sqrt = (3*Va)/(math.pi*drydenParameters.Lv)
        w_sqrt = (3*Va)/(math.pi*drydenParameters.Lw)

        sigv = drydenParameters.sigmav*(math.sqrt(v_sqrt))
        sigw = drydenParameters.sigmaw*(math.sqrt(w_sqrt))

        self.Phi_u = [[exp_u]]
        self.Gamma_u = [[(drydenParameters.Lu/Va)*(1 - exp_u)]]
        self.H_u = [[drydenParameters.sigmau*math.sqrt(u_sqrt)]]

        self.Phi_v = MatrixMath.matrixScalarMultiply(math.exp(-exp_v), PV_matrix)
        self.Gamma_v = MatrixMath.matrixScalarMultiply(math.exp(-exp_v), GV_matrix)
        self.H_v = MatrixMath.matrixScalarMultiply(sigv, HV_matrix)

        self.Phi_w = MatrixMath.matrixScalarMultiply(math.exp(-exp_w), Pw_matrix)
        self.Gamma_w = MatrixMath.matrixScalarMultiply(math.exp(-exp_w), Gw_matrix)
        self.H_w = MatrixMath.matrixScalarMultiply(sigw, Hw_matrix)

        return

    def getDrydenTransferFns(self):
        return self.Phi_u, self.Gamma_u, self.H_u, self.Phi_v, self.Gamma_v, self.H_v, self.Phi_w, self.Gamma_w, self.H_w

    def getWind(self):
        return self.Wind

    def setWind(self, windState):
        self.Wind = windState

    def Update(self, uu = None, uv = None, uw = None):
        if uu == None:
            uu = random.gauss(0,1)
        if uv == None:
            uv = random.gauss(0,1)
        if uw == None:
            uw = random.gauss(0,1)

        xu = MatrixMath.matrixMultiply(self.Phi_u, self.x_u)
        gu = MatrixMath.matrixScalarMultiply(uu, self.Gamma_u)
        xu_plus = MatrixMath.matrixAdd(xu, gu)

        Wu = MatrixMath.matrixMultiply(self.H_u, xu_plus)
        self.Wind.Wu = Wu[0][0]
        self.x_u = xu_plus

        xv = MatrixMath.matrixMultiply(self.Phi_v, self.x_v)
        gv = MatrixMath.matrixScalarMultiply(uv, self.Gamma_v)
        xv_plus = MatrixMath.matrixAdd(xv, gv)

        Wv = MatrixMath.matrixMultiply(self.H_v, xv_plus)
        self.Wind.Wv = Wv[0][0]
        self.x_v = xv_plus

        xw = MatrixMath.matrixMultiply(self.Phi_w, self.x_w)
        gw = MatrixMath.matrixScalarMultiply(uw, self.Gamma_w)
        xw_plus = MatrixMath.matrixAdd(xw, gw)

        Ww = MatrixMath.matrixMultiply(self.H_w, xw_plus)
        self.Wind.Ww = Ww[0][0]
        self.x_w = xw_plus

        return

    def reset(self):
        self.Wind = States.windState()
        self.x_u = [[0.0]]
        self.x_v = [[0.0], [0.0]]
        self.x_w = [[0.0], [0.0]]
