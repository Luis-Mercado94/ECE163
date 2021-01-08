import math
import random
from ..Containers import States
from ..Utilities import MatrixMath
from ..Constants import VehiclePhysicalConstants as VPC

class WindModel():
    def __init__(self, dT=VPC.dT, Va=VPC.InitialSpeed, drydenParameters=VPC.DrydenNoWind):
        # initialize wind states, gust model
        self.Wind = States.windState()

        self.x_u = [[0.0]]
        self.x_v = [[0.0],
                    [0.0]]
        self.x_w = [[0.0],
                    [0.0]]

        self.Phi_u = [[1.0]]
        self.Gamma_u = [[0.0]]
        self.H_u = [[1.0]]

        self.phi_v = [[1.0, 0.0],
                      [0.0, 1.0]]
        self.Gamma_v = [[0.0],
                        [0.0]]
        self.H_v = [[1.0, 1.0]]

        self.Phi_w = [[1.0, 0.0],
                      [0.0, 1.0]]
        self.Gamma_w = [[0.0],
                        [0.0]]
        self.H_w = [[1.0, 1.0]]

        # determine discrete transfer functions for gust models
        self.drydenParam = drydenParameters
        self.CreateDrydenTransferFns(dT, Va, self.drydenParam)
        return

    def CreateDrydenTransferFns(self, dT, Va, drydenParameters):
        # if there is no wind set dryden transfer functions to 1/0
        if drydenParameters == VPC.DrydenNoWind:
            self.Phi_u = [[1.0]]
            self.Gamma_u = [[0.0]]
            self.H_u = [[1.0]]

            self.Phi_v = [[1.0, 0.0],
                          [0.0, 1.0]]
            self.Gamma_v = [[0.0],
                            [0.0]]
            self.H_v = [[1.0, 1.0]]

            self.Phi_v = [[1.0, 0.0],
                         [0.0, 1.0]]
            self.Gamma_v = [[0.0],
                            [0.0]]
            self.H_v = [[1.0, 1.0]]
            return

        # return if Va is close to zero
        if math.isclose(Va, 0.0):
            return

        # calculate dryden transfer functions
        self.Phi_u = [[math.exp((-Va * dT) / drydenParameters.Lu)]]
        self.Gamma_u = [[(drydenParameters.Lu / Va) * (1 - math.exp((-Va * dT) / drydenParameters.Lu))]]
        self.H_u = [[drydenParameters.sigmau * math.sqrt((2 * Va) / (math.pi * drydenParameters.Lu))]]

        phi_vM = [[1 - ((Va * dT) / drydenParameters.Lv), -((Va / drydenParameters.Lv)**2) * dT],
                  [dT, 1 + ((Va * dT) / drydenParameters.Lv)]]
        self.Phi_v = MatrixMath.matrixScalarMultiply(math.exp((-Va * dT) / drydenParameters.Lv), phi_vM)
        gamma_vM = [[dT],
                    [((((drydenParameters.Lv / Va)**2) * ((math.exp((Va * dT) / drydenParameters.Lv)) - 1))
                      - ((drydenParameters.Lv * dT) / Va))]]
        self.Gamma_v = MatrixMath.matrixScalarMultiply(math.exp((-Va * dT) / drydenParameters.Lv), gamma_vM)
        alphaHv = [[1, Va / (math.sqrt(3) * drydenParameters.Lv)]]
        self.H_v = MatrixMath.matrixScalarMultiply(drydenParameters.sigmav * math.sqrt((3 * Va) / (math.pi * drydenParameters.Lv)), alphaHv)

        phi_wM = [[1 - ((Va * dT) / drydenParameters.Lw), -((Va / drydenParameters.Lw)**2) * dT],
                  [dT, 1 + ((Va * dT) / drydenParameters.Lw)]]
        self.Phi_w = MatrixMath.matrixScalarMultiply(math.exp((-Va * dT) / drydenParameters.Lw), phi_wM)
        gamma_wM = [[dT],
                    [((((drydenParameters.Lw / Va)**2) * ((math.exp((Va * dT) / drydenParameters.Lw)) - 1))
                      - ((drydenParameters.Lw * dT) / Va))]]
        self.Gamma_w = MatrixMath.matrixScalarMultiply(math.exp((-Va * dT) / drydenParameters.Lw), gamma_wM)
        alphaHw = [[1, Va / (math.sqrt(3) * drydenParameters.Lw)]]
        self.H_w = MatrixMath.matrixScalarMultiply(drydenParameters.sigmaw * math.sqrt((3 * Va) / (math.pi * drydenParameters.Lw)), alphaHw)
        return

    def Update(self, uu=None, uv=None, uw=None):
        # use random value for noise if noise inputs are not given
        if uu is None:
            rando_u = random.gauss(0, 1)
        else:
            rando_u = uu
        if uv is None:
            rando_v = random.gauss(0, 1)
        else:
            rando_v = uv
        if uw is None:
            rando_w = random.gauss(0, 1)
        else:
            rando_w = uw

        # update wind gust
        self.x_u = MatrixMath.matrixAdd(MatrixMath.matrixMultiply(self.x_u, self.Phi_u), MatrixMath.matrixScalarMultiply(rando_u, self.Gamma_u))
        wu = MatrixMath.matrixMultiply(self.H_u, self.x_u)
        self.Wind.Wu = wu[0][0]

        self.x_v = MatrixMath.matrixAdd(MatrixMath.matrixMultiply(self.Phi_v, self.x_v), MatrixMath.matrixScalarMultiply(rando_v, self.Gamma_v))
        wv = MatrixMath.matrixMultiply(self.H_v, self.x_v)
        self.Wind.Wv = wv[0][0]

        self.x_w = MatrixMath.matrixAdd(MatrixMath.matrixMultiply(self.Phi_w, self.x_w), MatrixMath.matrixScalarMultiply(rando_w, self.Gamma_w))
        ww = MatrixMath.matrixMultiply(self.H_w, self.x_w)
        self.Wind.Ww = ww[0][0]
        return

    def getWind(self):
        # return wind state
        return self.Wind

    def reset(self):
        # reset wind model and coloring filters
        self.Wind = States.windState()

        self.x_u = [[0.0]]
        self.x_v = [[0.0],
                    [0.0]]
        self.x_w = [[0.0],
                    [0.0]]
        return

    def setWind(self, windState):
        # set wind state to given state
        self.Wind = windState
        return

    def getDrydenTransferFns(self):
        # retrieve dryden transfer functions
        Phi_u = self.Phi_u
        Gamma_u = self.Gamma_u
        H_u = self.H_u

        Phi_v = self.Phi_v
        Gamma_v = self.Gamma_v
        H_v = self.H_v

        Phi_w = self.Phi_w
        Gamma_w = self.Gamma_w
        H_w = self.H_w
        return Phi_u, Gamma_u, H_u, Phi_v, Gamma_v, H_v, Phi_w, Gamma_w, H_w


