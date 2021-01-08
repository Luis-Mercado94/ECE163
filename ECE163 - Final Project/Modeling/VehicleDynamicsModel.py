import math
from ..Containers import States
from ..Containers import Inputs
from ..Utilities import MatrixMath
from ..Utilities import Rotations
from ..Constants import VehiclePhysicalConstants as VPC


class VehicleDynamicsModel():
    def __init__(self, dT=VPC.dT):
        # keep track of internal state
        self.state = States.vehicleState()
        self.dot = States.vehicleState()
        self.dT = dT
        return

    def ForwardEuler(self, forcesMoments):
        # find derivative and integrate
        dot = self.derivative(self.state, forcesMoments)
        state = self.IntegrateState(self.dT, self.state, dot)
        return state

    def IntegrateState(self, dT, state, dot):
        newState = States.vehicleState()

        # integrate R
        Rexp = self.Rexp(dT, state, dot)
        newState.R = MatrixMath.matrixMultiply(Rexp, state.R)

        # integrate euler angles
        newState.yaw, newState.pitch, newState.roll = Rotations.dcm2Euler(newState.R)

        # integrate body rates
        newState.p = state.p + (dot.p * dT)
        newState.q = state.q + (dot.q * dT)
        newState.r = state.r + (dot.r * dT)

        # integrate position
        newState.pn = state.pn + (dot.pn * dT)
        newState.pe = state.pe + (dot.pe * dT)
        newState.pd = state.pd + (dot.pd * dT)

        # integrate velocities
        newState.u = state.u + (dot.u * dT)
        newState.v = state.v + (dot.v * dT)
        newState.w = state.w + (dot.w * dT)

        # integrate airspeed and chi
        newState.Va = state.Va
        newState.alpha = state.alpha
        newState.beta = state.beta
        newState.chi = math.atan2(dot.pe, dot.pn)
        return newState

    def Rexp(self, dT, state, dot):
        # declare I and current/ acceleration vectors
        I = [[1, 0, 0],
             [0, 1, 0],
             [0, 0, 1]]
        current = [[state.p],
                   [state.q],
                   [state.r]]
        accel = [[dot.p],
                 [dot.q],
                 [dot.r]]

        # find w, wx, wx2, and magnitude of w
        w = MatrixMath.matrixAdd(current, MatrixMath.matrixScalarMultiply(dT/2, accel))
        wx = MatrixMath.matrixSkew(w[0][0], w[1][0], w[2][0])
        wx2 = MatrixMath.matrixMultiply(wx, wx)
        wmag = math.hypot(w[0][0], w[1][0], w[2][0])

        # check if approx for sin and cos term needed
        if wmag <= 0.2:
            sinc = ((dT) - (((dT ** 3) * (wmag ** 2))/ 6) + (((dT ** 5) * (wmag ** 4)) / 120))
            cosc = (((dT ** 2) / 2) + (((dT ** 4) * (wmag ** 2)) / 24) + (((dT ** 6) * (wmag ** 4)) / 720))
        else:
            sinc = (math.sin(wmag * dT)/wmag)
            cosc = ((1 - math.cos(wmag * dT))/(wmag ** 2))

        # add components together
        negsinc = MatrixMath.matrixScalarMultiply(-1, MatrixMath.matrixScalarMultiply(sinc, wx))
        Re1 = MatrixMath.matrixAdd(negsinc, MatrixMath.matrixScalarMultiply(cosc, wx2))
        Re = MatrixMath.matrixAdd(I, Re1)
        return Re

    def Update(self, forcesMoments):
        # integrate states and update internal state
        self.state = self.ForwardEuler(forcesMoments)
        return

    def derivative(self, state, forcesMoments):
        dot = States.vehicleState()
        R = Rotations.euler2DCM(state.yaw, state.pitch, state.roll)

        # derivative of euler angles
        bodyRate = [[state.p],
                    [state.q],
                    [state.r]]
        tEuler = [[1, math.sin(state.roll)*math.tan(state.pitch), math.cos(state.roll)*math.tan(state.pitch)],
               [0, math.cos(state.roll), -math.sin(state.roll)],
               [0, math.sin(state.roll)/math.cos(state.pitch), math.cos(state.roll)/math.cos(state.pitch)]]

        dotEuler = MatrixMath.matrixMultiply(tEuler, bodyRate)
        dot.roll = dotEuler[0][0]
        dot.pitch = dotEuler[1][0]
        dot.yaw = dotEuler[2][0]

        # derivative of positions
        vel = [[state.u],
               [state.v],
               [state.w]]

        dotPos = MatrixMath.matrixMultiply(MatrixMath.matrixTranspose(R), vel)
        dot.pn = dotPos[0][0]
        dot.pe = dotPos[1][0]
        dot.pd = dotPos[2][0]

        # derivative of body rates
        Gamma3 = VPC.Jzz/VPC.Jdet
        Gamma4 = VPC.Jxz/VPC.Jdet
        Gamma5 = (VPC.Jzz - VPC.Jxx)/VPC.Jyy
        Gamma6 = VPC.Jxz/VPC.Jyy
        Gamma8 = VPC.Jxx/VPC.Jdet
        p1 = [[(VPC.Gamma1*state.p*state.q) - (VPC.Gamma2*state.q*state.r)],
              [(Gamma5*state.p*state.r) - (Gamma6*((state.p**2) - (state.r**2)))],
              [(VPC.Gamma7*state.p*state.q) - (VPC.Gamma1*state.q*state.r)]]
        p2 = [[(Gamma3*forcesMoments.Mx) + (Gamma4*forcesMoments.Mz)],
              [forcesMoments.My/VPC.Jyy],
              [(Gamma4*forcesMoments.Mx) + (Gamma8*forcesMoments.Mz)]]

        dotBRate = MatrixMath.matrixAdd(p1, p2)
        dot.p = dotBRate[0][0]
        dot.q = dotBRate[1][0]
        dot.r = dotBRate[2][0]

        # derivative of velocities
        xVel = [[(state.r * state.v) - (state.q * state.w)],
                [(state.p * state.w) - (state.r * state.u)],
                [(state.q * state.u) - (state.p * state.v)]]
        f = [[forcesMoments.Fx],
             [forcesMoments.Fy],
             [forcesMoments.Fz]]
        fm = MatrixMath.matrixScalarMultiply(1 / VPC.mass, f)

        dotVel = MatrixMath.matrixAdd(xVel, fm)
        dot.u = dotVel[0][0]
        dot.v = dotVel[1][0]
        dot.w = dotVel[2][0]

        # derivative of R
        w = MatrixMath.matrixSkew(state.p, state.q, state.r)
        dot.R = MatrixMath.matrixScalarMultiply(-1, MatrixMath.matrixMultiply(w, R))
        return dot

    def getVehicleState(self):
        # return states
        return self.state

    def reset(self):
        # reset states and derivatives
        self.state = States.vehicleState()
        self.dot = States.vehicleState()
        return

    def resetVehicleState(self):
        # reset state
        self.state = States.vehicleState()
        return self.state

    def setVehicleState(self, state):
        # set vehicle state
        self.state = state
        return






