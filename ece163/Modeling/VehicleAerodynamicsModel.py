"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create VehicleAerodynamicsModel.py.
"""

import math


from ..Containers import States
from ..Containers import Inputs
from ..Modeling import VehicleDynamicsModel
from ..Modeling import WindModel
from ..Utilities import MatrixMath
from ..Constants import VehiclePhysicalConstants as VPC


# import ece163.Containers.States as States
# import ece163.Containers.Inputs as Inputs
# import ece163.Modeling.VehicleDynamicsModel as VehicleDynamicsModel
# import ece163.Modeling.WindModel as WindModel
# import ece163.Utilities.MatrixMath as MatrixMath
# import ece163.Constants.VehiclePhysicalConstants as VPC

class VehicleAerodynamicsModel:
    def __init__(self, initialSpeed = 25.0, initialHeight = -100.0):
        self.vehicleDynamics = VehicleDynamicsModel.VehicleDynamicsModel()
        self.wind = States.windState()
        self.windModel = WindModel.WindModel()
        self.vehicleDynamics.forceMoments = Inputs.forcesMoments()
        self.vehicleDynamics.state.u = initialSpeed
        self.vehicleDynamics.state.pd = initialHeight
        self.vehicleDynamics.state.pe = VPC.InitialEastPosition
        self.vehicleDynamics.state.pn = VPC.InitialNorthPosition

    def getVehicleState(self):
        return self.vehicleDynamics.state

    def setVehicleState(self, state):  # state is everything we need
        self.vehicleDynamics.state = state  # to predict where it is going
        return

    def gravityForces(self,state):
        self.gravity = Inputs.forcesMoments()
        Gravity_vector = [[0],
                          [0],
                          [VPC.g0]]
        x = MatrixMath.matrixMultiply(state.R, Gravity_vector)
        F = MatrixMath.matrixScalarMultiply(VPC.mass, x)
        self.gravity.Fx = F[0][0]
        self.gravity.Fy = F[1][0]
        self.gravity.Fz = F[2][0]

        return self.gravity

    def CalculateCoeff_alpha(self, alpha):
        e1 = math.e**(-VPC.M*(alpha-VPC.alpha0))
        e2 = math.e**(VPC.M*(alpha + VPC.alpha0))
        sigma = (1+e1+e2)/((1+e1)*(1+e2))

        CL_alpha = (1-sigma)*(VPC.CL0+(VPC.CLalpha*alpha)) + (sigma)*(2*math.sin(alpha)*math.cos(alpha))
        CD_alpha = (1-sigma)*(VPC.CDp + ((CL_alpha * alpha)**2/(math.pi*VPC.AR*VPC.e))) + (sigma*(2*(math.sin(alpha))**2))
        CM_alpha = VPC.CM0 + VPC.CMalpha*alpha

        return CL_alpha, CD_alpha, CM_alpha

    def aeroForces(self, state):
        self.aero = Inputs.forcesMoments()

        if (state.Va == 0):
            q_bar = 1
            b_bar = 1
        else:
            q_bar = (VPC.c*state.q)/(2*state.Va)
            b_bar = (VPC.b/(2*state.Va))

        p_bar = b_bar*state.p
        r_bar = b_bar*state.r
        CL = self.CalculateCoeff_alpha(state.alpha)
        F_lift = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S))*(CL[0] + VPC.CLq*q_bar)
        F_Drag = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S))*(CL[1] + VPC.CDq*q_bar)

        Force_Matrix = [[-F_Drag],
                        [-F_lift]]

        angle_matrix = [[math.cos(state.alpha), -math.sin(state.alpha)],
                        [math.sin(state.alpha), math.cos(state.alpha)]]

        Foce_body = MatrixMath.matrixMultiply(angle_matrix, Force_Matrix)

        m = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S)*(VPC.c))*(CL[2] + VPC.CMq*q_bar)

        l = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S)*(VPC.b))*((VPC.Clbeta*state.beta) + (VPC.Clp*p_bar) + (VPC.Clr*r_bar))
        n = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S)*(VPC.b))*((VPC.Cnbeta*state.beta) + (VPC.Cnp*p_bar) + (VPC.Cnr*r_bar))

        Fy = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S))*((VPC.CYbeta*state.beta) + (VPC.CYp*q_bar) + (VPC.CYr*r_bar))

        self.aero.Fx = Foce_body[0][0]
        self.aero.Fy = Fy
        self.aero.Fz = Foce_body[1][0]
        self.aero.Mx = l
        self.aero.My = m
        self.aero.Mz = n

        return self.aero

    def controlForces(self, state, controls):
        ctrl = Inputs.forcesMoments()

        F_Lift = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S))*((VPC.CLdeltaE*controls.Elevator))
        F_Drag = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S))*((VPC.CDdeltaE)*(controls.Elevator))
        fy = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S))*(VPC.CYdeltaA*controls.Aileron + VPC.CYdeltaR*controls.Rudder)

        Mx = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S)*(VPC.b))*(VPC.CldeltaA*controls.Aileron + VPC.CldeltaR*controls.Rudder)
        My = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S)*(VPC.c))*(VPC.CMdeltaE*controls.Elevator)
        Mz = ((1/2)*(VPC.rho)*(state.Va**2)*(VPC.S)*(VPC.b))*(VPC.CndeltaA*controls.Aileron + VPC.CndeltaR*controls.Rudder)

        Force_Matrix = [[-F_Drag],
                        [-F_Lift]]

        angle_matrix = [[math.cos(state.alpha), -math.sin(state.alpha)],
                        [math.sin(state.alpha), math.cos(state.alpha)]]

        Foce_body = MatrixMath.matrixMultiply(angle_matrix, Force_Matrix)

        prop_force = self.CalculatePropForces(state.Va, controls.Throttle)

        ctrl.Fx = Foce_body[0][0] + prop_force[0]
        ctrl.Fz = Foce_body[1][0]
        ctrl.Fy = fy
        ctrl.Mx = Mx + prop_force[1]
        ctrl.My = My
        ctrl.Mz = Mz

        return ctrl

    def CalculateAirspeed(self, state, wind):

        if ((wind.Wn == 0) and (wind.We == 0) and (wind.Wd == 0)):
            azimuth_w = 0
            elevation_w = 0
        else:
            azimuth_w = math.atan2(wind.We, wind.Wn)
            elevation_w = -math.asin(wind.Wd/math.hypot(wind.Wn, wind.We, wind.Wd))

        AE_matrix = [[math.cos(azimuth_w)*math.cos(elevation_w), math.sin(azimuth_w)*math.cos(elevation_w), -math.sin(elevation_w)],
                     [-math.sin(azimuth_w), math.cos(azimuth_w), 0],
                     [math.cos(azimuth_w)*math.sin(elevation_w), math.sin(azimuth_w)*math.sin(elevation_w), math.cos(elevation_w)]]

        wind_ned = [[wind.Wn], [wind.We], [wind.Wd]]
        wind_uvw = [[wind.Wu], [wind.Wv], [wind.Ww]]

        win_mat = MatrixMath.matrixMultiply(MatrixMath.matrixTranspose(AE_matrix), wind_uvw)

        w1_matrix = MatrixMath.matrixAdd(wind_ned, win_mat)

        Wi_matrix = MatrixMath.matrixMultiply(state.R, w1_matrix)

        ur = state.u - Wi_matrix[0][0]
        vr = state.v - Wi_matrix[1][0]
        wr = state.w - Wi_matrix[2][0]

        Va = math.hypot(ur, vr, wr)
        alpha = math.atan2(wr,ur)
        beta = math.asin(vr/Va)

        return Va, alpha, beta

    def CalculatePropForces(self, Va, Throttle):
        Kt = ((60)/(2*math.pi*VPC.KV))
        Vin = (VPC.V_max)*(Throttle)

        a = ((VPC.rho)*(VPC.D_prop**5)*(VPC.C_Q0))/(4*(math.pi**2))
        b = (((VPC.rho)*(VPC.D_prop**4)*(Va)*(VPC.C_Q1))/(2*math.pi) + ((Kt*Kt)/VPC.R_motor))
        c = ((VPC.rho)*(VPC.D_prop**3)*(Va**2)*(VPC.C_Q2) - (Kt)*(Vin/VPC.R_motor) + Kt*VPC.i0)

        try:
            Om = ((-b + math.sqrt(b**2 - 4*a*c))/(2*a))
        except:
            Om = 100

        J = ((2*math.pi*Va)/(Om*VPC.D_prop))

        CT = ((VPC.C_T0) + (VPC.C_T1*J) + (VPC.C_T2*(J**2)))
        CQ = ((VPC.C_Q0) + (VPC.C_Q1*J) + (VPC.C_Q2*(J**2)))

        Fx_prop = ((VPC.rho)*(Om**2)*(VPC.D_prop**4)*CT)/(4*(math.pi**2))
        Mx_prop = -((VPC.rho)*(Om**2)*(VPC.D_prop**5)*CQ)/(4*(math.pi**2))

        result = [Fx_prop, Mx_prop]

        return result

    def setWindModel(self, Wn=0, We=0, Wd=0, drydenParameters=VPC.DrydenNoWind):
        self.windModel.Wind.Wn = Wn
        self.windModel.Wind.We = We
        self.windModel.Wind.Wd = Wd
        self.windModel.drydenParameters = drydenParameters

    def getWindState(self):
        return self.windModel.Wind

    def updateForces(self, state, wind, controls):
        airSpeed = self.CalculateAirspeed(state, wind)

        state.Va = airSpeed[0]
        state.alpha = airSpeed[1]
        state.beta = airSpeed[2]

        G_force = self.gravityForces(state)
        Aero_force = self.aeroForces(state)
        cntrl_force = self.controlForces(state, controls)

        self.vehicleDynamics.forceMoments.Fx = G_force.Fx + Aero_force.Fx + cntrl_force.Fx
        self.vehicleDynamics.forceMoments.Fy = G_force.Fy + Aero_force.Fy + cntrl_force.Fy
        self.vehicleDynamics.forceMoments.Fz = G_force.Fz + Aero_force.Fz + cntrl_force.Fz
        self.vehicleDynamics.forceMoments.Mx = G_force.Mx + Aero_force.Mx + cntrl_force.Mx
        self.vehicleDynamics.forceMoments.My = G_force.My + Aero_force.My + cntrl_force.My
        self.vehicleDynamics.forceMoments.Mz = G_force.Mz + Aero_force.Mz + cntrl_force.Mz

        return self.vehicleDynamics.forceMoments

    def Update(self, controls):
        forces = self.updateForces(self.vehicleDynamics.state, self.windModel.Wind, controls)
        VehicleDynamicsModel.VehicleDynamicsModel.Update(self.vehicleDynamics, forces)

    def reset(self):
        VehicleAerodynamicsModel.__init__(self, initialSpeed=VPC.InitialSpeed, initialHeight=VPC.InitialDownPosition)