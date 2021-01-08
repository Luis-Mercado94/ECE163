"""
Luis Mercado, Pauline Mae F. Cuaresma, and Sabrina Fong
ECE 163 Final Project
Due 12/17/20

Added to this module for the final project is VOupdate, distance, avoidance. VOupdate is used to determine whether the
UAV will crash into the obstacle. Distance determines the distance between the UAV and obstacle. Avoidance determines
the direction of the avoidance (up or down).
Lines 260, 269 - 332 contain the additions to this module for the final project.
"""

import math
import random

"""
from ..Containers import States
from ..Containers import Inputs
from ..Modeling import VehicleDynamicsModel
from ..Modeling import WindModel
from ..Utilities import MatrixMath
from ..Constants import VehiclePhysicalConstants as VPC
"""

""" Imports """
import ece163.Containers.States as States
import ece163.Containers.Inputs as Inputs
import ece163.Modeling.VehicleDynamicsModel as VehicleDynamicsModel
import ece163.Modeling.WindModel as WindModel
import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Constants.VehiclePhysicalConstants as VPC

class VehicleAerodynamicsModel:
    def __init__(self,initialSpeed = 25.0, initialHeight = -100.0):
        self.vehicleDynamics = VehicleDynamicsModel.VehicleDynamicsModel()   # contains states
        # self.wind = States.windState()                # current wind states
        self.windModel = WindModel.WindModel()
        self.vehicleDynamics.state.u = initialSpeed     # initial speed along x-axis; from rudder
        self.vehicleDynamics.state.pd = initialHeight
        self.vehicleDynamics.state.pn = VPC.InitialNorthPosition
        self.vehicleDynamics.state.pe = VPC.InitialEastPosition
        """ FINAL PROJECT ADDITIONS """
        self.uav2 = States.vehicleState(99.2,0,-100,0,0,0,0,0,0,0,0,0)
        self.crash = False
        self.preDir = None
        self.Dir = None
        return


    def getVehicleState(self):
        return self.vehicleDynamics.state


    def setVehicleState(self, state):
        self.vehicleDynamics.state = state
        return


    def gravityForces(self, state):
        self.gravForce = Inputs.forcesMoments()
        """ note: Fz = m*R_ib*grav_vector """
        g_vector = [[0], [0], [VPC.g0 *VPC.mass]]                    # gravity vector
        force_mat = MatrixMath.matrixMultiply(state.R, g_vector)     # updating Fz

        """ Extracting forces """
        self.gravForce.Fx = force_mat[0][0]
        self.gravForce.Fy = force_mat[1][0]
        self.gravForce.Fz = force_mat[2][0]
        return self.gravForce


    def CalculateCoeff_alpha(self, alpha):
        """ Sigma Eqn """
        c1 = math.e ** (-VPC.M * (alpha - VPC.alpha0))
        c2 = math.e ** (VPC.M * (alpha + VPC.alpha0))
        sigma = (1 + c1 + c2) / ((1 + c1) * (1 + c2))

        """ Coefficients """
        CL_alpha = ((1 - sigma) * (VPC.CL0 + (VPC.CLalpha * alpha))) + ((sigma) * (2 * math.sin(alpha) * math.cos(alpha)))
        CD_alpha = ((1 - sigma) * (VPC.CDp + ((CL_alpha * alpha) ** 2 / (math.pi * VPC.AR * VPC.e)))) + (sigma * (2 * (math.sin(alpha)) ** 2))
        CM_alpha = VPC.CM0 + (VPC.CMalpha * alpha)
        return CL_alpha, CD_alpha, CM_alpha


    def CalculateAirspeed(self, state, wind):
        if((wind.Wn == 0) and (wind.We == 0) and (wind.Wd == 0)):
            Xw = 0
            Yw = 0
        else:
            Xw = math.atan2(wind.We, wind.Wn)  # "chi_wind" simlar to alpha
            Yw = -math.asin(wind.Wd / math.hypot(wind.Wn, wind.We, wind.Wd))


        """ matrices """
        W_ned = [[wind.Wn], [wind.We], [wind.Wd]]  # in inertial
        W_uvw = [[wind.Wu], [wind.Wv], [wind.Ww]]  # in Xw, "chi_wind", frame

        R_XYw = [[math.cos(Xw) * math.cos(Yw), math.sin(Xw) * math.cos(Yw), -math.sin(Yw)],
                [-math.sin(Xw), math.cos(Xw), 0],
                [math.cos(Xw) * math.sin(Yw), math.sin(Xw) * math.sin(Yw), math.cos(Yw)]]


        """ matrices and matrix manipulations """
        W_ig =MatrixMath.matrixMultiply(MatrixMath.matrixTranspose(R_XYw), W_uvw) # = (R_Xw)^T * Wg

        W_b = MatrixMath.matrixMultiply(state.R, MatrixMath.matrixAdd(W_ned, W_ig)) #wind total in body frame -- used to calculate alpha, beta

        """ Get Va, alpha, beta """
        states = [[state.u], [state.v], [state.w]]
        Va_mat = MatrixMath.matrixSubtract(states, W_b)

        """ alpha, beta, Va """
        alpha = math.atan2(Va_mat[2][0], Va_mat[0][0])  # angle of attack
        Va = math.hypot(Va_mat[0][0], Va_mat[1][0], Va_mat[2][0])  # Airspeed
        if (math.isclose(Va,0)):    #consider corner case
            beta = 0
        else:
            beta = math.asin(Va_mat[1][0] / Va)
        return Va, alpha, beta



    def aeroForces(self, state):
        self.aeroForce = Inputs.forcesMoments()

        constant = 0.5 * VPC.rho * (state.Va ** 2) * VPC.S  # for fx, fz, My
        trig_mat = [[math.cos(state.alpha), - math.sin(state.alpha)],
                    [math.sin(state.alpha), math.cos(state.alpha)]]

        CalcCoeffs = VehicleAerodynamicsModel.CalculateCoeff_alpha(self, state.alpha)  # get coeffs

        """ For corner case if Va = 0"""
        if (state.Va == 0):
            q_bar = 1
            b_norm = 1
        else:
            q_bar = (VPC.c * state.q)/(2 * state.Va) # to normalize q
            b_norm = VPC.b/ (2 * state.Va) # normalize with span

        """ Forces & Moments """
        F_lift = constant * (CalcCoeffs[0] + (VPC.CLq * q_bar))
        F_drag = constant * (CalcCoeffs[1] + (VPC.CDq * q_bar))
        M_g = (constant * VPC.c) * (CalcCoeffs[2] + (VPC.CMq * q_bar))

        """ For moments """
        self.aeroForce.Mx = (constant * VPC.b) * ((VPC.Clbeta * state.beta) + (VPC.Clp * state.p * b_norm) + (VPC.Clr * state.r * b_norm) )
        self.aeroForce.My = M_g
        self.aeroForce.Mz = (constant * VPC.b) * ((VPC.Cnbeta * state.beta) + (VPC.Cnp * state.p * b_norm) + (VPC.Cnr * state.r * b_norm) )

        """ Force matrix """
        forces = [[-F_drag],
                  [-F_lift]]

        F_bod = MatrixMath.matrixMultiply(trig_mat, forces)
        self.aeroForce.Fx = F_bod[0][0]
        self.aeroForce.Fy = constant * ((VPC.CYbeta * state.beta) + (VPC.CYp * b_norm * state.p) + (VPC.CYr * b_norm * state.r))
        self.aeroForce.Fz = F_bod[1][0]
        return self.aeroForce



    def controlForces(self, state, control):
        self.ctrlForce = Inputs.forcesMoments()
        constant = 0.5 * VPC.rho * state.Va ** 2 * VPC.S
        trig_mat = [[math.cos(state.alpha), - math.sin(state.alpha)],
                    [math.sin(state.alpha), math.cos(state.alpha)]]

        "Obtain Propeller forces"
        propForces = VehicleAerodynamicsModel.CalculatePropForces(self, state.Va, control.Throttle)

        """ For forces """
        F_lift = constant * (VPC.CLdeltaE * control.Elevator)
        F_drag = constant * (VPC.CDdeltaE * control.Elevator)
        forces = [[-F_drag],
                  [-F_lift]]
        F_bod = MatrixMath.matrixMultiply(trig_mat, forces)

        self.ctrlForce.Fx = F_bod[0][0] + propForces[0]
        self.ctrlForce.Fy = (constant * ((VPC.CYdeltaA * control.Aileron) + (VPC.CYdeltaR * control.Rudder)))
        self.ctrlForce.Fz = F_bod[1][0]

        """ For moments """
        self.ctrlForce.Mx = (constant * VPC.b) * ((VPC.CldeltaA * control.Aileron) + (VPC.CldeltaR * control.Rudder)) + propForces[1]
        self.ctrlForce.My = (constant * VPC.c) * (VPC.CMdeltaE * control.Elevator)
        self.ctrlForce.Mz = (constant * VPC.b) * ((VPC.CndeltaA * control.Aileron) + (VPC.CndeltaR * control.Rudder))
        return self.ctrlForce



    def CalculatePropForces(self, Va, Throttle):
        """ KT = KE equation from Prop Cheat Sheet """
        KT = 60/(2*math.pi*VPC.KV)

        """ V_in eqn from Prop Cheat Sheet """
        V_in = VPC.V_max * Throttle

        """ a,b,c constants for Omega & Omega """
        a = (VPC.rho * VPC.D_prop**5 * VPC.C_Q0)/ (4*math.pi**2)
        b = (VPC.rho * VPC.D_prop**4 * Va * VPC.C_Q1)/(2*math.pi) + (KT**2)/VPC.R_motor
        c = (VPC.rho * VPC.D_prop**3 * Va**2 * VPC.C_Q2) - (KT * (V_in/VPC.R_motor)) + (KT * VPC.i0)
        quad = b**2 - (4*a*c)
        if (quad < 0):
            Omega = 100
        else:
            Omega = (-b + math.sqrt(quad)) / (2*a) #alternative if we're taking sqrt of negative

        """ Calculating J, the advance ratio """
        J = (2 * math.pi * Va)/ (Omega * VPC.D_prop)

        """ C_T & and C_Q calcs """
        C_T = VPC.C_T0 + (VPC.C_T1 * J) + (VPC.C_T2 * J**2)
        C_Q = VPC.C_Q0 + (VPC.C_Q1 * J) + (VPC.C_Q2 * J**2)

        """ Fx_prop & Mx_prop """
        Fx = (VPC.rho * Omega**2 * VPC.D_prop**4 * C_T) / (4* math.pi ** 2)
        Mx = (VPC.rho * Omega**2 * VPC.D_prop**5 * C_Q) / (4* math.pi ** 2)
        return [Fx, -Mx]



    def setWindModel(self, Wn=0.0, We=0.0, Wd=0.0, drydenParameters = Inputs.drydenParameters()):
        # self.windModel = States.windState(Wn, We, Wd, 0, 0, 0)  #note: for now no wind gusts added
        self.windModel.Wind.Wn = Wn
        self.windModel.Wind.We = We
        self.windModel.Wind.Wd = Wd
        self.windModel.drydenParameters = drydenParameters
        return


    def getWindState(self):
        return self.windModel.Wind


    def updateForces(self, state, wind, controls):
        airSp = self.CalculateAirspeed(state, wind) #outputs Va, alpha, beta

        """ Extract alpha, beta, and Va. Then update states """
        state.Va = airSp[0]
        state.alpha = airSp[1]
        state.beta = airSp[2]

        """ Calling required functions to obtain all forces and moments """
        grav = self.gravityForces(state)
        aero = self.aeroForces(state)
        ctrl = self.controlForces(state, controls)

        """ Total forces & moments """
        forcesMoments = Inputs.forcesMoments()
        forcesMoments.Fx = grav.Fx + aero.Fx + ctrl.Fx
        forcesMoments.Fy = grav.Fy + aero.Fy + ctrl.Fy
        forcesMoments.Fz = grav.Fz + aero.Fz + ctrl.Fz
        forcesMoments.Mx = grav.Mx + aero.Mx + ctrl.Mx
        forcesMoments.My = grav.My + aero.My + ctrl.My
        forcesMoments.Mz = grav.Mz + aero.Mz + ctrl.Mz
        return forcesMoments


    def Update(self,controls):
        self.windModel.Update()
        self.vehicleDynamics.Update(self.updateForces(self.vehicleDynamics.state, self.windModel.Wind, controls))
        self.VOupdate(self.uav2)
        return


    def reset(self):
        VehicleAerodynamicsModel.__init__(self, initialSpeed = VPC.InitialSpeed, initialHeight = VPC.InitialDownPosition)
        return


    """ IMPLEMENTATION OF FINAL PROJECT """

    def distance(self, uav2):
        # call update/ make sure the vehicle state is the most up to date
        x = self.vehicleDynamics.state.pn - uav2.pn
        y = self.vehicleDynamics.state.pe - uav2.pe
        z = self.vehicleDynamics.state.pd - uav2.pd

        # check if UAV passed obstacle; if passed, return value greater than d_avo to indicate it is a safe distance
        # away from obstacle
        if x < 0:
            distance = math.hypot(x, y, z)
        else:
            safe = 51 # larger than d_avo
            distance = safe
        return distance


    def VOupdate(self, uav2):
        """
        calculating VO cone parameters and determining whether there will be a collision to avoid it
        :param uav2:
        :return:
        """
        d_oi = self.distance(uav2) #sdistance of uav to obstacle
        if d_oi == 0: # to avoid divide by zero error
            d_oi = 1
        r_pz = VPC.b # setting radius of obstacle to the wing span of UAV
        d_avo = 50 # distance of when to start avoiding [m]

        # find the difference in positions between the UAV and obstacle
        z_dif = abs(uav2.pd - self.vehicleDynamics.state.pd)
        y_dif = abs(uav2.pe - self.vehicleDynamics.state.pe)
        x_dif = abs(uav2.pn - self.vehicleDynamics.state.pn)

        # determine angles of the velocity obstacle cone
        theta_oi = math.asin(z_dif/d_oi)
        phi_oi = math.atan2(y_dif, x_dif)

        d_vo = (d_oi**2 - r_pz**2)/d_oi # length of VO cone
        r_vo = r_pz*((math.sqrt(abs(d_oi**2 - r_pz**2)))/d_oi) # radius of effective base of VO cone
        a_vo = math.atan2(r_vo, d_vo) # opening angle of VO cone

        mx = [[math.cos(theta_oi)*math.cos(phi_oi)],
              [math.cos(theta_oi)*math.sin(phi_oi)],
              [math.sin(theta_oi)]]

        # speed of UAV
        vo = [[self.vehicleDynamics.state.u, self.vehicleDynamics.state.v, self.vehicleDynamics.state.w]]

        D_vo = MatrixMath.matrixScalarMultiply(d_vo, mx) # VO cone symmetric axis vector
        A_vo = [[0,0,0]] # Position of the apex. it is all zero since the obstacle is static
        dif_mat = MatrixMath.matrixSubtract(vo, A_vo)
        tran = MatrixMath.matrixTranspose(dif_mat)

        check1 = MatrixMath.matrixDotProduct(tran, D_vo)
        check2 = abs(self.vehicleDynamics.state.Va)*d_vo
        # check angle boundaries and if the UAV is at a critical avoidance distance away from the obstacle
        if (d_oi < d_avo) and (check1[0][0]/check2 > math.cos(a_vo)):
            self.crash = True
        else:
            self.crash = False
        return self.crash

    def avoidance(self):
        """
           avoidance() takes in no parameters and returns the 0 or 1 by random. 0 refers to the UAV going down. 1 refers
           to the UAV going up.
        """
        # 1 = up, 0 = down
        w = random.randint(0, 1)
        return w

