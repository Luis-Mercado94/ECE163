"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create VehicleDynamicModel.py.
"""

import math


import ece163.Containers.States as States
import ece163.Utilities.MatrixMath as MatrixMath
import ece163.Utilities.Rotations as Rotations
import ece163.Constants.VehiclePhysicalConstants as VPC
from ..Containers import States
from ..Utilities import MatrixMath
from ..Utilities import Rotations
from ..Constants import VehiclePhysicalConstants as VPC


class VehicleDynamicsModel:
	def __init__(self, dT = VPC.dT):
		self.state = States.vehicleState()
		self.dot = States.vehicleState()
		self.dT = dT

	def getVehicleState(self):
		return self.state

	def setVehicleState(self, state): # state is everything we need
		self.state = state		# to predict where it is going
		return

	def resetVehicleState(self):
		resetState = States.vehicleState()
		self.setVehicleState(resetState)
		return
		# self.state.pn = 0.0
		# self.state.pe = 0.0
		# self.state.pd = 0.0
		# # velocities
		# self.state.u = 0.0
		# self.state.v = 0.0
		# self.state.w = 0.0
		# # Euler Angles
		# self.state.yaw = 0.0
		# self.state.pitch = 0.0
		# self.state.roll = 0.0
		#
		# self.state.R = Rotations.euler2DCM(self.state.yaw, self.state.pitch, self.state.roll)
		#
		# # body rates
		# self.state.p = 0.0
		# self.state.q = 0.0
		# self.state.r = 0.0
		return

	def derivative(self, state, forcesMoments):

		dot = self.dot

		state.R = Rotations.euler2DCM(state.yaw, state.pitch, state.roll)

		Wx = MatrixMath.matrixSkew(state.p,state.q,state.r)

		dot.R = MatrixMath.matrixScalarMultiply(-1, MatrixMath.matrixMultiply(Wx, state.R))

		R_tran = MatrixMath.matrixTranspose(state.R)
		vel = [[state.u],[state.v],[state.w]]
		ned = MatrixMath.matrixMultiply(R_tran,vel)

		dot.pn = ned[0][0]
		dot.pe = ned[1][0]
		dot.pd = ned [2][0]

		bodymatrix = [[1, (math.sin(state.roll))*(math.tan(state.pitch)), (math.cos(state.roll))*(math.tan(state.pitch))],
					  [0, math.cos(state.roll), -(math.sin(state.roll))],
					  [0, (math.sin(state.roll))/(math.cos(state.pitch)), (math.cos(state.roll))/(math.cos(state.pitch))]]

		forcematrix = [[forcesMoments.Fx],
					   [forcesMoments.Fy],
					   [forcesMoments.Fz]]

		ratematrix = [[state.p],
					  [state.q],
					  [state.r]]

		ypr = MatrixMath.matrixMultiply(bodymatrix, ratematrix)

		dot.roll = ypr[0][0]
		dot.pitch = ypr[1][0]
		dot.yaw = ypr[2][0]

		uvw = [[(state.r*state.v) - (state.q*state.w)],
			   [(state.p*state.w) - (state.r*state.u)],
			   [(state.q*state.u) - (state.p*state.v)]]

		massForce = MatrixMath.matrixScalarMultiply(1/VPC.mass, forcematrix)

		Duvw = MatrixMath.matrixAdd(uvw, massForce)

		dot.u = Duvw[0][0]
		dot.v = Duvw[1][0]
		dot.w = Duvw[2][0]

		Gamma3 = VPC.Jzz/VPC.Jdet
		Gamma4 = VPC.Jxz/VPC.Jdet
		Gamma5 = (VPC.Jzz-VPC.Jxx)/VPC.Jyy
		Gamma6 = VPC.Jxz/VPC.Jyy
		Gamma8 = VPC.Jxx/VPC.Jdet

		# m1 = [[VPC.Jxz*state.p*state.q + (VPC.Jyy-VPC.Jzz)*state.q*state.r],
		# 	  [VPC.Jxz*(state.r**2 - state.p**2) + (VPC.Jzz - VPC.Jxx)*state.p*state.r],
		# 	  [(VPC.Jxx - VPC.Jyy)*state.p*state.q - VPC.Jxz*state.q*state.r]]
		#
		# m2 = [[forcesMoments.Mx],
		# 	  [forcesMoments.My],
		# 	  [forcesMoments.Mz]]
		#
		# m3 = MatrixMath.matrixAdd(m1, m2)
		#
		# ERmatrix = MatrixMath.matrixMultiply(VPC.JinvBody, m3)

		m1 = [[VPC.Gamma1*state.p*state.q - VPC.Gamma2*state.q*state.r],
			  [Gamma5*state.p*state.r - Gamma6*(state.p**2 - state.r**2)],
			  [VPC.Gamma7*state.p*state.q - VPC.Gamma1*state.q*state.r]]

		m2 = [[Gamma3*forcesMoments.Mx + Gamma4*forcesMoments.Mz],
			  [forcesMoments.My*(1/VPC.Jyy)],
			  [Gamma4*forcesMoments.Mx + Gamma8*forcesMoments.Mz]]

		ERmatrix = MatrixMath.matrixAdd(m1, m2)

		dot.p = ERmatrix[0][0]
		dot.q = ERmatrix[1][0]
		dot.r = ERmatrix[2][0]

		return dot

	def reset(self):
		VehicleDynamicsModel.__init__(self)

		return

	def Rexp(self, dT, state, dot):

		pk = state.p + (dot.p*(dT/2))
		qk = state.q + (dot.q*(dT/2))
		rk = state.r + (dot.r*(dT/2))

		w = math.hypot(pk, qk, rk)
		wx = MatrixMath.matrixSkew(pk, qk, rk)

		I = [[1,0,0],
			 [0,1,0],
			 [0,0,1]]

		if w <= 0.2:
			t1 = dT - ((dT ** 3) * (w ** 2)) / 6 + ((dT ** 5) * w ** 4) / 120
			t2 = (dT ** 2) / 2 - ((dT ** 4) * (w ** 4)) / 24 + ((dT ** 6) * (w ** 4)) / 720
			p1 = MatrixMath.matrixScalarMultiply(t1, wx)
			p2 = MatrixMath.matrixScalarMultiply(t2, MatrixMath.matrixMultiply(wx, wx))
		else:
			p1 = MatrixMath.matrixScalarMultiply((math.sin(w * dT) / w), wx)
			p2 = MatrixMath.matrixScalarMultiply((1 - math.cos(w * dT)) / w**2, MatrixMath.matrixMultiply(wx, wx))

		p3 = MatrixMath.matrixSubtract(I, p1)

		Rexp = MatrixMath.matrixAdd(p3, p2)

		return Rexp

	def IntegrateState(self, dT, state, dot):
		# newState = States.vehicleState()
		Rexp = VehicleDynamicsModel.Rexp(self, dT, state, dot)

		state.R = MatrixMath.matrixMultiply(Rexp, state.R)

		DCM = Rotations.dcm2Euler(state.R)
		# state.dcm = DCM

		state.yaw = DCM[0]
		state.pitch = DCM[1]
		state.roll = DCM[2]

		state.pe = state.pe + (dot.pe * dT)
		state.pn = state.pn + (dot.pn * dT)
		state.pd = state.pd + (dot.pd * dT)
		state.u = state.u + (dot.u * dT)
		state.v = state.v + (dot.v * dT)
		state.w = state.w + (dot.w * dT)
		state.p = state.p + (dot.p * dT)
		state.q = state.q + (dot.q * dT)
		state.r = state.r + (dot.r * dT)
		state.chi = math.atan2(dot.pe, dot.pn)

		newState = state

		return newState

	def ForwardEuler(self, forcesMoments):
		state = self.getVehicleState()
		dot = VehicleDynamicsModel.derivative(self, state, forcesMoments)
		
		newState = VehicleDynamicsModel.IntegrateState(self, self.dT, state, dot)

		return newState

	def Update(self, forcesMoments):
		self.state = self.ForwardEuler(forcesMoments)
		return
