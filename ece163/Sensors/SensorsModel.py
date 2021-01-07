"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create SensorModel.py that implement the functions to simulate actual sensor outputs
"""

import math
import random
from ece163.Modeling import VehicleAerodynamicsModel
from ece163.Utilities import MatrixMath
from ..Containers import Sensors
from ..Constants import VehiclePhysicalConstants as VPC
from ..Constants import VehicleSensorConstants as VSC

class GaussMarkov():
    def __init__(self, dT = VPC.dT, tau = 1e6, eta = 0.0):
        self.dT = dT
        self.tau = 1e6
        self.tau = tau
        self.eta = 0.0
        self.eta = eta
        self.v = 0
        return

    def reset(self):
        self.dT = VPC.dT
        self.tau = 1e6
        self.eta = 0.0
        self.v = 0

        return

    def update(self, vnoise = None):
        if vnoise == None or vnoise == 0:
            vnoise = random.gauss(0, self.eta)
        self.v = (self.v * (math.exp(-self.dT/self.tau))) + vnoise
        return self.v

class GaussMarkovXYZ():
    def __init__(self, dT = VPC.dT, tauX = 1e6, etaX = 0.0, tauY = None, etaY = None,
                 tauZ = None, etaZ = None):
        self.gaussX = GaussMarkov(dT, tauX, etaX)
        #if blocks to ceck if Y and Z are none to set them to the X value
        if tauY == None:
            tauY = tauX
        if etaY == None:
            etaY = etaX
        if tauZ == None:
            tauZ = tauX
        if etaZ == None:
            etaZ = etaX

        self.gaussY = GaussMarkov(dT, tauY, etaY)
        self.gaussZ = GaussMarkov(dT, tauZ, etaZ)
        return

    def reset(self):
        self.gaussX.reset()
        self.gaussY.reset()
        self.gaussZ.reset()

        return

    def update(self, vXnoise = None, vYnoise = None, vZnoise = None):

        gmX = self.gaussX.update(vXnoise)
        gmY = self.gaussY.update(vYnoise)
        gmZ = self.gaussZ.update(vZnoise)

        return gmX, gmY, gmZ

class SensorsModel():
    def __init__(self, areoModel = VehicleAerodynamicsModel.VehicleAerodynamicsModel(), taugyro = VSC.gyro_tau,
                 etagyro = VSC.gyro_eta, tauGPS = VSC.GPS_tau, etaGPSHorizontal = VSC.GPS_etaHorizontal,
                 etaGPSVertical = VSC.GPS_etaVertical, gpsUpdateHz = VSC.GPS_rate):
        self.sTrue = Sensors.vehicleSensors()
        self.sNoisy = Sensors.vehicleSensors()
        self.biases = self.initializeBiases()
        self.sigma = self.initializeSigmas()
        self.gyroGM = GaussMarkovXYZ(VPC.dT, taugyro, etagyro)
        self.gpsGM = GaussMarkovXYZ(1/gpsUpdateHz, tauGPS, etaGPSHorizontal, tauGPS, etaGPSHorizontal,
                                    tauGPS, etaGPSVertical)
        self.dt = VPC.dT
        self.UpdateTicks = 0
        self.gpsTicksUpdate = (1/VPC.dT)*VSC.GPS_rate
        self.VAM = areoModel

        return
    def getSensorsNoisy(self):
        return self.sNoisy

    def getSensorsTrue(self):
        return self.sTrue

    def initializeBiases(self, gyroBias=VSC.gyro_bias, accelBias=VSC.accel_bias, magBias=VSC.accel_bias,
                         baroBias=VSC.baro_bias, pitotBias=VSC.pitot_bias):
        sensorBiases = Sensors.vehicleSensors()

        sensorBiases.gyro_x = gyroBias * random.uniform(-1, 1)
        sensorBiases.gyro_y = gyroBias * random.uniform(-1, 1)
        sensorBiases.gyro_z = gyroBias * random.uniform(-1, 1)

        sensorBiases.accel_x = accelBias * random.uniform(-1, 1)
        sensorBiases.accel_y = accelBias * random.uniform(-1, 1)
        sensorBiases.accel_z = accelBias * random.uniform(-1, 1)

        sensorBiases.mag_x = magBias * random.uniform(-1, 1)
        sensorBiases.mag_y = magBias * random.uniform(-1, 1)
        sensorBiases.mag_z = magBias * random.uniform(-1, 1)

        sensorBiases.baro = baroBias * random.uniform(-1, 1)
        sensorBiases.pitot = pitotBias * random.uniform(-1, 1)

        sensorBiases.gps_n = 0.0
        sensorBiases.gps_e = 0.0
        sensorBiases.gps_alt = 0.0
        sensorBiases.gps_sog = 0.0
        sensorBiases.gps_cog = 0.0

        return sensorBiases

    def initializeSigmas(self, gyroSigma=VSC.gyro_sigma, accelSigma=VSC.accel_sigma, magSigma=VSC.mag_sigma,
                         baroSigma=VSC.baro_sigma, pitotSigma=VSC.pitot_sigma,
                         gpsSigmaHorizontal=VSC.GPS_sigmaHorizontal, gpsSigmaVertical=VSC.GPS_sigmaVertical,
                         gpsSigmaSOG=VSC.GPS_sigmaSOG, gpsSigmaCOG=VSC.GPS_sigmaCOG):

        senorSigma = Sensors.vehicleSensors()

        senorSigma.gyro_x = gyroSigma
        senorSigma.gyro_y = gyroSigma
        senorSigma.gyro_z = gyroSigma

        senorSigma.accel_x = accelSigma
        senorSigma.accel_y = accelSigma
        senorSigma.accel_z = accelSigma

        senorSigma.mag_x = magSigma
        senorSigma.mag_y = magSigma
        senorSigma.mag_z = magSigma

        senorSigma.baro = baroSigma
        senorSigma.pitot = pitotSigma

        senorSigma.gps_n = gpsSigmaHorizontal
        senorSigma.gps_e = gpsSigmaHorizontal
        senorSigma.gps_alt = gpsSigmaVertical
        senorSigma.gps_sog = gpsSigmaSOG
        senorSigma.gps_cog = gpsSigmaCOG

        return senorSigma

    def reset(self):
        self.__init__()
        return

    def update(self):
        self.sTrue = self.UpdateSensorsTrue(self.sTrue, self.VAM.vehicleDynamics.state, self.VAM.vehicleDynamics.dot)
        self.sNoisy = self.updateSensorsNoisy(self.sTrue, self.sNoisy, self.biases, self.sigma)
        self.UpdateTicks += 1

        return

    def updateAccelsTrue(self, state, dot):

        accel_x = dot.u + state.q*state.w - state.r*state.v + VPC.g0*math.sin(state.pitch)
        accel_y = dot.v + (state.r*state.u) - (state.p*state.w) - (VPC.g0*math.cos(state.pitch)*math.sin(state.roll))
        accel_z = dot.w + (state.p*state.v) - (state.q*state.u) - (VPC.g0*math.cos(state.pitch)*math.cos(state.roll))

        return accel_x, accel_y, accel_z

    def updateGyrosTrue(self, state):

        gyro_x = state.p
        gyro_y = state.q
        gyro_z = state.r

        return gyro_x,gyro_y, gyro_z

    def updateMagsTrue(self, state):

        mags = MatrixMath.matrixMultiply(state.R, VSC.magfield)

        return mags[0][0], mags[1][0], mags[2][0]

    def updatePressureSensorsTrue(self, state):

        baro = -(VPC.rho*VPC.g0*-state.pd) + VSC.Pground
        pitot = (1/2)*(VPC.rho)*(state.Va**2)

        return baro, pitot

    def updateGPSTrue(self, state, dot):

        gps_n = state.pn
        gps_e = state.pe
        gps_alt = -state.pd
        gps_sog = math.hypot(state.u, state.v, state.w)
        gps_cog = math.atan2(dot.pe, dot.pn)

        return gps_n, gps_e, gps_alt, gps_sog, gps_cog

    def updateSensorsNoisy(self, trueSensors=Sensors.vehicleSensors(), noisySensors=Sensors.vehicleSensors(),
                           sensorBiases=Sensors.vehicleSensors(), sensorSigmas=Sensors.vehicleSensors()):
        sn = Sensors.vehicleSensors()

        gbX, gbY, gbZ = self.gyroGM.update() #Gauss Markov process
        #Gyro Noisy functions
        sn.gyro_x = trueSensors.gyro_x + sensorBiases.gyro_x + gbX + random.gauss(0, sensorSigmas.gyro_x)
        sn.gyro_y = trueSensors.gyro_y + sensorBiases.gyro_y + gbY + random.gauss(0, sensorSigmas.gyro_y)
        sn.gyro_z = trueSensors.gyro_z + sensorBiases.gyro_z + gbZ + random.gauss(0, sensorSigmas.gyro_z)
        #Accelerometer noisy functions
        sn.accel_x = trueSensors.accel_x + sensorBiases.accel_x + random.gauss(0, sensorSigmas.accel_x)
        sn.accel_y = trueSensors.accel_y + sensorBiases.accel_y + random.gauss(0, sensorSigmas.accel_y)
        sn.accel_z = trueSensors.accel_z + sensorBiases.accel_z + random.gauss(0, sensorSigmas.accel_z)
        #Magnetometer noisy functions
        sn.mag_x = trueSensors.mag_x + sensorBiases.mag_x + random.gauss(0, sensorSigmas.mag_x)
        sn.mag_y = trueSensors.mag_y + sensorBiases.mag_y + random.gauss(0, sensorSigmas.mag_y)
        sn.mag_z = trueSensors.mag_z + sensorBiases.mag_z + random.gauss(0, sensorSigmas.mag_z)
        #Pressure Noisy functions
        sn.baro = trueSensors.baro + sensorBiases.baro + random.gauss(0, sensorSigmas.baro)
        sn.pitot = trueSensors.pitot + sensorBiases.pitot + random.gauss(0, sensorSigmas.pitot)

        if (self.UpdateTicks % self.gpsTicksUpdate == 0): #Check to see if ready to update
            gpsBN, gpsBE, gpsAlt = self.gpsGM.update()
            #GPS noisy functions
            sn.gps_n = trueSensors.gps_n + gpsBN + random.gauss(0, sensorSigmas.gps_n)
            sn.gps_e = trueSensors.gps_e + gpsBE + random.gauss(0, sensorSigmas.gps_e)
            sn.gps_alt = trueSensors.gps_alt + gpsAlt + random.gauss(0, sensorSigmas.gps_alt)
            sn.gps_sog = trueSensors.gps_sog + random.gauss(0, sensorSigmas.gps_sog*100)

            if (math.isclose(trueSensors.gps_cog, 0)):
                sn.gps_cog = trueSensors.gps_cog + random.gauss(0, sensorSigmas.gps_cog*100)
            else:
                sn.gps_cog = trueSensors.gps_cog + random.gauss(0, (sensorSigmas.gps_cog*VPC.InitialSpeed)/trueSensors.gps_sog)
        else: #Keep GPS values the same since not ready to update
            sn.gps_n = noisySensors.gps_n
            sn.gps_e = noisySensors.gps_e
            sn.gps_alt = noisySensors.gps_alt
            sn.gps_sog = noisySensors.gps_sog
            sn.gps_cog = noisySensors.gps_cog

        return sn

    def UpdateSensorsTrue(self, prevTrueSensors, state, dot):
        st = Sensors.vehicleSensors()
        #Setting all true values
        st.gyro_x, st.gyro_y, st.gyro_z = self.updateGyrosTrue(state)
        st.accel_x, st.accel_y, st.accel_z = self.updateAccelsTrue(state, dot)
        st.mag_x, st.mag_y, st.mag_z = self.updateMagsTrue(state)
        st.baro, st.pitot = self.updatePressureSensorsTrue(state)

        if (self.UpdateTicks % self.gpsTicksUpdate == 0):# Check if ready to update
            st.gps_n, st.gps_e, st.gps_alt, st.gps_sog, st.gps_cog = self.updateGPSTrue(state, dot)
        else:#Else keep everything the same
            st.gps_n = prevTrueSensors.gps_n
            st.gps_e = prevTrueSensors.gps_e
            st.gps_alt = prevTrueSensors.gps_alt
            st.gps_sog = prevTrueSensors.gps_sog
            st.gps_cog = prevTrueSensors.gps_cog

        return st
