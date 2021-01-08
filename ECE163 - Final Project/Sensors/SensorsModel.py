import math
import random
from ..Modeling import VehicleAerodynamicsModel
from ..Utilities import MatrixMath
from ..Containers import Sensors
from ..Constants import VehiclePhysicalConstants as VPC
from ..Constants import VehicleSensorConstants as VSC


class GaussMarkov:
    def __init__(self, dT=VPC.dT, tau=1e6, eta=0.0):
        # initialize internal variables
        self.dT = dT
        self.tau = tau
        self.eta = eta

        # save previous tau and eta in case values equal to None
        self.tau_prev = 0.0
        self.eta_prev = 0.0

        self.v = 0.0
        return

    def reset(self):
        # reset noise model
        self.v = 0.0
        return

    def update(self, vnoise=None):
        # assign tau and eta to previous value if they have None value
        if self.tau is None and self.eta is None:
            self.tau = self.tau_prev
            self.eta = self.eta_prev

        # generate noise if none passed into function
        if vnoise is None:
            vnoise = random.gauss(0, self.eta)

        # compute noise model
        v = (math.exp(-self.dT / self.tau) * self.v) + vnoise

        # update internal values
        self.v = v
        self.tau_prev = self.tau
        self.eta_prev = self.eta
        return v


class GaussMarkovXYZ:
    def __init__(self, dT=VPC.dT, tauX=1e6, etaX=0.0, tauY=None, etaY=None, tauZ=None, etaZ=None):
        # initialize internal variables
        self.dT = dT
        self.tauX = tauX
        self.etaX = etaX

        # if Y and Z values equal to None, set them to X values
        if tauY is None and etaY is None:
            self.tauY = tauX
            self.etaY = etaX
        else:
            self.tauY = tauY
            self.etaY = etaY

        if tauZ is None and etaZ is None:
            self.tauZ = tauX
            self.etaZ = etaX
        else:
            self.tauZ = tauZ
            self.etaZ = etaZ

        # initialize noise models for X, Y, and Z
        self.gmX = GaussMarkov(dT, self.tauX, self.etaX)
        self.gmY = GaussMarkov(dT, self.tauY, self.etaY)
        self.gmZ = GaussMarkov(dT, self.tauZ, self.etaZ)

        self.gmX.v = 0.0
        self.gmY.v = 0.0
        self.gmZ.v = 0.0
        return

    def reset(self):
        # reset noise model
        self.gmX.v = 0.0
        self.gmY.v = 0.0
        self.gmZ.v = 0.0
        return

    def update(self, vXnoise=None, vYnoise=None, vZnoise=None):
        # generate and update noise model
        vX = self.gmX.update(vXnoise)

        vY = self.gmY.update(vYnoise)

        vZ = self.gmZ.update(vZnoise)
        return vX, vY, vZ


class SensorsModel:
    def __init__(self, aeroModel=VehicleAerodynamicsModel.VehicleAerodynamicsModel(), taugyro=VSC.gyro_tau,
                 etagyro=VSC.gyro_eta, tauGPS=VSC.GPS_tau, etaGPSHorizontal=VSC.GPS_etaHorizontal,
                 etaGPSVertical=VSC.GPS_etaVertical, gpsUpdateHz=VSC.GPS_rate):
        # initialize internal variables
        self.aeroModel = aeroModel

        self.taugyro = taugyro
        self.etagyro = etagyro

        self.tauGPS = tauGPS
        self.etaGPSVertical = etaGPSVertical
        self.etaGPSHorizontal = etaGPSHorizontal

        self.gpsUpdateHz = gpsUpdateHz

        self.true = Sensors.vehicleSensors()
        self.noisy = Sensors.vehicleSensors()
        self.biases = self.initializeBiases()
        self.sigmas = self.initializeSigmas()

        self.gyroGM = GaussMarkovXYZ(VPC.dT, taugyro, etagyro)
        self.gpsGM = GaussMarkovXYZ(1/gpsUpdateHz, tauGPS, etaGPSHorizontal, tauGPS, etaGPSVertical)

        self.dT = aeroModel.vehicleDynamics.dT
        self.state = aeroModel.vehicleDynamics.state
        self.dot = aeroModel.vehicleDynamics.dot

        # update counters
        self.updateTicks = 0
        self.updateGPSticks = int(1 / (VPC.dT * gpsUpdateHz))
        return

    def getSensorsNoisy(self):
        # retrieve noisy sensor values
        noisy = self.noisy
        return noisy

    def getSensorsTrue(self):
        # retrieve true sensor values
        true = self.true
        return true

    def initializeBiases(self, gyroBias=VSC.gyro_bias, accelBias=VSC.accel_bias, magBias=VSC.mag_bias,
                         baroBias=VSC.baro_bias, pitotBias=VSC.pitot_bias):
        # compute initial biases for sensors
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
        # compute initial sigmas for sensors
        sensorSigmas = Sensors.vehicleSensors()

        sensorSigmas.gyro_x = gyroSigma
        sensorSigmas.gyro_y = gyroSigma
        sensorSigmas.gyro_z = gyroSigma

        sensorSigmas.accel_x = accelSigma
        sensorSigmas.accel_y = accelSigma
        sensorSigmas.accel_z = accelSigma

        sensorSigmas.mag_x = magSigma
        sensorSigmas.mag_y = magSigma
        sensorSigmas.mag_z = magSigma

        sensorSigmas.baro = baroSigma

        sensorSigmas.pitot = pitotSigma

        sensorSigmas.gps_n = gpsSigmaHorizontal
        sensorSigmas.gps_e = gpsSigmaHorizontal
        sensorSigmas.gps_alt = gpsSigmaVertical
        sensorSigmas.gps_sog = gpsSigmaSOG
        sensorSigmas.gps_cog = gpsSigmaCOG
        return sensorSigmas

    def reset(self):
        # reset internal values
        self.true = Sensors.vehicleSensors()
        self.prevTrueSensors = Sensors.vehicleSensors()
        self.noisy = Sensors.vehicleSensors()
        self.biases = self.initializeBiases()
        self.sigmas = self.initializeSigmas()
        self.gyroGM = GaussMarkovXYZ(VPC.dT, self.taugyro, self.etagyro)
        self.gpsGM = GaussMarkovXYZ(1 / self.gpsUpdateHz, self.tauGPS, self.etaGPSHorizontal, self.tauGPS, self.etaGPSVertical)
        self.dT = self.aeroModel.vehicleDynamics.dT
        self.state = self.aeroModel.vehicleDynamics.state
        self.dot = self.aeroModel.vehicleDynamics.dot

        self.updateTicks = 0
        return

    def update(self):
        # update internal true and noisy values
        self.true = self.updateSensorsTrue(self.true, self.state, self.dot)
        self.noisy = self.updateSensorsNoisy(self.true, self.noisy, self.biases, self.sigmas)

        # update ticks
        self.updateTicks += 1
        return

    def updateAccelsTrue(self, state, dot):
        # update true acceleration values
        accel_x = dot.u + (state.q * state.w) - (state.r * state.v) + (VPC.g0 * math.sin(state.pitch))
        accel_y = (dot.v + (state.r * state.u) - (state.p * state.w) - (VPC.g0 * math.cos(state.pitch) *
                                                                        math.sin(state.roll)))
        accel_z = (dot.w + (state.v * state.p) - (state.q * state.u) - (VPC.g0 * math.cos(state.pitch) *
                                                                        math.cos(state.roll)))
        return accel_x, accel_y, accel_z

    def updateGPSTrue(self, state, dot):
        # update true GPS values
        gps_n = state.pn
        gps_e = state.pe
        gps_alt = -state.pd
        gps_sog = math.hypot(state.u, state.v, state.w)
        gps_cog = math.atan2(dot.pe, dot.pn)
        return gps_n, gps_e, gps_alt, gps_sog, gps_cog

    def updateGyrosTrue(self, state):
        # update true gyro values
        gyro_x = state.p
        gyro_y = state.q
        gyro_z = state.r
        return gyro_x, gyro_y, gyro_z

    def updateMagsTrue(self, state):
        # update true magnetic field values
        mag = MatrixMath.matrixMultiply(state.R, VSC.magfield)
        mag_x = mag[0][0]
        mag_y = mag[1][0]
        mag_z = mag[2][0]
        return mag_x, mag_y, mag_z

    def updatePressureSensorsTrue(self, state):
        # update true pressure values
        baro = (VPC.rho * VPC.g0 * state.pd) + VSC.Pground
        pitot = (VPC.rho * (state.Va**2)) / 2
        return baro, pitot

    def updateSensorsNoisy(self, trueSensors=Sensors.vehicleSensors(), noisySensors=Sensors.vehicleSensors(),
                           sensorBiases=Sensors.vehicleSensors(), sensorSigmas=Sensors.vehicleSensors()):
        # update noisy sensor values
        noisy = Sensors.vehicleSensors()

        gx, gy, gz = self.gyroGM.update()
        noisy.gyro_x = trueSensors.gyro_x + sensorBiases.gyro_x + gx + random.gauss(0, sensorSigmas.gyro_x)
        noisy.gyro_y = trueSensors.gyro_y + sensorBiases.gyro_y + gy + random.gauss(0, sensorSigmas.gyro_y)
        noisy.gyro_z = trueSensors.gyro_z + sensorBiases.gyro_z + gz + random.gauss(0, sensorSigmas.gyro_z)

        noisy.accel_x = trueSensors.accel_x + sensorBiases.accel_x + random.gauss(0, sensorSigmas.accel_x)
        noisy.accel_y = trueSensors.accel_y + sensorBiases.accel_y + random.gauss(0, sensorSigmas.accel_y)
        noisy.accel_z = trueSensors.accel_z + sensorBiases.accel_z + random.gauss(0, sensorSigmas.accel_z)

        noisy.mag_x = trueSensors.mag_x + sensorBiases.mag_x + random.gauss(0, sensorSigmas.mag_x)
        noisy.mag_y = trueSensors.mag_y + sensorBiases.mag_y + random.gauss(0, sensorSigmas.mag_y)
        noisy.mag_z = trueSensors.mag_z + sensorBiases.mag_z + random.gauss(0, sensorSigmas.mag_z)

        noisy.baro = trueSensors.baro + sensorBiases.baro + random.gauss(0, sensorSigmas.baro)

        noisy.pitot = trueSensors.pitot + sensorBiases.pitot + random.gauss(0, sensorSigmas.pitot)

        # update noisy GPS values only if ready to, otherwise HOLD and use previous values
        if self.updateTicks % self.updateGPSticks == 0:
            gpsN, gpsE, gpsALT = self.gpsGM.update()
            noisy.gps_n = trueSensors.gps_n + gpsN + random.gauss(0, sensorSigmas.gps_n)
            noisy.gps_e = trueSensors.gps_e + gpsE + random.gauss(0, sensorSigmas.gps_e)
            noisy.gps_alt = trueSensors.gps_alt + gpsALT + random.gauss(0, sensorSigmas.gps_alt)

            noisy.gps_sog = trueSensors.gps_sog + random.gauss(0, sensorSigmas.gps_sog)

            if math.isclose(noisy.gps_sog, 0):
                noisy.gps_cog = trueSensors.gps_cog + random.gauss(0, sensorSigmas.gps_cog * 100)
            else:
                noisy.gps_cog = trueSensors.gps_cog + random.gauss(0, ((sensorSigmas.gps_cog * VPC.InitialSpeed) /
                                                                       trueSensors.gps_sog))
        else:
            noisy.gps_n = noisySensors.gps_n
            noisy.gps_e = noisySensors.gps_e
            noisy.gps_alt = noisySensors.gps_alt
            noisy.gps_sog = noisySensors.gps_sog
            noisy.gps_cog = noisySensors.gps_cog
        return noisy

    def updateSensorsTrue(self, prevTrueSensors, state, dot):
        # update true sensor values
        true = Sensors.vehicleSensors()
        true.gyro_x, true.gyro_y, true.gyro_z = self.updateGyrosTrue(state)
        true.accel_x, true.accel_y, true.accel_z = self.updateAccelsTrue(state, dot)
        true.mag_x, true.mag_y, true.mag_z = self.updateMagsTrue(state)
        true.baro, true.pitot = self.updatePressureSensorsTrue(state)

        # update noisy GPS values only if ready to, otherwise HOLD and use previous values
        if self.updateTicks % self.updateGPSticks == 0:
            true.gps_n, true.gps_e, true.gps_alt, true.gps_sog, true.gps_cog = self.updateGPSTrue(state, dot)
        else:
            true.gps_n = prevTrueSensors.gps_n
            true.gps_e = prevTrueSensors.gps_e
            true.gps_alt = prevTrueSensors.gps_alt
            true.gps_sog = prevTrueSensors.gps_sog
            true.gps_cog = prevTrueSensors.gps_cog
        return true
