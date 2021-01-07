"""
Luis Mercado
lurmerca 1658336
ECE 163

Objective: Create all rotation functions
"""

import math
from . import MatrixMath

def euler2DCM(yaw, pitch, roll):
    firstRot = [[math.cos(yaw), math.sin(yaw), 0], #first rotation matrix for yaw angle
                [-(math.sin(yaw)), math.cos(yaw), 0],
                [0,0,1]]

    secondRot = [[math.cos(pitch), 0, -(math.sin(pitch))], #second rotation matrix for pitch angle
                [0, 1, 0],
                [math.sin(pitch), 0, math.cos(pitch)]]

    thirdRot = [[1, 0, 0], #third rotation matrix for roll angle
                [0, math.cos(roll), math.sin(roll)],
                [0, -(math.sin(roll)), math.cos(roll)]]

    DCM1 = MatrixMath.matrixMultiply(thirdRot, secondRot) #first half of the final matrix multiplication
    FullDCM = MatrixMath.matrixMultiply(DCM1, firstRot) #second half of big matric multiplication

    return FullDCM #returning final DCM

def dcm2Euler(EA):

    if (EA[0][2] > 1): #Grim lock check
        EA[0][2] = 1
    elif (EA[0][2] < -1):
        EA[0][2] = -1

    pitch = -(math.asin(EA[0][2]))
    roll = math.atan2(EA[1][2], EA[2][2]) #DCM to euler conversitions
    yaw = math.atan2(EA[0][1], EA[0][0])

    result = [yaw, pitch, roll]

    return result

def ned2enu (NEDpoints):
    R = [[0,1,0], #matrix to switch north and east, and convert down to up
        [1,0,0],
        [0,0,-1]]

    result = MatrixMath.matrixMultiply(NEDpoints, R) #mulitply to convert

    return result
