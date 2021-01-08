import math
from . import MatrixMath


def dcm2Euler(DCM):
    if DCM[0][2] > 1:  # force value to +/-1 for pitch
        DCM[0][2] = 1
    elif DCM[0][2] < -1:
        DCM[0][2] = -1
    yaw = math.atan2(DCM[0][1], DCM[0][0]) # find euler angles
    pitch = -math.asin(DCM[0][2])
    row = math.atan2(DCM[1][2], DCM[2][2])
    return yaw, pitch, row


def euler2DCM(yaw, pitch, roll):
    ry = [[math.cos(yaw), math.sin(yaw), 0], # matrix for yaw
          [-math.sin(yaw), math.cos(yaw), 0],
          [0, 0, 1]]
    rp = [[math.cos(pitch), 0, -math.sin(pitch)], # matrix for pitch
          [0, 1, 0],
          [math.sin(pitch), 0, math.cos(pitch)]]
    rr = [[1, 0, 0], # matrix for roll
          [0, math.cos(roll), math.sin(roll)],
          [0, -math.sin(roll), math.cos(roll)]]
    dcm = MatrixMath.matrixMultiply(rp, ry) # combine matrices to get DCM
    dcm = MatrixMath.matrixMultiply(rr, dcm)
    return dcm


def ned2enu(points):
    result = points
    for i in result: # switch x and y columns
        temp = i[1]
        i[1] = i[0]
        i[0] = temp
    for i in result: # invert z
        i[2] *= -1
    return result
