import numpy as np

from . import matrixUtil
from . import ellipsoildUtil

# global variable: 
GRS80 = ellipsoildUtil.Ellipsoild("GRS80")

##################################################################
# base GNSS Center of Wuhan University,  professor niu xiaojin
#
# reference coordinate system: 
#   Navigation coordinate system, NORTH -   EAST -  DOWN
#   Body coordinate system,       FORWARD - RIGHT - DOWN

# yaw/pitch/roll：
# yaw:   ψ (Psi),  -180 ~ 180（0~360） degree
# pitch: θ (Theta),-90  ~ 90 
# roll:  φ (Phi),  -180 ~ 180 

#
##################################################################

# transfer euler angle to directinon consin matrix(DCM)
# 
# input: 3*1 vector,euler angle, roll-pitch-yaw, in radian,
#         the order of axis rotation is: Z(yaw)->Y(pitch)->X(roll),      
#         which rotates from referece coordinate system 2 body coordinate system,
#         right hand rule, North->East->Down
# output: 3*3 matrix,dirction consin matrix(DCM),from body to reference
def eurlerAngle2Dcm(eulerAngle):
    roll = eulerAngle[0]
    pitch = eulerAngle[1]
    yaw = eulerAngle[2]

    DCM = np.eye(3)
    # column 1  
    DCM[0,0] =  np.cos(pitch) * np.cos(yaw)
    DCM[1,0] =  np.cos(pitch) * np.sin(yaw)
    DCM[2,0] = -np.sin(pitch)   

    # column 2
    DCM[0,1] = -np.cos(roll) * np.sin(yaw) + np.sin(roll) * np.sin(pitch) * np.cos(yaw)
    DCM[1,1] =  np.cos(roll) * np.cos(yaw) + np.sin(roll) * np.sin(pitch) * np.sin(yaw)
    DCM[2,1] =  np.sin(roll) * np.cos(pitch)   

    # column 3
    DCM[0,2] =  np.sin(roll) * np.sin(yaw) + np.cos(roll) * np.sin(pitch) * np.cos(yaw)
    DCM[1,2] = -np.sin(roll) * np.cos(yaw) + np.cos(roll) * np.sin(pitch) * np.sin(yaw)
    DCM[2,2] =  np.cos(roll) * np.cos(pitch)

    return matrixUtil.gramschitd(DCM)

# transfer directinon consin matrix(DCM) to euler angle
# 
# input: 3*3 matrix,dirction consin matrix(DCM)
# output: 3*1 vector,euler angle, roll-pitch-yaw, in radian
# def dcm2EulerAngle(dcm):
#     dcm = dcm.reshape((3,3))    
    
#     # -sin(pitch) != 1 or -1
#     if dcm[2,0] - 1 < 1e-7 or dcm[2,0] + 1 < 1e-7:
#         pitch_0 = -np.arcsin(dcm[2,0])
#         pitch_1 = np.pi - pitch_0
        
#         roll_0 = np.arctan2(dcm[2,1]/np.cos(pitch_0), dcm[2,2]/np.cos(pitch_0))
#         roll_1 = np.arctan2(dcm[2,1]/np.cos(pitch_1), dcm[2,2]/np.cos(pitch_1))
        
#         yaw_0 = np.arctan2(dcm[1,0]/np.cos(pitch_0), dcm[0,0]/np.cos(pitch_0))
#         yaw_1 = np.arctan2(dcm[1,0]/np.cos(pitch_1), dcm[0,0]/np.cos(pitch_1))
        
#         return np.array([[roll_0,pitch_0,yaw_0],[roll_1,pitch_1,yaw_1]])
        
#     else:
#         yaw = 0
#         # -sin(pitch) == -1
#         if dcm[2,0] + 1 > 1e-7:
#             pitch = np.pi/2
#             roll = yaw + np.arctan2(dcm[0,1],dcm[0,2])
            
#             return np.array([roll,pitch,yaw])
#         else:
#             pitch = -np.pi/2
#             roll = -yaw + np.arctan2(-dcm[0,1],-dcm[0,2])
            
#             return np.array([roll,pitch,yaw])

def dcm2EulerAngle2(dcm):
    pitch = np.arctan2(-dcm[2,0],np.sqrt(np.power(dcm[2,1],2) + np.power(dcm[2,2],2)))
    
    roll = np.arctan2(dcm[2,1],dcm[2,2])
    
    yaw = np.arctan2(dcm[1,0],dcm[0,0])
    
    return np.array([roll,pitch,yaw])       

    
# transfer quaternion 2 directinon consin matrix(DCM)
# 
# input:4*1 vector,real numbler,quaternion
# output: 3*3 matrix
def quaternion2Dcm(quaternion):
    q0 = quaternion[0]    
    q1 = quaternion[1]
    q2 = quaternion[2]
    q3 = quaternion[3]
       
    DCM = np.ndarray((3,3))    
    # column 1  
    DCM[0,0] = q0**2 + q1**2 -q2**2 - q3**2 
    DCM[1,0] = 2 * (q1*q2 + q0*q3)
    DCM[2,0] = 2 * (q1*q3 - q0*q2) 

    # column 2
    DCM[0,1] = 2 * (q1*q2 - q0*q3)
    DCM[1,1] = q0**2 - q1**2 + q2**2 - q3**2 
    DCM[2,1] = 2 * (q2*q3 + q0*q1)

    # column 3
    DCM[0,2] = 2 * (q1*q3 + q0*q2) 
    DCM[1,2] = 2 * (q2*q3 - q0*q1) 
    DCM[2,2] = q0**2 - q1**2 -q2**2 + q3**2  

    return matrixUtil.gramschitd(DCM)

def dcm2Quaternion(dcm):
    eulerAngle = dcm2EulerAngle2(dcm)
    quaternion = eulerAngle2Quaternion(eulerAngle)
    return quaternion    

# transfer euler angle to quaternion
# 
# input: 3*1 vector,euler angle, roll-pitch-yaw, in radian,
#        the order of axis rotation is: Z(yaw)->Y(pitch)->X(roll),      
#        which rotates from referece coordinate system 2 body coordinate system,
#        right hand rule, North->East->Down
# output: 4*1 matrix, quaternion
def eulerAngle2Quaternion(eulerAngle):
    roll = eulerAngle[0]
    pitch = eulerAngle[1]
    yaw = eulerAngle[2]
    
    quaternion = np.array([0.0,0,0,0])
    quaternion[0] = np.cos(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    quaternion[1] = np.sin(roll/2) * np.cos(pitch/2) * np.cos(yaw/2) - np.cos(roll/2) * np.sin(pitch/2) * np.sin(yaw/2)
    quaternion[2] = np.cos(roll/2) * np.sin(pitch/2) * np.cos(yaw/2) + np.sin(roll/2) * np.cos(pitch/2) * np.sin(yaw/2)
    quaternion[3] = np.cos(roll/2) * np.cos(pitch/2) * np.sin(yaw/2) - np.sin(roll/2) * np.sin(pitch/2) * np.cos(yaw/2)

    return quaternion

def quaternion2EulerAngle(quaternion):
    dcm = quaternion2Dcm(quaternion)
    eulerAngle = dcm2EulerAngle2(dcm)    
    return eulerAngle
    

# transform the rotation-vector to DCM
# 
# input: 3*1 vector, x-axis/y-axis/z-axis
# output: 3*3, DCM
def rotationVector2Dcm(x,y,z):    
    rotationVector = np.array([x,y,z]).reshape((3,1))
    
    model_rv = np.sqrt(x*x + y*y + z*z)
    asm_rv = matrixUtil.getASM(rotationVector)
    
    part1 = np.eye(3)
    part2 = np.sin(model_rv) / model_rv * asm_rv
    part3 = (1- np.cos(model_rv)) / model_rv**2 * asm_rv @ asm_rv
    
    dcm = part1 + part2 + part3   
    return dcm

# input: roatation vector, np.array([x,y,z])
def rotationVector2Quaternion(rotationVector):
    # rotationVector = np.array([x,y,z]).reshape((3,1))        
    model_rv = np.sqrt(rotationVector[0]*rotationVector[0] + rotationVector[1]*rotationVector[1] + rotationVector[2]*rotationVector[2])
    quternion = np.array([0.0,0,0,0])
    # quaternionReal
    quternion[0] = np.cos(0.5 * model_rv) 
    # quaternionImage
    quternion[1:4] = np.sin(0.5 * model_rv) / (0.5 * model_rv) * 0.5 * rotationVector    
    
    return quternion

# compute the angle velocity of the special latitue in navigation system (NORTH - EAST - DOWN)
#  input: latitue
# output: 3*1 vector
def computeOmega_ie_n(latitude):
    omega_ie = GRS80.wie
    omega_ie_n = np.array([omega_ie * np.cos(latitude), 0, -omega_ie * np.sin(latitude)])
    return omega_ie_n

def computeOmega_en_n(latitude,height,velocityNorth,velocityEast,ellipsoild = GRS80):
    # a = ellipsoild.a
    # f = ellipsoild.f
    Rn = ellipsoild.computeRn(latitude)
    Rm = ellipsoild.computeRm(latitude)
    omega_en_n = np.array([velocityEast / (Rn + height), 
                           -velocityNorth / (Rm + height), 
                           -velocityEast * np.tan(latitude) / (Rn + height)])
    return omega_en_n

def computerOmega_in_n(latitude,height,velocityNorth,velocityEast,ellipsoild):
    w_en_n = computeOmega_en_n(latitude,height,velocityNorth,velocityEast,ellipsoild)
    w_ie_n = computeOmega_ie_n(latitude)    
    w_in_n = w_en_n + w_ie_n
    
    return w_in_n
    

# computer DSM by double vectors
# 
#  input: v1r/v2r/v1b/v2b, 3*1 vector  
#      
# output: DSM
def doubleVectorsAttitude(v1r,v2r,v1b,v2b):
    vr = np.zeros((3,3))    

    vr[:,0] = matrixUtil.vectorNormalization(v1r)
    vr[:,1] = matrixUtil.vectorNormalization(np.cross(v1r,v2r))
    vr[:,2] = matrixUtil.vectorNormalization(np.cross(np.cross(v1r,v2r),v1r))
    # vr = matrixUtil.normalization(vr)

    vb = np.zeros((3,3))
    vb[:,0] = matrixUtil.vectorNormalization(v1b)
    vb[:,1] = matrixUtil.vectorNormalization(np.cross(v1b,v2b))
    vb[:,2] = matrixUtil.vectorNormalization(np.cross(np.cross(v1b,v2b),v1b))
    # vb = matrixUtil.normalization(vb)

    return vr @ (vb.T)

# Static analytic coarse alignment
# navigation coordinate: north/east/down
# 
#  input: g_b,      3*1 vector
#         wib_b,    3*1 vector
#         latitude, rad
#      
# output: DSM
def staticCoarseAlignment(f_b,wib_b,latitude,height = 0,ellipsoild = GRS80):
    gravity = ellipsoild.computeGravity(latitude,height)    
    g_n = np.array([0, 0, gravity])
    
    wie = ellipsoild.wie
    wie_n = np.array([wie * np.cos(latitude), 0, -wie * np.sin(latitude)])

    g_b = -f_b
    wie_b = wib_b

    return doubleVectorsAttitude(g_n,wie_n,g_b,wie_b)






    
