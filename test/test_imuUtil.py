import sys

sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np

from src import imuUtil
from src import unitTools

def test_eurlerAngle2Dcm():
    roll = unitTools.degree2radian(60)
    pitch = unitTools.degree2radian(45)
    yaw = unitTools.degree2radian(30)
    eulerAngle = np.array([roll,pitch,yaw])

    dcm = imuUtil.eurlerAngle2Dcm(eulerAngle)
    dcm_check = np.array([[0.6123725,0.2803301,0.7391989],[0.3535534,0.7391989,-0.5732233],[-0.7071068,0.6123725,0.3535534]])
    # print(dcm)
    assert(np.max(np.abs(np.subtract(dcm,dcm_check))) < 0.001)

def test_eulerAngle2Quanternion():
    roll = unitTools.degree2radian(60)
    pitch = unitTools.degree2radian(45)
    yaw = unitTools.degree2radian(30)
    eulerAngle = np.array([roll,pitch,yaw])

    quaternion = imuUtil.eulerAngle2Quaternion(eulerAngle)
    quaternion_check = np.array([0.8223632,0.3604234, 0.4396797, 0.02226])
    # print(quaternion)
    assert(np.max(np.abs(np.subtract(quaternion,quaternion_check))) < 0.001)

def test_dcmFromQuaternion():
    quaternion = np.array([[0.8223632],[0.3604234], [0.4396797], [0.02226]])
    dcmFromQuaternion = imuUtil.quaternion2Dcm(quaternion)
    
    dcm_check = np.array([[0.6123725,0.2803301,0.7391989],[0.3535534,0.7391989,-0.5732233],[-0.7071068,0.6123725,0.3535534]])
    assert(np.max(np.abs(np.subtract(dcmFromQuaternion,dcm_check))) < 0.001)

def test_eulerAngleFromDcm():
    roll = unitTools.degree2radian(60)
    pitch = unitTools.degree2radian(45)
    yaw = unitTools.degree2radian(30)
    eulerAngle = np.array([roll,pitch,yaw])
    dcm = imuUtil.eurlerAngle2Dcm(eulerAngle)
    
    eulerAngleFromDcm = imuUtil.dcm2EulerAngle2(dcm)
    eulerAngle_check = np.array([1.04719755,0.78539816,0.52359878])
    # print(imuUtil.dcm2EulerAngle(dcm))
    # assert(np.max(np.abs(np.subtract(eulerAngleFromDcm[0],eulerAngle_check))) < 0.001)

def test_dcmFromRotationVector():
    dcmFromRotationVector = imuUtil.rotationVector2Dcm(0.7668133,0.9354339,0.047359)
    # print(dcmFromRotationVector)
    dcm_check = np.array([[0.6123725,0.2803301,0.7391989],[0.3535534,0.7391989,-0.5732233],[-0.7071068,0.6123725,0.3535534]])
    assert(np.max(np.abs(np.subtract(dcmFromRotationVector,dcm_check))) < 0.001)

def test_rotationVector2quaternion():
    rv = np.array([0.7668133,0.9354339,0.047359])
    q = imuUtil.rotationVector2Quaternion(rv)
    quaternion_check = np.array([0.8223632,0.3604234, 0.4396797, 0.02226])
    assert(np.max(np.abs(np.subtract(q,quaternion_check))) < 0.001)    

def test_staticCoarseAlignment():
    # 100C data
    g_i_imu = np.array([-0.104341081,0.118158883,-9.78100349])
    w_i_imu = np.array([0.000540942,0.000305581,-0.000382525])
    
    eulerAngle_b_imu = np.array([0,0,unitTools.degree2radian(180)])
    dcm_b_imu = imuUtil.eurlerAngle2Dcm(eulerAngle_b_imu)
    
    g_b = dcm_b_imu @ g_i_imu
    wib_b = dcm_b_imu @ w_i_imu
    
    # CPT data
    # g_b = np.array([-0.215111879,-0.033849041,9.805542419])
    # wib_b = np.array([-4.40541E-05,-5.22588E-05,5.10151E-05])

    # dcm = imuUtil.staticCoarseAlignment(g_b,wib_b,unitTools.degree2radian(31.1646379336))      
    # eulerAngle = unitTools.radian2degree(imuUtil.dcm2EulerAngle2(dcm))   
    
    # eulerAngleTest = np.array([0.779669204,0.706718828,151.4885478]) 
    # assert(np.max(np.abs(np.subtract(eulerAngle,eulerAngleTest))) < 2)
 
    # print(eulerAngle)

