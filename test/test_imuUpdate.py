
import sys

sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np


from src import imuUpdater
from src import unitTools
from src import imuUtil
from src import imuFactory
from src import ellipsoildUtil
   
# def test_update_100C():    
#     refernceEllipsoild = ellipsoildUtil.Ellipsoild()
    
#     # 100C para
#     imuParameter = imuFactory.ImuPara(1.0E-8,2.0E-8,1)
#     # euler angle of CPT to body 
#     eulerAngle_b_imu = unitTools.degree2radian(np.array([0,0,180]))   
    
#     #inatial value
#     positionK_2 = np. array([unitTools.degree2radian(121.4479958),unitTools.degree2radian(31.16558855),-8.9804])
#     positionK_1 = np. array([unitTools.degree2radian(121.4479958),unitTools.degree2radian(31.16558856),-8.98])
#     velocityK_2 = np.array([0.0003,-0.0009,	-0.0001])
#     velocityK_1 = np.array([0.0015,-0.0001,	-0.0014])  
#     eulerAngle_k_1 = unitTools.degree2radian(np.array([0.711088399,0.714204131,151.472879])) 
      
#     # file dir
#     rawImuFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220706_100C\\motion_20220706145455_100C_AGB7707_rawimu.ASC"
    
#     updater = imuUpdater.ImuPVAUpdate(positionK_2,positionK_1,
#                  velocityK_2,velocityK_1,
#                  eulerAngle_k_1,
#                  imuParameter,
#                  eulerAngle_b_imu,
#                  rawImuFileDir,def test_novatelUpdate():

#                  refernceEllipsoild)
#     # output file dir
#     outputFilleDir = "E:\\project_201903\\ecarx-tools-202106\\imuTodef test_updatePosition():
# ols\\data\\20220706_100C\\motion_20220706145455_100C_AGB7707_rawimu_compute.txt"
#     updater.result2File(outputFilleDir) 
def test_update():
    ellipsoild = ellipsoildUtil.Ellipsoild()
    
    # CPT para
    imuPara = imuFactory.ImuPara(0.1/(3600.0*256.0),0.05/np.power(2,15),100,
                                      3600,3600,3600,3600,
                                      unitTools.degree2radian(0.012) / np.power(3600,0.5),100 * 10E-9 * 9.8 / np.power(3600,0.5),
                                      unitTools.degree2radian(0.5)/3600,1250 * 10E-9 * 9.8,
                                      100 * 10E-9,100 * 10E-9)
    # euler angle of CPT to body     
    eulerAngle_imu_b = unitTools.degree2radian(np.array([180,0,90]))  
    # CPT inatial value
    positionK_2 = np. array([unitTools.degree2radian(31.16422804478),unitTools.degree2radian(121.44929669384),-5.8393])
    positionK_1 = np. array([unitTools.degree2radian(31.16422821342),unitTools.degree2radian(121.44929703574),-5.8391])    
    velocityK_2 = np.array([1.8311,3.3165,-0.0072])
    velocityK_1 = np.array([1.8324,3.3113,-0.0071])    
    eulerAngle_k_1 = unitTools.degree2radian(np.array([-0.017828612,0.226908041,60.60174824]))
    # initial error state       
    pk_1 = np.diag([0.02*0.02,0.02*0.02,0.05*0.05,
                0.01*0.01,0.01*0.01,0.01*0.01,
                unitTools.degree2radian(0.007*0.007),unitTools.degree2radian(0.007*0.007),unitTools.degree2radian(0.01*0.01),
                imuPara.sigma_gb,imuPara.sigma_gb,imuPara.sigma_gb,                  
                imuPara.sigma_ab,imuPara.sigma_ab,imuPara.sigma_ab,                    
                imuPara.sigma_gs,imuPara.sigma_gs,imuPara.sigma_gs,                   
                imuPara.sigma_as,imuPara.sigma_as,imuPara.sigma_as                   
                ]) 
    # file dir
    rawImuFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_RAWIMU.ASC"
    gpsPosFile = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_bestgnsspos.ASC"
    gpsVelFile = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_bestgnssvel.ASC"
    outputFilleDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_compute.txt"
      
    novatelUpdate = imuUpdater.NovatelRawImuAscFileUpdater(positionK_2, positionK_1, velocityK_2, velocityK_1, eulerAngle_k_1, eulerAngle_imu_b, imuPara, pk_1, ellipsoild,
                 rawImuFileDir)
    
    resultList = novatelUpdate.update()   
    
    novatelUpdate.result2File(outputFilleDir,resultList)