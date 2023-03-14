import sys

sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np

from src import imuGpsFusionModel
from src import imuUtil
from src import ellipsoildUtil
from src import unitTools
from src import imuFactory
from src import imuUpdater

class TestGnssPVObsEquation():
    def setup(self):        
        # known parameter
        numberMeasureVariable = 6
        numberStateVariable = 21
        position = np.array([unitTools.degree2radian(121.44929703574),unitTools.degree2radian(31.16422821342),-5.8391])  
        velocity = np.array([1.8324,3.3113,-0.0071])
        eulerAngle_imu_b = unitTools.degree2radian(np.array([180,0,90])) 
        cbn = imuUtil.eurlerAngle2Dcm(eulerAngle_imu_b)
        leverArm = [0.6,1.1,1.2]
        gpsPos = np.array([unitTools.degree2radian(121.44929703575),unitTools.degree2radian(31.16422821346),-5.8399]) 
        gpsVel = np.array([1.8321,3.3112,-0.0072])
        rmsGpsPos = [3,3,3]
        rmsGpsVel = [0.05,0.05,0.05]
        wib_b = [1,2,3]

        ellipsoild = ellipsoildUtil.Ellipsoild()

        self.obsEquation = imuGpsFusionModel.GnssPVObsEquation(numberMeasureVariable,numberStateVariable,position,velocity,cbn,leverArm,gpsPos,gpsVel,rmsGpsPos,rmsGpsVel,wib_b,ellipsoild) 

    def test_setR(self):        
        self.obsEquation.setGnssPVObsEquationPara()        
        print("hello")
 
class TestImuErrorStateEquation():
    def setup(self):       
       # CPT para
        imuPara = imuFactory.ImuPara(0.1/(3600.0*256.0),0.05/np.power(2,15),100,
                                        3600,3600,3600,3600,
                                        unitTools.degree2radian(0.012) / np.power(3600,0.5),100 * 10E-9 * 9.8 / np.power(3600,0.5),
                                        unitTools.degree2radian(0.5)/3600,1250 * 10E-9 * 9.8,
                                        100 * 10E-9,100 * 10E-9)
        
        ellipsoild = ellipsoildUtil.Ellipsoild()
        
        # euler angle of CPT to body     
        eulerAngle_imu_b = unitTools.degree2radian(np.array([180,0,90]))  
        # CPT inatial value
        positionK_2 = np. array([unitTools.degree2radian(31.16422804478),unitTools.degree2radian(121.44929669384),-5.8393])
        positionK_1 = np. array([unitTools.degree2radian(31.16422821342),unitTools.degree2radian(121.44929703574),-5.8391])    
        velocityK_2 = np.array([1.8311,3.3165,-0.0072])
        velocityK_1 = np.array([1.8324,3.3113,-0.0071])    
        eulerAngle_k_1 = unitTools.degree2radian(np.array([-0.017828612,0.226908041,60.60174824])) 
        
              
                   
        resultList = []
        
        # initial error state       
        # pk_1 = np.diag([1.0,1.0,1.5,
        #             0.01,0.01,0.01,
        #             unitTools.degree2radian(0.007),unitTools.degree2radian(0.007),unitTools.degree2radian(0.01),
        #             imuPara.sigma_gb,imuPara.sigma_gb,imuPara.sigma_gb,                  
        #             imuPara.sigma_ab,imuPara.sigma_ab,imuPara.sigma_ab,                    
        #             imuPara.sigma_gs,imuPara.sigma_gs,imuPara.sigma_gs,                   
        #             imuPara.sigma_as,imuPara.sigma_as,imuPara.sigma_as                   
        #             ])
        pk_1 = np.diag([0.02,0.02,0.05,
                    0.01,0.01,0.01,
                    unitTools.degree2radian(0.007),unitTools.degree2radian(0.007),unitTools.degree2radian(0.01),
                    imuPara.sigma_gb,imuPara.sigma_gb,imuPara.sigma_gb,                  
                    imuPara.sigma_ab,imuPara.sigma_ab,imuPara.sigma_ab,                    
                    imuPara.sigma_gs,imuPara.sigma_gs,imuPara.sigma_gs,                   
                    imuPara.sigma_as,imuPara.sigma_as,imuPara.sigma_as                  
                    ])
        # pk_1 = np.eye(21) * 0.0
        
        # file dir
        rawImuFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_RAWIMU.ASC"
               
        self.novatelUpdate = imuUpdater.NovatelRawImuAscFileUpdater(positionK_2, positionK_1, velocityK_2, velocityK_1, eulerAngle_k_1, eulerAngle_imu_b, imuPara, pk_1, ellipsoild,
                 rawImuFileDir)       
        
    def test_kmUpdate(self):
         # leval arm        
        levalArm = [0.0,-0.6,-0.012]        
        rmsGpsPos = np.array([1.0,1.0,1.5])
        rmsGpsVel = np.array([0.007,0.008,0.02]) 
        
         # file dir
        # rawImuFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_RAWIMU.ASC"
        gpsPosFile = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_bestgnsspos.ASC"
        gpsVelFile = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_bestgnssvel.ASC"
        
        outputFilleDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_compute_2.txt"         
        
        resultList = self.novatelUpdate.kfUpdate(gpsPosFile,gpsVelFile,levalArm,rmsGpsPos,rmsGpsVel)
        self.novatelUpdate.result2File(outputFilleDir,resultList)
  