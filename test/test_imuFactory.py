import sys

sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np

from src import imuFactory
from src import unitTools

class TestImuFactory:
    def setup(self):
        self.imuFac = imuFactory.ImuFactory()       

    def test_read(self):
        
        # CPT para
        imuParameter = imuFactory.ImuPara(0.1/(3600.0*256.0),0.05/np.power(2,15),100,
                                        3600,3600,3600,3600,
                                        unitTools.degree2radian(0.012) / np.power(3600,0.5),100 * 10E-6 / np.power(3600,0.5),
                                        unitTools.degree2radian(0.5)/3600,1250 * 10E-9,
                                        100 * 10E-9,100 * 10E-9) 
        rawImuFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_RAWIMU.ASC"
        self.imuFac.getRawImuDataFromNovateAscFile(rawImuFileDir,imuParameter)
        print("hello")
        
    def test_bestgnsspos(self):
        inputFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_bestgnsspos.ASC"
        self.imuFac.getBestGnssPosDataFromNovateAscFile(inputFileDir)
        print("hello")
    
    def test_bestgnssVel(self):
        inputFileDir = "E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220714_CPT\\20220714145029_CPT_BF05633_MOTION_bestgnssvel.ASC"
        self.imuFac.getBestGnssVelDataFromNovateAscFile(inputFileDir)
        print("hello")
     
     