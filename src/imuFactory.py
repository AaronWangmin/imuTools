# import imuData 

from ssl import VerifyMode

import numpy as np
import re

from . import unitTools

class ImuPara:    
    def __init__(self,gyroScaleFactor,acceScaleFactor,outputRate,
                 Tgb,Tab,Tgs,Tas,   
                 arw,vrw,                                              
                 sigma_gb,sigma_ab, 
                 sigma_gs,sigma_as) -> None:
        
        self.gyroScaleFactor = gyroScaleFactor
        self.acceScaleFactor = acceScaleFactor
        self.outputRate = outputRate
        # correlative time, default value: 3600 s
        self.Tgb = Tgb
        self.Tab = Tab
        self.Tgs = Tgs
        self.Tas = Tas
        # angular random walk, eg: 0.012 deg/sqrt(hour) 
        self.arw = arw 
        # volocity random walk, eg: 100 ug/sqrt(hour)        
        self.vrw = vrw
        # gyro bias stability(std),eg: 0.5 deg /hour
        self.sigma_gb = sigma_gb
        # acce bias stability(std),eg: 1250 ug /hour
        self.sigma_ab = sigma_ab
         # gyro  scale factor repeatability,eg: 100 ppm
        self.sigma_gs = sigma_gs
        # acce scale factor repeatability,eg: 100 ppm
        self.sigma_as = sigma_as       

class ImuFactory:
    def __init__(self) -> None:
        self.imuRecords = [] 
        self.gpsPosRecords = []  
        self.gpsVelRecords = []    

    # get dvelocity and dAngle from novatel rawimuxa format file, *.asc
    def getRawImuDataFromNovateAscFile(self,fileDir,imuPara):
        f = open(fileDir, "r")        

        for strLine in f.readlines(): 
            
            startIndex = 0
            value = strLine.split(",") 
            time   =  float(value[startIndex + 12])
            zAccel =  float(value[startIndex + 14]) * imuPara.acceScaleFactor 
            yAccel = -float(value[startIndex + 15]) * imuPara.acceScaleFactor 
            xAccel =  float(value[startIndex + 16]) * imuPara.acceScaleFactor 
            zGyro  =  float(value[startIndex + 17]) * imuPara.gyroScaleFactor 
            yGyro  = -float(value[startIndex + 18]) * imuPara.gyroScaleFactor 
            xGyro  =  float(value[startIndex + 19].split("*")[0]) * imuPara.gyroScaleFactor 
            
            imuRecord = [time,xAccel,yAccel,zAccel,xGyro,yGyro,zGyro]
            
            self.imuRecords.append(imuRecord)
                       
        f.close()
        
    # get gps position from novatel bestgnsspos format file, *.asc
    def getBestGnssPosDataFromNovateAscFile(self,fileDir):
        f = open(fileDir, "r")        

        for strLine in f.readlines(): 
            
            startIndex = 0            
            # value = strLine.split(",") 
            value = re.split('[;,]',strLine) 
            solutionStatus = value[10]
            # extract valid gpspose    
            if(solutionStatus == "SOL_COMPUTED" ):
                time   =  float(value[startIndex + 6])
                # rad
                latitude        = unitTools.degree2radian(float(value[startIndex + 12]))
                longitude       = unitTools.degree2radian(float(value[startIndex + 13])) 
                # m
                height          = -float(value[startIndex + 14])            
                sigmaLatitude   = float(value[startIndex + 17])  
                sigmaLongitude  = float(value[startIndex + 18]) 
                sigmaHeight     = float(value[startIndex + 19])  
                
                gpsPosRecord = [time,latitude,longitude,height,sigmaLongitude,sigmaLatitude,sigmaHeight]
                
                self.gpsPosRecords.append(gpsPosRecord)
                       
        f.close()
    
    # get gps velocity from novatel bestgnssvel format file, *.asc
    def getBestGnssVelDataFromNovateAscFile(self,fileDir):
        f = open(fileDir, "r")        

        for strLine in f.readlines(): 
            
            startIndex = 0
            # value = strLine.split(",")
            value = re.split('[;,]',strLine) 
            solutionStatus = value[10]
            # extract valid gpspose    
            if(solutionStatus == "SOL_COMPUTED" ):             
                time          =  float(value[startIndex + 6])
                horizontalVel =  float(value[startIndex + 14])
                azimuth       =  unitTools.degree2radian(float(value[startIndex + 15]))
                verticalVel   = -float(value[startIndex + 16])       
                
                # compute, m/s
                northVel      = horizontalVel * np.cos(azimuth)      
                eastVel       = horizontalVel * np.sin(azimuth)
                
                gpsVelRecord = [time,northVel,eastVel,verticalVel]
                
                self.gpsVelRecords.append(gpsVelRecord)
                       
        f.close()
    
    # get corrected delta-velocity and delta-dAngle by bias 
    def correctImuRecord(self):
        
        pass
        
            
           


#--------------------------------test-------------
# novatel 100C-IMU para
# novatel100CPara = ImuPara(1.0E-8,2.0E-8,200)  

# fd = ImuFactory()
# fd.RawImuData\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220706_100C\\driver_20220706145455_RZC_SHANGHAI_AGB7707_rawimuxa.ASC",novatel100CPara)
# print(len(fd.imuRecords))
# print(len(fd.ZAccelList))
