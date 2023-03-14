from asyncio.windows_events import NULL
import numpy as np

from . import imuUtil
from . import imuFactory
from . import imuGpsFusionModel
from . import kalmanFilter
from . import ellipsoildUtil
from . import matrixUtil
from . import unitTools

# globle variable 
GRS80 = ellipsoildUtil.Ellipsoild()

def recursionLine(xk_2,xk_1):
    x = xk_1 + (xk_1 - xk_2) / 2
    return x

# recursion the values of wie_n,wen_n,v_n,g_n at time k-1/2,by the values at time k-1 and k-2
def recursionLineValues(positionK_1,positionK_2,velocityK_1,velocityK_2,ellipsoild):    
    positionRecursion = recursionLine(positionK_2,positionK_1)
    velocityRecursion = recursionLine(velocityK_2,velocityK_1)
    
    latitudeRecursion = positionRecursion[0]
    heightRecursion = positionRecursion[2]
    wie_n_recursion = imuUtil.computeOmega_ie_n(latitudeRecursion)
    wen_n_recursion = imuUtil.computeOmega_en_n(latitudeRecursion,heightRecursion,velocityRecursion[0],velocityRecursion[1],ellipsoild)  
    gravity_recursion = ellipsoild.computeGravity(latitudeRecursion,heightRecursion)   
    
    return  positionRecursion,velocityRecursion,wie_n_recursion,wen_n_recursion,gravity_recursion

class ImuAttitudeUpdater():    
    # input:       
    #        tk,                time at tk  
    #        positionK,         3*1 vector,the rescursion-position by positionK_1 and positionK_2, longitude/latitude/height,
    #        velocityK,         3*1 vector,the rescursion-velocity by velocityK_1 and velocityK_2, north/east/height  
    #        tk_1,              time at tk-1    
    #        cbn_k_1,           3*3 matrix, the known-attitude at tk-1 moment
    #        dAngleK_1,         3*1 vector, the output of imu at tk-1 moment
    #        dAngleK,
    #        ellipsoild,
    # output: 3*1 vector,       the updated velocity 
    def __init__(self,tk,positionK,velocityK,tk_1, cbn_k_1, dAngleK_1, dAngleK, ellipsoild=GRS80) -> None:
        self.tk = tk
        self.positionK = positionK
        self.velocityK = velocityK
        self.tk_1 = tk_1
        self.cbn_k_1 = cbn_k_1
        self.dAngleK_1 = dAngleK_1
        self.dAngleK = dAngleK
        self.ellipsoild =  ellipsoild   
    # update Body-coordinate-system by rotationVector (FORWARD - RIGHT - DOWN)
    # 
    #  input: dAngleCurrent, 3*1 vector,the angle-increment of current time
    #         dAnglePrevious,3*1 vector,the angle-increment of previous time 
    # output: Quaternion
    def updateQuaternion_bk_1_bk(self):
        # asm = matrixUtil.getASM(dAnglePrevious)
        rotationVector = self.dAngleK + 1/12 * np.cross(self.dAngleK_1, self.dAngleK)
        quarternion_bk_1_bk = imuUtil.rotationVector2Quaternion(rotationVector)  
        return quarternion_bk_1_bk

    # update Navigation-coordinate-system by rotationVector (FORWARD - RIGHT - DOWN)
    # 
    #  input: latitude,
    #         height,
    #         velocityNorth,
    #         velocityEast,
    #         dTime,          
    # output: Quaternion
    def updateQuarternion_nk_nk_1(self):
        # (positionRecursion,velocityRecursion,wie_n_recursion,wen_n_recursion,gravity_recursion) = self.recursionLineValues()
        omega_ie_n = imuUtil.computeOmega_ie_n(self.positionK[0])
        omega_en_n = imuUtil.computeOmega_en_n(self.positionK[0],self.positionK[2],self.velocityK[0],self.velocityK[1])
        rotationVector = (omega_ie_n + omega_en_n) * (self.tk - self.tk_1)  
        quarternion_nk_nk_1 = imuUtil.rotationVector2Quaternion(rotationVector) 
        return quarternion_nk_nk_1

    # update current Quaternion
    def updateQuarternion_nb_k(self):
        quarternion_bk_1_bk = self.updateQuaternion_bk_1_bk()
        quarternion_nk_nk_1 = self.updateQuarternion_nk_nk_1()
        quarternion_nb_tk_1 = self.cbn_k_1
        
        quarternion_nb_k = matrixUtil.quarternionMultiply(quarternion_nk_nk_1,quarternion_nb_tk_1)
        quarternion_nb_k = matrixUtil.quarternionMultiply(quarternion_nb_k,quarternion_bk_1_bk)
        quarternion_nb_k = matrixUtil.vectorNormalization(quarternion_nb_k)        
        
        return quarternion_nb_k
    
class ImuVelocityUpdater():
    
    def __init__(self,
                 tk,velocityRecursion,wie_n_recursion,wen_n_recursion,gravityRecursion,
                 dAngleK,dVelocityK,
                 tk_1,velocityK_1, cbn_k_1, 
                 dAngleK_1,dVelocityK_1, ellipsoild) -> None:
        self.tk = tk       
        self.veloctityRecursion = velocityRecursion
        self.gravityRecursion = gravityRecursion
        self.wie_n_recursion = wie_n_recursion
        self.wen_n_recursion = wen_n_recursion                
        self.dAngleK = dAngleK
        self.dVelocityK = dVelocityK
        self.tk_1 = tk_1
        self.cbn_k_1 = cbn_k_1
        self.velocityK_1 = velocityK_1
        self.dAngleK_1 = dAngleK_1
        self.dVelocityK_1 = dVelocityK_1
        self.ellipsoild = ellipsoild     
            
    # computer the first part of velocity by harmful acceleration
    # input: wie_n_recursion,
    #        wen_n_recursion,
    #        gravity_recursion    #        
    #        tk,tk_1, float, the time of current and previous
    # output: 
    def computeHarmfulDeltaVelocity(self):        
        deltaVelocityByHarmfulAcce = (self.gravityRecursion - np.cross(2 * self.wie_n_recursion + self.wen_n_recursion,self.veloctityRecursion)) \
            * (self.tk -self.tk_1)
            
        return deltaVelocityByHarmfulAcce
    
    # computer the rotation-compensation term
    # input: dAngle_k,      3*1 vector, the delta-angle of output of imu 
    #        dVelocity_k,   3*1 vector, the delta-velocity of output of imu                 
    #        
    # output: 3*1 vector 
    def computerRotationCompensation(self):        
        rc = 1/2 * np.cross(self.dAngleK,self.dVelocityK)
        return rc
    
    # computer the sculling-compensation term
    # input: dAngle_k,          3*1 vector, the delta-angle of output of imu at time of k
    #        dAngle_k_1         at time of k
    #        dVelocity_k,       3*1 vector, the delta-velocity of output of imu at time of k
    #        dVecocity_k_1      at time of k    #                
    #        
    # output: 3*1 vector 
    def computerScullingCompensation(self):        
        sc = 1/12 * (np.cross(self.dAngleK_1,self.dVelocityK) + np.cross(self.dVelocityK_1,self.dAngleK))
        return sc
    
    
    def computeUsefulDeltaVelocity(self):       
        xi = (self.wie_n_recursion + self.wen_n_recursion) * (self.tk - self.tk_1)
        
        rotationCompensation = self.computerRotationCompensation()
        scullingCompensation = self.computerScullingCompensation()
        dV_fk_bk_1 = self.dVelocityK_1 + rotationCompensation + scullingCompensation       
               
        deltaVelocityByUsefulAcce = (np.eye(3) - 1/2 * matrixUtil.getASM(xi)) @ self.cbn_k_1 @ dV_fk_bk_1
        
        return deltaVelocityByUsefulAcce
    
    # update the veloctiy
    # input: positionK_2,       3*1 vector, the known-position at tk-2 moment
    #        positionK_1        
    #        velocityK_2,       3*1 vector, the known-velocity at tk-2 moment
    #        velocityK_1,       
    #        cbn_k_1,           3*3 matrix, the known-attitude at tk-1 moment
    #
    #        dAngleK_1,         3*1 vector, the output of imu at tk-1 moment
    #        dAngleK,
    #        dVelocityK_1,     
    #        dVelocityK,
    #        tk_1,tk            time at tk-1,tk
    #
    # output: 3*1 vector,       the updated velocity 
    def velocityUpdate(self):         
        dv_gcor_k_n = self.computeHarmfulDeltaVelocity()
        
        dv_fk_n = self.computeUsefulDeltaVelocity()
        
        vk_n = self.velocityK_1 + dv_gcor_k_n + dv_fk_n        
        
        return vk_n

class ImuPositionUpdate():    
    def __init__(self,tk,velocityK,tk_1,positionK_1,velocityK_1,rmRecursion,ellipsoild) -> None:
        self.tk = tk
        self.velocityK = velocityK
        self.tk_1 = tk_1
        self.positionK_1 = positionK_1
        self.velocityK_1 = velocityK_1
        self.rmRecursion = rmRecursion
        # self.Rn = Rn
        self.ellipsoild = ellipsoild        
     
    # update position
    # navigation coordinate: north/east/down
    # 
    #  input: PositionPrevious,      3*1 vector, latitude/longitude/height
    #         velocityPrevious,      3*1 vector
    #         velocityCurrent,       3*1 vector
    #         timePerioud
    #         timeCurrent   
    # output: position, 3*1 vector
    def updatePosition(self):        
        velocityAverage = (self.velocityK_1 + self.velocityK) / 2      
        dTime = self.tk - self.tk_1
        
        # update height  
        heightCurrent = self.positionK_1[2] - velocityAverage[2] * dTime
        
        # update latitude    
        # Rm = self.ellipsoild.computeRm(self.positionK_1[1])
        heightAverage = (heightCurrent + self.positionK_1[2]) /2
        latitudeCurrent = self.positionK_1[0] + velocityAverage[0] / (self.rmRecursion + heightAverage) * dTime
        
        # update longitude
        latitudeAverage = (latitudeCurrent + self.positionK_1[0]) / 2
        Rn = self.ellipsoild.computeRn(latitudeCurrent)
        logitudeCurrent = self.positionK_1[1] + velocityAverage[1] / ((Rn + heightAverage) * np.cos(latitudeAverage)) *dTime
        
        position = np.array([latitudeCurrent,logitudeCurrent,heightCurrent])
        return position

# input:positionK_2, positionK_1, 3*1 vector, known value
#       velocityK_2, velocityK_1
#       eulerAngle_k_1
#       imuParameter
#       eulerAngle_imu_b
#       rawImuFileDir,
#       rawImuRecord, 6*1 vector: [time,xAccel,yAccel,zAccel,xGyro,yGyro,zGyro]          
#       ellipsoild
# output: (time,cbn,position,velocity)     
class ImuPVAUpdate:
    def __init__(self,
                 positionK_2,positionK_1,
                 velocityK_2,velocityK_1,
                 eulerAngle_k_1,               
                 eulerAngle_imu_b,
                 imuPara,
                 pk_1,             
                 ellipsoild,
                 levalArm = [0,0,0]) -> None:
        self.positionK_2 = positionK_2
        self.positionK_1 = positionK_1
        self.velocityK_2 = velocityK_2
        self.velocityK_1 = velocityK_1
        self.eulerAngle_k_1 = eulerAngle_k_1     
        self.eulerAngle_imu_b = eulerAngle_imu_b
        self.imuPara = imuPara
        # covariance matrix of error state variable 
        self.pk_1 = pk_1
        self.ellipsoild = ellipsoild        
        
        # error state equation
        self.errorStateEquation = NULL
        
        # add error equation
    
    # update one pva-record by one imu record    
    def pvaUpdate(self,imuRecordK_1,imuRecordK):             
        quaternionK_1 = imuUtil.eulerAngle2Quaternion(self.eulerAngle_k_1)        
        # euler angle of imu to body             
        dcm_imu_b = imuUtil.eurlerAngle2Dcm(self.eulerAngle_imu_b)                                    
        
        # imu measurement value at tk-1
        tk_1         = imuRecordK_1[0]
        dVelocityK_1 = dcm_imu_b @ imuRecordK_1[1:4]
        dAngleK_1    = dcm_imu_b @ imuRecordK_1[4:7]
        
        # imu measurement value at tk
        tk         = imuRecordK[0]
        dVelocityK = dcm_imu_b @ imuRecordK[1:4]      
        dAngleK    = dcm_imu_b @ imuRecordK[4:7] 
        
        # rescursion value by the value at time k_2 and K_1       
        (positionRecursion,velocityRecursion,wie_n_recursion,wen_n_recursion,gravity_recursion) = \
            recursionLineValues(self.positionK_1,self.positionK_2,self.velocityK_1,self.velocityK_2,self.ellipsoild)
                        
        # attitude update
        imuAttitudeUpdater = ImuAttitudeUpdater(tk,positionRecursion,velocityRecursion,tk_1,quaternionK_1,dAngleK_1,dAngleK,self.ellipsoild)
        attitudeCurrent = imuAttitudeUpdater.updateQuarternion_nb_k()
        cbn_k_1 = imuUtil.quaternion2Dcm(quaternionK_1)
        
        # velocity update
        imuVelocityUpdater =  ImuVelocityUpdater(tk,velocityRecursion,wie_n_recursion,wen_n_recursion,gravity_recursion,
                                                    dAngleK,dVelocityK,
                                                    tk_1,self.velocityK_1,cbn_k_1,
                                                    dAngleK_1,dVelocityK_1,self.ellipsoild)
        velocityCurrent = imuVelocityUpdater.velocityUpdate()
        
        # postition update
        rmRecursion = self.ellipsoild.computeRm(positionRecursion[0])
        imuPositionUpdater = ImuPositionUpdate(tk,velocityCurrent,tk_1,self.positionK_1,self.velocityK_1,rmRecursion,self.ellipsoild)
        positionCurrent = imuPositionUpdater.updatePosition()
        
        # state error estimate        
        self.errorStateEquation = imuGpsFusionModel.ImuErrorStateEquation(21,18,self.positionK_1,self.velocityK_1,cbn_k_1,imuUtil.quaternion2Dcm(attitudeCurrent),
                                                                     dVelocityK,dAngleK,tk,tk_1,self.pk_1,self.imuPara,self.ellipsoild)
        self.errorStateEquation.setErrorStateEquation()
        (xk,pmCurrent) = self.errorStateEquation.timeUpdate()                        
        
        # rescursion
        self.positionK_2 = self.positionK_1
        self.positionK_1 = positionCurrent
        self.velocityK_2 = self.velocityK_1
        self.velocityK_1 = velocityCurrent
        quaternionK_1 = attitudeCurrent
        
        self.pk_1 = pmCurrent        
    
        return imuRecordK[0],positionCurrent,velocityCurrent,attitudeCurrent,pmCurrent[np.arange(21),np.arange(21)]    

class NovatelRawImuAscFileUpdater(ImuPVAUpdate):        
    def __init__(self, positionK_2, positionK_1, velocityK_2, velocityK_1, eulerAngle_k_1, eulerAngle_imu_b, imuPara, pk_1, ellipsoild,
                 rawImuAscFileDir) -> None:
        super().__init__(positionK_2, positionK_1, velocityK_2, velocityK_1, eulerAngle_k_1, eulerAngle_imu_b, imuPara, pk_1, ellipsoild)      
        self.rawImuAscFileDir = rawImuAscFileDir  
        
    def update(self):
        # get imu data from  rawimu file
        fd = imuFactory.ImuFactory()
        fd.getRawImuDataFromNovateAscFile(self.rawImuAscFileDir,self.imuPara)    
        imuRecords = fd.imuRecords      
               
        resultList = []             
        for i in range(2,len(imuRecords)-5):                        
            (time,position,velocity,quaternion,variance) = self.pvaUpdate(imuRecords[i-1],imuRecords[i])         
            
            result = np.concatenate((position,velocity,quaternion,variance))  
            result = np.insert(result,0,time)
            resultList.append(result)                              
        
        return resultList
    
    def kfUpdate(self,gpsPosFile,gpsVelFile,levalArm = [0,0,0],rmsGpsPos = [0,0,0],rmsGpsVel = [0,0,0]):         
        # get imu data from  rawimu file
        fd = imuFactory.ImuFactory()
        fd.getRawImuDataFromNovateAscFile(self.rawImuAscFileDir,self.imuPara)  
        fd.getBestGnssPosDataFromNovateAscFile(gpsPosFile)
        fd.getBestGnssVelDataFromNovateAscFile(gpsVelFile)  
        imuRecords = fd.imuRecords   
        gpsPosRecords = fd.gpsPosRecords
        gpsVelRecords = fd.gpsVelRecords  
        
        # initial value                    
        # errorStateEquation = NULL
        
        # index of gpsPosRecords
        m = 0
        # gps sample rate 
        gpsRate = 100
                               
        # inital value
        (time,position,velocity,attitude,variance) = (0,0,0,0,0)
        resultList = []        
        # xk_1 = np.zeros(21)        
        for i in range(2,len(imuRecords)-5,gpsRate):            
            
            for j in range(gpsRate): 
                if(i < len(imuRecords)-5):                 
                    (time,position,velocity,attitude,variance) = self.pvaUpdate(imuRecords[i-1],imuRecords[i]) 
                    self.pk_1 = np.diag(variance)                   
                     
                    result = np.concatenate((position,velocity,attitude,variance))  
                    result = np.insert(result,0,time)
                    resultList.append(result)
                    i += 1   
                             
            # kf measure update(gpsPosition,gpsVelocity)            
            if(m < len(gpsPosRecords)):               
                # measure equation                 
                gpsPos    = np.array(gpsPosRecords[m][1:4])
                rmsGpsPos = np.array(gpsPosRecords[m][4:])
                gpsVel    = np.array(gpsVelRecords[m][1:])
                m += 1
                               
                dAngleK = imuRecords[i][4:]
                # kalman filter
                gnssPvObsEquation = imuGpsFusionModel.GnssPVObsEquation(6,21,position,velocity,imuUtil.quaternion2Dcm(attitude),levalArm,
                                                                        gpsPos,gpsVel,rmsGpsPos,rmsGpsVel,dAngleK,self.ellipsoild)
                gnssPvObsEquation.setGnssPVObsEquationPara()           
                (xk,pk) = gnssPvObsEquation.measureUpdate(self.pk_1)                
                               
                # update position
                resultList[-1][1:4] = gnssPvObsEquation.position                
                # update velocity
                resultList[-1][4:7] = gnssPvObsEquation.velocity
                # update attitude
                resultList[-1][7:11] = imuUtil.dcm2Quaternion(gnssPvObsEquation.cbn)
                
                resultList[-1][11:] = pk[np.arange(21),np.arange(21)] 
                
                # recursion              
                # self.pk_1 = pk                
                # pk_1 设置为初值
                self.pk_1 = np.diag([0.02,0.02,0.05,
                    0.01,0.01,0.01,
                    unitTools.degree2radian(0.007),unitTools.degree2radian(0.007),unitTools.degree2radian(0.01),
                    self.imuPara.sigma_gb,self.imuPara.sigma_gb,self.imuPara.sigma_gb,                  
                    self.imuPara.sigma_ab,self.imuPara.sigma_ab,self.imuPara.sigma_ab,                    
                    self.imuPara.sigma_gs,self.imuPara.sigma_gs,self.imuPara.sigma_gs,                   
                    self.imuPara.sigma_as,self.imuPara.sigma_as,self.imuPara.sigma_as                  
                    ])                  
                
                self.positionK_1 = resultList[-1][1:4]
                self.velocityK_1 = resultList[-1][4:7]
                self.eulerAngle_k_1 = imuUtil.quaternion2EulerAngle(resultList[-1][7:11])          
            
        return resultList   
        
    def result2File(self,outputFilleDir,resultList):         
        fo = open(outputFilleDir, "w+")
        for r in range(0,len(resultList),100):
            record = resultList[r]
            msg = ""
            
            # time 
            msg += str(record[0]) + ","           
            latitude = unitTools.radian2degree(record[1])
            logitude = unitTools.radian2degree(record[2])
            # logitude,latitude
            msg += str(latitude) + "," + str(logitude) + "," 
            
            # other fileds           
            for f in range(3,len(record)):
                msg += str(record[f]) + ","
                
            fo.write(msg + "\n")           
               
        fo.close() 
    
        
   
    
        