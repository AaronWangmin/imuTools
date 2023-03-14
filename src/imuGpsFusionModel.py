import numpy as np

from . import ellipsoildUtil
from . import matrixUtil
from . import imuUtil

# imu error state vector: 21*1 vector, [dr,dv,phi,bg,ba,sg,sa]
#   
#   dr                : 3*1 vector,N-E-D position error in navigation system
#   dv                : 3*1 vector,N-E-D velocity error in navigation system
#   phi               : 3*1 vector,attitude error in navigation system#
#   bg / ba / sg / sa : 3*1 vecotr,imu  groy bias / acce bias / groy scale error / acce scale error
# 

class ImuErrorStateEquation:
    
    # postion / velcocity / attitude: 3*1 vector,in navigation system
    #   position: longitude / latitude / height
    #   velocity: North / East / Down
    #   attitude: 3*3,dcm
    # fb / wib                      : 3*1 vector, imu measure data
    # Tgb,Tab,Tgs,Tas               : 3*1 vector
    def __init__(self,numberErrorStateVarible,numberWhiteNoiseVariable,
                 positionK_1,velocityK_1,cbn_k_1,cbn_k,fb,wib,tk,tk_1,
                 pk_1,                
                 imuPara,
                 ellipsoild) -> None:
        
        self.numberErrorStateVarible = numberErrorStateVarible
        self.numberWhiteNoiseVariable = numberWhiteNoiseVariable
        
        self.position = positionK_1
        self.velocity = velocityK_1
        self.cbn_k_1 = cbn_k_1
        self.cbn_k = cbn_k
        self.fb = fb
        self.wib = wib
        self.tk = tk
        self.tk_1 = tk_1
         
        # covariance of error state variable
        self.pk_1 = pk_1 
        self.imuPara = imuPara       
        self.ellipsoild = ellipsoild
        
        # 临时计算参数
        self.latitude = self.position[0]
        self.height = self.position[2]
        self.velocityN = self.velocity[0]
        self.velocityE = self.velocity[1]
        self.velocityD = self.velocity[2]
        self.wie = self.ellipsoild.wie
        
        self.rm = self.ellipsoild.computeRm(self.latitude)
        self.rn = self.ellipsoild.computeRn(self.latitude)
        self.g = self.ellipsoild.computeGravity(self.latitude,self.height)[2]
        
        # 计算并设置各系数矩阵
        # phi matrix, n*n
        self.phi = np.eye(self.numberErrorStateVarible)
        # gamma matrix, n*s
        # self.gamma     = np.zeros((self.numberErrorStateVarible,self.numberWhiteNoiseVariable))  
        self.gamma     = np.eye(self.numberErrorStateVarible)      
        # covariance matrix of state-variable, n*n
        self.pk_1         = pk_1
        # covariance matrix of white-noise, n*n   
        self.Q      = np.eye(self.numberErrorStateVarible) * 0.0    
        
        # kf time update: 
        # predictX: 21 * 1, x/y/z (m),/vx/vy/xz (m/s),phix/phiy/phiz (rad),
        self.predictX = np.zeros(self.numberErrorStateVarible)
        self.predictP = np.eye(self.numberErrorStateVarible) * 0.0              
    
    # deltaX(k) = phi(k/k-1) * deltaX(k_1) + gamma * w(k_1) 
    # deltaX(k) = phi(k/k-1) * deltaX(k_1) + w(k_1) 
    def setErrorStateEquation(self):
        self.setPhiMatrix()
        # self.setGammaMatrix()
        # self.setPMatrix()
        self.setQMatrix()         
    
    # derivatives(deltaX(t)) = F(t) * delta(t) + G(t) * w(t)
    # set F-matrix   
    def comupterFMatrix(self):
        Frr = self.computerFrr()
        Fvr = self.computerFvr()
        Fvv = self.computerFvv()
        Fphir = self.computerFphir()
        Fphiv = self.computerFphiv()
        
        F = np.zeros((self.numberErrorStateVarible,self.numberErrorStateVarible))
        F[0:3,0:3] = Frr
        F[0:3,3:6] = np.eye(3)
        
        F[3:6,0:3] = Fvr
        F[3:6,3:6] = Fvv
        F[3:6,6:9] = matrixUtil.getASM(self.cbn_k @ self.fb)
        F[3:6,12:15] = self.cbn_k
        F[3:6,18:21] = self.cbn_k @ np.diag(self.fb)
        
        F[6:9,0:3] = Fphir
        F[6:9,3:6] = Fphiv
        win_n = imuUtil.computerOmega_in_n(self.latitude,self.height,self.velocityN,self.velocityE,self.ellipsoild)
        F[6:9,6:9] = -matrixUtil.getASM(win_n)
        F[6:9,9:12] = -self.cbn_k
        # ???? wib
        F[6:9,15:18] = - self.cbn_k @ np.diag(self.wib)
        
        F[9:12,9:12]   = (-1 / self.imuPara.Tgb) * np.eye(3)         
        F[12:15,12:15] = (-1 / self.imuPara.Tab) * np.eye(3)
        F[15:18,15:18] = (-1 / self.imuPara.Tgs) * np.eye(3)
        F[15:18,15:18] = (-1 / self.imuPara.Tas) * np.eye(3)
        
        return F   
    
    # position error in navigation coordinate, m
    def computerFrr(self):
        Frr = np.zeros((3,3))
        
        Frr[0,0] = - self.velocityD / (self.rm + self.height)
        Frr[0,2] =   self.velocityN / (self.rm + self.height)
        
        Frr[1,0] =   self.velocityE * np.tan(self.latitude) / (self.rn + self.height)
        Frr[1,1] = -(self.velocityD + self.velocityN * np.tan(self.latitude)) / (self.rn + self.height)
        Frr[1,2] =   self.velocityE / (self.rn + self.height)
        
        return Frr    
    
    def computerFvr(self):
        Fvr = np.zeros((3,3))
        
        Fvr[0,0] = - 2 * self.velocityE * self.wie * np.cos(self.latitude) /  (self.rm + self.height) \
            - np.power(self.velocityE / np.cos(self.latitude),2) / ((self.rm + self.height) * (self.rn + self.height))
        Fvr[0,2] =  self.velocityN * self.velocityD / np.power(self.rm + self.height,2) \
            - np.power(self.velocityE,2) * np.tan(self.latitude) / np.power(self.rn + self.height,2)
            
        Fvr[1,0] =   2 * self.wie * (self.velocityN * np.cos(self.latitude) - self.velocityD * np.sin(self.latitude)) / (self.rm + self.height) \
            + self.velocityN * self.velocityE * np.power(1 / np.cos(self.latitude),2) / ((self.rm + self.height) * (self.rn + self.height))
        Fvr[1,2] =  (self.velocityE * self.velocityD + self.velocityN * self.velocityE * np.tan(self.latitude)) \
            / np.power(self.rn + self.height,2)   
            
        Fvr[2,0] =   2 * self.wie * self.velocityE * np.sin(self.latitude) / (self.rm + self.height)
        Fvr[2,2] = - np.power(self.velocityE,2) / np.power(self.rn + self.height,2) \
            - np.power(self.velocityN,2) / np.power(self.rm + self.height,2) \
            + 2 * self.g / (np.power(self.rm * self.rn,0.5) + self.height)
            
        return Fvr
    
    def computerFvv(self):
        Fvv = np.zeros((3,3))
        
        Fvv[0,0] =  self.velocityD / (self.rm + self.height)
        Fvv[0,1] = -2 * (self.wie * np.sin(self.latitude) + self.velocityE * np.tan(self.latitude) / (self.rn + self.height))    
        Fvv[0,2] = self.velocityN / (self.rm + self.height)
        
        Fvv[1,0] =  2 * self.wie * np.sin(self.latitude) + self.velocityE * np.tan(self.latitude) / (self.rn + self.height)
        Fvv[1,1] =  (self.velocityD + self.velocityN * np.tan(self.latitude)) / (self.rn + self.height)
        Fvv[1,2] =  2 * self.wie * np.cos(self.latitude) + self.velocityE / (self.rn + self.height)
        
        Fvv[2,0] = -2 * self.velocityN / (self.rm + self.height)
        Fvv[2,1] = -2 * (self.wie * np.cos(self.latitude) + self.velocityE / (self.rn + self.height))
        
        return Fvv
    
    def computerFphir(self):
        Fphir = np.zeros((3,3))
        
        Fphir[0,0] = -self.wie * np.sin(self.latitude) / (self.rm + self.height)
        Fphir[0,2] =  self.velocityE / np.power(self.rn + self.height,2)
        
        Fphir[1,2] = -self.velocityN / np.power(self.rm + self.height,2)
        
        Fphir[2,0] = -self.wie * np.cos(self.latitude) / (self.rm + self.height) \
            - self.velocityE * np.power(1/np.cos(self.latitude),2) / ((self.rm + self.height) * (self.rn + self.height))
        Fphir[2,2] = -self.velocityE * np.tan(self.latitude) / np.power(self.rn + self.height,2)

        return Fphir
    
    def computerFphiv(self):
        Fphiv = np.zeros((3,3))
        
        Fphiv[0,1] =  1 / (self.rn + self.height)
        
        Fphiv[1,0] = -1 / (self.rm + self.height)
        
        Fphiv[2,1] = - np.tan(self.latitude) / (self.rm + self.height)
        
        return Fphiv
       
    # deltaX(k) = phi(k/k-1) * deltaX(k_1) + w(k_1)
    # phi(k/k-1)= I + F(k-1) * deltaT
    def setPhiMatrix(self):
         FM = self.comupterFMatrix()
         self.phi = np.eye(self.numberErrorStateVarible) + FM * (self.tk - self.tk_1)
         
    def setGammaMatrix(self):
        # self.gamma[3:21,0:18] = np.eye(self.numberWhiteNoiseVariable)  
        pass           
  
    # input:    P,errorStateCovarianceMatrix
    #           Q,whiteNoiseCovarianceMatrix
    def setPMatrix(self):
        pass    
    
    def computerGMatrix(self,cbn):
        GM = np.zeros((self.numberErrorStateVarible,self.numberWhiteNoiseVariable))
        GM[3:6  ,0:3] = cbn
        GM[6:9  ,3:6] = cbn
                
        GM[9:,6:] = np.eye(12)
        return GM
    
    # ??? VRW,ARW, sigma
    def computerPowerSpectralDensityQMatrix(self):
        psdm = np.eye(self.numberWhiteNoiseVariable)
        psdm[0:3,0:3] = np.power(self.imuPara.vrw,2) * np.eye(3)
        psdm[3:6,3:6] = np.power(self.imuPara.arw,2) * np.eye(3)
        psdm[6:9,6:9]     = (2 * np.power(self.imuPara.sigma_gb,2) / self.imuPara.Tgb) * np.eye(3)
        psdm[9:12,9:12]   = (2 * np.power(self.imuPara.sigma_ab,2) / self.imuPara.Tab) * np.eye(3)
        psdm[12:15,12:15] = (2 * np.power(self.imuPara.sigma_gs,2) / self.imuPara.Tgs) * np.eye(3)
        psdm[15:18,15:18] = (2 * np.power(self.imuPara.sigma_as,2) / self.imuPara.Tas) * np.eye(3)
        return psdm
    
    def setQMatrix(self):
        # phiMatrix = self.setPhiMatrix()
        GMatrix_k_1 = self.computerGMatrix(self.cbn_k_1)
        GMatrix_k = self.computerGMatrix(self.cbn_k)
        psdm = self.computerPowerSpectralDensityQMatrix()
        
        QMatrix = 1 / 2 * (self.phi @ GMatrix_k_1 @ psdm @ GMatrix_k_1.T @ self.phi.T \
            + GMatrix_k @ psdm @ GMatrix_k.T) * (self.tk - self.tk_1)
        
        self.Q = QMatrix       
    
    # predict
    def timeUpdate(self,xk_1 = np.zeros(21)):
        xk_predict = self.phi @ xk_1
        
        pk_predict = self.phi @ self.pk_1 @ (self.phi).T + self.Q       
        
        return xk_predict,pk_predict       
  

# deltaZ = H * deltaX + w        
class GnssPVObsEquation:
    # input:
    #       position/velocity/cbn   : imu computer value              
    #       leverArm                : 3*1 vector, imu to gnss antenna, x/y/down
    #       gpsPose/gpsVel          : gps value
    #       wib_b                   :
    #       ellipsoild              :
    def __init__(self,numberMeasureVariable,numberStateVariable,position,velocity,cbn,leverArm,gpsPos,gpsVel,rmsGpsPos,rmsGpsVel,wib_b,ellipsoild) -> None:
        self.numberMeasureVariable = numberMeasureVariable
        self.numberStateVariable = numberStateVariable
        
        self.position = position
        self.velocity = velocity        
        self.cbn = cbn
        
        self.leverArm = leverArm
        self.gpsPose = gpsPos
        # m
        self.gpsVel = gpsVel        
        self.rmsGpsPos = rmsGpsPos
        self.rmsGpsVel = rmsGpsVel        
        
        self.ellipsoild = ellipsoild
        self.win_n = imuUtil.computerOmega_in_n(self.position[1],self.position[2],self.velocity[0],self.velocity[1],self.ellipsoild)
        self.wib_b = wib_b
        
        # Z vector: 
        self.ZVector = []        
        # H matrix: 6*21
        self.HMatrix = np.zeros((numberMeasureVariable,numberStateVariable))
        # R matrix:
        self.RMatrix = np.zeros((numberMeasureVariable,numberMeasureVariable))
        
        # temp 
        self.DMatrix = np.zeros((3,3))
    
    def setGnssPVObsEquationPara(self):
        self.xyz2LonLatHeight()
        
        self.computeZVector()
        self.computeHr()
        self.computeHv()
        self.setRMatrix() 
        
    # Gnss observation vector    #  
    def computeZVector(self):
        # gnss pos-obs vector, Logitude/Latitude/Height       
        gpsPosByImu = self.position + self.DMatrix @ self.cbn @ self.leverArm               
        # difference between gpsPosByImu and gpsPose, N/E/D,(m)
        ZPos = np.linalg.inv(self.DMatrix) @ (gpsPosByImu - self.gpsPose)        
        
        # gnss vel-obs vector, N/E/D (m/s)
        gpsVelByImu = self.velocity - matrixUtil.getASM(self.win_n) @ self.cbn @ self.leverArm \
            - self.cbn @ (matrixUtil.getASM(self.leverArm) @ self.wib_b)
        ZVel = gpsVelByImu - self.gpsVel
        
        self.ZVector = np.concatenate((ZPos,ZVel))  
   
    # Gnss pos obs coefficient-matrix    
    def computeHr(self):        
        self.HMatrix[0:3,0:3] = np.eye(3) 
        
        self.HMatrix[0:3,6:9] = (self.cbn @ matrixUtil.getASM(self.leverArm)) 
    
    # Gnss velocity obs coefficient-matrix           
    def computeHv(self):
        self.HMatrix[3:6,3:6]   = np.eye(3)
        
        self.HMatrix[3:6,6:9]   = - matrixUtil.getASM(self.win_n) @ (self.cbn @ matrixUtil.getASM(self.leverArm)) \
            - self.cbn @ matrixUtil.getASM(np.cross(self.leverArm,self.wib_b))
            
        self.HMatrix[3:6,9:12]  = - self.cbn @ matrixUtil.getASM(self.leverArm)
        
        self.HMatrix[3:6,15:18] = - self.cbn @ matrixUtil.getASM(self.leverArm) @ np.diag(self.wib_b)
       
    # observation white noise matrix
    # input:    rmsPos: 3*1 vector, m
    #           rmsVel: 3*1 vector, m/s
    #       
    def setRMatrix(self):        
        r = np.concatenate((self.rmsGpsPos,self.rmsGpsVel)).reshape(6,-1)
        R = r @ r.T
        self.RMatrix = R
    
    # transfer meter(NED) to LonLatHeight(navigation)    
    def xyz2LonLatHeight(self):
        rm = self.ellipsoild.computeRm(self.position[1])
        rn = self.ellipsoild.computeRn(self.position[1])
        DMatrix = np.diag([1 / (rm + self.position[2]),
                     1 / ((rn + self.position[2]) * np.cos(self.position[1])),
                     -1])                           
       
        self.DMatrix = DMatrix
        
    def measureUpdate(self,pk_predict,xk_predict = np.zeros(21)):
        K = pk_predict @ (self.HMatrix).T @ np.linalg.inv(self.HMatrix @ pk_predict @ (self.HMatrix).T + self.RMatrix)
        
        xk = xk_predict + K @ (self.ZVector - self.HMatrix @ xk_predict)
        
        n = self.numberStateVariable
        pk = (np.eye(n) - K @ self.HMatrix) @ pk_predict @ (np.eye(n) - K @ self.HMatrix).T + K @ self.RMatrix @ (K).T        
        
        # update pva
        self.position = self.position - self.DMatrix @ xk[0:3]
        self.velocity = self.velocity - xk[3:6]
        self.cbn = np.linalg.inv(np.eye(3) - matrixUtil.getASM(xk[6:9])) @ self.cbn 
        
        # update bg/ba/sg/sa        
                
        return xk,pk
        
    
    
    
        
        
            
        
    
    