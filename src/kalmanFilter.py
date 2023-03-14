import numpy as np

class KalmanFilter:
    # xK = phi * xK_1 + gamma * wK_1
    # zK = Hk  * xK   + vk
    #
    # xK_1:     n*1, state vector
    # pK_1:     n*n, 
    #
    # phi:      n*n, stateTransformMatrix
    # gamma:    n*s,
    # Q:        s*s, system noise variance matrix
    #
    # Z:        m*1, measure vector
    # H:        m*n, measure matrix
    # R:        m*m, measure noise variance matrix
    def __init__(self,xK_1,pK_1,
                 phi,gamma,QK_1,
                 Z,H,R) -> None:
        self.x_1 = xK_1
        self.p_1 = pK_1
        
        self.phi = phi
        self.gamma = gamma
        self.QK_1 = QK_1
        self.Z = Z
        self.H = H
        self.R = R
        
    def timeUpdate(self):
        xk_predict = self.phi @ self.x_1
        
        pk_predict = self.phi @ self.p_1 @ (self.phi).T + self.gamma @ self.QK_1 @ (self.gamma).T
        
        return xk_predict,pk_predict
    
    def measureUpdate(self,xk_predict,pk_predict):
        # (xk_predict,pk_predict) = self.timeUpdate
        
        K = pk_predict @ (self.H).T @ np.linalg.inv(self.H @ pk_predict @ (self.H).T + self.R)
        
        xk = xk_predict + K @ (self.Z - self.H @ xk_predict)
        
        n = len(self.x_1)
        pk = (np.eye(n) - K @ self.H) @ pk_predict @ (np.eye(n) - K @ self.H).T + K @ self.R @ (K).T
        
        return xk,pk
    
           
        
        
          