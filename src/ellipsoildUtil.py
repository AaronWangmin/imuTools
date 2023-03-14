import numpy as np

class Ellipsoild:
    def __init__(self,
                 name = "GRS80",
                 a = 6378137.0,
                 f = 1/298.257222101,
                 wie = 7.292115E-5,
                 GM = 3.986005E+14,
                 gamma_a = 9.7803267715,
                 gamma_b = 9.8321863685) -> None:
        self.name = name
        self.a = a
        self.f = f
        self.wie = wie
        self.GM = GM
        self.gamma_a = gamma_a
        self.gamma_b = gamma_b
        
    # compute semi-minor axis of ellipsoid
    def computeB(self):
        b = self.a * (1 - self.f)
        return b

    # compute the First Eccentricity of ellipsoid
    def computeE(self):
        b = self.computeB()
        e = np.power(self.a*self.a - b*b,0.5) / self.a
        return e

    # compute Radius of meridian curvature,子午圈曲率半径
    def computeRm(self,latitude):
        e = self.computeE()
        Rm = self.a * (1 - e *e) / np.power(1 - e * e * np.sin(latitude) * np.sin(latitude),3/2)
        return Rm

    def computeRn(self,latitude):
        e = self.computeE()
        Rn = self.a / np.power(1 - e * e * np.sin(latitude) * np.sin(latitude),1/2)
        return Rn
    
    # compute the gravity
    # output: 3*1 vector, N/E/D 
    def computeGravity(self,latitude,height = 0.0):
        b = self.computeB()
        m = np.power((self.wie),2)  * np.power((self.a),2) * b / self.GM
        gamma_phi = (self.a * self.gamma_a * np.power(np.cos(latitude),2)  + b * self.gamma_b * np.power(np.sin(latitude),2)) / \
            np.power(np.power(self.a,2) * np.power(np.cos(latitude),2) + b * b * np.power(np.sin(latitude),2),0.5)
        gravity = gamma_phi * (1 - 2/self.a * (1 + self.f + m - 2 * self.f * np.power(np.sin(latitude),2)) * height \
            + 3/np.power(self.a,2) * height * height)
        
        return np.array([0.0,0.0,gravity])    