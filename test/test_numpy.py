
import numpy as np

a = np.array([[1,2,3],[4,5,6]])
print(a)
print(a.T)

b = np.array([10,20,30])
print(b)
print(b.T)
print(np.shape(b.T))

c = np.array([[1],[2],[3]])
print(c)

print(7.292115 * np.power(10.0,-5))


def test_write():
    positionResult = []
    
    for i in range(10):
        positionCurrent = [1,2,3]
        positionResult.append(positionCurrent)

    fo = open("E:\\project_201903\\ecarx-tools-202106\\imuTools\\data\\20220706_100C\\driver_20220706145455_RZC_SHANGHAI_AGB7707_imucompute.txt", "w+")
    for info in positionResult:
        fo.write(str(info[0]) + "," + str(info[1]) + "," + str(info[2]) + "," + "\n")
    fo.close()   
    
    
def test_forInRange():
    for i in range(2,10):
        print(i)
        
def testEye():
    identify = np.eye(3)
    print(identify)
    
def testOnes():
    z1 = np.zeros(3)
    z2 = np.zeros((2,3))
    
    print(z2)
    
def testSize():
    a = [1,2,4]
    b = len(a)
    print("hello")
    
def testNumpyArrya():
    z = []
    
    pos = [1,2,3]
    
    vel = [4,5,6]
    
    # z.append(pos)
    # z.append(vel)
    
    z = np.concatenate((pos,vel))
    
    print(z)
    
def testMultiply():   
    rmsPos = [1,2,3]  
    rmsVel = [4,5,6]    
    r = np.concatenate((rmsPos,rmsVel)).reshape(6,-1)
    # rt = r.T
    # r_col = np.ndarray((6,1))  
    # r_col = r
    R = r @ r.T
    return R
    
def test_arrayInsert():
    a = np.array([1,2,3])
    a = [1,2,3]
    # b = np.array([4,5,6])
    a = np.insert(a,0,5)
    print(a)
    
def test_lastOne():
    a = np.array([1,2,3])    
    b = a[-2]
    print(b)
    
def test_diagA():
    a = [1,2,3]
    b = np.diag([4,5,6])
    print(b)
    
def test_lastArray():
    a = np.array([[1,2,3],[4,5,6],[7,8,9]])
    b = a[-2][1:]

    print("hello")
            

    
    
    
    