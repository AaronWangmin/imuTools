import sys
sys.path.insert(0,"E:\\project_201903\\ecarx-tools-202106\\imuTools")

import numpy as np

from src import matrixUtil

def test_gramschitd():    
    gramschitd_before = np.array([[1,2,-1],[-1,3,1],[4,-1,0]]).T
    gramschitd_after = matrixUtil.gramschitd(gramschitd_before)
    gramschitd_check = np.array([[0.408,0.816,-0.408],[-0.577,0.577,0.577],[0.707,0,0.707]]).T
    assert(np.max(np.abs(np.subtract(gramschitd_check,gramschitd_after))) < 0.001)
    
def test_vectorrmalization():
    a = np.array([1,2,3,4])
    a = matrixUtil.vectorNormalization(a)
    a_test = [0.182574186,0.365148372,0.547722558,0.730296743]
    assert(np.max(np.abs(np.subtract(a,a_test))) < 0.001) 

def test_getASM():
    asm = matrixUtil.getASM(np.array([1,2,3]))
    asm_test = np.matrix([[0,-3,2],
                        [3,0,-1],
                        [-2,1,0]])
    assert(np.max(np.abs(np.subtract(asm,asm_test))) < 0.001)
    
def test_quarternionMultiply():
    q1 = np.array([1,2,3,4])
    q2 = np.array([4,5,6,7])
    
    q = matrixUtil.quarternionMultiply(q1,q2)
    q_test = np.array([-52,10,24,20])
    assert(np.max(np.abs(np.subtract(q,q_test))) < 0.001)
    

   
    