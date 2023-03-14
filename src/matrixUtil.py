import numpy as np

# compute the vector-projection of alpha on beta
# 
def projection(alpha,beta):
    return (np.inner(alpha,beta) / np.inner(beta,beta)) * beta

# normalization column matrix
# 
def normalization(inputMatrix):
    (n,m) = np.shape(inputMatrix)
    for i in range(m):
        # normL2 = np.sqrt(inputMatrix[:,i] @ inputMatrix[:,i])
        normL2 = np.linalg.norm(inputMatrix[:,i])        
        inputMatrix[:,i] = inputMatrix[:,i] * (1/normL2)        
        i += 1
    
    return inputMatrix

# normalization a vector    
def vectorNormalization(a):
    normL2 = np.sqrt(a @ a)
    return a / normL2

# Gramschitd unit-orthgonal
#  
def gramschitd(inputMatrix):   
    (n,m) = np.shape(inputMatrix)
    beta = np.ndarray((n,m))
        
    for i in range(m):
        beta[:,i] = inputMatrix[:,i]
        
        # from 2th colunm compute
        if i > 0:                     
            for j in range(i):                            
                beta[:,i] = beta[:,i] - projection(inputMatrix[:,i],beta[:,j])               
    
    return normalization(beta)   

# get a antisysmmetric matrix of a 3*1 vector
# 
# input: 3*1 vector,
# output:3*3 matrix
def getASM(xyzVector):
    x = xyzVector[0]
    y = xyzVector[1]
    z = xyzVector[2]

    asm = np.zeros([3,3])
    # column 1
    asm[1,0] = z
    asm[2,0] = -y
    
    # column 2
    asm[0,1] = -z
    asm[2,1] = x

    # column 3
    asm[0,2] = y
    asm[1,2] = -x

    return asm        

def getVectorFromAsm(asm):
    xyzVector = np.array([-asm[1,2], asm[0,2],-asm[0,1]])
    return xyzVector

# power of a antisysmmetric matrix
# 
# input: 3*3 vector,
# output:3*3 matrix
def pow_asm(asm,n):
    xyzVector = getVectorFromAsm(asm)
    norm_xyzVector = np.linalg.norm(xyzVector)
    if n % 2 == 1:  
        return np.power(-1,(n-1)/2) * np.power(norm_xyzVector,n-1) * (asm)
    else:
        return np.power(-1,(n-2)/2) * np.power(norm_xyzVector,n-2) * (asm @ asm)

# exponents of a antisysmmetric matrix
# 
# input: 3*3 vector,
# output:3*3 matrix
def exponents_asm(asm):
    xyzVector = getVectorFromAsm(asm)
    norm_xyzVector = np.linalg.norm(xyzVector)
    
    exponentsAsm = np.ones((3,3)) + np.sin(norm_xyzVector) / norm_xyzVector * asm + \
        (1 - np.cos(norm_xyzVector)/np.power(norm_xyzVector,2)) * (asm @ asm)
    
    return exponentsAsm
    
  
# quarternion multipy
# input:
def quarternionMultiply(q1,q2):
    M_q1 = np.zeros([4,4])
    M_q1[0,:] = - q1
    M_q1[:,0] = q1
    M_q1[1:,1:] = getASM(q1[1:])
    M_q1[1,1] = M_q1[2,2] = M_q1[3,3] = q1[0]
        
    return M_q1 @ q2


    


