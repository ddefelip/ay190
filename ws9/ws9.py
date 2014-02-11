#!/usr/bin/env python

import numpy as np
import scipy.sparse.linalg as sp
import scipy.linalg as sl
import matplotlib.pyplot as pl
import time

# load the data files

LSE1_m = np.loadtxt('./LSE1_m.dat')
LSE2_m = np.loadtxt('./LSE2_m.dat')
LSE3_m = np.loadtxt('./LSE3_m.dat')
LSE4_m = np.loadtxt('./LSE4_m.dat')
LSE5_m = np.loadtxt('./LSE5_m.dat')

LSE1_bvec = np.loadtxt('./LSE1_bvec.dat')
LSE2_bvec = np.loadtxt('./LSE2_bvec.dat')
LSE3_bvec = np.loadtxt('./LSE3_bvec.dat')
LSE4_bvec = np.loadtxt('./LSE4_bvec.dat')
LSE5_bvec = np.loadtxt('./LSE5_bvec.dat')

# convenient way to store these different sizes matrices 
# in one object

m = [LSE1_m, LSE2_m, LSE3_m, LSE4_m, LSE5_m]
bvec = [LSE1_bvec, LSE2_bvec, LSE3_bvec, LSE4_bvec, LSE5_bvec]

print 'loaded all files'

# vectors of solutions which we will fill up later

LSE1_x = np.zeros(LSE1_bvec.size)
LSE2_x = np.zeros(LSE2_bvec.size)
LSE3_x = np.zeros(LSE3_bvec.size)
LSE4_x = np.zeros(LSE4_bvec.size)
LSE5_x = np.zeros(LSE5_bvec.size)

x = [LSE1_x, LSE2_x, LSE3_x, LSE4_x, LSE5_x]

# see what the dimension are, and confirm that the determinants 
# are nonzero. The the "slogdet" functino actually returns the
# natural log of the determinant. I used this because the determinants
# were so big that they were causing the normal determinant function
# to return "inf" for the larger matrices.

for i in range(len(m)):
    print 'LSE'+str(i+1)+'_bvec has size',bvec[i].size
    print 'LSE'+str(i+1)+'_m has dimensions',m[i].shape
    print 'determinant of LSE'+str(i+1)+'_m is',np.linalg.linalg.slogdet(m[i])[1]



#### GAUSS ELIMINATION FUNCTIONS###

def gauss(A_matrix, b_vector):
    
    # make sure the matrix is square
    assert(A_matrix.shape[0] == A_matrix.shape[1])
    
    # number of rows (or columns)
    n = A_matrix.shape[0]
    
    for k in np.arange(0, n):
        for j in np.arange(k+1, n):
            # equation to replace j-th equation with on k-th iteration
            # make sure to first change the value in the bvector, 
            # since that depends on unchanged values in the A matrix.
            b_vector[j] = -A_matrix[j,k]/A_matrix[k,k] * b_vector[k] + b_vector[j]
            A_matrix[j] = -A_matrix[j,k]/A_matrix[k,k] * A_matrix[k] + A_matrix[j]

    

def backsubstitute(A_matrix, x_vector, b_vector):
    
    assert(A_matrix.shape[0] == A_matrix.shape[1])
    
    for j in range(x_vector.size)[::-1]:      
        for k in range(j, x_vector.size)[::-1]:
            # the [::-1] returns the range in reverse order
            # which is the order I want it since I am
            # backsubstituting everything
            
            if k == j:
                # the end of the backsubstituting chain for a 
                # row j, based on III.10.5 in the notes
                
                x_vector[j] += b_vector[j]
                x_vector[j] *= (1.0 / A_matrix[j][j])
                
            else:
                
                x_vector[j] -= ( A_matrix[j][k] * x_vector[k] )
    

# Test Gaussian Elimination functions on 3x3 matrix

testA = np.array([[2., 5., -3.], [1., -7., 4.], [-6., -2., 2.]])
testb = np.array([3., -1., -4.])
testx = np.zeros(3)

gauss(testA, testb)
backsubstitute(testA, testx, testb)

print 'actual answer is x = [1, 2, 3]'
print 'result from test:',testx



# RUN NumPy solver on the 5 matrices,
# timing how long it takes for each one 
print 'NumPy\'s "solve" function'
for i in range(5):
    
    start = time.time()
    
    np.linalg.solve(m[i], bvec[i]) 
    
    stop = time.time()
    print 'solved LSE matrix',i+1,'in',stop-start,'seconds'  
    
    
# RUN SciPy solver in the same way   
print 'SciPy\'s "solve" function'
for i in range(5):
    
    start = time.time()
    
    sl.solve(m[i], bvec[i]) 
    
    stop = time.time()
    print 'solved LSE matrix',i+1,'in',stop-start,'seconds'    


# RUN Gaussian Elimination 
print 'Gaussian Elimination'
for i in range(5):
    
    start = time.time()
    
    gauss(m[i], bvec[i])
    backsubstitute(m[i], x[i], bvec[i])
    
    stop = time.time()
    print 'solved LSE matrix',i+1,'in',stop-start,'seconds'

    