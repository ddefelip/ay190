#!/usr/bin/env python

import numpy as np

def dft(x):
    
    i = np.complex(0, 1)  # lets "i" be the imaginary number "i"
    w = np.exp(-2*np.pi*i / x.size)
    m = np.zeros((x.size, x.size), dtype=np.complex)
        # initializes transformation matrix
    
    for j in np.arange(x.size):
        for k in np.arange(x.size):
            m[j][k] = w**(j*k)
            
    return np.dot(m, x) # matrix multiplication
    
    

    