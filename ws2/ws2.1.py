#!/usr/bin/env python

import numpy as np

# initialize array of xn's
x = np.zeros(16, dtype=np.float32)
x[0] = np.float32(1)
x[1] = np.float32(1.0/3.0)

thirteenthirds = np.float32(13.0/3.0)
fourthirds = np.float32(4.0/3.0)

def f(a, b):
    return np.float32(thirteenthirds*a) - np.float32(fourthirds*b)

#implement recursion
for i in np.arange(2, 16):
    x[i] = np.float32(f(x[i-1], x[i-2]))
    
print 'using float32 numbers'
for i in np.arange(x.size):
    print x[i]

print 'using float64 numbers'
print (1.0/3.0)**15

#errors
print 'abs error is ', (1.0/3.0)**15 - x[15]
print 'rel error is ', ((1.0/3.0)**15 - x[15]) / (1.0/3.0)**15

    


