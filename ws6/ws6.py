#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as pl
from timeit import timeit
from dft import dft
import pylab

pl.clf()

# make sure dft(x) and np.fft.fft(x) output the same thing.
# when I run this section, they do.

test = pylab.randn(5)
myfunc = dft(test)
npfunc = np.fft.fft(test)

print 'dft gives:'
for item in myfunc:
    print item
    
print 'np.fft.fft gives:'
for item in npfunc:
    print item



# compare times for dft and ff

N = np.arange(10, 101)  #sizes of the array to be transformed
t1 = np.zeros(N.size)
t2 = np.zeros(N.size)

# code from worksheet to use "timeit"
for n in N:
    t1[n-10] = timeit("dft(x)", number=100, setup="from dft import dft; import pylab; x=pylab.randn(%d)" % n)
    t2[n-10] = timeit("np.fft.fft(x)", number=100, setup="import numpy as np; import pylab; x=pylab.randn(%d)" % n)
    


# part b: plot time for dft and fft vs n
mydft, = pl.plot(N,t1,'r',linewidth=2)  
fft, = pl.plot(N,t2,'g',linewidth=2)

pl.xlim(min(N),max(N))
pl.ylim(min(t2)/1.05, max(t1)*1.05)

pl.title('Speed of Discrete Fourier Transform vs Size of array')
pl.xlabel("N")
pl.ylabel("t (seconds)")

pl.legend((mydft,fft), ("My dft function","NumPy's fft function"), loc=(0.03,0.80), frameon=False)
pl.plot()

pl.savefig('ws6-b.pdf')



# plot time on a log scale to highlight time differences better
pl.clf()

mydft, = pl.plot(N,np.log10(t1),'r',linewidth=2)  
fft, = pl.plot(N,np.log10(t2),'g',linewidth=2)
pl.plot()  
    
pl.xlim(min(N),max(N))
pl.ylim(min(np.log10(t2))*1.05, max(np.log10(t1))*1.05)

pl.title('Speed of Discrete Fourier Transform vs Size of array')
pl.xlabel("N")
pl.ylabel("log10(t)")

pl.legend((mydft,fft), ("My dft function","NumPy's fft function"), loc=(0.50,0.40), frameon=False)
pl.plot()

pl.savefig("ws6-c.pdf")
