import numpy as np
import scipy as sp
import scipy.integrate
import scipy.misc
#import gmpy2 as gp2
from mpmath import *
#from decimal import *
mp.prec = 150

a = mp.mpc('0.0000')
I = mp.mpc('1')*1j
#b = mp.mpf('1.6732058845458435945971158485')
#a = mp.mpf('')
#b = mp.mpf('')
print((a - I)/(a + I))
print(mp.mpf(0))
#nprint(b/a,n=80,strip_zeros=False)
#print(b/a+0.1)
#nprint(b/a+0.1,n=80,strip_zeros=False)

#gp2.set_context(gp2.context())
#gp2.get_context()
#a = gp2.mpfr('4.5820048684854234019414943230')
#b = gp2.mpfr('1.6732058845458435945971158485')
#a = gp2.mpfc('3.001000000000000000000000000000')
#b = gp2.mpfr('1.000000000000000000000000000000')
#print(a.digits())
#c = gp2.div(b,a)
#print("a = ")
#print(a)
#print(a.digits())
#print("b = ")
#print(b)
#print("b/a = ")
#print(c)
#print(c.precision)
#print(c.as_mantissa_exp())
#print(c.digits())
#d = gp2.div(a,b)*c
#print(d)
#print(d)
#print(d.precision)
#print(d.as_mantissa_exp())
#print(d.digits())
#getcontext().prec = 112
#print(Decimal(34972394723047274923749724097230)/Decimal(12743977492374973497479239749211))





