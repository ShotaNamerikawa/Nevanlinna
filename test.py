from Nevanlinna import Nevanlinna, Nevanlinna_optimize,Integral
import temperature_spir
import numpy as np

wmax = 2 # eV
T = 116 # K
k_num = 25 # The number of Hardy basis
temp = temperature_spir.Temperature(T,wmax)
sample = temp.omegaF_sp[int(1/2*temp.omegaF_sp.shape[0]):]
#sample = np.pi*(2*np.arange(36)+1)/100
def green(x):
    return 1/(x + 0.1j) 
sigma = 0.5
mu = 0
'''
def green(x):
    return -1j/np.sqrt(2*np.pi)/sigma*np.exp(-1/2*(x - mu)**2/sigma**2) 
'''
nev_func = -green(sample*1j)
# Theta_(M+1) = 0
p_num = 32000
points = 2/(p_num)*np.arange(wmax*p_num) - wmax
test = Nevanlinna(sample,nev_func,energy_range=[-wmax,wmax],delta=0.001,k_num = k_num)
integ = Integral(points,test.spectral_func)
# Theta_(M+1) = 0
print("random coefficients")
#   fig,ax =plt.subplots()
value = []
coefficient = np.random.rand(k_num*4)*np.sqrt(np.pi)/(2*k_num)
for x in points:
    value.append(test.spectral_func(x,coefficient))


