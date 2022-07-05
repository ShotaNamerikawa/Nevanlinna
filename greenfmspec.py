import numpy as np
from scipy import integrate
import sys
sys.path.append("/home/shota/wannier_utils/src")
from wannier_utils.temperature_spir import Temperature

def boson_sp(x,para):
    return 2*np.sinh(x)/(np.cosh(x))**3

def fermion_sp(x,para):
    return 1/(np.cosh(2*x))**2

def gausiann(x,para):
    return 1/(np.sqrt(2*np.pi)*para[1])*np.exp(-1/(2*para[1]**2)*(x-para[0])**2)

class Convert:
    '''
    Args:
        spectrum(function)

    '''
    def __init__(self,spectrum,temp,omega_max,file,range=[-10,10],statistics="F"):
        self.spectrum = spectrum
        self.range = np.array(range)
        self.temp = temp
        self.omega_max = omega_max
        self.T = Temperature(temp,self.omega_max)
        if statistics == "F":
            self.points = self.T.omegaF_sp
        else:
            self.points = self.T.omegaB_sp
        self.integral()
        self.write(file)


    def integrand_real(self,x,para):
        return -x/(para[0]**2 + x**2)*self.spectrum(x,para[1:])

    def integrand_imag(self,x,para):
        return -para[0]/(para[0]**2 + x**2)*self.spectrum(x,para[1:])

    def integral(self):
        self.green = np.zeros([self.points.shape[0],2],dtype=np.float64)
        for x in np.arange(self.points.shape[0]):
            self.green[x,0] = integrate.quad(self.integrand_real,self.range[0],self.range[1],([self.points[x],0,0.01]))[0] 
        for x in np.arange(self.points.shape[0]):
            self.green[x,1] = integrate.quad(self.integrand_imag,self.range[0],self.range[1],([self.points[x],0,0.01]))[0] 
    
    def write(self,file):
        print(file)
        np.savetxt(file,np.hstack([self.points[:,None],self.green]))   
    
if __name__ == "__main__":
    test = Convert(gausiann,110,10,"green_gauss_delta.dat")