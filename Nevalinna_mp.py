import numpy as np
import scipy as sp
import scipy.integrate
import scipy.misc
import sys


def derivative(func,x,para,n=1,dx = 1e-6):
    def wrapper(y):
        return func(y,para)
    return sp.misc.derivative(wrapper, x, n=n, dx=dx)


class Nevalinna:
    '''class to calculate analytic continuation
    refer to PBL 126, 056402
    Nevalinna function: x + yi (y > 0) ->  x' + y'i (y'>= 0)
    
    Args:
        nev_func:Nevalinna function (- retarded Green's function)'s value at sample_points
        points:sample points of Matsubara frequency on upper half plane 
        contra: contractive(C+ -> D) function made from nev_func
    '''
    def __init__(self,sample_points,nev_func,energy_range = [-10,10],delta = 0.1,lamb = 1e-4):
        self.nev_func = nev_func
        self.points = sample_points
        self.p_num = self.points.shape[0]
        self.contra = self.mebius(self.nev_func)
        self.phi = np.zeros([self.p_num],dtype=np.complex128)
        self.phi_n_gamma()
        self.energy_range = energy_range
        self.delta = delta
        self.lamb = lamb

    def mebius(self,z):
        return (z-1j)/(z+1j)
        
    def last_contra(self,z):
        return 0

    def phi_n_gamma(self):
        contras = np.zeros([self.p_num,self.p_num],dtype=np.complex128)
        contras[0,:] = self.contra
        for x in np.arange(self.p_num-1):
            self.phi[x] = contras[x,x]
            a = (self.points[x+1:] - self.points[x])/(self.points[x+1:] - self.points[x].conj())
            c = self.phi[x]*a  
            contras[x+1,x+1:] = -(self.phi[x] - contras[x,x+1:])/(a-c*contras[x,x+1:])

    def interpolate(self,z):
        value = self.last_contra(z)
        for x in np.arange(self.p_num-1,-1,-1):
            a = (z - self.points[x])/(z-self.points[x].conj())
            b = self.phi[x]
            c = a*b
            value = (a*value + b)/(c*value + 1)
        value = 1j*(1 + value)/(1 - value)
        return value

    def interpolate_make(self,z,para,func):
        value = func(z,para)
        for x in np.arange(self.p_num-1,-1,-1):
            a = (z - self.points[x])/(z-self.points[x].conj())
            b = self.phi[x]
            c = a*b
            value = (a*value + b)/(c*value + 1)
        value = 1j*(1 + value)/(1 - value)
        return value

    def hardy_func(self, z, coeff):
        '''
        Args:
            coeff: coeff[:coeff.shape[0]/4] correspond to Re(a_k)
                   coeff[coeff.shape[0]/4:2/4*coeff.shape[0]] correspond to Im(a_k)
                   coeff[2/4*coeff.shape[0]:3/4*coeff.shape[0]] correspond to Re(b_k)
                   coeff[3/4*coeff.shape[0]:coeff.shape[0]] correspond to Im(b_k)
        '''
        k_num = int(coeff.shape[0]/4)
        return np.sum((coeff[:k_num]+1j*coeff[k_num:2*k_num])*Hardy.basis(z, np.arange(0,k_num))\
               + (coeff[2*k_num:3*k_num] + 1j*coeff[3*k_num:4*k_num])*Hardy.basis(z, np.arange(0,k_num)).conj())

    def spectral_func(self,x,para):
        '''Return spectral function value on x + 1j*self.delta
        '''
        return 1/np.pi*self.interpolate_make(x+1j*self.delta,para,self.hardy_func).imag

    def second_derivative(self,x,para):
        return 1/np.pi*derivative(self.spectral_func,x,para,n=2)

    def square_second(self,x,para):
        return self.second_derivative(x,para)**2

    def spectral_integral(self,para):
        print(para)
        return np.abs(1 - sp.integrate.quad(self.spectral_func,self.energy_range[0],\
               self.energy_range[1],args=para)[0])**2\
              +self.lamb*sp.integrate.quad(self.square_second,self.energy_range[0],\
               self.energy_range[1],args=para)[0]

    def optimize(self, k_num = 2, first_guess="random",coefficient = None):
        if type(coefficient) == np.ndarray:
            coefficient = coefficient
        elif first_guess == "random":
            coefficient = np.zeros([4*k_num])
            coefficient = np.random.rand(k_num*4)*np.sqrt(np.pi)/(2*k_num)
        else:
            print("error! Pass coefficient or set first_guess = \"random\"")
            exit
        self.hardy_coeff = sp.optimize.minimize(self.spectral_integral,coefficient,options={'verbose':2},jac="2-point")

class Hardy:
    def __init__(self):
        pass
    @classmethod
    def basis(cls,z,k):
        return 1/(np.sqrt(np.pi)*(z + 1j))*((z - 1j)/(z + 1j))**k


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import temperature_spir

    wmax = 6
    T = 10
    temp = temperature_spir.Temperature(T,wmax)
    sample1 = temp.omegaF_sp
    sample = sample1[int(len(sample1)/2):]*1j
    print(sample)
    #sample = np.pi*1j*np.arange(1,20)/100
    #print(sample)
    def green(x):
        return 1/(x + 0.2 + 0.2j) + 1/(x-0.3 + 0.1j)
    nev_func = -green(sample)
    test = Nevalinna(sample,nev_func,energy_range=[-5,5],delta=0.01)

    k_num = 25
    '''
    coefficient = np.array([1.33512512e-02,4.80458274e-02,7.04075187e-02,2.94869186e-02
,4.88623168e-02,7.66514719e-02,3.25837545e-02,1.48880043e-02
,8.77477746e-02,6.24962947e-02,2.76095069e-02,2.67440589e-05
,9.34805964e-02,4.58229347e-03,4.41375675e-02,3.35952680e-02
,6.62940652e-02,6.70646037e-02,5.29958011e-02,2.25296633e-02
,3.49488139e-02,6.16449270e-02,5.13069225e-02,1.16835728e-02
,6.74102318e-02,7.12937902e-02,4.27429042e-02,4.49172574e-02
,7.09440422e-02,6.00994591e-02,2.97764565e-02,6.25629929e-02
,4.64994115e-02,6.24328615e-02,2.14973116e-03,3.07507013e-02
,6.95656643e-02,9.36340974e-02,5.64026084e-02,2.53032708e-02])
'''
    coefficient = np.random.rand(4*k_num)*np.sqrt(np.pi)/(2*k_num)
#    print(sp.integrate.quad(lambda x:test.spectral_func(x,coefficient),-4,4))
#    coefficient[0:25] = 1
#    coefficient[50:75] = 1
#    print(test.spectral_integral(coefficient))
#    test.optimize(k_num=k_num)
#    coefficient = test.coefficient
    #coefficient = np.zeros(4)
    coefficients = [np.zeros([4*k_num]), np.random.rand(4*k_num)*np.sqrt(np.pi)/(2*k_num), np.array([1.33512512e-02,4.80458274e-02,7.04075187e-02,2.94869186e-02
,4.88623168e-02,7.66514719e-02,3.25837545e-02,1.48880043e-02
,8.77477746e-02,6.24962947e-02,2.76095069e-02,2.67440589e-05
,9.34805964e-02,4.58229347e-03,4.41375675e-02,3.35952680e-02
,6.62940652e-02,6.70646037e-02,5.29958011e-02,2.25296633e-02
,3.49488139e-02,6.16449270e-02,5.13069225e-02,1.16835728e-02
,6.74102318e-02,7.12937902e-02,4.27429042e-02,4.49172574e-02
,7.09440422e-02,6.00994591e-02,2.97764565e-02,6.25629929e-02
,4.64994115e-02,6.24328615e-02,2.14973116e-03,3.07507013e-02
,6.95656643e-02,9.36340974e-02,5.64026084e-02,2.53032708e-02])]
 #   coefficient = np.random.rand(4*k_num)*np.sqrt(np.pi)/(2*k_num)]
    label = ["\Theta_(R+1) = 0","Random coefficient Hardy","Optimized Hardy"]
    #fig,ax = plt.subplots()
    points = 0.001*np.arange(8000) - 4
    value = []   
    color = ["r","g","black"]
    for y in range(2):
        fig,ax = plt.subplots() 
        value = []
        for x in points:
            value.append(test.spectral_func(x,coefficients[y]))
        ax.plot(points,value,label=label[y],lw=1.0,c=color[y])
        ax.plot(points,-green(points).imag/np.pi,label="correct")
    #ax.plot(points,value,label="Nevalinna")
        ax.legend()
        plt.show()



    ''' Whether Hardy = 0 coresponds to Theta_M+1 = 0
    for x in np.arange(10):
        print("Theta_M+1 = 0")
        print(test.interpolate(x))
        print("Theta_M+1 = Hardy")
        print(test.interpolate_make(x,np.array([0,0]),test.hardy_func))
    '''

    '''comparision plot block
    fig,ax = plt.subplots()
    points = 0.001*np.arange(4000) -2
    value = []
    for x in points:
        value.append(test.interpolate(x).imag)
    ax.plot(points,-green(points).imag,label="correct")
    ax.plot(points,value,label="Nevalinna")
    ax.legend()
    plt.show()
    '''