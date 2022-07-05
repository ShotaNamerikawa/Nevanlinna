'''Nevanlinna Module.
This Module contains Classes and a function to calculate 
analytic continuation of Nevanlina functions(C_+ -> C^-_+) which contain
Green functions on C_+

Main class:Nevanlinna
Sub class:Hardy
Function:derivative

examples:
    
'''
import numpy as np
import scipy as sp
import scipy.integrate
import scipy.misc


def derivative(func,x,*args,n=1,dx = 1e-6,**keyargs):
    '''This function returns n-th order derivative
    This function returns n-th order derivative 
    at x point of functions which have one argument 
    and several parameters.

    Args:
        func:A function to be differentiated.
        x:The point where derivative is calculated.
        para:all parameter values func includes.
        n:The class number of derivative. Default is 1.
        dx:Distance between points to be calculate derivative.
    '''
    def wrapper(y):
        return func(y,*args,**keyargs)
    return sp.misc.derivative(wrapper, x, n=n, dx=dx)

class Hardy:
    '''Class to calculate Hardy basis.
    '''
    def __init__(self):
        pass

    @classmethod
    def basis(cls,z,k):
        '''Calculate k-th Hardy basis value at z
        '''
        return 1/(np.sqrt(np.pi)*(z + 1j))*((z - 1j)/(z + 1j))**k

    @classmethod
    def diff(cls,z,k):
        return 1/np.sqrt(np.pi)*1/(z + 1j)**(k + 2)*(z - 1j)**(k - 1)*()

class Nevanlinna:
    '''class to calculate analytic continuation of Nevanlinna function
    Refer to PRL 126, 056402
    Nevanlinna function: C_+ -> C^-_+ : x + yi (y > 0) ->  x' + y'i (y'>= 0)
    This class calculates analytic continuation of Nevanlinna function on C_+
    from values of it at Matsubara frequency points on C_+

    Attributes:
        nev_func:A Nevanlinna function (=a negative Green function)'s value at
                 sample points at Matsubara frequency on C_+
        points:Sample points of Matsubara frequency on C_+ 
        contra:A contractive(C_+ -> D) function made from nev_func
        energy_range: The range of path along which Functional integral is done
        delta:Functional integral is done along path parallel to and delta eV 
              away from real axis.
        lamb:Parameter which controls functionals
        phi:Values of n-th contractive function at the n-th Matsubara
            frequency sample points
        k_num:The number of Hardy basis to expand Theta_(M+1)

    '''
    def __init__(self,sample_points,nev_func,energy_range = [-10,10],delta = 0.001,lamb = 1e-4,k_num = 25,dtype=np.complex256):
        self.dtype = dtype
        self.nev_func = nev_func.astype(self.dtype)[::-1]
        self.points = sample_points.astype(self.dtype)[::-1]*1j
        self.p_num = self.points.shape[0]
        self.contra = self.mebius(self.nev_func)
        self.pick_check()
#        if self.pick_check() == False:
#            print("Pick matrix is not positive semi definite! ")
#        else:
#            print("Pick matrix is semi definite.")
        
        self.phi = np.zeros([self.p_num],dtype=self.dtype)
        self.phi_n_gamma()
        #print(format(self.phi[self.p_num-1],".34e"))
        #self.phi_n_gamma_2()
        self.energy_range = energy_range
        self.delta = delta
        self.lamb = lamb
        self.k_num = k_num       

    def pick_check(self):
        '''Check whether Pick matrix is positive semi-definite or not.
        '''
        pick = (1 - self.contra[:,None]*self.contra.conj()[None,:])\
               /(1 - self.mebius(self.points)[:,None]*self.mebius(self.points).conj()[None,:])
        print(pick[np.arange(self.p_num),np.arange(self.p_num)])
        #pick = pick #+ np.eye(self.p_num)*10**(-10)
        v = sp.linalg.eigvalsh(pick)
        #print(v)
        #sp.linalg.cholesky(pick)
        return np.all(v >= 0)

    def mebius(self,z):
        '''Mebius transformation of z
        '''
        return (z-1j)/(z+1j)
        
    def last_contra(self,z):
        '''Theta_(M+1) function's value at z
        '''
        return 0

    def phi_n_gamma(self):
        '''The method to calculate phi_n
        '''
        contras = np.zeros([self.p_num,self.p_num],dtype=self.dtype)
        contras[0,:] = self.contra
        for x in np.arange(self.p_num-1):
            self.phi[x] = contras[x,x]
            a = (self.points[x+1:] - self.points[x])/(self.points[x+1:] - self.points[x].conj())
            c = self.phi[x].conj()*a  
            contras[x+1,x+1:] = -(self.phi[x] - contras[x,x+1:])/(a-c*contras[x,x+1:])
        self.phi[self.p_num-1] = contras[self.p_num-1,self.p_num-1]

    def phi_n_gamma_2(self):
        phis =np.zeros([self.p_num],dtype=self.dtype)
        phis[0] = self.contra[0]
        for x in np.arange(self.p_num - 1):
            abcds = np.eye(2,dtype=self.dtype)
            for y in np.arange(x + 1):
                a = (self.points[x+1] - self.points[y])/(self.points[x+1] - self.points[y].conj())
                b = phis[y]
                c = a*b.conj()
                d = 1
                matrix = np.array([[a,b],[c,d]],dtype=self.dtype)
                abcds = abcds@matrix
            phis[x+1] = (-abcds[1,1]*self.contra[x+1]+abcds[0,1])/(abcds[1,0]*self.contra[x+1]-abcds[0,0])
        #print(format(phis[-1],".34e"))
        self.phi = phis    


    def abcd_cal(self,z):
        a = (z - self.points)/(z - self.points.conj())
        b = self.phi
        c = a*b.conj()
        matrix = np.array([[a,b],[c,np.ones(a.shape[0])]],dtype=self.dtype)
        matrix_product = matrix[:,:,0]
        for x in np.arange(self.p_num-1):
            matrix_product = matrix_product@matrix[:,:,x+1]
        return matrix_product.ravel()

    def interpolate(self,z):
        '''Interpolated values of a Nevanlinna function at z
        '''
        value = self.last_contra(z)
        abcd = self.abcd_cal(z)
        value = (abcd[0]*value + abcd[1])/(abcd[2]*value + abcd[3]) 
        value = 1j*(1 + value)/(1 - value)
        return value

    def interpolate_make(self,z,para,func,pre_abcd = None):
        '''value at z of an interpolated function which has parameters
        When Theta_(M+1) includes parameters, use this method.

            Args:
                z:points at which interpolated value is calculated
                para:Parameters contained by func
                func:Theta_(M+1) function
        '''
        func_value = func(z,para)
        if type(pre_abcd) == np.ndarray:
            abcd = pre_abcd
        else:
            abcd = self.abcd_cal(z)
        disk_value = (abcd[0]*func_value + abcd[1])\
                          *1/(abcd[2]*func_value + abcd[3]) 
        value = 1j*(1 + disk_value)/(1 - disk_value)
        return value

    def hardy_func(self, z, coeff):
        '''The value at z of a function spanned by Hardy basis

        Args:
            coeff: Coefficients correspond to Hardy basis and 
                   their conjugate
                   coeff[:coeff.shape[0]/4] correspond to Re(a_k)
                   coeff[coeff.shape[0]/4:2/4*coeff.shape[0]] 
                  correspond to Im(a_k)
                   coeff[2/4*coeff.shape[0]:3/4*coeff.shape[0]] 
                  correspond to Re(b_k)
                   coeff[3/4*coeff.shape[0]:coeff.shape[0]] 
                  correspond to Im(b_k)
        '''
        k_num = int(coeff.shape[0]/4)
        return np.sum((coeff[:k_num]+1j*coeff[k_num:2*k_num])*Hardy.basis(z, np.arange(0,k_num))\
               + (coeff[2*k_num:3*k_num] + 1j*coeff[3*k_num:4*k_num])*Hardy.basis(z, np.arange(0,k_num)).conj())

    def spectral_func(self,x,para,pre_abcd = None):
        '''Return spectral function value on x + 1j*self.delta
            Args:
                x:Real values on which spectral values are calculated
                para:Parameters included in a Nevanlinna function
        '''
        if type(pre_abcd) == np.ndarray:
            value = 1/np.pi*self.interpolate_make(x+1j*self.delta,para,self.hardy_func,pre_abcd = pre_abcd).imag 
            return value
        else:
            return 1/np.pi*self.interpolate_make(x+1j*self.delta,para,self.hardy_func).imag

    def write_parameters(self,prefix):
        np.savetxt("{}_M.dat".format(prefix),np.hstack([self.points.imag[:,None],self.nev_func[:,None]]))
        if hasattr(self,"hardy_coeff"):
            with open("{}_Hardy.dat","w") as f:
                f.write("#This file contains coefficients of Hardy basis\n")
                f.write("#corresponding to Matubara frequency in {}_M.dat\n".format(prefix))
                f.write("#delta is {} eV".format(str(self.delta)))
                f.write("\n#Correspondance between row and coefficients\n\n")
                f.write("#NUM=Total coefficients number = 4*n\n")
                f.write("#sum_(k=0)^(n)[a_k h_k + b_k h_k^*]\n\n")
                f.write("#0 - 1/4*NUM: Re(a_k)\n")
                f.write("#1/4*NUM - 2/4*NUM: Im(a_k)\n")
                f.write("#0 - 1/4*NUM: Re(b_k)\n")
                f.write("#1/4*NUM - 2/4*NUM: Im(b_k)\n\n")
            with open("{}_Hardy.dat","a") as f:
                np.savetxt(f,self.hardy_coeff)
        #with open("{}_Matsubara.dat".format(prefix),"w") as f:
            
class Integral:
    def __init__(self,points,integrand):
        self.points = points
        self.integrand = integrand

    def integral(self,*args,**keyargs):
        value = self.value_cal(*args,**keyargs)
        return scipy.integrate.simps(value,self.points)
    
    def value_cal(self,*args,**keyargs):
        return self.integrand(self.points,*args,**keyargs)

class Nevanlinna_optimize(Nevanlinna):
    def __init__(self,real_points,*args,**keyargs):
        super().__init__(*args,**keyargs)        
        self.real_points = real_points
        self.abcd = self.abcd_cal(self.real_points+1j*self.delta)
        self.spectral = Integral(self.real_points,self.spectral_func)
        self.spectral_sec = Integral(self.real_points,self.square_second)
        self.spectral_diff = Integral(self.real_points,self.spectral_func_differential)
        self.spectral_diff_sec = Integral(self.real_points,self.square_second_differential)

    def abcd_cal(self,z):
        a = (z[:,None] - self.points[None,:])/(z[:,None] - self.points.conj()[None,:])
        b = np.tile(self.phi[None,:],(z.shape[0],1))
        c = a*b.conj()
        matrix = np.array([[a,b],[c,np.tile(np.ones(self.p_num)[None,:],(z.shape[0],1))]],dtype=self.dtype)
        matrix_product = matrix[:,:,:,0] 
        for x in np.arange(self.p_num-1):
            matrix_product= np.einsum('ijk,jlk -> ilk',matrix_product,matrix[:,:,:,x+1])
        return matrix_product.reshape([matrix_product.shape[0]*matrix_product.shape[1],matrix_product.shape[2]])

    def jac(self,z,para,func,pre_abcd):
        '''This returns gradient of func as to para
        '''
        az = pre_abcd[0]
        bz = pre_abcd[1]
        cz = pre_abcd[2]
        dz = pre_abcd[3]    
        func_value = func(z,para)
        disk_value = (az*func_value + bz)\
                    *1/(cz*func_value + dz) 
        func_derivative = np.zeros([4*self.k_num,z.shape[0]],dtype=np.complex256)
        value = Hardy.basis(z[None,:],np.arange(self.k_num)[:,None])
        func_derivative[:self.k_num] = value
        func_derivative[self.k_num:2*self.k_num] = 1j*value
        func_derivative[2*self.k_num:3*self.k_num] = value.conj()
        func_derivative[3*self.k_num:] = 1j*value.conj()
        value =  2j/(1-disk_value)**2*(az*dz - bz*cz)\
               *1/(cz*disk_value+dz)**2
        value = value[None,:]*func_derivative
        #print("value is")
        #print(value.shape)
        #print("value end")
        return value

    def hardy_func(self, z, coeff):
        '''The value at z of a function spanned by Hardy basis

        Args:
            coeff: Coefficients correspond to Hardy basis and 
                   their conjugate
                   coeff[:coeff.shape[0]/4] correspond to Re(a_k)
                   coeff[coeff.shape[0]/4:2/4*coeff.shape[0]] 
                  correspond to Im(a_k)
                   coeff[2/4*coeff.shape[0]:3/4*coeff.shape[0]] 
                  correspond to Re(b_k)
                   coeff[3/4*coeff.shape[0]:coeff.shape[0]] 
                  correspond to Im(b_k)
        '''
        k_num = int(coeff.shape[0]/4)
        value =  np.sum((coeff[:k_num]+1j*coeff[k_num:2*k_num])[None,:]*Hardy.basis(z[:,None], np.arange(0,k_num)[None,:])\
               + (coeff[2*k_num:3*k_num] + 1j*coeff[3*k_num:4*k_num])[None,:]*Hardy.basis(z[:,None], np.arange(0,k_num)[None,:]).conj(),axis=1)
        return value

    def spectral_func_differential(self,x,para,pre_abcd):
        #print("spectral func differential")
        return 1/np.pi*self.jac(x+1j*self.delta,para,self.hardy_func,pre_abcd).imag

    def second_derivative(self,x,para):
        '''Second derivative at x
        '''
        value = 1/np.pi*derivative(self.spectral_func,x,para,n=2,**{"pre_abcd":self.abcd})
        return value

    def square_second(self,x,para):
        '''Squared Second derivative at x
        '''
        return self.second_derivative(x,para)**2

    def second_derivative_differential(self,x,para,pre_abcd):
        '''Second derivative differential at x
        '''
        return 1/np.pi*derivative(self.spectral_func_differential,x,para,pre_abcd,n=2)

    def square_second_differential(self,x,para,pre_abcd):
        '''Squared Second derivative at x
        '''
        #print("square second differential")
        #print(para.shape)
        return self.second_derivative_differential(x,para,pre_abcd)**2

    def spectral_integral(self,para):
        '''A functional which is used to make spectral function proper
        '''
        print(para)
        return np.abs(1 - sp.integrate.quad(self.spectral_func,self.energy_range[0],\
               self.energy_range[1],args=para)[0])**2\
              +self.lamb*sp.integrate.quad(self.square_second,self.energy_range[0],\
               self.energy_range[1],args=para)[0]

    def spectral_integral_sample(self,para):
        print("x")
        print(para)
        print("spectral part")
        value = np.abs(1 - self.spectral.integral(para,**{"pre_abcd":self.abcd}))**2
        print(value)
        value2 =self.lamb*self.spectral_sec.integral(para)
        print("lambda part")
        print("value2")
        print(value2)
        value += value2
        print("sum")
        print(value)
        return value

    def spectral_integral_differential(self,para):
        '''gradient as to para of spectral integral 
        '''
        value =(-1)*sp.integrate.quad_vec(self.spectral_func_differential,self.energy_range[0],\
               self.energy_range[1],args=(para,))[0]
        #print("value shape is")
        #print(value.shape)
        value*= 2*(1 - sp.integrate.quad(self.spectral_func,self.energy_range[0],\
               self.energy_range[1],args=para)[0])
        value+=self.lamb*sp.integrate.quad_vec(self.square_second_differential,self.energy_range[0],\
               self.energy_range[1],args=(para,))[0]
        return value

    def spectral_integral_differential_sample(self,para):
        '''gradient as to para of spectral integral 
        '''
        value = 2*(1 - self.spectral.integral(para,**{"pre_abcd":self.abcd}))\
                *(-1)*self.spectral_diff.integral(para,self.abcd)\
                + self.lamb*self.spectral_diff_sec.integral(para,self.abcd)
        return value

    def optimize(self, first_guess="random",coefficient = None):
        '''A method to make "good" Nevanlinna functions
        This function returns parameter values of Nevanlinna function
        at which spectral functional has an extremum.
        This method uses the CG method to calculate extremum.

            Args:
                first_guess:A flag to set initial guess of coefficients
                            of Hardy basis
                            If it is not "random", coefficients have to
                            be passed into this method.
                coefficient:Initial guess of coefficients of Hardy basis  
        '''
        if type(coefficient) == np.ndarray:
            coefficient = coefficient
        elif first_guess == "random":
            #coefficient = np.random.rand(self.k_num*4)*np.sqrt(np.pi)/(2*np.sqrt(2)*self.k_num)
            coefficient = np.zeros([self.k_num*4])
        else:
            print("error! Pass coefficient or set first_guess = \"random\"")
            exit  
        self.result = sp.optimize.minimize(self.spectral_integral_sample,\
                           coefficient,\
                           jac=self.spectral_integral_differential_sample)
        '''
        self.result = sp.optimize.basinhopping(self.spectral_integral_sample,\
                           coefficient,\
                           minimizer_kwargs={"jac":self.spectral_integral_differential_sample})
        '''
        self.hardy_coeff = self.result['x']

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import temperature_spir
    import numpy as np
    from Nevanlinna import Nevanlinna
    wmax = 10 # eV
    T = 116 # K
    k_num = 25 # The number of Hardy basis
    temp = temperature_spir.Temperature(T,wmax)
    sample = temp.omegaF_sp[int(1/2*temp.omegaF_sp.shape[0]):].astype(np.float128)
    #sample = np.pi*(2*np.arange(36)+1)/100
    # define green function
    #def green(x):
    #    return 1/(x + 0.001j) #+ 1/(x - 2 + 0.5j) #+ 1/(x+3 + 1j)
    # end
    sigma = 0.5
    mu = 0
    def green(x):
        return 1/np.sqrt(2*np.pi)/sigma*np.exp(1/2*x**2/sigma**2) 
    nev_func = -green(sample*1j)
    print(nev_func)
    p_num = 6000
    points = 2/(p_num)*np.arange(wmax*p_num) - wmax
  #  np.savetxt("ifile6.txt",np.hstack([sample[:,None],-nev_func.real[:,None],-nev_func.imag[:,None]]))
    '''
    fig,ax =plt.subplots()
    value64 = []
    test64 = Nevanlinna(sample,nev_func,energy_range=[-wmax,wmax],delta=0.001,k_num = k_num,dtype=np.complex128)
    for x in points:
        value64.append(1/np.pi*test64.interpolate(x+test64.delta*1j).imag)
    ax.plot(points,-green(points+1j*test64.delta).imag/np.pi,label="correct")
    ax.plot(points,value64,label="float64",lw=1.0)
    ax.set_xlim([-6,6])
    ax.set_xlabel(r"$\omega$ [eV]")
    ax.set_ylabel(r"A($\omega$)")
    ax.legend()
  #  plt.savefig("float64_compare3.png")
    plt.show()
    '''
    fig,ax =plt.subplots()
    value128 = []
    test128 = Nevanlinna(sample,nev_func,energy_range=[-wmax,wmax],delta=0.001,k_num = k_num)
    for x in points:
        value128.append(1/np.pi*test128.interpolate(x+test128.delta*1j).imag)
    ax.plot(points,-green(points+1j*test128.delta).imag/np.pi,label="correct")
    ax.plot(points,value128,label="float128",lw=1.0)
    ax.set_xlim([-6,6])
    ax.set_xlabel(r"$\omega$ [eV]")
    ax.set_ylabel(r"A($\omega$)")
    ax.legend()
    #plt.savefig("float128_compare3.png")
    plt.show()

    #cpp_spectrum = np.loadtxt("spectrum_2.dat")
    #fig,ax =plt.subplots()
    #ax.plot(points,-green(points+1j*test128.delta).imag/np.pi,label="correct")
    #ax.plot(cpp_spectrum[:,0],cpp_spectrum[:,1],label="cpp_nevanlinna",lw=1.0)
    #ax.set_xlim([-6,6])
    #ax.set_xlabel(r"$\omega$ [eV]")
    #ax.set_ylabel(r"A($\omega$)")
    #ax.legend()
    #plt.savefig("gmp_used_2.png")
    #plt.show()

    #
#    print(nev_func)
#    a = np.array([1],dtype=np.float128)
#    for x in range(17):
#        a*= 10**(-100)
#    print(a)
#    for x in [4,5,6]:
#        print(format((nev_func[-1*x]*a)[0],".34e"))
    #np.savetxt("ifile.txt",np.hstack([sample[:,None],np.zeros(sample.shape[0])[:,None],nev_func[:,None]]),fmt="%.34e")
    #print("green function")
    #print(nev_func)
    # Theta_(M+1) = 0
    # Theta_(M+1) = 0
    #print("random coefficients")
    #value = []
#    coefficient = np.random.rand(k_num*4).astype(np.float128)*np.sqrt(np.pi)/(2*k_num)
    
#    test = Nevanlinna_optimize(points,sample,nev_func,energy_range=[-wmax,wmax],delta=0.001,k_num = k_num)
    #integ = Integral(points,test.spectral_func)
    #print("integtral value is")
    #print(integ.integral(coefficient,**{"pre_abcd":test.abcd}))
#    test.optimize(coefficient=coefficient)
    #print("random coefficients")
    #print(coefficient)
    #print("optimized coefficients")
    #print(test.hardy_coeff)
#    fig,ax = plt.subplots() 
#    value = []
    #ax.plot(points,-green(points+test.delta*1j).imag/np.pi,label="correct")
    #ax.plot(points,value,lw=1.0,label="random")
    #ax.plot(points,value0,lw=1.0,label="zeros")
    #value = test.spectral_func(points,test.hardy_coeff)
#    ax.plot(points,value,lw=1.0,label="optimized")
#    ax.legend()
#    plt.show()
#    print(test.result)
    
    '''
    wmax = 4 # eV
    T = 70 # K
    k_num = 10 # The number of Hardy basis
    color = ["r","g","black"]
    temp = temperature_spir.Temperature(T,wmax)
    sample = temp.omegaF_sp[int(1/2*temp.omegaF_sp.shape[0]):]
    print(sample)
    print(sample > 0)
    '''
    '''
    def green(x):
        return 1/(x+1 + 0.1j) #+ 1/2*1/(x-2 + 0.2j) 
    nev_func = -green(sample*1j)   
    test = Nevanlinna(sample,nev_func,energy_range=[-4,4],delta=0.01)
    '''
    '''
    value = []
    for x in np.arange(len(sample)):
        value.append(test.interpolate(sample[x]*1j))
    print(nev_func)
    print(value)
    '''