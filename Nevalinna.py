import numpy as np

class Nevalinna:
    '''class to calculate analytic continuation
    refer to PBL 126, 056402
    Nevalinna function: x + yi (y > 0) ->  x' + y'i (y'>= 0)
    
    Args:
        nev_func:Nevalinna function (- retarded Green's function)'s value at sample_points
        points:sample points of Matsubara frequency on upper half plane 
        contra: contractive(C+ -> D) function made from nev_func
    '''
    def __init__(self,sample_points,nev_func):
        self.nev_func = nev_func
        self.points = sample_points
        self.p_num = self.points.shape[0]
        self.contra = self.mebius(self.nev_func)
        self.phi = np.zeros([self.p_num],dtype=np.complex128)
        self.phi_n_gamma()

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

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    sample = np.pi*1j*np.arange(1,20)/100
    def green(x):
        return 1/(x + 0.5j) + 1/(x -1 + 0.2j)
    nev_func = -green(sample)
    test = Nevalinna(sample,nev_func)
    fig,ax = plt.subplots()
    points = 0.001*np.arange(4000) -2
    value = []
    for x in points:
        value.append(test.interpolate(x).imag)
    ax.plot(points,-green(points).imag,label="correct")
    ax.plot(points,value,label="Nevalinna")
    ax.legend()
    plt.show()