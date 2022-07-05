import numpy as np
import scipy as sp
import scipy.integrate

class Multi_quad:
    '''class to integral vector function
    function(x,para)   
    '''
    def __init__(self,func):
        self.func = func
        #self.quad = np.vectorize(self.single_quad,excluded=[2])

    def single_quad(self,down,up,args):
        print(self.func)
        sp.integrate.quad(self.func,down,up,args=args)

if __name__ == "__main__":
    a = Multi_quad(lambda x,a:a[0]*x + a[1])
    a.single_quad(0,1,np.array([1,1]))
    #print(a.quad(np.array([0]),np.array([1]),np.array([1,1])))