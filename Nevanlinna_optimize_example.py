if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import temperature_spir
    import numpy as np
    from Nevanlinna import Nevanlinna
    wmax = 4 # eV
    T = 116 # K
    k_num = 25 # The number of Hardy basis
    temp = temperature_spir.Temperature(T,wmax)
    sample = temp.omegaF_sp[int(1/2*temp.omegaF_sp.shape[0]):]
    #sample = np.pi*(2*np.arange(36)+1)/100
    def green(x):
        return 1/(x + 1j)# +1/(x-0.1 + 1j) + 1/(x-0.4 + 1j) + 1/( x- 2 + 0.3j) + 1/(x+2 + 0.2j)# + 1/(x-3 +0.01j)
    sigma = 0.5
    mu = 0
    '''
    def green(x):
        return -1j/np.sqrt(2*np.pi)/sigma*np.exp(-1/2*(x - mu)**2/sigma**2) 
    '''
    nev_func = -green(sample*1j)
    print(nev_func)
    test = Nevanlinna(sample,nev_func,energy_range=[-wmax,wmax],delta=0.01,k_num = k_num)
    # Theta_(M+1) = 0
    points = 0.001*np.arange(wmax*2000) - wmax

    print("zero function")
    fig,ax =plt.subplots()
    value = []
    for x in points:
        value.append(test.spectral_func(x,np.zeros([4*k_num])))
    ax.plot(points,-green(points).imag/np.pi,label="correct")
    ax.plot(points,value,lw=1.0)
    ax.legend()
    plt.show()

    print("random coefficients")
    fig,ax =plt.subplots()
    value = []
    for x in points:
        value.append(test.spectral_func(x,np.random.rand(k_num*4)*np.sqrt(np.pi)/(2*k_num)))
    ax.plot(points,-green(points).imag/np.pi,label="correct")
    ax.plot(points,value,lw=1.0)
    ax.legend()
    plt.show()


 #   test.optimize()
    print("optimized coefficients")
    fig,ax = plt.subplots() 
    value = []
    for x in points:
        value.append(test.spectral_func(x,np.array([ 0.00296191,  0.006132   , 0.03278565 , 0.01775154 ,0.02600615,0.00350407
,0.00725673,0.01153532,0.0162966,0.02231397,0.00496188,0.00684689
,0.00162981,-0.00180955,-0.0149339,0.00620279,-0.00486561,-0.03064669
,-0.03692244,-0.04268056,-0.05808026,-0.07916564,-0.06860665,-0.08858662
,-0.11562656,0.01086565,0.01406621,0.00415616,0.00548112,0.03218414
,0.01050323,0.02906836,0.03238034,0.02582782,0.02789163,0.00120861
,-0.00809432,-0.00264358,0.01791576,0.00667019,0.0017338,0.00129529
,-0.03880184,-0.03904825,-0.03331913,-0.04335583,-0.05247995,-0.07135753
,-0.08700835,-0.09715323,0.00965046,0.02117042,0.00621813,0.00344794
,0.00433004,0.0030193,0.02529826,0.01621484,0.02280941,0.01149651
,0.01166262,0.02176173,0.00800317,-0.00765209,-0.01124727,-0.00301909
,-0.03307326,-0.02460192,-0.0456735,-0.03488437,-0.01711854,-0.07140144
,-0.10522444,-0.04907437,-0.08538139,0.0137159,0.01748945,0.02875527
,0.01401661,0.00686403,0.00151671,0.027729,0.01653994,0.0332112
,0.03069185,-0.00725201,-0.00457973,-0.00176671,-0.00237834,0.01682081
,0.01057036,-0.01867512,-0.0640402,0.00286598,0.00017601,-0.08319444
,-0.05618301,-0.03788125,-0.15767226,-0.13076316])))
    ax.plot(points,-green(points+test.delta*1j).imag/np.pi,label="correct")
    ax.plot(points,value,lw=1.0)
    plt.show()
