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
    test = Nevanlinna(sample,nev_func,energy_range=[-5,5],delta=0.01)

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
    #ax.plot(points,value,label="Nevanlinna")
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
    ax.plot(points,value,label="Nevanlinna")
    ax.legend()
    plt.show()
    '''