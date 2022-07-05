if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import temperature_spir
    import numpy as np
    from Nevanlinna import Nevanlinna
    wmax = 2 # eV
    T = 70 # K
    k_num = 10 # The number of Hardy basis
    color = ["r","g","black"]
    temp = temperature_spir.Temperature(T,wmax)
    sample = temp.omegaF_sp[int(1/2*temp.omegaF_sp.shape[0]):]
    def green(x):
        return 1/(x + 0.1j) 
    nev_func = -green(sample*1j)   
    test = Nevanlinna(sample,nev_func,energy_range=[-4,4],delta=0.001)
    coefficients = [np.zeros([4*k_num]), np.random.rand(4*k_num)*np.sqrt(np.pi)/(2*k_num), np.array([0.061524  ,0.04909083,0.06532009,0.04096199,0.0049477,0.07045386
,0.08328574,0.0286432,0.07831385,0.04628823,0.03556203,0.1110606
,0.12542929,0.09560942,0.04028036,0.04673346,0.00090295,0.0403369
,0.03134908,0.00645551,0.05447442,0.04051745,0.04143725,0.0563686
,0.09000172,0.02996258,0.0362983,0.06416604,0.05804415,0.01939144
,0.0344877, 0.02517031,0.14658374,0.02014969,0.04733996,0.02145694
,0.13702928,0.0366991, 0.04329489,0.0698173 ])]

    # Theta_(M+1) = 0
    points = 0.001*np.arange(2*wmax*1000) - wmax
    label = [r"$\Theta$_(R+1) = 0","Random coefficient Hardy","Optimized Hardy"] 
    fig,ax = plt.subplots() 
    value = []


    for x in points:
        value.append(test.spectral_func(x,np.zeros([4])))
    ax.plot(points,value,label=label[0],lw=1.0,c=color[0])
    ax.plot(points,-green(points+test.delta*1j).imag/np.pi,label="correct")
    ax.legend()
    plt.show()

    for y in range(2):
        fig,ax = plt.subplots() 
        value = []
        for x in points:
            value.append(test.spectral_func(x,coefficients[y+1]))
        ax.plot(points,value,label=label[y+1],lw=1.0,c=color[y+1])
        ax.plot(points,-green(points+test.delta*1j).imag/np.pi,label="correct")

    #ax.plot(points,value,label="Nevanlinna")
        ax.legend()
        plt.show()