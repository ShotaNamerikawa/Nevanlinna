import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

fig,ax = plt.subplots()
for i in np.arange(10):
    spectrum = np.loadtxt("spectrum_up_Fe_{}_Co_{}.dat".format(str(100 - 10*i),str(10*i)))
    ax.plot(spectrum[:,0],spectrum[:,1],c=cm.jet_r(i/10))
    spectrum = np.loadtxt("spectrum_down_Fe_{}_Co_{}.dat".format(str(100 - 10*i),str(10*i)))
    ax.plot(spectrum[:,0],-spectrum[:,1],c=cm.jet_r(i/10))
ax.set_xlabel("energy [eV]")
ax.set_ylabel("DOS [1/eV/unit cell]")
ax.set_xlim([-2,2])

#fig.colorbar(aspect=40, pad=0.08, orientation="vertical")
#plt.colorbar()
plt.show()
