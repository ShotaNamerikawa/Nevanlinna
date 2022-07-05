import sys
sys.path.append("/home/shota/wannier_utils/src/wannier_utils")
import hamiltonian
import numpy as np
import numpy.linalg

hr = hamiltonian.HamR(file_hr = "/home/shota/Fe/Fe_hr.dat")
hk = hamiltonian.HamK(hr,np.array([0.0,0.0,0.0]))
print(np.linalg.inv(hk.mat)[1,1])