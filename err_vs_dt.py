import matplotlib.pyplot as plt
import numpy as np

metodo=['file\velocity_verlet', 'file\rk4']
dts=[0.001, 0.002, 0.004, 0.008, 0.016]
omega=1

for m in metodo:
    err=[]
    for dt in dts:
        data=np.loadtxt(f"{m}_{dt}.dat")
        x0=data[0,1]
        v0=data[0,2]
        err.append(data[-1, 4]/(np.sqrt(x0*x0+(v0/omega)**2)))
    plt.figure()
    plt.loglog(dts, err, '-o', label=metodo)
    plt.xlabel('dt[s]')
    plt.ylabel('errore/A')
    plt.show()
    
