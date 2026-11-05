import matplotlib.pyplot as plt
import numpy as np

metodo=['velocity_verlet', 'rk4']
dts=[0.001, 0.002, 0.004, 0.008, 0.016]
omega=1
plt.figure()
for m  in metodo:
    err=[]
    for dt in dts:
        data=np.loadtxt(f"file/{m}_{dt}.dat")
        x0=data[0,1]
        v0=data[0,2]
        err.append(data[-1, 4]/(np.sqrt(x0*x0+(v0/omega)**2)))
    plt.plot(np.log(dts),np.log(err), '-o', label=m)

plt.xlabel('dt[s]')
plt.ylabel('errore/A')
plt.title('errore in funzione di Dt')
plt.axline((0,0), slope=2, label='retta di pendenza 2')
plt.axline((0,0), slope=3, color='red',label='retta di pendenza 3')
plt.xlim(np.log(dts[0]), np.log(dts[-1]))
plt.legend()
plt.savefig(f"grafici/errore_confronto.png", dpi=600)
plt.show()
    
