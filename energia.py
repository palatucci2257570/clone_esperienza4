import matplotlib.pyplot as plt
import numpy as np

metodo=['velocity_verlet', 'rk4']
omega=1
plt.figure()
for m in metodo:
    data=np.loadtxt(f"file/{m}1_{0.2}.dat")
    t=data[:,0]
    E=data[:,3]
    plt.plot(t ,(E-E[0])/E[0], '-', label=m)
plt.xlim(0,200)
plt.ylim(-0.007)
plt.xlabel('t[s]')
plt.ylabel('(E-E0)/E0')
plt.axline((0,0), slope=0, color='black', label='retta E0')
plt.title('energia in funzione del tempo')
plt.legend()
plt.savefig(f"energia_vs_t.png", dpi=600)
plt.show()
    
