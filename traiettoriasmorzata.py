import matplotlib.pyplot as plt
import numpy as np

gamma=[0, 0.16, 0.63, 3.15, 6.26, 15.65]
plt.figure()
for dt in gamma:
    data=np.loadtxt(f"file/gamma{dt}.dat")
    t = data[:,0]
    theta = data[:,1] 
    plt.plot(t, theta, '-', label=f'gamma={dt} [1/s]')
plt.xlabel('t[s]')
plt.ylabel('theta [rad]')
plt.title('traiettoria gamma')
plt.legend()
plt.savefig(f"grafici/traiettoriagamma.png", dpi=600)
plt.show()
    
plt.figure()
for dt in gamma:
    data=np.loadtxt(f"file/gamma{dt}.dat")
    t = data[:,0]
    E = data[:,3] 
    plt.plot(t, E, '-', label=f'gamma={dt} [1/s]')
plt.xlabel('t[s]')
plt.ylabel('energia [joule]')
plt.title('energia gamma')
plt.legend()
plt.savefig(f"grafici/energiagamma.png", dpi=600)
plt.show()

plt.figure()
for dt in gamma:
    data=np.loadtxt(f"file/gamma{dt}.dat")
    v = data[:,2]
    theta = data[:,1] 
    plt.plot(theta,v, '-', label=f'gamma={dt} [1/s]')
plt.xlabel('theta[rad]')
plt.ylabel('omega[rad/s]')
plt.title('fasi gamma')
plt.legend()
plt.savefig(f"grafici/fasigamma.png", dpi=600)
plt.show()
