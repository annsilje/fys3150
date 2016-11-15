import matplotlib.pyplot as plt
import numpy as np
import sys

data40 = np.genfromtxt('../output/40x40.txt', unpack=True)
data60 = np.genfromtxt('../output/60x60.txt', unpack=True)
data100 = np.genfromtxt('../output/100x100.txt', unpack=True)
data140 = np.genfromtxt('../output/140x140.txt', unpack=True)

ylabels=[r'Mean energy per spin $\langle E/N\rangle$', 
         r'Mean magnetization per spin $\langle M/N\rangle$',
         r'Mean magnetization per spin $\langle |M|/N\rangle$', 
         r'Specific heat $C_V$',
         r'Susceptibility $\chi$']
         
filenames=['mean_energy.pdf', 'mean_magnetization.pdf', 'mean_abs_magnetization.pdf', 'heat_capacity.pdf', 'susceptibility.pdf']
T_C = [data40[0][np.argmax(data40[5])], data60[0][np.argmax(data60[5])], data100[0][np.argmax(data100[5])], data140[0][np.argmax(data140[5])]]
L = [40, 60, 100, 140]
a = []
print('T_C(40x40)=%f' % (T_C[0]))
print('T_C(60x60)=%f' % (T_C[1]))
print('T_C(100x100)=%f' % (T_C[2]))
print('T_C(140x140)=%f' % (T_C[3]))

for i in range(len(T_C)):
    for j in range(i+1, len(T_C)):
        a.append((T_C[i] - T_C[j])/(1/L[i] - 1/L[j]))  
        
a = np.mean(np.array(a))
print("Average a: %f" % a)

T_C_Inf = []

for i in range(len(T_C)):
    T_C_Inf.append(T_C[i] - a/L[i])
    
T_C_Inf = np.mean(np.array(T_C_Inf))
print("Average T_C(L=Inf): %f" % T_C_Inf)

for i, label in enumerate(ylabels):
    plt.figure()
    plt.plot(data40[0], data40[i+1], linestyle='--', marker='o')
    plt.plot(data60[0], data60[i+1], linestyle='--', marker='s')
    plt.plot(data100[0], data100[i+1], linestyle='--', marker='8')
    plt.plot(data140[0], data140[i+1], linestyle='--', marker='v')
    plt.ylabel(label)
    plt.xlabel('Temperature')
    plt.legend(['40x40', '60x60', '100x100', '140x140'], bbox_to_anchor=(0.25, 1), loc=1, borderaxespad=0.)
    plt.tight_layout()
    plt.savefig('../fig/%s' % (filenames[i]))
    plt.show()
    
    
