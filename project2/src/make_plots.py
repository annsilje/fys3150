import matplotlib.pyplot as plt
import numpy as np
import glob, os
from collections import OrderedDict


def latex_float(float_number, digits=2):
    float_str = '%.*e' % (digits, float(float_number))
    base, exp = float_str.split("e")
    if int(exp)==0:
        return '$%s$' % (base)
    elif base.startswith('1.') and float(base).is_integer():
        return '$10^{%d}$' % int(exp)
    return '$%s \\times 10^{%d}$' % (base,int(exp))
    


def read_data_files():
    max_rho_max = 0
    data = dict()

    files = sorted(glob.glob('../output/*particle*.txt'))
    for filename in files:
        rho_max, n = np.genfromtxt(filename, max_rows=1)
        u1, u2, u3 = np.genfromtxt(filename, unpack=True, skip_header=1)
        if rho_max > max_rho_max: max_rho_max = rho_max 
    
        title = os.path.splitext(os.path.basename(filename))[0].replace('_', ' ')
        pot_str = title[:-15] if 'One' not in title else title
        freq_str = title[-4:] if 'One' not in title else ' '
        pot = data.setdefault(pot_str,OrderedDict()) 
        freq = pot.setdefault(freq_str, {})
        
        freq['u1'] = u1
        freq['u2'] = u2
        freq['u3'] = u3
        freq['rho_max'] = rho_max
        freq['n'] = n

    return data, max_rho_max
    
def make_plots(data):    
    for pot_str, potential in data.items():
        fig = plt.figure() if 'One' in pot_str else plt.figure(figsize=(10,12))
            
        for i, (freq_str, freq) in enumerate(potential.items()):
            h = freq['rho_max']/(freq['n']+1)
            rho = np.linspace(h, freq['rho_max']-h, freq['n'])
            
            #Plot Eigenvectors
            if freq_str != ' ': 
                plt.subplot(3,2, i+1)
                plt.title(r'$\omega_r=%s$' % (freq_str))
            plt.plot(rho, freq['u1']**2)
            plt.plot(rho, freq['u2']**2)
            plt.plot(rho, freq['u3']**2)
            plt.xlim([0, freq['rho_max']])
            
            plt.legend([r'$u_1$', r'$u_2$', r'$u_3$'])
            plt.xlabel(r'$\rho$')
            plt.ylabel(r'$u(r)^2$')
    
        plt.suptitle(pot_str, fontsize=18)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig('../fig/%s.pdf' % (pot_str.replace(' ', '_')))
    

def pot2txt(potential):
    if potential=='0': return 'One electron'
    if potential=='1': return 'No interaction'
    if potential=='2': return 'Interaction'


def make_runtime_table():
    iterations, freq, pot, jacobi, *arma = np.genfromtxt('../output/exec_time.txt', unpack=True)
    print("\\begin{tabularx}{\\textwidth}{l r r r r r r }")
    print("\\hline")
    print("Case & $\omega_r$ & Step Length & $\\rho_{max}$ & Iterations & Jacobi & eig\_sym\\\\")
    print("\\hline\\hline")
    
    with open('../output/exec_time.txt') as runtimes:
        for line in runtimes:
            iterations, freq, pot, h, rho_max, jacobi, *arma = line.split()
            if arma:
                print('%s & %4.2f & %s & %2.0f & %d & %7.2fs & %7.4fs  \\\\' % (
                        pot2txt(pot), float(freq), latex_float(h), float(rho_max), 
                        int(iterations), float(jacobi), float(arma[0])))
            else:
                print('%s & %4.2f & %s & %2.0f & %d & %7.2fs &  \\\\' % (
                        pot2txt(pot), float(freq), latex_float(h), float(rho_max), 
                        int(iterations), float(jacobi)))

    print("\\hline");
    print("\\end{tabularx}")
    

def plot_repulsion(rho_max):
    plt.figure()
    rho = np.linspace(0, rho_max, 200)
    plt.plot(rho[1:], 1/rho[1:])
    plt.title('Repulsive Coloumb interaction')
    plt.xlabel(r'$\rho$')
    plt.ylabel(r'$\frac{1}{\rho}$')
    plt.ylim([0,1])
    plt.tight_layout()
    plt.savefig('../fig/Repulsion.pdf')


def main():
    data, max_rho_max = read_data_files()
    make_plots(data)
    make_runtime_table()
    plot_repulsion(max_rho_max)
    
    plt.show()
    
    
if __name__ == '__main__':
    main()


