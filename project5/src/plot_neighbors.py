import numpy as np
import matplotlib.pyplot as plt
import glob
import scipy.optimize

def plot(agents, savings):
    plt.figure()
    symbol = ['o','s','<','d','v']
    filenames1 = sorted(glob.glob('../output/wealth_pdf_savings_%s_alpha_*_gamma_0.00_agents_%s.txt' % (savings, agents)))
    filenames2 = sorted(glob.glob('../output/cycle_data_savings_%s_alpha_*_gamma_0.00_agents_%s.txt' % (savings, agents)))
    for i, (filename1, filename2) in enumerate(zip(filenames1,filenames2)):

        m, counts = np.genfromtxt(filename1, unpack=True)
        _,_, _, accepted, max_interactions = np.genfromtxt(filename2, unpack=True)
        
        alpha = filename1[40:44]
        counts = np.trim_zeros(counts, 'b')/int(agents)
        m = m[:len(counts)] + 0.025
        
#        plt.loglog(m, counts, symbol[i] ,label=(r'$\alpha=%s$' % filename1[40:44]), markeredgewidth=0.0)
        plt.loglog(m, counts, label=(r'$\alpha=%s$' % filename1[40:44]))
        count = 0
        n = len(counts) - 1 
        while count < 0.20 :
            count += counts[n]
            n -= 1
        counts = counts
        
        f2 = lambda x, a: (x**(-1 -a))
        p, s = scipy.optimize.curve_fit(f2, m[-n:], counts[-n:])
        
        print(r'%s & %s & %s & %4.2f & %6.4f & %2.0f\%%\\' % (agents, savings, alpha, p[0], s[0][0], np.mean(accepted)/1e7*100))
        
    plt.xlabel(r'Money')
    plt.ylabel(r'Probability distribution')   
    plt.xlim(1e-2,1e3)
    plt.ylim(1e-7, 1e0) 
    plt.legend(numpoints=1)
    plt.savefig('../fig/neighbor_pdf_savings_%s_agents_%s.pdf' % (savings, agents))


agents = ['500', '1000']
savings = ['0.00', '0.50']


print(r"\begin{tabularx}{\textwidth}{X X X r r r}")
print(r"\hline")
print(r"$N$ & $\lambda$ & $\alpha$ & $\nu$ & $\sigma_{\nu}^2$ & Accepted\\")
print(r"\hline\hline")

for agent in agents:
    for save in savings:
        plot(agent, save)
print(r"\hline");
print(r"\end{tabularx}")
plt.show()
