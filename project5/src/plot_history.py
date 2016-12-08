import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize
import glob

def plot(agents, alpha, savings):
    plt.figure()

    filenames1 = sorted(glob.glob('../output/wealth_pdf_savings_%s_alpha_%s_gamma_*_agents_%s.txt' % (savings, alpha, agents)))
    filenames2 = sorted(glob.glob('../output/cycle_data_savings_%s_alpha_%s_gamma_*_agents_%s.txt' % (savings, alpha, agents)))
    symbol = ['o','s','<','d','v']
    for i, (filename1,filename2) in enumerate(zip(filenames1, filenames2)):
        gamma = filename1[51:55]
        m, counts = np.genfromtxt(filename1, unpack=True)
        _,_, _, accepted, max_interactions = np.genfromtxt(filename2, unpack=True, skip_header=1)
        
        counts = np.trim_zeros(counts, 'b')/int(agents)
        m = m[:len(counts)] + 0.025
        
        plt.loglog(m, counts, label=(r'$\gamma=%s$' % filename1[51:55]))
        
        count = 0
        n = len(counts) - 1 
        while count < 0.20 :
            count += counts[n]
            n -= 1
        counts = counts
        
        f2 = lambda x, a: (x**(-1 -a))
        p, s = scipy.optimize.curve_fit(f2, m[-n:], counts[-n:])
        
        print(r'%s & %s & %s & %s & %d & %4.2f & %6.4f & %2.0f\%%\\' % (agents, savings, alpha, gamma, np.max(max_interactions), p[0], s[0][0], np.mean(accepted)/1e7*100))
    plt.xlabel(r'Money')
    plt.ylabel(r'Probability distribution')    
    plt.xlim(10**-1.8,10**2.2)
    plt.ylim(10**-6.5, 10**-0.8)
    plt.legend(numpoints=1)
    plt.savefig('../fig/history_pdf_savings_%s_alpha_%s_agents_%s.pdf' % (savings, alpha, agents))


       


agents = ['1000']
savings = ['0.00', '0.50']
alpha = ['1.00', '2.00']

print(r"\begin{tabularx}{\textwidth}{X X X X X X X r}")
print(r"\hline")
print(r"$N$ & $\lambda$ & $\alpha$ & $\gamma$ & $\max(c_{ab})$ & $\nu$ & $\sigma_{\nu}^2$ & Accepted \\")
print(r"\hline\hline")

for agent in agents:
    for save in savings:
        for a in alpha:
            plot(agent, a, save)
            
            
print(r"\hline");
print(r"\end{tabularx}")

plt.show()
