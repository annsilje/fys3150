import numpy as np
import matplotlib.pyplot as plt
from matplotlib import markers
import glob
import scipy.special
import scipy.optimize




def alphas(filenames):
    print(r"\begin{tabularx}{\textwidth}{X X X X X r}")
    print(r"\hline")
    print(r"$\lambda$ & $\alpha$ & $\gamma$ & $N$ & $\nu$ & $\sigma_{\nu}^2$\\")
    print(r"\hline\hline")
    
    for i, filename in enumerate(filenames):
        agents = int(filename[63:-4])
        alpha = float(filename[40:44])
        gamma = float(filename[51:55])
        savings = float(filename[29:33])
        
        m, counts = np.genfromtxt(filename, unpack=True)
        counts = np.trim_zeros(counts, 'b')/agents
        m = m[:len(counts)] + 0.025
        
        count = 0
        n = len(counts) - 1 
        while count < 0.20 :
            count += counts[n]
            n -= 1
        
        f2 = lambda x, a: (x**(-1 -a))
        p, s = scipy.optimize.curve_fit(f2, m[-n:], counts[-n:])
        
        print(r'%4.2f & %4.2f & %4.2f & %d & %4.2f & %6.4f \\' % (savings, alpha, gamma, agents, p[0], s[0][0]))
        
        plt.figure(i)
        plt.title(r'$\lambda=%4.2f$ , $\alpha=%4.2f$, $\gamma=%4.2f$, $N=%d$' % (savings, alpha, gamma, agents))    
        plt.plot(m[-n:], counts[-n:], label='pdf')
        plt.plot(m[-n:], m[-n:]**(-1 - p[0]), label=r'fitted, $\alpha=%4.2f$' % (p[0]))
        plt.xlabel(r'Money')
        plt.ylabel(r'Probability distribution')
        plt.legend()

    print(r"\hline");
    print(r"\end{tabularx}")
    
    plt.tight_layout()
    plt.show() 
    
    

filenames = sorted(glob.glob('../output/wealth_pdf_savings_*_alpha_*_gamma_*_agents_*.txt'))    
alphas(filenames)
