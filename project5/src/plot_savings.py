import numpy as np
import matplotlib.pyplot as plt
from matplotlib import markers
import glob
import scipy.special
import scipy.optimize


def pdfs(filenames):
    data_point = ['ro','gs','b<','kd']
    theoretical = ['r', 'g', 'b', 'k']
    print(r"\begin{tabularx}{\textwidth}{X X X X}")
    print(r"\hline")
    print(r"$N$ & $\lambda$ & $\nu$ & $\sigma_{\nu}^2$\\")
    print(r"\hline\hline")       
            
    for i, filename in enumerate(filenames):

        agents = int(filename[63:66])
        savings = float(filename[29:33])
        m, counts = np.genfromtxt(filename, unpack=True)
        counts = np.trim_zeros(counts, 'b')/agents
        m = m[:len(counts)] + 0.025
        
        count = 0
        j = len(counts) - 1 
        while count < 0.20 :
            count += counts[j]
            j -= 1
        
        f2 = lambda x, a: (x**(-1-a))
        p, s = scipy.optimize.curve_fit(f2, m[-j:], counts[-j:])
        
        print(r'%d & %4.2f & %4.2f & %6.4f \\' % (agents, savings, p[0], s[0][0]))
        
        n = 1 + 3*savings/(1-savings)
        pdf = n**n/scipy.special.gamma(n)*m**(n-1)*np.exp(-n*m)*0.05
        
        plt.figure(1)
        plt.plot(m, counts, data_point[i], label=(r'$\lambda=%4.2f$' % savings))
        plt.plot(m, pdf, theoretical[i])
        plt.figure(2)
        plt.loglog(m, counts, data_point[i], label=(r'$\lambda=%4.2f$' % savings))  
        
        plt.figure(i+3)
        plt.plot(m[-j:], counts[-j:], 'o', label='data')
        plt.plot(m[-j:], m[-j:]**(-1 - p[0]), label=r'$\nu=%4.2f$' % (p[0]))
        plt.xlabel('Money')
        plt.ylabel('Probability Distribution')
        plt.margins(x=0.05, y=0.05)
        plt.legend(numpoints=1)
        plt.tight_layout()
        plt.savefig('../fig/Pareto_exponent_savings_%4.2f.pdf' % savings)  
     
    print(r"\hline");
    print(r"\end{tabularx}")
 
        
    plt.figure(1)    
    plt.legend(numpoints=1)
    plt.xlim([0, 2.5])
    plt.ylabel('Probability distribution')
    plt.xlabel('Money')
    plt.tight_layout()
    plt.savefig('../fig/savings_pdf.pdf')
    
    plt.figure(2)
    plt.legend(numpoints=1)
    plt.ylabel('Probability distribution')
    plt.xlabel('Money')
    plt.tight_layout()
    plt.savefig('../fig/savings_pdf_loglog.pdf')
    plt.show()


filenames = sorted(glob.glob('../output/wealth_pdf_savings_*_alpha_0.00_gamma_0.00_agents_500.txt'))    
pdfs(filenames)
