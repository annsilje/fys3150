import matplotlib.pyplot as plt
import numpy as np
import scipy.special
import sys

if len(sys.argv) < 5:
    print("Missing input!")
    print("python plot_hist.py <l> <a> <g> <N>")
    sys.exit(1)

savings = sys.argv[1]
alpha = sys.argv[2]
gamma = sys.argv[3]
agents = sys.argv[4]

filename = '../output/cycle_data_savings_%s_alpha_%s_gamma_%s_agents_%s.txt' % (savings, alpha, gamma, agents)
cycle, mu, sigma2, *_ = np.genfromtxt(filename, unpack=True)

plt.figure()
plt.plot(cycle, mu)
plt.ylabel(r'$\mu_m$')
plt.xlabel('Monte Carlo Cycle')
plt.margins(x=0.1,y=0.1)
plt.tight_layout()
plt.savefig('../fig/mu_savings_%s_alpha_%s_gamma_%s_agents_%s.pdf' % (savings, alpha, gamma, agents))

plt.figure()
plt.plot(cycle, sigma2)
plt.ylabel(r'$\sigma^2_m$')
plt.xlabel('Monte Carlo Cycle')
plt.margins(x=0.1,y=0.1)
plt.tight_layout()
plt.savefig('../fig/sigma_savings_%s_alpha_%s_gamma_%s_agents_%s.pdf' % (savings, alpha, gamma, agents))


filename = '../output/wealth_pdf_savings_%s_alpha_%s_gamma_%s_agents_%s.txt' % (savings, alpha, gamma, agents)

m, counts = np.genfromtxt(filename, unpack=True)
counts = np.trim_zeros(counts, 'b')
m = m[:len(counts)]
n = 1 + 3*float(savings)/(1-float(savings))
pdf = n**n/scipy.special.gamma(n)*m**(n-1)*np.exp(-n*m)

plt.figure()
plt.plot(m, pdf)
plt.hist(m, weights=counts, bins=m, normed=True)
plt.ylabel('Normalized distribution of agents')
plt.xlabel('Money')
plt.tight_layout()
plt.savefig('../fig/histogram_savings_%s_alpha_%s_gamma_%s_agents_%s.pdf' % (savings, alpha, gamma, agents))

m = m + 0.025

plt.figure()
plt.semilogy(m, counts/int(agents), 'o')
plt.tight_layout()
plt.ylabel('Probability distribution')
plt.xlabel('Money')
plt.tight_layout()
plt.savefig('../fig/pdf_log_savings_%s_alpha_%s_gamma_%s_agents_%s.pdf' % (savings, alpha, gamma, agents))

plt.show()
