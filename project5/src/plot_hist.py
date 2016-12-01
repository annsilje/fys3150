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

filename = '../output/expectation_values_savings_%s_alpha_%s_gamma_%s_agents_%s.txt' % (savings, alpha, gamma, agents)
cycle, mu, sigma2 = np.genfromtxt(filename, unpack=True)

plt.figure()
plt.plot(cycle, mu)
plt.ylabel(r'$\mu_m$')
plt.xlabel('Monte Carlo Cycle')
plt.tight_layout()

plt.figure()
plt.plot(cycle, sigma2)
plt.ylabel(r'$\sigma^2_m$')
plt.xlabel('Monte Carlo Cycle')
plt.tight_layout()



filename = '../output/wealth_pdf_savings_%s_alpha_%s_gamma_%s_agents_%s.txt' % (savings, alpha, gamma, agents)

m, counts = np.genfromtxt(filename, unpack=True)
counts = np.trim_zeros(counts, 'b')
m = m[:len(counts)]
n = 1 + 3*float(savings)/(1-float(savings))
x = m + 0.025

pdf = n**n/scipy.special.gamma(n)*(x)**(n-1)*np.exp(-n*x)

plt.figure()
plt.title(r'$n=%d$' % n)
plt.plot(m, pdf)
plt.hist(m, weights=counts, bins=m, normed=True)
plt.tight_layout()

plt.figure()
plt.semilogy(m, counts)
plt.tight_layout()


plt.show()
