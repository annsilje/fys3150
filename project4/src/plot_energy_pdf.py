import matplotlib.pyplot as plt
import numpy as np
import sys, glob


if len(sys.argv) < 3:
    print('Missing input')
    print('Usage: python make_plots.py <L> <T> [random]')
    sys.exit(1)


l = sys.argv[1]
T = sys.argv[2]
random = '_start_random' if len(sys.argv) > 3 else '_start_all_up'

filename = '../output/%sx%s_temp%s%s_energy_pdf.txt' % (l, l, T, random)

data = np.genfromtxt(filename, unpack=True, skip_header=1)
with open(filename) as pdf_file:
    header = pdf_file.readline().split()
    
std = float(header[-1])

data_1 = np.trim_zeros(data[1], 'b')
len_data_1 = len(data_1)
data_1 = np.trim_zeros(data_1, 'f')
offset = len_data_1 - len(data_1)  
mean = np.dot(data[0][offset:offset+len(data_1)], data_1)/np.sum(data_1)
    
plt.plot(data[0][offset:offset+len(data_1)], data_1/np.sum(data_1))
plt.axvline(mean+std, color='r')
plt.axvline(mean-std, color='r')
plt.ylabel('Probability distribution')
plt.xlabel('Energy')
plt.margins(x=0.01, y=0.01)
plt.tight_layout()
plt.savefig('../fig/energy_pdf_temp%s.pdf' % T)
plt.show()
