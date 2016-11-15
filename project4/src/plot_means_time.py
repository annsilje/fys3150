import matplotlib.pyplot as plt
import numpy as np
import sys, glob


def plot_time(T, ylims):

    filenames = ['../output/20x20_temp%s_start_random.txt' % (T), '../output/20x20_temp%s_start_all_up.txt' % (T)]
    legends = ['Random', 'All up']
    plotnames = ['accepted', 'energy', 'abs_magnetization']

    style = ['-', '--', '-.', ':']
        
    for j, filename in enumerate(filenames): 
        data = np.genfromtxt(filename, unpack=True)
        ylabels=[r'Percentage of accepted spin flips', r'Mean energy per spin $<E/N>$', 
              r'Mean magnetization per spin $<|M|/N>$']
        
        for i, label in enumerate(ylabels):
            plt.figure(i)
            plt.plot(data[0], data[i+1], style[j])
            plt.ylabel(label)
            plt.xlabel('Monte Carlo Cycle')
            plt.legend(legends)
            plt.margins(x=0.05)
            plt.ylim(ylims[i])
            plt.gca().get_yaxis().get_major_formatter().set_useOffset(False)
            plt.tight_layout()
        
    for i, name in enumerate(plotnames):
        plt.figure(i)
        plt.savefig('../fig/20x20_%s%s.pdf' % (name, T))
    plt.show()


ylims = [[0.0005, 0.0015],[-1.999, -1.993],[0.98, 1.02]]
plot_time('1.0', ylims)

ylims = [[0.266, 0.273],[-1.245, -1.230],[0.42, 0.50]]
plot_time('2.4', ylims)

