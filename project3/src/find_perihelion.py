import matplotlib.pyplot as plt
import numpy as np
import argparse
from scipy.signal import argrelextrema


def read_data(read_all):

    if not read_all:
        min_x, min_y = np.genfromtxt('../output/Mercury_Perihelion_Location.txt', unpack=True)
        return min_x, min_y
            
    filename = '../output/Mercury_Modified_Newtonian_Gravity.txt'
    chunk = 1000000
    more_data = True

    min_x = list()
    min_y = list()
    lines_read = 1 
    while more_data: 
        
        data = np.genfromtxt(filename, unpack=True, skip_header=lines_read, max_rows=chunk)
        lines = len(data[0]) if len(data) != 0 else 0
        lines_read = lines_read + lines
        print("Read %d lines from %s" % (lines_read, filename))
        
        if lines > 0:    
            x = data[3] - data[0]
            y = data[4] - data[1]
            r = np.sqrt(x**2 + y**2)
            idx = argrelextrema(r, np.greater)
            min_x = min_x + x[idx].tolist()
            min_y = min_y + y[idx].tolist()
            
        more_data = False if lines < chunk else True



    min_x = np.array(min_x)
    min_y = np.array(min_y)
    np.savetxt('../output/Mercury_Perihelion_Location.txt', np.transpose([min_x, min_y]))
    return min_x, min_y


def plot_perihelion_precession(min_x, min_y):
    radians = np.arctan(min_y/min_x)
    arcsecs = radians*(360*3600/(2*np.pi))

    plt.plot(arcsecs)
    plt.ylabel(r'Perihelion precession $\theta_p$ (arc seconds)')
    plt.xlabel('Revolution')
    #plt.title('Angular change in position of Mercury\'s perihelion')
    plt.savefig('../fig/Perihelion_Precession.pdf')
    plt.show()

    print("Value of the first perihelion angle: %4.2f arcseconds" % arcsecs[0])
    print("Value of the last perihelion angle: %4.2f arcseconds" % arcsecs[-1])
    print("Difference: %4.2f arcseconds" %(arcsecs[-1] - arcsecs[0]))

def setup_parser():
    parser = argparse.ArgumentParser(description='Find the perihelion coordinates for each revolution and save them to file. Plots the angle.')
    parser.add_argument('--read', action="store_true", help='Reads the original (x,y,z) values and creates a new file with only the perihelion coordinates')
    return parser;

def main():
    parser = setup_parser()
    args = parser.parse_args()
    min_x, min_y = read_data(args.read)
    plot_perihelion_precession(min_x, min_y)
    
if __name__ == '__main__':
    main()
