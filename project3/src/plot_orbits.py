import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from cycler import cycler
import sys, os, argparse


def latex_float(float_number, digits=2):
    float_str = '%.*e' % (digits, float(float_number))
    base, exp = float_str.split("e")
    if int(exp)==0:
        return '$%s$' % (base)
    elif base.startswith('1.') and float(base).is_integer():
        return '$10^{%d}$' % int(exp)
    return '$%s \\times 10^{%d}$' % (base,int(exp))


def setup_parser():
    parser = argparse.ArgumentParser(description='Plot orbits')
    parser.add_argument('filename', metavar='FILE', type=str, help='file with planet positions')
    parser.add_argument('--dim2', action="store_true", help='plot in 2d by ignoring the z-coordinate')
    parser.add_argument('--planet', action="store", dest="planet", type=str, help='Name of planet to plot alone with the and the Sun')
    return parser;

def read_data(filename):
    header = np.genfromtxt(filename, dtype=str, unpack=True, max_rows=1)
    data = np.genfromtxt(filename, unpack=True, skip_header=1)
    return header.tolist(), data


def plot_orbits(header, data, filename, planet, dim2=False):
    cols = len(data)

    filename = os.path.splitext(os.path.basename(filename))[0]
    title = filename.replace('_',' ')

    iterations = int(header.pop())
    step_size = float(header.pop())
    no_bodies = int(cols/3)      

    fig = plt.figure(figsize=(10,8))
    cm = plt.get_cmap('nipy_spectral')
    ax = fig.gca() if dim2 else fig.gca(projection='3d')
    ax.set_prop_cycle(cycler('color', [cm(i/no_bodies) for i in range(no_bodies)]))

    legends = list()
    for i in range(no_bodies):
        if planet and header[i] != 'Sun' and planet!=header[i]:
            continue
        legends.append(header[i])
        if dim2:
            if header[i] == 'Sun':
                ax.plot(data[i*3], data[i*3+1],'.', linewidth=3, rasterized=True)
            else:
                ax.plot(data[i*3], data[i*3+1], rasterized=True)
        else:
            if header[i] == 'Sun':
                ax.plot(data[i*3], data[i*3+1], data[i*3+2],'.' , linewidth=3, rasterized=True)
            else:
                ax.plot(data[i*3], data[i*3+1], data[i*3+2], rasterized=True)
    
    #ax.set_title('%s \n $%4.2f  years, \\Delta t=$%s $years$' % (title, years, latex_float(step_size,1))).set_y(1.05)
    plt.legend(legends, numpoints=1, fontsize='xx-large')    
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    planet = '_' + planet if planet else ''
    if not dim2:
        ax.set_zlabel('z')
    plt.tight_layout()
    if dim2:
        plt.savefig('../fig/%s%s_2D.pdf' % (filename, planet))
    else:
        plt.savefig('../fig/%s%s.pdf' % (filename, planet))
    plt.show()


def main():
    parser = setup_parser()
    args = parser.parse_args()
    filename = args.filename
    print(args.planet)
    header, data = read_data(filename)
    plot_orbits(header, data, filename, args.planet, args.dim2)
    
        
if __name__ == '__main__':
    main()
