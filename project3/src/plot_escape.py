import matplotlib.pyplot as plt
import numpy as np
import glob


def main():
    files = sorted(glob.glob('../output/Two_Bodies_Escape*.txt'))
    legends = [r'Sun', r'$v_{E}=\sqrt{8\pi^2}$', r'$v_{E}=8$', r'$v_{E}=8.8$', r'$v_{E}=8.88$']
    plt.plot([0],[0],'o')
    for filename in files:
        data = np.genfromtxt(filename, unpack=True, skip_header=1)
        x, y = data[3], data[4]
        plt.plot(x,y)

    plt.legend(legends, numpoints=1)         
    plt.xlabel('x')
    plt.ylabel('y')
    plt.tight_layout()
    plt.savefig('../fig/Two_Bodies_Escape_2D.pdf')
    plt.show()
    
        
if __name__ == '__main__':
    main()
