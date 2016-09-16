import matplotlib.pyplot as plt
import numpy as np
import glob

def latex_float(float_number, digits=2):
    float_str = '%.*e' % (digits, float(float_number))
    base, exp = float_str.split("e")
    if int(exp)==0:
        return '$%s$' % (base)
    elif base.startswith('1.') and float(base).is_integer():
        return '$10^{%d}$' % int(exp)
    return '$%s \\times 10^{%d}$' % (base,int(exp))
    

def make_timetable():
    #Print recorded runtime as latex table
    print("\\begin{tabular}{r|r|r|r|r|r}")
    print("n & Generic & Custom & G/C & LUD & LUD/C \\\\")
    print("\\hline")

    with open('../output/runtimes.txt') as runtimes:
        for line in runtimes:
            n, gen, spec, diff, *lud = line.split()
            if lud:
                print('%s & %ss & %ss & %s & %ss & %s \\\\' % (
                        latex_float(n), latex_float(gen), 
                        latex_float(spec), latex_float(diff),
                        latex_float(lud[0]), latex_float(lud[1])))
            else:
                print('%s & %ss & %ss & %ss & & \\\\' % (
                        latex_float(n), latex_float(gen),
                        latex_float(spec), latex_float(diff)))
    print("\\hline");
    print("\\end{tabular}")


def make_solution_plots_and_table():
    files = sorted(glob.glob('../output/*_steps.txt'), reverse=True)

    cols = 2
    rows = int(np.ceil(len(files)/cols))
    fig1 = plt.figure(1, figsize=(10,12))
    fig2 = plt.figure(2, figsize=(10,12))
    fig3 = plt.figure(3, figsize=(10,6))
    fig4 = plt.figure(4, figsize=(10,12))

    step_length = []
    max_rel_error = []
    
    #Read files and compare solutions
    for i, filename in enumerate(files):
        x, v_generic, v_custom, *v_lud = np.genfromtxt(filename, unpack=True)
        analytical = 1 - (1 - np.exp(-10)) * x - np.exp(-10*x)
        rel_error = np.abs((v_generic - analytical)/analytical)
        n = len(x) + 1
    
        plt.figure(1)
        plt.subplot(rows, cols, i+1)
        plt.title(r'n=%s' % (latex_float(n)))
        plt.plot(x, v_generic)
        plt.plot(x, analytical)
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.legend(['Generic', 'Analytical']) 

        plt.figure(2)
        plt.subplot(rows, cols, i+1)
        plt.title(r'n=%s' % (latex_float(n)))
        plt.plot(x, v_custom)
        plt.plot(x, analytical)
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.legend(['Custom', 'Analytical']) 
                
        if v_lud:
            plt.figure(3)
            plt.subplot(2, cols, i+1)
            plt.title(r'n=%s' % (latex_float(n)))
            plt.plot(x, v_lud[0])
            plt.plot(x, analytical)
            plt.xlabel('x')
            plt.ylabel('u(x)')
            plt.legend(['LU-decomposition', 'Analytical']) 

        step_length.append(1/(n+1))
        max_rel_error.append(np.max(rel_error))
        
        plt.figure(4)
        plt.subplot(rows, cols, i+1)
        plt.title(r'n=%s' % (latex_float(n)))
        plt.plot(x, rel_error)
        plt.margins(x=0.05, y=0.05)
        plt.xlabel('x')
        plt.ylabel('Relative Error')
    
    plt.figure(1)    
    plt.tight_layout()
    plt.savefig('../fig/generic_vs_analytical.pdf')
    
    plt.figure(2)    
    plt.tight_layout()
    plt.savefig('../fig/custom_vs_analytical.pdf')

    plt.figure(3)    
    plt.tight_layout()
    plt.savefig('../fig/lud_vs_analytical.pdf')
    
    plt.figure(4)
    plt.tight_layout()
    plt.savefig('../fig/rel_error.pdf')

    plt.show()

    #Maximum error table
    print("\\begin{tabular}{r|r|r}")
    print("Step size (h) & Max relative error($\\epsilon$) & Max $log(\\epsilon)$\\\\")
    print("\\hline")
    for h, err in zip(step_length, max_rel_error):
        print("%s & %s & %s\\\\" % (latex_float(h, 0), latex_float(err, 8), latex_float(np.log(err), 2)))

    print("\\hline");
    print("\\end{tabular}")


make_timetable()
make_solution_plots_and_table()



        
