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
    print("\\begin{tabular}{r|r|r|r}")
    print("Steps (n) & Generic & Specific & Difference\\\\")
    print("\\hline")

    with open('runtimes.txt') as runtimes:
        for line in runtimes:
            line = line.split()
            print('%s & %ss & %ss & %ss \\\\' % (latex_float(line[0]), latex_float(line[1], 4),
                latex_float(line[2], 4), latex_float(line[3], 4)))
    print("\\hline");
    print("\\end{tabular}")

def make_solution_plots_and_table():
    files = sorted(glob.glob('*_steps.txt'), reverse=True)

    cols = 2
    rows = int(np.ceil(len(files)/cols))
    fig = plt.figure(figsize=(10,12))

    step_length = []
    max_rel_error = []
    
    

    for i, filename in enumerate(files):
        x, v_generic, v_specific = np.genfromtxt(filename, unpack=True)
        analytical = 1 - (1 - np.exp(-10)) * x - np.exp(-10*x)
        n = len(x) + 1
    
        rel_error = np.abs((analytical - v_generic)/analytical)
    
        plt.subplot(rows, cols, i+1)
        plt.title('n=%d' % (n))
        plt.plot(x, v_generic)
        plt.plot(x, analytical)
        plt.xlabel('x')
        plt.ylabel('u(x)')
        plt.legend(['numerical', 'analytical']) 

        step_length.append(1/(n+1))
        max_rel_error.append(np.max(rel_error))
        
    plt.tight_layout()
    plt.savefig('numerical_vs_analytical.pdf')

    fig = plt.figure()
    plt.plot(np.log(step_length), np.log(max_rel_error))
    plt.xlabel('Step size ($log(h)$)')
    plt.ylabel(r'Maximum relative error ($log(\epsilon)$)')
    plt.tight_layout()
    plt.savefig('error.pdf')
    plt.show()


    print("\\begin{tabular}{r|r}")
    print("Step size (h) & Maximum relative error($\\epsilon$)\\\\")
    print("\\hline")
    for h, err in zip(step_length, max_rel_error):
        print("%s & %s \\\\" % (latex_float(h, 0), latex_float(err, 8)))

    print("\\hline");
    print("\\end{tabular}")


make_timetable()
make_solution_plots_and_table()



        
