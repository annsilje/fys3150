//System includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>

//Custom includes
#include "lib.h"

using namespace std;

void tridiagonal_solver(double* a, double* b, double* c, double* f, double* v, int n){
    double* d = new double[n];
    double* w = new double[n];
    int i;
        
    //Setup diagonal elements and column vector after gauss elimination 
    d[0] = b[0];
    w[0] = f[0];

    for(i = 1; i < n; i++){
        d[i] = b[i] - a[i] * c[i-1] / d[i-1]; // 3 flop
        w[i] = f[i] - a[i] * w[i-1] / d[i-1];  // 3 flop
    }
    
    //Backward substitution
    v[n-1] = w[n-1] / d[n-1];
    for(i = n-2; i >= 0; i--) v[i] = (w[i] - c[i] * v[i+1]) / d[i]; // 3 flop
    
    delete []d;
    delete []w;
    
}

void custom_solver(double* f, double* v, int n){
    double* d = new double[n];
    double* w = new double[n];
    int i;
    
    //Setup diagonal elements and column vector after gauss elimination 
    d[0] = 2;
    w[0] = f[0];
    
    for(i = 1; i < n; i++){
        d[i] = (i + 2) / (double)(i+1); // 1 flop 
        w[i] = f[i] + w[i-1] / d[i-1];  // 2 flop   
    }
    
    //Backward substitution
    v[n-1] = w[n-1] / d[n-1];
    for(i = n-2; i >= 0; i--) v[i] = (w[i] + v[i+1]) / d[i]; // 2 flop
    
    delete []d;
    delete []w;
}

void run_program(ofstream& fs_clock, int n){
    //step length
    double h = 1/double(n+1);
    //The three-diagonal matrix elements
    double* a = new double[n];
    double* b = new double[n];
    double* c = new double[n];
    //Grid points and column vector of equation system
    double* x = new double[n];
    double* f = new double[n];
    //The unknown variable
    double* v_generic = new double[n];
    double* v_custom = new double[n];
    //Iterator
    int i;    

    //Setup tridiagonal elements, gridpoints and column vector f
    for(i = 0; i < n; i++){
        a[i] = -1; //First element is not used
        b[i] = 2;
        c[i] = -1; //Last element is not used
        x[i] = i*h;
        f[i] = 100 * exp(-10 * x[i]);
    } 
    
    
    //Solve equation system with generic solver
    clock_t start = clock();
    tridiagonal_solver(a, b, c, f, v_generic, n);
    for(i = 0; i < n; i++) v_generic[i] *= pow(h,2);
    double runtime_generic = (clock() - start)/(double)CLOCKS_PER_SEC;
    fs_clock << n << " " << runtime_generic << " ";
    
    //Solve equation system with special solver
    start = clock();
    custom_solver(f, v_custom, n);
    for(i = 0; i < n; i++) v_custom[i] *= pow(h,2);
    double runtime_custom = (clock() - start)/(double)CLOCKS_PER_SEC;
    fs_clock << runtime_custom << " "; 
    fs_clock << runtime_generic / runtime_custom << " ";
    
    if (n<=1000){ //Skip if n is too large due to high memory usage and time consumption
        //Solve equation system with library functions
        double** a_matrix = (double**) matrix(n, n, sizeof(double));
        int j;
        for (i=0;i<n;i++){
            for (j=0; j<n; j++){
                a_matrix[i][j] = 0;
            }
            a_matrix[i][i] = b[i];
            if (i > 0) a_matrix[i-1][i] = c[i];
            if (i > 0) a_matrix[i][i-1] = a[i];
        }
        double flag; 
        int* permutations = new int[n];
        start = clock();
        ludcmp(a_matrix, n, permutations, &flag); //the origianl a_matrix is lost
        lubksb(a_matrix, n, permutations, f); //the original f is lost
        for(i = 0; i < n; i++) f[i] *= pow(h,2);
        double runtime_lud = (clock() - start)/(double)CLOCKS_PER_SEC;
        fs_clock << runtime_lud << " " ;
        fs_clock << runtime_lud / runtime_custom << " ";
        free_matrix((void**)a_matrix);
        delete []permutations;
    }
    fs_clock << endl;
    
    //Save results
    ofstream fs;
    stringstream ss;
    ss << "../output/" << n << "_steps.txt"; 
    fs.open(ss.str().c_str());
    fs << setiosflags(ios::showpoint) << setprecision(16);   
    if (n<=1000){
        for(i=1;i<n;i++) fs << setw(22) << x[i] << setw(22) << v_generic[i] 
                            << setw(22) << v_custom[i] 
                            << setw(22) << f[i] << endl;
    }
    else {
        for(i=1;i<n;i++) fs << setw(22) << x[i] << setw(22) << v_generic[i] 
                            << setw(22) << v_custom[i] << endl;
                            
    }
    fs.close();
    
    //Free memory
    delete []a;
    delete []b;
    delete []c;
    delete []x;
    delete []f;
    delete []v_generic;
    delete []v_custom;
}


int main(int argc, char* argv[]){
    if(argc < 2){
        cout << "Missing input argument: n\n";
        cout << "Usage: " << argv[0] << " <n> [n] [n] ...\n";
        cout << "\t<n> - number of steps in first iteration of algorithm\n";
        cout << "\t[n] - number of steps in next iteration of algorithm\n";
        cout << "Example: " << argv[0] << " 10\n";
        return 1;
    }

    ofstream fs_clock;
    fs_clock.open("../output/runtimes.txt");   
    fs_clock << setiosflags(ios::showpoint | ios::scientific) << setprecision(10);

    int i;
    for (i=1; i<argc;i++) run_program(fs_clock, atoi(argv[i]));
    
    fs_clock.close();   
    return 0;
}
