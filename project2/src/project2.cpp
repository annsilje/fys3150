#include <armadillo>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <iomanip>

using namespace std;
using namespace arma;

//globals
const double zero = 1.0e-10;
enum potential {ONE, TWO_NO_INTERACTION, TWO_INTERACTION};

/* Searches the upper triangular part of matrix A of dimensions nxn for
*  the larges absolute value.
*
*  Returns the largest absolute value and the indices k and l of this element are set
*/
double max_off_diagonal(mat& A, int& k, int& l, int n){
    double max_value = 0;
    
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            double a_ij = fabs(A(i,j));
            if (a_ij > max_value){
                max_value = a_ij;
                k = i;
                l = j;
            }
        }
    }
    return max_value;
}

/*
*/
void jacobi_transform(mat& A, mat& R, int k, int l, int n){
    //determine sin_x and cos_x
    double sin_x, cos_x;
    if(fabs(A(k,l)) > zero){
        double t, tau, root;
        tau = (A(l,l) - A(k,k))/(2*A(k,l));
        root = sqrt(1.0 + pow(tau,2));
        //comment if
        t = (tau >= 0) ? 1.0/(tau + root) : -1.0/(root - tau);
        cos_x = 1/sqrt(1+pow(t,2));
        sin_x = cos_x*t;
    }else{
        cos_x = 1.0;
        sin_x = 0.0;
    }

    //Calculate the new transformed matrix elements 
    double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
    a_kk = A(k,k);
    a_ll = A(l,l);
    A(k,k) = pow(cos_x,2)*a_kk - 2.0*cos_x*sin_x*A(k,l) + pow(sin_x,2)*a_ll;
    A(l,l) = pow(sin_x,2)*a_kk + 2.0*cos_x*sin_x*A(k,l) + pow(cos_x,2)*a_ll;
    A(k,l) = 0.0; 
    A(l,k) = 0.0;
    
    for(int i = 0; i < n; i++){
        if(i != k && i != l){
            a_ik = A(i,k);
            a_il = A(i,l);
            A(i,k) = cos_x*a_ik - sin_x*a_il;
            A(k,i) = A(i,k);
            A(i,l) = cos_x*a_il + sin_x*a_ik;
            A(l,i) = A(i,l);
        }
        //Calculate new eigenvectors
        r_ik = R(i,k);
        r_il = R(i,l);
        R(i,k) = cos_x*r_ik - sin_x*r_il;
        R(i,l) = cos_x*r_il + sin_x*r_ik;
    }
    return;
}


int jacobi(mat& A, mat& R, int n, int max_iter){
    int iterations = 0;
    int k, l;
    double max_off_diag = max_off_diagonal(A, k, l, n);

    while(max_off_diag > zero && iterations < max_iter ){
        jacobi_transform(A, R, k, l, n);
        max_off_diag = max_off_diagonal(A, k, l, n);
        iterations++;
    }
    
    return iterations;
}


void setup_matrix(mat& A, double h, int n, potential potential, double omega_r=1.0){
    
    A.fill(0.0);
    //Setup constant elements
    double e = -1.0 / pow(h,2);
    double d_const = 2 / pow(h,2);
    double* rho = new double[n];
    
    for(int i = 0; i < n; i++){
        if(i > 0) A(i,i-1) = e;
        if(i < n-1) A(i,i+1) = e;
        A(i,i) = d_const;
        rho[i] = (i+1)*h;
    }
    //Add potential function
    switch (potential){
        case ONE: {
            for(int i = 0; i < n; i++) A(i,i) += pow(rho[i],2);
            break;
        }
        case TWO_NO_INTERACTION: {
            for(int i = 0; i < n; i++){
                A(i,i) += pow(omega_r, 2)*pow(rho[i], 2);
            }
            break;
        }

        case TWO_INTERACTION: {
            for(int i = 0; i < n; i++){
                A(i,i) += pow(omega_r, 2)*pow(rho[i], 2) + 1/rho[i];
            }
            break;
        }
    }
    delete []rho;
}



string create_filename(potential potential, double omega_r=1.0){
    if(potential == ONE){
        return "../output/One_particle.txt";
    }
    stringstream ss;
    if(potential == TWO_NO_INTERACTION){
        ss << setiosflags(ios::showpoint|ios::fixed) << setprecision(2); 
        ss << "../output/Two_particles_no_interaction_frequency_" << setfill('0') << setw(4)  << omega_r << ".txt";
    }
    if(potential == TWO_INTERACTION){
        ss << setiosflags(ios::showpoint|ios::fixed) << setprecision(2); 
        ss << "../output/Two_particles_with_interaction_frequency_" << setfill('0') << setw(4)  << omega_r << ".txt";
    }

    return ss.str(); 

}

double get_rho_max(potential potential, double omega_r){
    //Decent values for rho_max is found by trial and error
    //These values give reasonable results for n=>400 for all cases
    if(potential == ONE) return 10;
    if(fabs(omega_r - 0.01) < zero) return 50;
    if(fabs(omega_r - 0.10) < zero) return 20;
    if(fabs(omega_r - 0.50) < zero) return 15;
    if(fabs(omega_r - 1.00) < zero) return 10;
    if(fabs(omega_r - 5.00) < zero) return 5;
    cout << "Warning: Using untested omega_r. Results may be bad" << endl;
    return 10; 
}



bool test_orthogonality(mat& R, int n){

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            vec w_i = R.col(i);
            vec w_j = R.col(j);
            double dot_product = dot(w_i, w_j);
            //cout << "Dot product: " << dot_product << " i=" << i << " j=" << j << endl;
            if(i != j && fabs(dot_product) > zero) return false;
            if(i == j && fabs(dot_product - 1) > zero) return false;
        }
    }
    return true;
}


bool test_max_off_diagonal(){
    //Set up test case
    int n = 5;
    mat A(n,n);
    for(int i = 0; i < n; i++){
        for(int j = i + 1; j < n; j++){
            A(i,j) = i + j;
            A(j,i) = A(i,j);
        }
        A(i,i) = 1; 
    }
    //Perform unit test
    int k, l;
    double max_value = max_off_diagonal(A, k, l, n);
    if (max_value == 2*n-3 && k == n-2 && l == n-1) return true;
    
    return false; 
}


int test_eigenvalues(vec& diag1, vec& diag2, int n){
    int counter = 0;
    for(int i=0; i < n; i++){
        if(fabs(diag1(i) - diag2(i)) < zero) counter++;
    }
    return counter;
}


void solve_problem(ofstream& fs_clock, int n, potential potential, bool do_tests, double omega_r=1.0){
    if(potential==ONE) cout << "Solving for one particle" << endl;
    if(potential==TWO_NO_INTERACTION) cout << "Solving for two particles without " 
                                           << "interaction, frequency=" << omega_r << endl;
    if(potential==TWO_INTERACTION) cout << "Solving for two particles with "
                                           << "interaction, frequency=" << omega_r << endl;
    
    //Init system                                   
    double rho_max = get_rho_max(potential, omega_r);
    double h = rho_max/(1+n);
    int max_iter = n*1000;
    mat A(n,n); 
    mat U(n,n, fill::eye); //Eigenvectors
    setup_matrix(A, h, n, potential, omega_r);
    
    //Solve system
    double start = clock();
    int iter = jacobi(A, U, n, max_iter);
    fs_clock << iter << " " << omega_r << " " << potential << " " << h << " " << rho_max << " "
             << (clock() - start)/(double) CLOCKS_PER_SEC << " ";
    
    //Output results
    vec diag = diagvec(A);
    uvec idx = sort_index(diag);
    diag = diag(idx);
    cout << setiosflags(ios::showpoint | ios::fixed) << setprecision(4);
    cout << "First three eigenvalues from Jacobi: " << diag(0) << " " << diag(1) << " " << diag(2) << endl;

    //Save results 
    ofstream fs;
    string filename = create_filename(potential, omega_r);
    fs.open(filename.c_str());
    fs << rho_max << " " << n << endl;
    fs << setiosflags(ios::showpoint) << setprecision(16);   
    for(int i = 0; i < n; i++){
        fs << U(i, idx(0)) << " " << U(i, idx(1)) << " " << U(i, idx(2)) << endl; 
    }
    fs.close();

    //Perform some tests
    if(do_tests){
        setup_matrix(A, h, n, potential, omega_r);
        //Find three lowest eigenvalues using arma::eig_sym
        double start = clock();
        vec diag2 = eig_sym(A); 
        fs_clock << " " << (clock() - start)/(double) CLOCKS_PER_SEC;
        cout << "First three eigenvalues from arma::eig_sym: " << diag2(0) << " " << diag2(1) << " " << diag2(2) << endl; 
        cout << "Is orthogonality of eigenvectors preserved?: " << test_orthogonality(U, n) << endl; 
        cout << "Unit test result of max_off_diagonal: " << test_max_off_diagonal() << endl;
        cout << "How many eigenvalues are equal for the two methods?: " << test_eigenvalues(diag, diag2, n) << endl;
    }
    fs_clock << endl;
    cout << endl;
    

}

void print_help(){
        cout << "Usage: project2 <n> [test]\n";
        cout << "\t<n> - number of steps for algorithm, recommended n=400\n";
        cout << "\t[test] - add an arbitrary argument to run tests during execution\n";
        cout << "Example: project2 400 something\n";
}



int main(int argc, char** argv){
    if(argc < 2){
        cout << "Missing input arguments \n";
        print_help();
        return 1;
    }
        
    int n = atoi(argv[1]);
    bool do_tests = (argc > 2) ? true : false;
    ofstream fs_clock;
    fs_clock.open("../output/exec_time.txt");
    fs_clock << setiosflags(ios::showpoint | ios::fixed) << setprecision(6);   

    solve_problem(fs_clock, n, ONE, do_tests);

    solve_problem(fs_clock, n, TWO_NO_INTERACTION, do_tests, 0.01);
    solve_problem(fs_clock, n, TWO_NO_INTERACTION, do_tests, 0.1);
    solve_problem(fs_clock, n, TWO_NO_INTERACTION, do_tests, 0.5);
    solve_problem(fs_clock, n, TWO_NO_INTERACTION, do_tests, 1.0);
    solve_problem(fs_clock, n, TWO_NO_INTERACTION, do_tests, 5.0);

    solve_problem(fs_clock, n, TWO_INTERACTION, do_tests, 0.01);
    solve_problem(fs_clock, n, TWO_INTERACTION, do_tests, 0.1);
    solve_problem(fs_clock, n, TWO_INTERACTION, do_tests, 0.5);
    solve_problem(fs_clock, n, TWO_INTERACTION, do_tests, 1.0);
    solve_problem(fs_clock, n, TWO_INTERACTION, do_tests, 5.0);

    fs_clock.close();    
    return 0;
}






