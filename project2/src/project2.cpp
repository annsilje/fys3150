#include <armadillo>
#include <iostream>
#include <cmath>
#include <ctime>
#include <iomanip>

using namespace std;
using namespace arma;

//globals
const double zero = 1.0e-10;

double find_max_off_diag(mat& A, int& k, int& l, int n, double prev_max){
    double max_value = prev_max;
    double start_max = max_value;
    int start_k = k;
    int start_l = l;
    //search column k
    for(int i = 0; i < start_k; i++){
        double a_ik = fabs(A(i,start_k));
        if(a_ik > max_value){
            max_value = a_ik;
            k = i;
            l = start_l;
        }
    }
    //search column l
    for(int i = 0; i < start_l; i++){
        double a_il = fabs(A(i,start_l));
        if(a_il > max_value){
            max_value = a_il;
            k = i;
            l = start_l;
        }
    }
    //search row k 
    for(int j = start_k + 1; j < n; j++){
        double a_kj = fabs(A(start_k,j));
        if(a_kj > max_value){
            max_value = a_kj;
            l = j;
            k = start_k;
        }
    }
    //search row l
    for(int j = start_l + 1; j < n; j++){
        double a_lj = fabs(A(start_l,j));
        if(a_lj > max_value){
            max_value = a_lj;
            l = j;
            k = start_k;
        }
    }
    if(max_value > start_max){
        return max_value;
    }
    max_value = 0;
    //brute force search the rest
    for(int i = 0; i < n; i++){
        if(!(i == k || i == l)){
            //cout << "i=" << i << endl;
            for(int j = i+1; j < n; j++){
                if (!(j == k || j == l)){
                    double a_ij = fabs(A(i,j));
                    //cout << "A(" << i << "," << j << ")=" << a_ij << ", max_value= " << max_value << endl;
                    if (a_ij > max_value){
                        max_value = a_ij;
                        k = i;
                        l = j;
                    }
                }else{
                    //cout << "Skipping column j=" << j << ", k=" << k << ", l=" << l << endl;
                }
            }
        }else{
           //cout << "Skipping row i=" << i << ", k=" << k << ", l=" << l << endl;
        }
    }
    return max_value;
}

double stupid_max(mat& A, int& k, int& l, int n){
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
} //end jacobi_transform


int jacobi(mat& A, mat& R, int n, int max_iter){
    int iterations = 0;
    int k = 0;
    int l = 1;
    double max_off_diag = fabs(A(k,l));

    while(max_off_diag > zero && iterations < max_iter ){
        jacobi_transform(A, R, k, l, n);
  //      max_off_diag = find_max_off_diag(A, k, l, n, max_off_diag);       
        max_off_diag = stupid_max(A, k, l, n);
        iterations++;
    }
    
    return iterations;
}


bool test_orthogonality(mat& R, int n){

    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            vec w_i = R.col(i);
            vec w_j = R.col(j);
            double dot_product = dot(w_i, w_j);
            if(i != j && fabs(dot_product) > zero) return false;
            if(i == j && fabs(dot_product - 1) > zero) return false;
        }
    }
    return true;
}


int main(int argc, char** argv){
    if(argc < 2){
        cout << "Missing input argument: n\n";
        cout << "Usage: " << argv[0] << " <n> [rho_max]\n";
        cout << "\t<n> - number of steps for algorithm\n";
        cout << "\t[rho_max] - upper boundary of rho, where u(rho_max) = 0";
        cout << "Example: " << argv[0] << " 10\n";
        return 1;
    }
    
    
    int n = atoi(argv[1]);
    int max_iter = 1000*n;
    double rho_max = (argc == 3) ? atof(argv[2]) : 8.0;
    double h = rho_max/(n+1);
    
    mat A(n,n, fill::zeros); 
    mat R(n,n, fill::eye);
    
    //Setup A
    double e = -1.0/pow(h,2);
    double d_const = 2/pow(h,2);
    for(int i = 0; i < n; i++){
        if(i>0) A(i, i-1) = e;
        if(i<n-1) A(i, i+1) = e;
        A(i,i) = d_const + pow(i*h,2);
    }
    //Copy A
    mat A_copy(A);
    double start = clock();
    int iter = jacobi(A, R, n, max_iter);
    cout << "Runtime Jacobi:" << (clock() - start)/(double) CLOCKS_PER_SEC << "seconds"<< endl;
    vec diag = sort(diagvec(A));
    cout << "finished after " << iter << " iterations (max=" << max_iter << ")" <<endl;
    cout << "diag_123: " << diag(0) << " " << diag(1) << " " << diag(2) << endl;
    bool ortho_test = test_orthogonality(R, n);
    cout << "Ortho test result: " << ortho_test << endl; 

    start = clock();
    vec diag2 = eig_sym(A_copy);
    cout << "Runtime arma::eig_sym:" << (clock() - start)/(double) CLOCKS_PER_SEC << "seconds"<< endl;
    cout << "diag_123: " << diag2(0) << " " << diag2(1) << " " << diag2(2) << endl; 
}






