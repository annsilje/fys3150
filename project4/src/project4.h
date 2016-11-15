#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <random>
#include <ctime>
#include "mpi.h"
#include "lib.h"


class Experiment{
private:
    int m_spins;
    long int m_cycles;   
    double m_temp;
    bool m_random_start;
    std::mt19937_64& m_generator;
    
    int m_rank;
    int m_no_rank;
    
    int** m_lattice;
    double m_relative_pdf[17];  
       
    int m_max_abs_energy;
    int* m_energy_pdf;
    long int m_accepted = 0;
    long int m_total_accepted = 0;
    
    //P(E)
    std::ofstream m_pdf_file; 
    std::string m_pdf_filename;    
    bool m_pdf_file_open = false;
    
    //means per cycle
    std::ofstream m_mean_file;
    std::string m_mean_filename;
    bool m_mean_file_open = false;
    
    double m_energy = 0;
    double m_magnetization = 0;
    
    double m_acc_E = 0;
    double m_acc_M = 0;
    double m_acc_E2 = 0;
    double m_acc_M2 = 0;
    double m_acc_abs_M = 0;
    
    double m_local_acc_E = 0;
    double m_local_acc_M = 0;
    double m_local_acc_E2 = 0;
    double m_local_acc_M2 = 0;
    double m_local_acc_abs_M = 0;
    
    double m_mean_E = 0;
    double m_mean_M = 0;
    double m_mean_abs_M = 0;
    double m_var_E = 0;
    double m_var_M = 0;
    double m_heat_capacity = 0;
    double m_susceptibility = 0;
    
    long int start_cycle();
    long int end_cycle();

    void write(int cycle);
    void reset_means();
    void report();
    void metropolis(bool write_to_file);
    void accumulate_averages();

public:
    Experiment(int spins, int cycles, double temp, bool random_start, std::mt19937_64& m_generator, int rank, int no_rank);
    virtual ~Experiment();
    
    double get_mean_energy(){return m_mean_E;}
    double get_mean_abs_magnetization(){return m_mean_abs_M;}
    double get_mean_magnetization(){return m_mean_M;}
    double get_heat_capacity(){return m_heat_capacity;}
    double get_susceptibility(){return m_susceptibility;}
    
    void set_temperature(double temp);
    
    void init();
    void run(bool write_to_file);

    void compute_averages(int cycles, bool local);
    void compute_energy_pdf();

};

inline int periodic(int i, int limit, int add){ 
    return (i+limit+add)%(limit);
}
