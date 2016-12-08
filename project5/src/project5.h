#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <ctime>
#include <cmath>
#include "lib.h"

inline bool is_close(double a, double b){
    return fabs(a - b) < 1e-8;
}

class Simulation{
private:
    int m_agents;
    int m_transactions;
    int m_interactions;
    int m_max_interactions;
    double m_savings;
    double m_alpha;
    double m_gamma;
    double* m_money;
    int** m_c;
    
public:
    Simulation(double savings, double start_money, int agents, double alpha, double gamma);
    virtual ~Simulation();
    
    double* get_money(){return m_money;}
    int get_interactions(){return m_interactions;}
    int get_max_interactions(){return m_max_interactions;}
    void run(std::mt19937_64& generator);
    
};

class Experiment{
private:
    int m_cycles;
    double m_start_money;
    int m_agents;
    double m_savings;
    double m_alpha;
    double m_gamma;
    double m_bin_width;
    int m_bins;
    
    long* m_acc_bins; 
    
    std::string m_pdf_file;
    std::ofstream m_exp_file;
    
public:

    Experiment(int cycles, int agents, double savings, double alpha, double gamma);
    virtual ~Experiment();
    
    void run(std::mt19937_64& generator);
    void write_pdf();
    void accumulate_bins(Simulation& simulation, int cycle);
    void write_cycle_data(Simulation& simulation, int cycle);
};
