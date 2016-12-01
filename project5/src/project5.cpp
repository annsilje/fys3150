#include "project5.h"

using namespace std;



/*
**********************************************************************************************************
* Simulation
**********************************************************************************************************
*/

Simulation::Simulation(double savings, double start_money, int agents, double alpha, double gamma): 
    m_agents(agents), m_transactions(10000000), m_savings(savings), m_alpha(alpha), m_gamma(gamma){
    
    m_money = new double[m_agents];
    m_interactions = (int**) matrix(m_agents, m_agents, sizeof(int));
    for(int i=0; i<m_agents; i++){
        m_money[i] = start_money;
        for(int j=0; j<m_agents;j++) m_interactions[i][j] = 0; 
    }       
}


Simulation::~Simulation(){
    delete [] m_money;
    free_matrix((void**) m_interactions);
}


void Simulation::run(mt19937_64& generator){
    uniform_int_distribution<int> pick_agent(0, m_agents-1 );
    uniform_real_distribution<double> pick_amount(0.0, 1.0);
    uniform_real_distribution<double> pick_interaction(0.0, 1.0);

    for(int i=0;i<m_transactions; i++){
        int a = pick_agent(generator);
        int b = pick_agent(generator);
        double epsilon = pick_amount(generator);
        
        double delta_m = fabs(m_money[a] - m_money[b]);
        double probability;
        if (delta_m < 1 || is_close(m_alpha, 0)) {
            probability = 1;
        }
        else if(is_close(m_gamma, 0)){
            probability = pow(delta_m, -m_alpha);
        } 
        else{
            int c = m_interactions[a][b];
            probability = pow(delta_m, -m_alpha)*pow(c + 1, m_gamma);
        }
        
        double interaction = pick_interaction(generator);        
        //double probability = pow(fabs(m_money[a] - m_money[b]), -m_alpha)*pow(c + 1, m_gamma);
        //if (probability > interaction || is_close(m_money[a], m_money[b]) || is_close(m_alpha, 0)){
        if(probability > interaction){
            double total_money = m_money[a] + m_money[b];
            m_money[a] = m_savings*m_money[a] + epsilon * (1 - m_savings) * total_money;
            m_money[b] = m_savings*m_money[b] + (1 - epsilon) * (1 - m_savings) * total_money;
            m_interactions[a][b]++;
            m_interactions[b][a]++;
        }
        
    }
}

/*
**********************************************************************************************************
* Experiment
**********************************************************************************************************
*/

Experiment::Experiment(int cycles, int agents, double savings, double alpha, double gamma): 
    m_cycles(cycles), m_start_money(1.0), m_agents(agents), m_savings(savings), m_alpha(alpha), m_gamma(gamma), m_bin_width(0.05){
    
    m_bins = (int) m_start_money*m_agents/m_bin_width;
    m_acc_bins = new int[m_bins]; 
    
    for(int i=0; i<m_bins;i++) m_acc_bins[i] = 0;
    
    stringstream ss;
    ss << setiosflags(ios::showpoint) << fixed << setprecision(2);
    ss << "../output/wealth_pdf_savings_" << savings << "_alpha_" << alpha << "_gamma_" << gamma << "_agents_" << agents << ".txt";
    m_pdf_file = ss.str();
    
    stringstream ss2;
    ss2 << setiosflags(ios::showpoint) << fixed << setprecision(2);
    ss2 << "../output/expectation_values_savings_" << savings << "_alpha_" << alpha << "_gamma_" << gamma << "_agents_" << agents << ".txt";
    m_exp_file.open(ss2.str().c_str());
    m_exp_file << setiosflags(ios::showpoint) << fixed << setprecision(8);
}


Experiment::~Experiment(){
    delete [] m_acc_bins;
    m_exp_file.close();
}


void Experiment::write(){
    ofstream hist_file;
    hist_file.open(m_pdf_file.c_str());
    
    hist_file << setiosflags(ios::showpoint);
    for (int i=0; i<m_bins; i++){
        hist_file << setw(15) << setprecision(8) << i*m_bin_width << " ";
        hist_file << setw(15) << setprecision(8) << m_acc_bins[i]/m_cycles << endl;
    }
    hist_file.close();
}


void Experiment::write_expectation_values(int cycle){
    double sigma = 0;
    double mu = 0;
    for(int i=0; i<m_bins; i++){
        double m = (i+0.5)*m_bin_width;
        double m2 = m*m;
        double pdf = m_acc_bins[i]/(double)cycle/(double)m_agents;
        sigma += m2*pdf;
        mu += m*pdf;
    }
    sigma = sigma/m_bins - pow((mu/m_bins),2);

    m_exp_file << setw(5) << cycle << " " << setw(10) << mu/m_bins << " " << setw(10) << sigma << endl;
}


void Experiment::accumulate_bins(Simulation& simulation, int cycle){
    double* money = simulation.get_money();    
    for(int i=0; i<m_agents; i++){
        int bin = (int) floor(money[i]/m_bin_width);
        m_acc_bins[bin]++;
    }
    if(cycle > 1 && cycle % 10 == 0){
        write_expectation_values(cycle);
    }
}


void Experiment::run(mt19937_64& generator){
    for (int i=0; i<m_cycles; i++){
        Simulation simulation(m_savings, m_start_money, m_agents, m_alpha, m_gamma);
        simulation.run(generator);
        accumulate_bins(simulation, i);
    }
}



/*
**********************************************************************************************************
* main
**********************************************************************************************************
*/
int main(int argc, char** argv){

    if (argc < 6){
        cout << "Missing input! Usage: project5 <n> <l> <N> <a> <g>" << endl;
        cout << "<n> - number of Monte Carlo cycles" << endl;
        cout << "<l> - amount of money to save in each transaction, range [0,1]" << endl;
        cout << "<N> - number of agents" << endl;
        cout << "<a> - alpha, neighbour iteraction exponent" << endl;
        cout << "<g> - gamma, historical interaction exponent" << endl;
        return 1;
    }

    //Init
    random_device seed;
    mt19937_64 generator(seed());
    clock_t start = clock();

    int cycles = atoi(argv[1]);
    double savings = atof(argv[2]);
    int agents = atoi(argv[3]);
    double alpha = atof(argv[4]);
    double gamma = atof(argv[5]);
        
    //Do experiment  
    Experiment experiment(cycles, agents, savings, alpha, gamma);
    experiment.run(generator);
    
    //Save distribution
    experiment.write();

    cout << "Finished after " << (clock() - start)/(double) CLOCKS_PER_SEC << "seconds" << endl; 
    
}
