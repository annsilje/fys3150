#include "project4.h"

using namespace std;


/*
******************************************************************************
* Experiment class
******************************************************************************
*/

Experiment::Experiment(int spins, int cycles, double temp, bool random_start, std::mt19937_64& generator, int rank, int no_rank): 
    m_spins(spins), m_cycles(cycles), m_temp(temp), m_random_start(random_start), m_generator(generator), m_rank(rank), m_no_rank(no_rank){
    
    string txt = random_start ? "_start_random" : "_start_all_up"; 
    stringstream ss;
    ss << setiosflags(ios::showpoint) << fixed;
    ss << "../output/" << spins << "x" << spins <<"_temp" << setw(3) << setprecision(1)  << temp << txt << ".txt";
    m_mean_filename = ss.str(); 
        
    stringstream ss2;
    ss2 << setiosflags(ios::showpoint) << fixed;
    ss2 << "../output/" << spins << "x" << spins << "_temp" << setw(3) << setprecision(1) << temp << txt << "_energy_pdf.txt"; 
    m_pdf_filename = ss2.str();
    
    m_lattice = (int**) matrix(spins, spins, sizeof(int));
    m_max_abs_energy = spins*spins*2;
    m_energy_pdf = new int[m_max_abs_energy*2+1];
    
    set_temperature(m_temp);
}

Experiment::~Experiment(){
    free_matrix((void**) m_lattice);
    delete []m_energy_pdf;
    
    if(m_mean_file_open){
        m_mean_file.close();
    }
    if(m_pdf_file_open){
        m_pdf_file.close();
    }
}

void Experiment::set_temperature(double temp)
{
    m_temp = temp;
    
    for (int dE=-8; dE<=8; dE++) m_relative_pdf[dE+8] = 0;
    for (int dE=-8; dE<=8; dE+=4) m_relative_pdf[dE+8] = exp(-dE/m_temp);
}

long int Experiment::start_cycle(){
    int steps = m_cycles/m_no_rank;
    return m_rank*steps;    
}

long int Experiment::end_cycle(){
    int steps = m_cycles/m_no_rank;
    return (m_rank+1)*steps;
}

void Experiment::init(){
    if(m_random_start){
        cout << "Random start" << endl;
        uniform_int_distribution<int> uniform_int(0, 1);
        for (int i = 0; i < m_spins; i++){
            for (int j = 0; j < m_spins; j++){
                m_lattice[i][j] = (uniform_int(m_generator)==1) ? 1 : -1;
                m_magnetization += (double) m_lattice[i][j];
            }
        }
        for (int i = 0; i < m_spins; i++){
            for (int j = 0; j < m_spins; j++){
                m_energy -= (double)m_lattice[i][j]*(m_lattice[periodic(i, m_spins, -1)][j] +
                                                     m_lattice[i][periodic(j, m_spins, -1)]);
            }
        }
    }
    else{
        cout << "All up start" << endl;
        for (int i=0; i<m_spins;i++){
            for (int j=0; j<m_spins;j++) m_lattice[i][j] = 1;    
        }
        m_magnetization = m_spins*m_spins;
        m_energy = -2*m_spins*m_spins;
        
    }
    cout << "Starting energy: " << m_energy << endl;
    cout << "Starting magnetization: " << m_magnetization << endl;
    cout << "Temperature: " << m_temp << endl;
    cout << "Lattice: " << m_spins << "x" << m_spins << endl;

}


void Experiment::metropolis(bool store_energy){
    uniform_int_distribution<int> uniform_int(0, m_spins - 1);
    uniform_real_distribution<double> uniform_real(0.0, 1.0);
    
    for (int i=0;i<m_spins; i++){
        for(int j=0;j<m_spins;j++){
            int x = uniform_int(m_generator);
            int y = uniform_int(m_generator);

            int dE = 2*m_lattice[y][x]*(m_lattice[y][periodic(x, m_spins, -1)] + 
                                        m_lattice[periodic(y, m_spins, -1)][x] + 
                                        m_lattice[y][periodic(x, m_spins, 1)] +
                                        m_lattice[periodic(y, m_spins, 1)][x]);
                                        
            if(uniform_real(m_generator) <= m_relative_pdf[dE+8]){
                m_lattice[y][x] *= -1; //flip spin
                m_magnetization += (double) 2*m_lattice[y][x];
                m_energy += (double) dE;
                m_accepted++;
            }
            
        }
    }
    if(store_energy) m_energy_pdf[(int)m_energy+m_max_abs_energy]++;
}

void Experiment::reset_means(){
    m_acc_E = 0;
    m_acc_M = 0;
    m_acc_E2 = 0;
    m_acc_M2 = 0;
    m_acc_abs_M = 0;
    
    m_local_acc_E = 0;
    m_local_acc_M = 0;
    m_local_acc_E2 = 0;
    m_local_acc_M2 = 0;
    m_local_acc_abs_M = 0;
}



void Experiment::run(bool write_to_file){
    reset_means();
    m_accepted = 0;
    m_total_accepted = 0;
    
    for (long int i=start_cycle()+1; i<=end_cycle()+1; i++){
        metropolis(false);
        accumulate_averages();
        if (i % 100 == 0 && write_to_file && m_no_rank == 1) write(i);
    }
    MPI_Reduce(&m_local_acc_E, &m_acc_E, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_local_acc_M, &m_acc_M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_local_acc_E2, &m_acc_E2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_local_acc_M2, &m_acc_M2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_local_acc_abs_M, &m_acc_abs_M, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&m_accepted, &m_total_accepted, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if (m_rank == 0) report();
}


void Experiment::compute_energy_pdf(){
    reset_means();
    int cycles = 10000;
    for (int i=0; i<cycles; i++) {
        metropolis(true);
        accumulate_averages();
    }
    
    compute_averages(cycles, true);
    
    if(!m_pdf_file_open){
        m_pdf_file.open(m_pdf_filename.c_str());
        m_pdf_file_open = true;
    }
    
    m_pdf_file << setiosflags(ios::showpoint);
    m_pdf_file << "#Computed variance: " << m_var_E << " Computed standard deviation: " << sqrt(m_var_E) << endl;
    for (int E= -m_max_abs_energy; E <= m_max_abs_energy; E+=4){
        m_pdf_file << setw(15) << setprecision(8) << (double)E/m_spins/m_spins;
        m_pdf_file << setw(8) << m_energy_pdf[E+m_max_abs_energy] << endl;
    }
}


void Experiment::compute_averages(int cycles, bool local){
    double mean_E2;
    double mean_M2;
    
    if (local){
        //use mean values accumulated in one process
        m_mean_E = m_local_acc_E/(double) cycles;
        m_mean_M = m_local_acc_M/(double) cycles;
        m_mean_abs_M = m_local_acc_abs_M/(double) cycles;
        
        mean_E2 = m_local_acc_E2/(double) cycles;
        mean_M2 = m_local_acc_M2/(double) cycles;
    } else {
        //use mean values accumulated for all processes
        m_mean_E = m_acc_E/(double) cycles;
        m_mean_M = m_acc_M/(double) cycles;
        m_mean_abs_M = m_acc_abs_M/(double) cycles;
        
        mean_E2 = m_acc_E2/(double) cycles;
        mean_M2 = m_acc_M2/(double) cycles;
    }
    
    m_mean_E = m_mean_E/m_spins/m_spins;
    m_mean_M = m_mean_M/m_spins/m_spins;
    m_mean_abs_M = m_mean_abs_M/m_spins/m_spins;
    
    mean_E2 = mean_E2/pow(m_spins, 4);
    mean_M2 = mean_M2/pow(m_spins, 4);
    
    m_var_E = (mean_E2 - m_mean_E*m_mean_E);
    m_var_M = (mean_M2 - m_mean_abs_M*m_mean_abs_M);
        
    m_heat_capacity = m_var_E/m_temp/m_temp;
    m_susceptibility = m_var_M/m_temp;
        
    
    
}

void Experiment::report(){

    compute_averages(m_cycles, false);
    
    cout << "Completed " << m_cycles << " Monte Carlo cycles with " << m_spins*m_spins << " spins. " << endl;
    cout << m_total_accepted << " of " << m_cycles*m_spins*m_spins << " spin flips were accepted" << endl;
    cout << setiosflags(ios::showpoint);
    cout << "Temperature:         " << setw(15) << setprecision(8) << m_temp << endl;
    cout << "<E> per spin:        " << setw(15) << setprecision(8) << m_mean_E << endl;
    cout << "Variance E per spin: " << setw(15) << setprecision(8) << m_var_E << endl; 
    cout << "Heat capacity:       " << setw(15) << setprecision(8) << m_heat_capacity << endl; 
    cout << "<M> per spin:        " << setw(15) << setprecision(8) << m_mean_M << endl;
    cout << "variance M per spin: " << setw(15) << setprecision(8) << m_var_M << endl; 
    cout << "Susceptibility:      " << setw(15) << setprecision(8) << m_susceptibility << endl; 
    cout << "<|M|> per spin:      " << setw(15) << setprecision(8) << m_mean_abs_M << endl;
    
}


void Experiment::accumulate_averages(){
    m_local_acc_E += m_energy;
    m_local_acc_M += m_magnetization;
    m_local_acc_E2 += pow(m_energy, 2);
    m_local_acc_M2 += pow(m_magnetization, 2);
    m_local_acc_abs_M += fabs(m_magnetization);
}


void Experiment::write(int cycle){
    if(!m_mean_file_open){
        m_mean_file.open(m_mean_filename.c_str()); 
        m_mean_file_open = true;    
    }

    compute_averages(cycle, true);
    
    m_mean_file << setiosflags(ios::showpoint);
    m_mean_file << setw(8) << cycle << " ";
    m_mean_file << setw(15) << setprecision(8) << m_accepted/(double)cycle/m_spins/m_spins << " ";
    m_mean_file << setw(15) << setprecision(8) << m_mean_E << " ";
    m_mean_file << setw(15) << setprecision(8) << m_mean_abs_M << " ";
    m_mean_file << endl;
}


/*
*******************************************************************************
* Start here
*******************************************************************************
*/

void do_experiment(int spins, mt19937_64& generator, int rank, int no_rank){
    stringstream ss;
    ss << "../output/" << spins << "x" << spins << ".txt";
        
    ofstream fs;
    if(rank == 0) fs.open(ss.str().c_str());
    
    
    double start_T = 2.1;
    double end_T = 2.35;
    double step_T = 0.01;
    long int cycles = 1000000;
    
    clock_t start = clock();
    
    Experiment experiment(spins, cycles, start_T, true, generator, rank, no_rank);    
    experiment.init();
    for (double T = start_T; T<=end_T; T+=step_T){
        experiment.set_temperature(T);
        experiment.run(false);
        
        if (rank == 0) {
            experiment.compute_averages(cycles, false);
        
            fs << setiosflags(ios::showpoint);
            fs << setw(15) << setprecision(8) << T;
            fs << setw(15) << setprecision(8) << experiment.get_mean_energy();
            fs << setw(15) << setprecision(8) << experiment.get_mean_magnetization();
            fs << setw(15) << setprecision(8) << experiment.get_mean_abs_magnetization();
            fs << setw(15) << setprecision(8) << experiment.get_heat_capacity();
            fs << setw(15) << setprecision(8) << experiment.get_susceptibility() << endl;
        }
        
    }
    
    cout << "Process " << rank << " finished with " << spins << "x" << spins << " lattice after " << (clock() - start)/(double) CLOCKS_PER_SEC << "seconds" << endl; 
    
    if (rank == 0) fs.close();
}


int main(int argc, char** argv){
    int no_rank;
    int rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &no_rank);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);     


    random_device seed;
    mt19937_64 generator(seed());
       
    if (argc == 1){
        //no arguments, run full scenario
        
        do_experiment(40, generator, rank, no_rank);
        do_experiment(60, generator, rank, no_rank);
        do_experiment(100, generator, rank, no_rank);
        do_experiment(140, generator, rank, no_rank);

    }
    else if (argc < 4){
        cout << "Missing input" << endl;
        cout << "Usage: project4 <L> <n> <T> [random]" << endl;
        cout << "<L> - lattice size, yielding LxL number of spins" << endl;
        cout << "<n> - number of monte carlo cycles" << endl;
        cout << "<T> - temperature" << endl;
        cout << "[random] - arbitrary argument. Starts the experiment with random spins, otherwise all spins start up." << endl;
        
    }
    else{
       int spins = atoi(argv[1]);
       long int cycles = atoi(argv[2]);
       double temp = atof(argv[3]);
       bool random = (argc > 4) ? true : false;
       clock_t start = clock();
             
       if (no_rank != 1){
            cout << "Warning! No output files are produced when using more than one process ." << endl;
       }
       
       Experiment experiment(spins, cycles, temp, random, generator, rank, no_rank);
       experiment.init();
       experiment.run(true);
       
       cout << "Process " << rank << " finished with " << spins << "x" << spins << " lattice after " << (clock() - start)/(double) CLOCKS_PER_SEC << "seconds" << endl;       
       
       if (no_rank == 1) experiment.compute_energy_pdf();
       
    }
    
    MPI_Finalize();

}
