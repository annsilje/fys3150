#include "system.h"
#include "Integrators/integrator.h"
#include "Potentials/potential.h"
#include "Potentials/newtoniangravity.h"
#include "Potentials/modifiednewtoniangravity.h"
#include "InitialConditions/initialcondition.h"
#include "particle.h"

#include <iostream>
#include <iomanip>
using namespace std;

System::System(string outfile){
	m_outFilename = outfile;
}

void System::computeForces() {
    resetAllForces();
    m_potential->resetPotentialEnergy();
    m_potential->computeForces(m_particles);
}

void System::resetAllForces() {
    for (int i=0; i<m_numberOfParticles; i++) {
        m_particles.at(i)->resetForces();
    }
}

void System::setPotential(Potential* potential) {
    m_potential = potential;
}

void System::setIntegrator(Integrator* integrator) {
    m_integrator = integrator;
}

void System::setInitialCondition(InitialCondition* initialCondition) {
    m_initialCondition = initialCondition;
    m_initialCondition->setupParticles(*this);
}

void System::setDt(double dt) {
    m_integrator->setDt(dt);
}

void System::integrate(int numberOfSteps) {
	double start = clock();
    m_integrateSteps = numberOfSteps;
    int printModulus = numberOfSteps/20;
    printIntegrateInfo(0, printModulus);
    for (int i=1; i<numberOfSteps+1; i++) {
        m_integrator->integrateOneStep(m_particles);
        printIntegrateInfo(i, printModulus);
        writePositionsToFile();
    }
    closeOutFile();

    cout << "Finished in " << (clock() - start)/(double) CLOCKS_PER_SEC << "seconds." << endl;
}

void System::addParticle(Particle* p) {
    m_particles.push_back(p);
    cout << p->getName() << " ";
    m_numberOfParticles += 1;
}

void System::computeKineticEnergy() {

	m_kineticEnergy = 0;
	for (int i=0; i < m_numberOfParticles; i++){
		Particle* p = m_particles.at(i);
		double m = p->getMass();
		double v2 = p->velocitySquared();
		m_kineticEnergy += 0.5*m*v2;
	}
}
void System::computeMomenta(){
	//angular momentum for all bodies around the center of mass at (0,0,0)
	m_linearMomentum = vec3(0,0,0);
	m_angularMomentum = vec3(0,0,0);
    for (int i=0; i < m_numberOfParticles; i++){
    	vec3 r = m_particles.at(i)->getPosition();
    	double m = m_particles.at(i)->getMass();
    	vec3 v = m_particles.at(i)->getVelocity();
    	vec3 p = m*v;
    	m_linearMomentum += p;
    	m_angularMomentum += r.cross(p);
    }
}

void System::printIntegrateInfo(int stepNumber, int printModulus) {
    if (stepNumber == 0) {
        cout << endl
             << " STARTING INTEGRATION "    << endl
             << "-------------------------" << endl
             << "  o Number of steps:     " << m_integrateSteps << endl
             << "  o Time step, dt:       " << m_integrator->getDt() << endl
			 << "  o Duration (years):    " << m_integrateSteps*m_integrator->getDt()<<endl
             << "  o Initial condition:   " << m_initialCondition->getName() << endl
             << "  o Number of particles: " << m_particles.size() << endl
             << "  o Potential in use:    " << m_potential->getName() << endl
             << "  o Integrator in use:   " << m_integrator->getName() << endl
			 << "  o Sun fixed to origin: " << m_integrator->getSunFixed() << endl
			 << "  o Ephemerides file:    " << m_initialCondition->getEphFile() << endl
             << endl;
    }
    else if (stepNumber % printModulus == 0 || stepNumber == 1) {

        computeKineticEnergy();
        computeMomenta();
        m_potentialEnergy   = m_potential->getPotentialEnergy();
        m_totalEnergy       = m_kineticEnergy + m_potentialEnergy;

        printf("Step: %10d E=%10.5e Ek=%10.5e Ep=%10.5e P=%10.5e L=%10.5e\n",
                       stepNumber, m_totalEnergy, m_kineticEnergy, m_potentialEnergy,
					   m_linearMomentum.length(), m_angularMomentum.length());
        fflush(stdout);

//        cout << fixed << setprecision(8);
//        cout << "Step: " << setw(10) << stepNumber <<
//        		" E=" << setw(10) << m_totalEnergy <<
//        		" Ek=" << setw(10) << m_kineticEnergy <<
//				" Ep=" << setw(10) << m_potentialEnergy <<
//				" P=" << setw(10) << m_linearMomentum.length() <<
//				" L=" << setw(10) <<m_angularMomentum.length() << endl;
    }
}

void System::removeLinearMomentum() {

    vec3 totalMomentum = vec3(0,0,0);
    int sun_index;
    for (int i=0; i< m_numberOfParticles; i++){
    	Particle* p = m_particles.at(i);
    	if (p->getName() != "Sun"){
    		totalMomentum += p->getMass()*p->getVelocity();
    	}
    	else{
    		sun_index = i;
    	}
    }
    Particle* sun = m_particles.at(sun_index);
    vec3& sun_velocity = sun->getVelocity();
    cout << endl;
    cout << "Total momentum of planets  " << totalMomentum << " solarmasses*AU/year"<< endl;
    cout << "Changing Sun veloctiy from " << sun->getVelocity()
    	 << " to " << -1.0*totalMomentum << endl;
    sun_velocity = -1.0*totalMomentum;

    totalMomentum = vec3(0,0,0);
    for (int i=0; i< m_numberOfParticles; i++){
    	Particle* p = m_particles.at(i);
      	totalMomentum += p->getMass()*p->getVelocity();
    }
    cout << "Total momentum of system   " << totalMomentum << " solarmasses*AU/year"<< endl;
}

void System::setFileWriting(bool writeToFile) {
    m_writeToFile = writeToFile;
}

void System::writePositionsToFile() {
	if (!m_writeToFile) return;

	if (m_outFileOpen == false){
		//On first write open the file and write a header
		m_outFile.open(m_outFilename, ios::out);
	    m_outFile << setiosflags(ios::scientific) << setprecision(10);
	    m_outFileOpen = true;

	    //Write planet names as header info, same order as data columns
		for (int i=0; i < m_numberOfParticles; i++){
			Particle* p = m_particles.at(i);
			m_outFile << p->getName() << " ";
		}
		//Add step size and number of steps to header also
		m_outFile << m_integrator->getDt() << " " << m_integrateSteps << endl;
	}

    for (int i=0; i < m_numberOfParticles; i++){
    	vec3 x = m_particles.at(i)->getPosition();
    	m_outFile << setw(12) << x[0] << " " <<
    				 setw(12) << x[1] << " " <<
					 setw(12) << x[2] << " ";
    }
    m_outFile << endl;

}

void System::closeOutFile() {
    if (m_writeToFile == true) {
        m_outFile.close();
        m_outFileOpen = false;
    }
}






