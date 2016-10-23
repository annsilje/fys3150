#include "examples.h"
#include "system.h"
#include "particle.h"
#include "Integrators/forwardeuler.h"
#include "Integrators/velocityverlet.h"
#include "Potentials/newtoniangravity.h"
#include "Potentials/modifiednewtoniangravity.h"
#include "Potentials/nopotential.h"
#include "InitialConditions/twobody.h"
#include "InitialConditions/threebody.h"
#include "InitialConditions/solarsystem.h"

using namespace std;

void Examples::twoBodyEuler(string ephFile, double dt, int steps) {
    double G = 4*M_PI*M_PI;
    string outfile = "../output/Two_Bodies_Forward_Euler.txt";

    System* twoBodySystem = new System(outfile);
    twoBodySystem->setIntegrator        (new ForwardEuler(twoBodySystem, dt, true));
    twoBodySystem->setPotential         (new NewtonianGravity(G));
    twoBodySystem->setInitialCondition  (new TwoBody(ephFile));
    twoBodySystem->setFileWriting       (true);
    //twoBodySystem->removeLinearMomentum ();
    twoBodySystem->integrate            (steps);
}


void Examples::twoBodyVerlet(string ephFile, double dt, int steps) {
    double G = 4*M_PI*M_PI;
    string outfile = "../output/Two_Bodies_Velocity_Verlet.txt";

    System* twoBodySystem = new System(outfile);
    twoBodySystem->setIntegrator        (new VelocityVerlet(twoBodySystem, dt, true));
    twoBodySystem->setPotential         (new NewtonianGravity(G));
    twoBodySystem->setInitialCondition  (new TwoBody(ephFile));
    twoBodySystem->setFileWriting       (true);
    //twoBodySystem->removeLinearMomentum ();
    twoBodySystem->integrate            (steps);
}

void Examples::twoBodyEscape(double dt, int steps){
    double G = 4*M_PI*M_PI;

    System* twoBodySystem = new System("../output/Two_Bodies_Escape.txt");
    twoBodySystem->setIntegrator        (new VelocityVerlet(twoBodySystem, dt, true));
    twoBodySystem->setPotential         (new NewtonianGravity(G));
    twoBodySystem->setInitialCondition  (new TwoBody("../input/init_sun_earth.txt"));
    twoBodySystem->setFileWriting       (true);
    //twoBodySystem->removeLinearMomentum ();
    twoBodySystem->integrate            (steps);

    System* twoBodySystem1 = new System("../output/Two_Bodies_Escape1.txt");
    twoBodySystem1->setIntegrator        (new VelocityVerlet(twoBodySystem1, dt, true));
    twoBodySystem1->setPotential         (new NewtonianGravity(G));
    twoBodySystem1->setInitialCondition  (new TwoBody("../input/init_sun_earth_bound1.txt"));
    twoBodySystem1->setFileWriting       (true);
    //twoBodySystem1->removeLinearMomentum ();
    twoBodySystem1->integrate            (steps);

    System* twoBodySystem2 = new System("../output/Two_Bodies_Escape2.txt");
    twoBodySystem2->setIntegrator        (new VelocityVerlet(twoBodySystem2, dt, true));
    twoBodySystem2->setPotential         (new NewtonianGravity(G));
    twoBodySystem2->setInitialCondition  (new TwoBody("../input/init_sun_earth_bound2.txt"));
    twoBodySystem2->setFileWriting       (true);
    //twoBodySystem2->removeLinearMomentum ();
    twoBodySystem2->integrate            (steps);

    System* twoBodySystem3 = new System("../output/Two_Bodies_Escape3.txt");
    twoBodySystem3->setIntegrator        (new VelocityVerlet(twoBodySystem3, dt, true));
    twoBodySystem3->setPotential         (new NewtonianGravity(G));
    twoBodySystem3->setInitialCondition  (new TwoBody("../input/init_sun_earth_bound3.txt"));
    twoBodySystem3->setFileWriting       (true);
    //twoBodySystem3->removeLinearMomentum ();
    twoBodySystem3->integrate            (steps);

    System* twoBodySystem4 = new System("../output/Two_Bodies_Escape4.txt");
    twoBodySystem4->setIntegrator        (new VelocityVerlet(twoBodySystem4, dt, true));
    twoBodySystem4->setPotential         (new NewtonianGravity(G));
    twoBodySystem4->setInitialCondition  (new TwoBody("../input/init_sun_earth_bound4.txt"));
    twoBodySystem4->setFileWriting       (true);
    //twoBodySystem4->removeLinearMomentum ();
    twoBodySystem4->integrate            (steps);
}

void Examples::newtonianMercury(double dt, int steps){
    double G = 4*M_PI*M_PI;
    string outfile = "../output/Mercury_Newtonian_Gravity.txt";

    System* twoBodySystem = new System(outfile);
    twoBodySystem->setIntegrator        (new VelocityVerlet(twoBodySystem, dt, true));
    twoBodySystem->setPotential         (new NewtonianGravity(G));
    twoBodySystem->setInitialCondition  (new TwoBody("../input/init_sun_mercury.txt", "Mercury"));
    twoBodySystem->setFileWriting       (true);
    //twoBodySystem->removeLinearMomentum ();
    twoBodySystem->integrate            (steps);
}

void Examples::modifiedNewtonianMercury(double dt, int steps){
    double G = 4*M_PI*M_PI;
    double c = 299792458/(149597871e3/(365.25*86400)); // AU/year
    string outfile = "../output/Mercury_Modified_Newtonian_Gravity.txt";

    System* twoBodySystem = new System(outfile);
    twoBodySystem->setIntegrator        (new VelocityVerlet(twoBodySystem, dt, true));
    twoBodySystem->setPotential         (new ModifiedNewtonianGravity(G, c));
    twoBodySystem->setInitialCondition  (new TwoBody("../input/init_sun_mercury.txt", "Mercury"));
    twoBodySystem->setFileWriting       (true);
    //twoBodySystem->removeLinearMomentum ();
    twoBodySystem->integrate            (steps);
}

void Examples::threeBodyProblem(string ephFile, double dt, int steps) {
	double G = 4*M_PI*M_PI;


    System* threeBodySystemFixed = new System("../output/Three_Bodies_Fixed.txt");
    threeBodySystemFixed->setIntegrator        (new VelocityVerlet(threeBodySystemFixed, dt, true));
    threeBodySystemFixed->setPotential         (new NewtonianGravity(G));
    threeBodySystemFixed->setInitialCondition  (new ThreeBody(ephFile));
    threeBodySystemFixed->setFileWriting       (true);
    //threeBodySystemFixed->removeLinearMomentum ();
    threeBodySystemFixed->integrate            (steps);


    System* threeBodySystem = new System("../output/Three_Bodies.txt");
    threeBodySystem->setIntegrator        (new VelocityVerlet(threeBodySystem, dt, false));
    threeBodySystem->setPotential         (new NewtonianGravity(G));
    threeBodySystem->setInitialCondition  (new ThreeBody(ephFile));
    threeBodySystem->setFileWriting       (true);
    threeBodySystem->removeLinearMomentum ();
    threeBodySystem->integrate            (steps);
}

void Examples::threeBodyBigJupiter(string ephFile, double dt, int steps){
	double G = 4*M_PI*M_PI;

	System* threeBodySystemFixed = new System("../output/Three_Bodies_Big_Jupiter_Fixed.txt");
	threeBodySystemFixed->setIntegrator        (new VelocityVerlet(threeBodySystemFixed, dt, true));
	threeBodySystemFixed->setPotential         (new NewtonianGravity(G));
	threeBodySystemFixed->setInitialCondition  (new ThreeBody(ephFile, 10));
	threeBodySystemFixed->setFileWriting       (true);
	//threeBodySystemFixed->removeLinearMomentum ();
	threeBodySystemFixed->integrate            (steps);

	System* threeBodySystem = new System("../output/Three_Bodies_Big_Jupiter.txt");
	threeBodySystem->setIntegrator        (new VelocityVerlet(threeBodySystem, dt, false));
	threeBodySystem->setPotential         (new NewtonianGravity(G));
	threeBodySystem->setInitialCondition  (new ThreeBody(ephFile, 10));
	threeBodySystem->setFileWriting       (true);
	threeBodySystem->removeLinearMomentum ();
	threeBodySystem->integrate            (steps);
}

void Examples::threeBodyHugeJupiter(string ephFile, double dt, int steps){
	double G = 4*M_PI*M_PI;

	System* threeBodySystemFixed = new System("../output/Three_Bodies_Huge_Jupiter_Fixed.txt");
	threeBodySystemFixed->setIntegrator        (new VelocityVerlet(threeBodySystemFixed, dt, true));
	threeBodySystemFixed->setPotential         (new NewtonianGravity(G));
	threeBodySystemFixed->setInitialCondition  (new ThreeBody(ephFile, 1000));
	threeBodySystemFixed->setFileWriting       (true);
	//threeBodySystemFixed->removeLinearMomentum ();
	threeBodySystemFixed->integrate            (steps);

	System* threeBodySystem = new System("../output/Three_Bodies_Huge_Jupiter.txt");
	threeBodySystem->setIntegrator        (new VelocityVerlet(threeBodySystem, dt, false));
	threeBodySystem->setPotential         (new NewtonianGravity(G));
	threeBodySystem->setInitialCondition  (new ThreeBody(ephFile, 1000));
	threeBodySystem->setFileWriting       (true);
	threeBodySystem->removeLinearMomentum ();
	threeBodySystem->integrate            (steps);
}


void Examples::solarSystemProblem(string ephFile, double dt, int steps) {
	double G = 4*M_PI*M_PI;
	string outfile = "../output/Solar_System.txt";

    System* solarSystem = new System(outfile);
    solarSystem->setIntegrator        (new VelocityVerlet(solarSystem, dt, false));
    solarSystem->setPotential         (new NewtonianGravity(G));
    solarSystem->setInitialCondition  (new SolarSystem(ephFile));
    solarSystem->setFileWriting       (true);
    solarSystem->removeLinearMomentum ();
    solarSystem->integrate            (steps);
}
