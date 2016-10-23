#include <iostream>
#include "examples.h"
#include <string>

using namespace std;

enum system{
	TWO_BODY_EULER = 0,
	TWO_BODY_VERLET = 1,
	TWO_BODY_ESCAPE = 2,
	THREE_BODY = 3,
	THREE_BODY_BIG_JUPITER = 4,
	THREE_BODY_HUGE_JUPITER = 5,
	NEWTONIAN_MERCURY = 6,
	MODIFIED_NEWTONIAN_MERCURY = 7,
	//THREE_BODY = 8,
	SOLAR_SYSTEM = 9
};

void print_help(){
	cout << "Usage: ./project3 <scenario> <dt> <steps> [eph_file]" << endl;
	cout << "<scenario> - integer representing the problem to solve. Available options are : " << endl;
	cout << "\t 0 - Earth-Sun problem with the Forward Euler algorithm, sun fixed" << endl;
	cout << "\t 1 - Earth-Sun problem with the Velocity Verlet algorithm, sun fixed" << endl;
	cout << "\t 2 - Planet-Sun with theoretical escape velocity, ignores [eph_file], sun fixed" << endl;
	cout << "\t 3 - Earth-Sun-Jupiter problem " << endl;
	cout << "\t 4 - Earth-Sun-Jupiter problem with 10xJupiter Mass" << endl;
	cout << "\t 5 - Earth-Sun-Jupiter problem with 1000xJupiter Mass" << endl;
	cout << "\t 6 - Mercury-Sun problem with the newtonian gravity, ignores [eph_file], sun fixed" << endl;
	cout << "\t 7 - Mercury-Sun problem with the modified newtonian gravity, ignores [eph_file], sun fixed" << endl;
	//cout << "\t 8 - Earth-Sun-Jupiter problem with the Velocity Verlet algorithm" << endl;
	cout << "\t 9 - Solar System problem with the Velocity Verlet algorithm" << endl;
	cout << "<dt> - time step in AU/year " << endl;
	cout << "<steps> - number of integration steps " << endl;
	cout << "[eph_file] - filename of ephemerides file, optional" << endl;
}

int main(int argc, char** argv) {
	if (argc < 4){
		print_help();
		return 1;
	}

	int scenario = atoi(argv[1]);
	double dt    = atof(argv[2]);
	int steps    = atoi(argv[3]);

	string filename;
	if (argc > 4){
		filename = string(argv[4]);
	}else {
		filename = "../input/init_solar_system_2015_12_01.txt";
	}

	switch(scenario){
		case TWO_BODY_EULER:
			Examples::twoBodyEuler(filename, dt, steps);
			break;
		case TWO_BODY_VERLET:
			Examples::twoBodyVerlet(filename, dt, steps);
			break;
		case TWO_BODY_ESCAPE:
			Examples::twoBodyEscape(dt, steps);
			break;
		case THREE_BODY:
			Examples::threeBodyProblem(filename, dt, steps);
			break;
		case THREE_BODY_BIG_JUPITER:
			Examples::threeBodyBigJupiter(filename, dt, steps);
			break;
		case THREE_BODY_HUGE_JUPITER:
			Examples::threeBodyHugeJupiter(filename, dt, steps);
			break;
		case NEWTONIAN_MERCURY:
			Examples::newtonianMercury(dt, steps);
			break;
		case MODIFIED_NEWTONIAN_MERCURY:
			Examples::modifiedNewtonianMercury(dt, steps);
			break;
//		case THREE_BODY:
//			Examples::threeBodyProblem(filename, dt, steps);
//			break;
		case SOLAR_SYSTEM:
			Examples::solarSystemProblem(filename, dt, steps);
			break;
		default:
			cout << "Error parsing input arguments. Unknown scenario" << endl;
			print_help();
			return 1;
	}

    return 0;
}
