#include "initialcondition.h"


InitialCondition::InitialCondition(std::string filename): m_filename(filename) {

}

void InitialCondition::readEphemerides(){
	std::cout << "Reading " << m_filename << std::endl;
	std::ifstream fs_eph;
	fs_eph.open(m_filename);
	std::string line;
	//Skip two header lines
	std::getline(fs_eph, line);
	std::getline(fs_eph, line);

	std::string name;
	double x, y, z;
	double dx, dy, dz;
	double mass;

	while(fs_eph >> name >> x >> y >> z >> dx >> dy >> dz >> mass){
		Particle* p = new Particle(name, vec3(x,y,z), vec3(dx,dy,dz), mass);
		m_input_data.push_back(p);
	}

	fs_eph.close();
}


std::string InitialCondition::getName() {
    return "Solar System";
}

InitialCondition::~InitialCondition(){
	for (int i = 0; i < m_input_data.size(); i++){
		delete m_input_data.at(i);
	}
}

