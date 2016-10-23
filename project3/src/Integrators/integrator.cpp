#include "integrator.h"
#include "../system.h"

Integrator::Integrator(System* system, double dt, bool sunFixed) {
    m_system = system;
    m_dt = dt;
    m_sun_fixed = sunFixed;
}

void Integrator::setDt(double dt) {
    m_dt = dt;
}

std::string Integrator::getName() {
    return "Unknown";
}
