#pragma once
#include <vector>
#include <string>
#include "../particle.h"

class Integrator {
protected:
    double          m_dt;
    bool            m_sun_fixed = false;
    class System*   m_system    = nullptr;

public:
    Integrator(class System* system, double dt, bool sunFixed);
    void setDt(double dt);
    double getDt() { return m_dt; }
    bool getSunFixed(){ return m_sun_fixed; }
    virtual std::string getName();
    virtual void integrateOneStep(std::vector<Particle*> particles) = 0;
};
