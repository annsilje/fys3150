#pragma once
#include <iostream>
#include <cmath>
#include <string>


class Examples {
public:
    static void twoBodyEuler(std::string ephFile, double dt, int steps);
    static void twoBodyVerlet(std::string ephFile, double dt, int steps);
    static void twoBodyEscape(double dt, int steps);
    static void threeBodyProblem(std::string ephFile, double dt, int steps);
    static void threeBodyBigJupiter(std::string ephFile, double dt, int steps);
    static void threeBodyHugeJupiter(std::string ephFile, double dt, int steps);
    static void newtonianMercury(double dt, int steps);
    static void modifiedNewtonianMercury(double dt, int steps);
    //static void threeBodyProblem(std::string ephFile, double dt, int steps);
    static void solarSystemProblem(std::string ephFile, double dt, int steps);
};
