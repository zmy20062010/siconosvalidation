    // User-defined main parameters
#include "SiconosKernel.hpp"

double l_1 = 1.0;//length of the  rod
double m_1 = 1.0; //mass of the  rod
double m = 2.0; //mass of the slider
double g = 10.0; //gravitational acceleration
double Height = 1.5; //height of the fixed point on the link with respect to the slider

double eN = 0.6; //normal restitution coefficient
double eT = 0.0; //tangential restitution coefficient
double mu = 0.2; //coefficient friction
double Fe = 15.0; //applied force on the slider

double theta = 0.5 + 1e-5;

// initial conditions
double theta_1_0 = 0.0;
double x_0 = 0.0;

double omega_1_0 = 0.0;
double v_0 = 0.0;