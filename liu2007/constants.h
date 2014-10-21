    // User-defined main parameters
#include "SiconosKernel.hpp"

double l_1 = 1.0;//length of the first link
double l_2 = 1.0;//length of the second link
double m_1 = 1.0; //mass of the first link
double m_2 = 1.0; //mass of the second link
double m = 2.0; //mass of the slider
double g = 10.0; //gravitational acceleration
double Height = 1.5; //height of the fixed point on the link with respect to the slider

double eN = 0.6; //normal restitution coefficient
double eT = 0.0; //tangential restitution coefficient
double mu = 0.2; //coefficient friction
double T_1 = 5.0;
double T_2 = 10.0;
double F = 15;

double theta = 0.5 + 1e-5;

// initial conditions
double x_1_0 = M_PI / 2.0;
double x_2_0 = M_PI / 2.0;
double x_0 = 0.0;

double v_1_0 = 0.0;
double v_2_0 = 0.0;
double v_0 = 0.0;