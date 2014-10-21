#ifdef _WIN32 
#define SICONOS_EXPORT extern "C" __declspec(dllexport) 
#else 
#define SICONOS_EXPORT extern "C" 
#endif  

#include "SiconosKernel.hpp"
#include "RuntimeException.hpp"

#include <math.h>
#include "constants.h"
#include <math.h>
#include <iostream>

SICONOS_EXPORT void mass(unsigned int sizeOfq, const double *q, double *mass, unsigned int sizeZ, double* z)
{
  //inertia matrix using (\theta)
  /* q -->  x_1,  x_2,  x    */
  /* v -->  v_1,  v_2,  v    */
  mass[0]  = (1.0 / 3 * m_1 + m_2) * l_1 * l_1;
  mass[1]  = 1.0 / 2 * m_2 * l_1 * l_2 * cos(q[0]-q[1]);
  mass[2]  = 0.0;

  mass[3]  = 1.0 / 2 * m_2 * l_1 * l_2 * cos(q[0]-q[1]);
  mass[4]  = 1.0 / 3 * m_2 * l_2 * l_2;  
  mass[5]  = 0.0;

  mass[6]  = 0.0;
  mass[7]  = 0.0;  
  mass[8]  = m;   
}

SICONOS_EXPORT void NNL(unsigned int sizeOfq, const double *q, const double *v, double *NNL, unsigned int sizeZ, double* z)
{
  /* q -->  x_1,  x_2,  x    */
  /* v -->  v_1,  v_2,  v    */
  NNL[0] = 1.0 / 2 * m_2 * l_1 * l_2 * sin(q[0] - q[1]) * v[1] * v[1];
  NNL[1] = 1.0 / 2 * m_2 * l_1 * l_2 * sin(q[0] - q[1]) * v[0] * v[0];
  NNL[2] = 0.0;
}

SICONOS_EXPORT void jacobNNLq(unsigned int sizeOfq, const double *q, const double *v, double *jacobq, unsigned int sizeOfZ, double* z)
{
  /* q -->  x_1,  x_2,  x    */
  /* v -->  v_1,  v_2,  v    */
  jacobq[0] = 1.0 / 2 * m_2 * l_1 * l_2 * cos(q[0] - q[1]) * v[1] * v[1];
  jacobq[1] = 1.0 / 2 * m_2 * l_1 * l_2 * cos(q[0] - q[1]) * v[0] * v[0];
  jacobq[2] = 0.0;

  jacobq[3] = - 1.0 / 2 * m_2 * l_1 * l_2 * cos(q[0] - q[1]) * v[1] * v[1];
  jacobq[4] = - 1.0 / 2 * m_2 * l_1 * l_2 * cos(q[0] - q[1]) * v[0] * v[0];
  jacobq[5] = 0.0;

  jacobq[6] = 0.0;
  jacobq[7] = 0.0;
  jacobq[8] = 0.0;  
}

SICONOS_EXPORT void jacobNNLv(unsigned int sizeOfq, const double *q, const  double *v, double *jacobv, unsigned int sizeOfZ, double* z)
{
  /* q -->  x_1,  x_2,  x    */
  /* v -->  v_1,  v_2,  v    */
  jacobv[0] = 0.0;
  jacobv[1] = m_2 * l_1 * l_2 * sin(q[0] - q[1]) * v[0];
  jacobv[2] = 0.0;

  jacobv[3] = m_2 * l_1 * l_2 * sin(q[0] - q[1]) * v[0];
  jacobv[4] = 0.0;
  jacobv[5] = 0.0;

  jacobv[6] = 0.0;
  jacobv[7] = 0.0;
  jacobv[8] = 0.0;  
}

SICONOS_EXPORT void FInt(double time, unsigned int sizeOfq, const double *q, const  double *v, double *FInt, unsigned int sizeOfZ, double* z)
{
  FInt[0] = 1.0 / 2 * m_1 * g * l_1 * sin(q[0]) + m_2 * g * l_2 * sin(q[0]);
  FInt[1] = 1.0 / 2 * m_2 * g * l_2 * sin(q[1]);
  FInt[2] = 0.0;
}

SICONOS_EXPORT void jacobFintq(double time, unsigned int sizeOfq, const double *q, const  double *v, double *jacobq, unsigned int sizeOfZ, double* z)
{
  jacobq[0] = 1.0 / 2 * m_1 * g * l_1 * cos(q[0]) + m_2 * g * l_2 * cos(q[0]);
  jacobq[1] = 0.0;
  jacobq[2] = 0.0;

  jacobq[3] = 0.0;
  jacobq[4] = 1.0 / 2 * m_2 * g * l_2 * cos(q[1]);
  jacobq[5] = 0.0;

  jacobq[6] = 0.0;
  jacobq[7] = 0.0;
  jacobq[8] = 0.0;
}

SICONOS_EXPORT void jacobFintv(double time, unsigned int sizeOfq, const double *q, const double *v, double *jacobv, unsigned int sizeOfZ, double* z)
{
  jacobv[0] = 0.0;
  jacobv[1] = 0.0;
  jacobv[2] = 0.0;

  jacobv[3] = 0.0;
  jacobv[4] = 0.0;
  jacobv[5] = 0.0;

  jacobv[6] = 0.0;
  jacobv[7] = 0.0;
  jacobv[8] = 0.0;
}

SICONOS_EXPORT void FExt(double time, unsigned int sizeOfq, double *FExt, unsigned int sizeOfZ, double* z)
{
  FExt[0] = T_1;
  FExt[1] = T_2;
  FExt[2] = F;
}

SICONOS_EXPORT void Gap(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* Gap, unsigned int sizeOfZ, double* z)
{
  Gap[0] = Height - l_1 * cos(q[0]) - l_2 * cos(q[1]); //normal gap
  Gap[1] = l_1 * sin(q[0]) + l_2 * sin(q[1]) - q[2];  //tangential gap
}

SICONOS_EXPORT void jacobGapq(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* jacobq, unsigned int sizeOfZ, double* z)
{
  jacobq[0] = l_1 * sin(q[0]);
  jacobq[1] = l_1 * cos(q[0]);

  jacobq[2] = l_2 * sin(q[1]);
  jacobq[3] = l_2 * cos(q[1]);

  jacobq[2] = 0.0;
  jacobq[3] = - 1.0;  
}

SICONOS_EXPORT void jacobGapt(unsigned int sizeOfq, const double* q, unsigned int sizeOfY, double* jacobt, unsigned int sizeOfZ, double* z)
{
  jacobt[0] = 0.0;
  jacobt[1] = 0.0;
}