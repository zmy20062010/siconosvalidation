// =============================== Robot arm sample (HuMAnsPa10) ===============================
//
// see modelRobot1.jpg for complete system view.
//
// Keywords: LagrangianDS, LagrangianLinear relation, Moreau TimeStepping, LCP.
//
// =============================================================================================
#include "SiconosKernel.hpp"
#include "RuntimeException.hpp"

#include <math.h>
#include "constants.h"
#include <math.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[])
{
  boost::timer time;
  time.restart();
  try
  {
    // ================= Creation of the model =======================
    // User-defined main parameters
    unsigned int nDof = 2;           // degrees of freedom for robot arm
    double t0 = 0;                   // initial computation time
    double T = 0.3;                  // final computation time
    double h = 1e-4;                 // time step
    // -------------------------
    // --- Dynamical Systems ---
    // -------------------------
    // ---> setting initial conditions
    /* q -->  theta_1,  x    */
    /* v -->  omega_1,  v    */
    SP::SiconosVector q0(new SiconosVector(nDof));
    SP::SiconosVector v0(new SiconosVector(nDof));
    q0->zero();
    v0->zero();

    (*q0)(0) = theta_1_0;
    (*q0)(1) = x_0;

    (*v0)(0) = omega_1_0;
    (*v0)(1) = v_0;

    // ---> designating the dynamical system
    cout << "Loading Dynamical System ..." << endl;
    SP::LagrangianDS dynamicalsystem(new LagrangianDS(q0, v0, "RodImpactPlugin:mass"));
    // external plug-in for the dynamical system in terms of different equation terms
    dynamicalsystem->setComputeNNLFunction("RodImpactPlugin", "NNL");
    dynamicalsystem->setComputeJacobianNNLqDotFunction("RodImpactPlugin", "jacobNNLv");
    dynamicalsystem->setComputeJacobianNNLqFunction("RodImpactPlugin", "jacobNNLq");

    dynamicalsystem->setComputeFIntFunction("RodImpactPlugin", "FInt");
    dynamicalsystem->setComputeJacobianFIntqDotFunction("RodImpactPlugin", "jacobFintv");
    dynamicalsystem->setComputeJacobianFIntqFunction("RodImpactPlugin", "jacobFintq");

    dynamicalsystem->setComputeFExtFunction("RodImpactPlugin", "FExt");

    // ----------------------------------
    // --- Interactions in the System ---
    // ----------------------------------
    // ---> building the interactions
    cout << "Loading Interactions ..." << endl;
    // non smooth laws
    SP::NonSmoothLaw nslaw(new NewtonImpactFrictionNSL(eN, eT, mu, 2));
    // relation of the interaction
    //SP::Relation relation(new LagrangianScleronomousR("BouncingManipulatorPlugin:h0", "BouncingManipulatorPlugin:G0"));
    SP::Relation relation(new LagrangianRheonomousR("RodImpactPlugin:Gap", "RodImpactPlugin:jacobGapq", "RodImpactPlugin:jacobGapt"));
    // interaction creation
    SP::Interaction interaction(new Interaction(2, nslaw, relation, 1));  
    // building the simulation model of the dynamical system
    SP::Model simulationmodel(new Model(t0, T));
    simulationmodel->nonSmoothDynamicalSystem()->insertDynamicalSystem(dynamicalsystem);
    simulationmodel->nonSmoothDynamicalSystem()->link(interaction, dynamicalsystem);

    // ---------------------------------------
    // --- Simulation Process of the System---
    // ---------------------------------------
    // --> setting simulation process
    cout << "Setting Simulation Parameters ..." << endl;
    // time discretisation into time stepping scheme
    SP::TimeDiscretisation t(new TimeDiscretisation(t0, h));
    SP::TimeStepping s(new TimeStepping(t));
    // one step integrators
    SP::OneStepIntegrator vOSI(new MoreauJeanOSI(dynamicalsystem, theta));
    s->insertIntegrator(vOSI);
    // one step non smooth problem
    SP::OneStepNSProblem osnspb(new FrictionContact(2, SICONOS_FRICTION_2D_ENUM));
    s->insertNonSmoothProblem(osnspb, SICONOS_OSNSP_TS_VELOCITY);

    cout << "=== End of model loading === " << endl;
    // =========================== End of model definition ===========================


    // ================================= Computation================================
    // --- Model initialization ---
    cout << "Initialization of the Model ..." << endl;
    simulationmodel->initialize(s);
    cout << "End of model initialization ..." << endl;

    int k = 0;
    int N = ceil((T - t0) / h); // Number of time steps
    boost::progress_display show_progress(N);

    // --- Get the values to be plotted ---
    // -> saved in a matrix dataPlot
    int outputSize = 7;
    SimpleMatrix dataPlot(N + 1, outputSize);
    // For the initial time step:
    SP::SiconosVector q = dynamicalsystem->q();
    SP::SiconosVector v = dynamicalsystem->velocity();

    dataPlot(k, 0) =  simulationmodel->t0();
    dataPlot(k, 1) = (*q)(0);
    dataPlot(k, 2) = (*q)(1);

    dataPlot(k, 3) = (*v)(0);
    dataPlot(k, 4) = (*v)(1);

    dataPlot(k, 5) = (*interaction->y(0))(0);
    dataPlot(k, 6) = (*interaction->lambda(1))(0); 

    cout << "====> Start computation ... " << endl << endl;

    while ((s->hasNextEvent()) && k <= N)
    {
      // get current time step
      k++;
      s->computeOneStep();

      dataPlot(k, 0) =  s->nextTime();
      dataPlot(k, 1) = (*q)(0);
      dataPlot(k, 2) = (*q)(1);

      dataPlot(k, 3) = (*v)(0);
      dataPlot(k, 4) = (*v)(1);

      dataPlot(k, 5) = (*interaction->y(0))(0);
      dataPlot(k, 6) = (*interaction->lambda(1))(0); 

      s->processEvents();
      ++show_progress;
    }

    cout << endl << "End of computation - Number of iterations done: " << k << endl;
    cout << "Computation Time " << time.elapsed()  << endl;
    // --- Output files ---
    cout << "Result File Outputing ..." << endl;
    dataPlot.resize(k, outputSize);
    ioMatrix::write("resultRodImpact.dat", "ascii" , dataPlot, "noDim");
  }
  catch (SiconosException e)
  {
    cout << e.report() << endl;
  }
  catch (...)
  {
    cout << "Exception caught in BouncingManipulator" << endl;
  }
}
