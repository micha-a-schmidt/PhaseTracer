/**
  Example program for PhaseTracer.

  We scan a parameter in a two field model and show the phases and phase transitions
  for each value of the parameter.
*/

#include <iostream>
#include <string>
#include <vector>

#include "SingletDM.hpp"
#include "transition_finder.hpp"
#include "phase_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"

int main() {

  LOGGER(debug);

  // Construct our model
  EffectivePotential::SingletDM model;
  // Comment out the parameter's yang was setting for this model just now
  // As it doens't seem to fit the FS model actually used.
  // Set the input parameters for SingletDM
  // const double mhiggs = 125.2;
  // const double mtop = 173.0;
  // const double VEV = 246.221;
  // double ms = mhiggs*0.5;
  // double lambda_hs = 0.01;
  // double muS = square(ms) - lambda_hs * square(VEV) / 2.;
  // double lambda_s = 2./square(mhiggs*VEV)*square(square(ms)-0.5*lambda_hs*square(VEV))+0.1;

  //PA: use example SKHA values for now
  const double mtop = 173.0;
  double muH2 = -8000;
  double lambda_hs = 0.1;
  double lambda_s =  0.1;
  double mus2 = 1000;
  
  std::vector<double> x(9);
  x[0] = mtop; // Qin = 173.0
  x[1] = mtop; // QEWSB = 173.0
  x[2] = muH2;  // HiggsIN
  x[3] = mus2; // muSInput
  x[4] = lambda_s; // LamSInput
  x[5] = lambda_hs;  // LamSHInput
  
  model.set_input(x);

  // Print masses of CP-even Higgs
  std::cout << "mhh = " << model.get_mh()[0] << std::endl;
  std::cout << "mss = " << model.get_mh()[1] << std::endl;
  
  return 0;


//  // Make PhaseFinder object and find the phases
//  PhaseTracer::PhaseFinder pf(model);
//  pf.set_seed(3);
//  pf.find_phases();
//  std::cout << pf;

//  // Make TransitionFinder object and find the transitions
//  // This takes the phase finder object as an argument. We must have already
//  // populated the phases by e.g. the find_phases method, as above
//  PhaseTracer::TransitionFinder tf(pf);
//  tf.find_transitions();
//  std::cout << tf;

//  // Print the data in a particular format for plotting
//  PhaseTracer::phase_plotter(tf, "SingletDM");
}
