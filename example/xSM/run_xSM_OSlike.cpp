/**
  Z2 real scalar singlet extension of
  the Standard Model 
*/

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "xSM_OSlike.hpp"

#include "phase_finder.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"
#include "phase_plotter.hpp"


int main(int argc, char* argv[]) {
  
  // Set level of screen  output
  LOGGER(info);
  
  // Construct our model
  EffectivePotential::xSM_OSlike model;
      
  // Make PhaseFinder object and find the phases
  PhaseTracer::PhaseFinder pf(model);    
  pf.find_phases();
  std::cout << pf;

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(pf);
  tf.find_transitions();
  std::cout << tf;

  // Print the data in a particular format for plotting
  PhaseTracer::phase_plotter(tf, "xSM");
  return 0;
}
