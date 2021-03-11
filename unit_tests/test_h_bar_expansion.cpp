#include "catch/catch.hpp"
#include "models/xSM_MSbar.hpp"
#include "models/SM_parameters.hpp"
#include "h_bar_expansion.hpp"
#include "transition_finder.hpp"
#include "logger.hpp"


TEST_CASE("h-bar expansion method", "[HBarExpansion]") {

  LOGGER(fatal);

  // Construct our model
  const double lambda_hs = 0.3;
  const double Q = 173.03;
  const double xi = 0.;
  const bool tree_level_tadpoles = true;
  const bool use_1L_EWSB_in_0L_mass = false;


  double ms = 0.5 * SM::mh;
  double lambda_s_min = 2. / square(SM::mh * SM::v) *
      square(square(ms) - 0.5 * lambda_hs * square(SM::v));
  double lambda_s = lambda_s_min + 0.1;
  
  // Construct our model
  auto model = EffectivePotential::xSM_MSbar::from_tadpoles(lambda_hs, lambda_s, ms, Q, xi, tree_level_tadpoles);
  model.set_daisy_method(EffectivePotential::DaisyMethod::None);
  model.set_use_1L_EWSB_in_0L_mass(use_1L_EWSB_in_0L_mass);

  // Make HBarExpansion object and find the phases
  PhaseTracer::HbarExpansion hb(model);
  hb.set_seed(0);
  Eigen::ArrayXd pseudo(2);
  pseudo << 0., model.get_v_tree_s();
  hb.add_pseudo_phase(pseudo);
  hb.find_phases();

  // Make TransitionFinder object and find the transitions
  PhaseTracer::TransitionFinder tf(hb);
  tf.set_TC_tol_rel(1e-12);
  tf.find_transitions();

  // Find gamma using high-temperature expansion
  PhaseTracer::HTExpansion ht(model);
  ht.set_seed(0);
  const double TC = tf.get_transitions()[0].TC;
  const auto ht_minima = ht.find_minima_at_t(TC);

  // Use minima with greatest Higgs
  double delta = 0.;
  for (const auto& m : ht_minima) {
    delta = std::max(delta, std::abs(m.x(0)));
  }

  const double gamma = delta / TC;

  CHECK(gamma == Approx(2.059489504));
  CHECK(TC == Approx(90.3332920749));
}
