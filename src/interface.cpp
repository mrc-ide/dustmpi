// No MPI code, lots of cpp11 code
#include <cpp11.hpp>

#include <dust/r/random.hpp>
#include "implementation.h"

[[cpp11::register]]
cpp11::sexp simulate_model(cpp11::list r_pars,
                           std::vector<double> initial_state,
                           int end_time,
                           int n_particles,
                           cpp11::sexp r_rng_ptr,
                           bool use_mpi) {
  const int freq = cpp11::as_cpp<int>(r_pars["freq"]);
  model::pars pars{cpp11::as_cpp<double>(r_pars["beta"]),
                   cpp11::as_cpp<double>(r_pars["gamma"]),
                   1.0 / static_cast<double>(freq),
                   freq};

  auto rng = dust::random::r::rng_pointer_get<rng_state_type>(r_rng_ptr);

  const auto result = run_simulation(pars, initial_state, end_time, n_particles, rng, use_mpi);

  cpp11::writable::doubles r_result(result.size());
  const auto n_state = initial_state.size();
  std::copy(result.begin(), result.end(), REAL(r_result));
  const auto dim = cpp11::writable::integers{(int)n_state, n_particles, end_time + 1};
  r_result.attr("dim") = dim;
  return r_result;
}
