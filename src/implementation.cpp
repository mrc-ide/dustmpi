// lots of cpp11 code, no mpi code
#include "implementation.h"

// Memory layout in the result is:
// [aaaaabbbbbcccccddddd]
// a = particle1, time1
// b = particle2, time1
// c = particle1, time2
// d = particle2, time2
std::vector<double> run_simulation(const model::pars& pars,
                                   const std::vector<double>& initial_state,
                                   int end_time,
                                   const int n_particles,
                                   dust::random::prng<rng_state_type> * rng) {
  auto pars_ptr = std::make_shared<model::pars>(pars);
  std::vector<model> models(n_particles, model(pars_ptr));

  const auto n_state = initial_state.size();

  // // Being generous here, but we only need this as n_state * n_threads
  // // as this is discarded.
  std::vector<double> state(n_state * n_particles);
  std::vector<double> state_next(n_state * n_particles);

  // // This one needs to be the whole shebang though
  std::vector<double> result(n_particles * n_state * (end_time + 1));

  // Can parallelise over i, but don't bother
  for (int i = 0; i < n_particles; ++i) {
    std::copy(initial_state.begin(), initial_state.end(),
              state.begin() + (i * n_state));
  }
  std::copy(state.begin(), state.end(), result.begin());

  // Can parallelise over i.
  for (int i = 0; i < n_particles; ++i) {
    auto& rng_state = rng->state(0);
    double * state_ptr =  state.data() + i * n_state;
    double * state_next_ptr = state_next.data() + i * n_state;
    for (int time = 0; time < end_time; ++time) {
      models[i].update(time,
                       state_ptr,
                       state_next_ptr,
                       rng_state);
      std::swap(state_ptr, state_next_ptr);

      const auto offset = (time + 1) * (n_state * n_particles) + i * n_state;
      std::copy_n(state_ptr, n_state, result.begin() + offset);
    }
  }

  return result;
}
