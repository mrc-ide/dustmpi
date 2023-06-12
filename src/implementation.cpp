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
                                   dust::random::prng<rng_state_type> * rng,
                                   const bool use_mpi) {

  auto pars_ptr = std::make_shared<model::pars>(pars);

  // For MPI, work out how much work this node should do.

  int rank = 0;
  std::vector<int> particles_per_node;   // How much work does each node do.
  std::vector<int> displacements;        // First particle each node runs.
  distribute_mpi_work(n_particles, use_mpi, rank,
                      particles_per_node, displacements);

  int node_particles = particles_per_node[rank];

  std::vector<model> models(node_particles, model(pars_ptr));

  const auto n_state = initial_state.size();
  std::vector<double> state(n_state * node_particles);
  std::vector<double> state_next(n_state * node_particles);

  // Storage for all the results.
  std::vector<double> all_results(n_particles * n_state * (end_time + 1));

  // Storage for just the results from this node
  std::vector<double> node_results(node_particles * n_state * (end_time + 1));

  for (int i = 0; i < node_particles; ++i) {
    std::copy(initial_state.begin(), initial_state.end(),
              state.begin() + (i * n_state));
  }
  std::copy(state.begin(), state.end(), node_results.begin());

  for (int i = 0; i < node_particles; ++i) {

    auto& rng_state = rng->state(0);
    double * state_ptr =  state.data() + i * n_state;
    double * state_next_ptr = state_next.data() + i * n_state;
    for (int time = 0; time < end_time; ++time) {
      models[i].update(time,
                       state_ptr,
                       state_next_ptr,
                       rng_state);
      std::swap(state_ptr, state_next_ptr);

      const auto offset = (time + 1) * (n_state * node_particles) + i * n_state;
      std::copy_n(state_ptr, n_state, node_results.begin() + offset);
    }
  }

  if (use_mpi) {
    collect_mpi_results(particles_per_node, displacements,
                        node_results, all_results);
    end_mpi();
  }

  return all_results;
}
