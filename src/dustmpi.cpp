#include "dustmpi.h"

void start_mpi() {
  MPI_Init(NULL, NULL);
}

int get_mpi_size() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

int get_mpi_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

void end_mpi() {
  MPI_Finalize();
}

void distribute_mpi_work(int n_particles, bool use_mpi, int &rank,
                         std::vector<int> &particles_per_node,
                         std::vector<int> &displacements) {
  if (use_mpi) {
    rank = get_mpi_rank();
    int mpi_size = get_mpi_size();

    // Distribute work between nodes - use ceiling of work / available nodes, to
    // avoid a straggler.

    int nodes_left = mpi_size;
    int cumul_work = 0;

    for (int node = 0; node < mpi_size; node++) {
      displacements.push_back(cumul_work);
      int work_for_this_node = (int) ceil(n_particles / (double) nodes_left);
      nodes_left--;
      particles_per_node.push_back(work_for_this_node);
      n_particles -= work_for_this_node;
      cumul_work += work_for_this_node;
    }

  } else {
    rank = 0;
    particles_per_node.push_back(n_particles);
    displacements.push_back(0);
  }
}

void collect_mpi_results(const std::vector<int> &particles_per_node,
                         const std::vector<int> &displacements,
                         const std::vector<double> &node_results,
                         std::vector<double> &all_results) {

  MPI_Gatherv(node_results.data(), node_results.size(), MPI_DOUBLE,
              all_results.data(), particles_per_node.data(),
              displacements.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);

}
