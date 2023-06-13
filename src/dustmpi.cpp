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

std::vector<double> combine_results(std::vector<double> res) {
  std::vector<double> all_res(res.size());
  MPI_Allreduce(res.data(), all_res.data(), res.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return all_res;
}
