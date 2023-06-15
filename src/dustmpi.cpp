#include "dustmpi.h"

[[cpp11::register]]
void start_mpi() {
  MPI_Init(NULL, NULL);
}

[[cpp11::register]]
int get_mpi_size() {
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  return size;
}

[[cpp11::register]]
int get_mpi_rank() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return rank;
}

[[cpp11::register]]
void end_mpi() {
  MPI_Finalize();
}

[[cpp11::register]]
double get_mpi_wtime() {
  return MPI_Wtime();
}


std::vector<double> combine_results(std::vector<double> res) {
  std::vector<double> all_res(res.size());
  MPI_Allreduce(res.data(), all_res.data(), res.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  return all_res;
}
