#pragma once

#include <mpi.h>
#include <vector>

void start_mpi();
int get_mpi_size();
int get_mpi_rank();
void end_mpi();
std::vector<double> combine_results(std::vector<double> results);
