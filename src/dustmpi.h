#pragma once

#include <mpi.h>
#include <vector>
#include <cmath>

void start_mpi();
int get_mpi_size();
int get_mpi_rank();

void distribute_mpi_work(int n_particles, bool use_mpi, int &rank,
                         std::vector<int> &particles_per_node,
                         std::vector<int> &displacements);

void collect_mpi_results(const std::vector<int> &particles_per_node,
                         const std::vector<int> &displacements,
                         const std::vector<double> &node_results,
                         std::vector<double> &all_results);
void end_mpi();
