simulate_model_to_file <- function(pars, y, end_time, n_particles,
                                   rng, output_file, use_mpi) {

  start_mpi()
  rank <- get_mpi_rank()

  res <- dustmpi:::simulate_model(pars, y, end_time, n_particles, rng, use_mpi)
  if (rank == 0) {
    saveRDS(res, output_file)
  }

  end_mpi()
}
