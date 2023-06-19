simulate_model_to_file <- function(pars, y, end_time, n_particles,
                                   rng, output_file, use_mpi) {

  start_mpi()
  rank <- get_mpi_rank()
  nodes <- get_mpi_size()

  start_t <- get_mpi_wtime()
  res <- dustmpi:::simulate_model(pars, y, end_time, n_particles, rng, use_mpi)
  end_t <- get_mpi_wtime()

  if (output_file != "") {
    saveRDS(res, paste0(output_file, ".", rank, collapse = ""))
  }

  end_mpi()
}
