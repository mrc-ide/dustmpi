simulate_model_to_file <- function(pars, y, end_time, n_particles,
                                   rng, output_file, use_mpi) {
  res <- simulate_model(pars, y, end_time, n_particles, rng, use_mpi)

  if ((!is.null(output_file)) && (length(res) > 0)) {
    saveRDS(res, output_file)
  }
  invisible()
}
