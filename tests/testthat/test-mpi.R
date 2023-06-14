test_mpi <- function(script, procs = 2) {
  cmdquote <- if (as.character(Sys.info()['sysname']) == "Windows") '"' else "'"

  suppressWarnings(
    system2("mpiexec", sprintf(" -n %d Rscript -e %s%s%s", procs,
                               cmdquote, script, cmdquote),
            stdout = TRUE, stderr = TRUE))

}

test_script <- function(n_particles, n_streams, out_file,
                        gamma, beta, freq, Si, Ii, Ri, time) {
  file_quote <- if (as.character(Sys.info()['sysname']) == "Windows") "'" else '"'

  gsub("\n", " ", sprintf(
  "dustmpi:::simulate_model_to_file(list(gamma = %s, beta = %s, freq = %s),
      c(%s, %s, %s, 0, 0), %s, %s,
      dust::dust_rng_pointer$new(42, n_streams = %s), %s%s%s, TRUE)",
      gamma, beta, freq, Si, Ii, Ri, time, n_streams, n_particles,
      file_quote, out_file, file_quote))
}

test_that("MPI Sim gives same result as single-processor", {
  pars <- list(gamma = 0.1, beta = 0.2, freq = 4)
  y <- c(1000, 100, 0, 0, 0)
  end_time <- 20
  rng1 <- dust::dust_rng_pointer$new(42, n_streams = 10)
  res1 <- simulate_model(pars, y, end_time, 10, rng1, FALSE)

  tmp <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)
  script <- test_script(10, 10, tmp, pars$gamma, pars$beta, pars$freq,
                        y[1], y[2], y[3], end_time)
  test_mpi(script, 10)
  res2 <- readRDS(tmp)

  expect_identical(dim(res1), dim(res2))
  expect_identical(res1, res2)
})

test_that("MPI Sim is faster than single-processor", {
  pars <- list(gamma = 0.1, beta = 0.2, freq = 4)
  y <- c(1000, 100, 0, 0, 0)
  end_time <- 5000
  rng1 <- dust::dust_rng_pointer$new(42, n_streams = 1000)
  t1 <- system.time(res1 <- simulate_model(pars, y, end_time, 1000, rng1, FALSE))

  script <- test_script(1000, 1000, NULL, pars$gamma, pars$beta, pars$freq, y[1], y[2], y[3], end_time)
  t2 <- system.time(test_mpi(script, 10))

  expect_true(t2[['elapsed']] < t1[['elapsed']])
})

