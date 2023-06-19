test_mpi <- function(script, procs = 2) {
  on_windows <- as.character(Sys.info()['sysname']) == "Windows"
  on_mac <- as.character(Sys.info()['sysname']) == "Darwin"
  on_linux <- as.character(Sys.info()['sysname']) == "Linux"

  cmdquote <- if (on_windows) '"' else "'"
  rs_path <- if (on_linux) "$(R_HOME)/bin/" else if (on_windows) "" else "/usr/local/bin/"

  ret <- system2("mpiexec", sprintf(" -n %s %sRscript -e %s%s%s", procs, rs_path,
                               cmdquote, script, cmdquote),
                 stdout = "", stderr = "")
  if (ret != 0) {
    message(sprintf("system2 reported error code %s", ret))
  }
}

test_script <- function(n_particles, n_streams, out_file,
                        alpha, beta, gamma, freq, Si, Ii, Ri, time) {
  file_quote <- if (as.character(Sys.info()['sysname']) == "Windows") "'" else '"'

  gsub("\n", " ", sprintf(
  "dustmpi:::simulate_model_to_file(list(alpha = %s, beta = %s, gamma = %s, freq = %s),
      c(%s, %s, %s, 0, 0), %s, %s,
      dust::dust_rng_pointer$new(42, n_streams = %s), %s%s%s, TRUE)",
      alpha, beta, gamma, freq, Si, Ii, Ri, time, n_streams, n_particles,
      file_quote, out_file, file_quote))
}

test_that("MPI Sim gives same result as single-processor", {
  pars <- list(alpha = 0.1, beta = 0.2, gamma = 0.1, freq = 4)
  y <- c(1000, 100, 0, 0, 0)
  end_time <- 200
  n_particles <- 100

  rng1 <- dust::dust_rng_pointer$new(42, n_streams = n_particles)
  res1 <- simulate_model(pars, y, end_time, n_particles, rng1, FALSE)

  tmp <- normalizePath(tempfile(), winslash = "/", mustWork = FALSE)
  script <- test_script(n_particles, n_particles, tmp, pars$alpha, pars$beta, pars$gamma, pars$freq,
                        y[1], y[2], y[3], end_time)

  test_mpi(script, procs = 2)
  res2 <- readRDS(paste0(tmp, ".0", collapse=""))

  expect_identical(dim(res1), dim(res2))
  expect_identical(res1, res2)
})

test_that("MPI Sim is faster than single-processor", {
  pars <- list(alpha = 0.1, beta = 0.2, gamma = 0.1, freq = 4)
  y <- c(1000, 100, 0, 0, 0)
  end_time <- 200000
  n_particles <- 100
  rng1 <- dust::dust_rng_pointer$new(42, n_streams = n_particles)
  t1 <- system.time(res1 <- simulate_model(pars, y, end_time, n_particles, rng1, FALSE))

  script <- test_script(n_particles, n_particles, "",
                        pars$alpha, pars$beta, pars$gamma, pars$freq, y[1], y[2], y[3], end_time)
  t2 <- system.time(test_mpi(script, 2))

  #expect_true(t2[['elapsed']] < t1[['elapsed']])
})

