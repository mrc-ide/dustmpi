test_that("Simulation code can run", {
  pars <- list(alpha = 0.1, beta = 0.2, gamma = 0.1, freq = 4)
  y <- c(1000, 100, 0, 0)
  end_time <- 20

  rng1 <- dust::dust_rng_pointer$new(42, n_streams = 1)
  res1 <- simulate_model(r_pars = pars, initial_state = y, end_time = end_time,
                         n_particles = 1, r_rng_ptr = rng1, use_mpi = FALSE)

  rng2 <- dust::dust_rng_pointer$new(42, n_streams = 2)
  res2 <- simulate_model(pars, y, end_time, 2, rng2, FALSE)

  expect_equal(drop(res1), res2[, 1, ])
  expect_true(all(res2[1, , -1] < 1000))
  expect_true(all(apply(res2[1:3, , ], 2:3, sum) == sum(y)))
})
