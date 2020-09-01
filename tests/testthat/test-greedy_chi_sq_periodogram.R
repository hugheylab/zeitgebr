context("greedy_chi_sq")

test_that("greedy-chi-sq periodogram works", {
  data(dams_sample)

  per <- dams_sample[,
                     greedy_chi_sq_periodogram(activity, sampling_rate = 1/60),
                     by = id]

})


