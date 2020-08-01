context("testing maxmin procedure on larger samples ")

# artificial data set
set.seed(0)
xy <- replicate(2, rnorm(150))

library(landmark)
library(testthat)

test_that("nothing fails unexpectedly", {
  expect_silent(landmarks_maxmin(xy))
  expect_silent(landmarks_maxmin(xy, num = 8))
  expect_silent(landmarks_maxmin(xy, dist_method = "maximum"))
  expect_silent(landmarks_maxmin(xy, dist_method = "manhattan"))
  expect_silent(landmarks_maxmin(xy, radius = 0.10))
  expect_silent(landmarks_maxmin(xy, radius = 0.10, frac = TRUE))
  expect_silent(landmarks_maxmin(xy, seed_index = 9L))
  expect_silent(landmarks_maxmin(xy, pick_method = "random"))
  expect_silent(landmarks_maxmin(xy, pick_method = "last"))
  expect_silent(landmarks_maxmin(dist(xy)))
})


lm <- landmark:::maxmin_pc(x = xy, eps = -1.0, n = 10, dist_f = identity, metric = 1, seed = 0, pick = 0, cover = FALSE)
#landmark:::maxmin_dist(x = dist(s), n_pts = nrow(s), eps = 0.1, n = 0, seed = 0, pick = 0, cover = FALSE)


landmarks_maxmin(xy, radius = 0.10, frac = TRUE)

landmark::landmarks_maxmin(xy, num = 5, seed_index = 8)
