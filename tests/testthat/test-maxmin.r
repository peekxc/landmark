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

test_that("maxmin creates a valid cover", {
  lm <- landmarks_maxmin(xy, cover = TRUE)
  all_points <- Reduce(union, lm$cover_set)
  expect_equal(sort(all_points), seq(nrow(xy)))

  lm <- landmarks_maxmin(xy, num = 20, cover = TRUE)
  all_points <- Reduce(union, lm$cover_set)
  expect_equal(sort(all_points), seq(nrow(xy)))

  lm <- landmarks_maxmin(xy, radius = quantile(dist(xy), probs = 0.10), cover = TRUE)
  all_points <- Reduce(union, lm$cover_set)
  expect_equal(sort(all_points), seq(nrow(xy)))
})

test_that("cover elements all have radius less than given threshold", {
  dist_xy <- as.matrix(dist(xy))
  radius <- quantile(dist_xy, probs = 0.10)
  lm <- landmarks_maxmin(xy, radius = radius, cover = TRUE)
  pts_per_cover <- sapply(lm$landmark, function(lm_idx){ sum(dist_xy[lm_idx,] <= radius) })
  expect_equal(pts_per_cover, sapply(lm$cover_set, length))
})

test_that("'dist' objects work equally well", {
  dist_xy <- dist(xy)
  radius <- quantile(dist_xy, probs = 0.10)
  expect_silent(landmarks_maxmin(dist_xy, radius = radius, cover = TRUE))

  ## Dist object version
  lm_dist <- landmarks_maxmin(dist_xy, radius = radius, cover = TRUE)

  ## Point cloud version
  lm_pc <- landmarks_maxmin(xy, radius = radius, cover = TRUE)

  expect_equal(lm_dist$landmark, lm_pc$landmark)
  expect_equal(lm_dist$cover_set, lm_pc$cover_set)
})

test_that("maxmin cover is correct", {
  dist_xy <- as.matrix(dist(xy))
  for (p in c(0.10, 0.15, 0.20, 0.25)){
    radius <- quantile(dist_xy, probs = p)
    lm <- landmarks_maxmin(xy, radius = radius, cover = TRUE)
    true_cover <- lapply(lm$landmark, function(lm_idx){ unname(sort(which(dist_xy[lm_idx,] <= radius))) })
    expect_equal(true_cover, lapply(lm$cover_set, sort))
  }
})

test_that("outputs match input parameters", {
  lm <- landmarks_maxmin(xy, num = 20, cover = TRUE)
  expect_equal(nrow(lm), 20L)
  lm <- landmarks_maxmin(xy, radius = 0.10, frac = TRUE, cover = FALSE)
  expect_true(is.vector(lm))
  lm <- landmarks_maxmin(xy, num = 25L, pick_method = "last", seed_index = 5L)
  expect_true(is.vector(lm))
  expect_equal(length(lm), 25L)
  expect_equal(lm[1], 5)

  expect_silent(landmarks_maxmin(xy, num = 25L, seed_index = "random"))
  lm <- landmarks_maxmin(xy, num = 25L, seed_index = "random")
  expect_equal(length(lm), 25L)
  expect_silent(landmarks_maxmin(xy, num = 25L, seed_index = "minmax"))
  lm <- landmarks_maxmin(xy, num = 25L, seed_index = "minmax")
  expect_equal(length(lm), 25L)
})

