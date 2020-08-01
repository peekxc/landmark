#' @name landmarks_maxmin
#' @title Ball-based Landmark Sets
#' @author Matt Piekenbrock
#' @author Jason Cory Brunson
#' @author Yara Skaf
#' @description Compute landmark sets based on fixed-radius balls.
#' @details This function uses the maxmin procedure to produce a set of evenly
#'   spaced landmark points from a data set. Maxmin is a simple greedy algorithm
#'   that is relatively efficient, but it has a tendency to pick out extremal
#'   points.
#'
#'   One, both, or neither of `num` and `radius` may be passed values. If
#'   neither is specified, then `num` is defaulted to `24L`. To generate a complete
#'   landmark set, use `radius = 0L`.
#' @references De Silva, Vin, and Gunnar E. Carlsson. "Topological estimation
#'   using witness complexes." SPBG 4 (2004): 157-166.
#' @references Dłotko, Paweł. "Ball Mapper: A Shape Summary for Topological Data
#'   Analysis." (2019). Web.
#' @param x a data matrix.
#' @param dist_method a character string specifying the distance metric to use;
#'   passed to `proxy::dist(method)`. Any distance measure in the \code{proxy}
#'   package is supported.
#' @param pick_method a character string specifying the method for selecting one
#'   among indistinguishable points, either `"first"` (the default), `"last"`,
#'   or `"random"`.
#' @param num a positive integer; the desired number of landmark points, or of
#'   sets in a ball cover.
#' @param radius a positive number; the desired radius of each landmark ball, or
#'   of each set in a ball cover.
#' @param frac logical; whether to treat `radius` as a fraction of the diameter
#'   of `x`.
#' @param seed_index an integer (the first landmark to seed the algorithm) or
#'   one of the character strings `"random"` (to select a seed uniformly at
#'   random) and `"minmax"` (to select a seed from the minmax set).
#' @param engine character string specifying the implementation to use; one of
#'   `"original"`, `"C++"`, or `"R"`. When not specified, the R engine is used.
#' @param cover logical; whether to return a data frame of landmark indices and
#'   cover sets (by member index) rather than only a vector of landmark indices.
#' @param extend_num,extend_radius length-two numeric vectors used to extend
#'   landmark parameters for cover set construction. See [extension()].
#' @example inst/examples/ex-landmarks-maxmin.r
#' @rdname landmarks_maxmin
#' @export
landmarks_maxmin <- function(
  x,
  dist_method = "euclidean",
  num = 24L, radius = -1.0, frac = FALSE,
  seed_index = 1L,
  pick_method = c("first", "last", "random"),
  engine = "C++",
  cover = FALSE,
  extend_num = extension(mult = 0, add = 0),
  extend_radius = extension(mult = 0, add = 0)
) {
  ## validate inputs
  stopifnot(is.matrix(x) || is(x, "dist"))

  ## Get number of points in the set
  n_pts <- as.integer(ifelse(is.matrix(x), NROW(x), attr(x, "Size")))
  stopifnot(n_pts > 0)

  ## Choose the tie-breaker
  if (missing(pick_method) || pick_method == "first"){
    pick_method <- 0L
  } else {
    pick_method <- match(pick_method, c("first", "last", "random")) - 1L
  }
  stopifnot(is.numeric(pick_method))

  ## Choose the distance method
  dist_method <- tolower(dist_method)
  stopifnot(dist_method %in%  c("maximum", tolower(proxy::pr_DB$get_entry_names())))
  metric <- c(1,2,3,3,0)[match(dist_method, c("euclidean", "manhattan", "maximum", "supremum"), nomatch = 5)]
  dist_f <- ifelse(metric > 0, identity, function(x,y) { as.numeric(proxy::dist(rbind(x,y), method = dist_method)) })

  # if neither parameter is specified, limit the set to 24 landmarks
  num <- ifelse(missing(num) && missing(radius), min(n_pts, 24L), as.integer(num))

  ## If fraction supplied, compute diameter and apply fraction to it to get radius
  if (frac) {
    stopifnot(!missing(radius) && radius > 0.0)
    ## pairwise distance + convex hull is actually faster to get diameter for point cloud
    diameter <- ifelse(is.matrix(x), max(dist(x[chull(x),,drop=FALSE])), max(x))
    radius <- max(1, radius * diameter)
  }

  # handle seed selection
  if (is.character(seed_index)) {
    seed_index <- switch (
      match.arg(seed_index, c("random", "minmax")),
      random = sample(nrow(x), size = 1L),
      minmax = {
        mm_idx <- minmax(x, dist_method = dist_method)
        mm_idx[[1L]]
      }
    )
  }
  stopifnot(is.numeric(seed_index), seed_index >= 1L, seed_index <= nrow(x))

  ## From here on, it's assumed all the parameters are fine, or handled on the C++ side
  if (is.matrix(x)){
    res <- maxmin_pc(x = x, eps = radius, n = num, dist_f = dist_f, metric = metric, seed = seed_index-1L, pick = pick_method, cover = cover)
  } else {
    res <- maxmin_dist(x = x, n_pts = n_pts, eps = radius, n = num, seed = seed_index-1L, pick = pick_method, cover = cover)
  }
  stopifnot(is.list(res))

  ## If cover requested, return a data.frame, otherwise just return the landmark indices
  if (cover){
    res <- data.frame(landmark = res$landmarks, cover_set = I(res$cover))
  } else {
    res <- as.vector(res$landmarks)
  }

  ## TODO: incorporate this extension
  # print warnings if a parameter was adjusted
  # ext_num <- num * (1 + extend_num[[1L]]) + extend_num[[2L]]
  # if (! is.null(num)) {
  #   if (NROW(res) > ext_num) {
  #     warning(sprintf("Required %d (> num = %d sets of radius %d.", NROW(res), ext_num, radius))
  #   } else if (NROW(res) < ext_num) {
  #     warning(sprintf("Only %d (< num = %d) distinct landmark points were found.", NROW(res), ext_num))
  #   }
  # }

  # return landmarks
  return(res)
}

#' @rdname landmarks_maxmin
#' @export
minmax <- function(x, y = NULL, dist_method = "euclidean") {
  ## Choose the distance method
  dist_method <- tolower(dist_method)
  stopifnot(is(x, "dist") || dist_method %in%  c("maximum", tolower(proxy::pr_DB$get_entry_names())))
  if (dist_method == "maximum"){ dist_method <- "supremum" } ## for compatibility w/ stats::dist

  # use distances from `x` if `y` is not specified
  if (missing(y) || is.null(y)) { y <- x }

  ## If dist object supplied, use lower-triangular indices per index
  if (is(x, "dist")){
    n <- attr(x, "Size")
    to_nat <- function(i,j,n){ ifelse(i < j, n*i - i*(i+1)/2 + j - i - 1, n*j - j*(j+1)/2 + i - j - 1) }
    max_dist <- sapply(seq(n), function(i){
      max(x[to_nat(i-1L, setdiff(seq(n), i)-1L, n)+1L])
    })
  } else {
    max_dist <- sapply(seq(n), function(i){ max(proxy::dist(x[i,,drop=FALSE], xy)) })
  }
  return(which(max_dist == min(max_dist)))
}

#' @rdname landmarks_maxmin
#' @export
maxmin <- function(x, y = NULL, dist_method = "euclidean") {
  ## Choose the distance method
  dist_method <- tolower(dist_method)
  stopifnot(is(x, "dist") || dist_method %in%  c("maximum", tolower(proxy::pr_DB$get_entry_names())))
  if (dist_method == "maximum"){ dist_method <- "supremum" } ## for compatibility w/ stats::dist

  # use distances from `x` if `y` is not specified
  if (missing(y) || is.null(y)) { y <- x }

  ## If dist object supplied, use lower-triangular indices per index
  if (is(x, "dist")){
    n <- attr(x, "Size")
    to_nat <- function(i,j,n){ ifelse(i < j, n*i - i*(i+1)/2 + j - i - 1, n*j - j*(j+1)/2 + i - j - 1) }
    min_dist <- sapply(seq(n), function(i){
      min(x[to_nat(i-1L, setdiff(seq(n), i)-1L, n)+1L])
    })
  } else {
    min_dist <- sapply(seq(n), function(i){ min(proxy::dist(x[i,,drop=FALSE], xy)) })
  }
  return(which(min_dist == max(min_dist)))
}

#' minmax_R <- function(
#'   x, y,
#'   dist_method = "euclidean"
#' ) {
#'
#'   # initialize minimum distance
#'   dist_min <- Inf
#'   # initialize minmax set
#'   mm_idx <- integer(0)
#'
#'   # across all points
#'   for (idx in seq(nrow(x))) {
#'
#'     # maximum distance from point
#'     dist_idx <- max(proxy::dist(
#'       x[idx, , drop = FALSE],
#'       y,
#'       method = dist_method)
#'     )
#'
#'     if (dist_idx == dist_min) {
#'       # if equal to reigning minimum distance, append to minmax set
#'       mm_idx <- c(mm_idx, idx)
#'     } else if (dist_idx < dist_min) {
#'       # if less than reigning minimum distance, replace and reinitialize
#'       dist_min <- dist_idx
#'       mm_idx <- c(idx)
#'     }
#'
#'   }
#'
#'   # return minmax subset
#'   mm_idx
#' }
#'
#' #' @rdname landmarks_maxmin
#' #' @export
#' maxmin <- function(
#'   x, y = NULL,
#'   dist_method = "euclidean"
#' ) {
#'
#'   # use distances from `x` if `y` is not specified
#'   if (is.null(y)) {
#'     y <- x
#'     self <- TRUE
#'   } else {
#'     self <- FALSE
#'   }
#'
#'   # update if/when C++ implementation is available
#'   maxmin_R(x = x, y = y, self = self, dist_method = dist_method)
#' }
#'
#' maxmin_R <- function(
#'   x, y, self,
#'   dist_method = "euclidean"
#' ) {
#'
#'   # initialize minimum distance
#'   dist_max <- 0
#'   # initialize maxmin set
#'   mm_idx <- integer(0)
#'
#'   # across all points
#'   for (idx in seq(nrow(x))) {
#'
#'     # minimum distance from point
#'     dist_idx <- min(proxy::dist(
#'       x[idx, , drop = FALSE],
#'       if (self) x[-idx, , drop = FALSE] else y,
#'       method = dist_method)
#'     )
#'
#'     if (dist_idx == dist_max) {
#'       # if equal to reigning maximum distance, append to maxmin set
#'       mm_idx <- c(mm_idx, idx)
#'     } else if (dist_idx > dist_max) {
#'       # if greater than reigning maximum distance, replace and reinitialize
#'       dist_max <- dist_idx
#'       mm_idx <- c(idx)
#'     }
#'
#'   }
#'
#'   # return minmax subset
#'   mm_idx
#' }
#'
#' #' @rdname landmarks_maxmin
#' #' @export
#' landmarks_maxmin <- function(
#'   x,
#'   dist_method = "euclidean", pick_method = "first",
#'   num = NULL, radius = NULL, frac = FALSE,
#'   seed_index = 1L,
#'   engine = NULL,
#'   cover = FALSE,
#'   extend_num = extension(mult = 0, add = 0),
#'   extend_radius = extension(mult = 0, add = 0)
#' ) {
#'   # validate inputs
#'   stopifnot(is.matrix(x))
#'   dist_method <- tolower(dist_method)
#'   stopifnot(dist_method %in% tolower(proxy::pr_DB$get_entry_names()))
#'   pick_method <- match.arg(pick_method, c("first", "last", "random"))
#'   if (is.null(engine)) engine <- "R"
#'   engine <- match.arg(engine, c("original", "C++", "R"))
#'   if (engine == "C++" && dist_method != "euclidean")
#'     warning("C++ engine is available only for Euclidean distances; ",
#'             "using R engine instead.")
#'
#'   # if neither parameter is specified, limit the set to 24 landmarks
#'   if (is.null(num) && is.null(radius)) {
#'     num <- min(nrow(unique(x)), 24L)
#'   }
#'   # apply `frac` to `radius`
#'   if (frac) {
#'     radius <- max(1, radius * nrow(x))
#'   }
#'   # validate parameters
#'   if (! is.null(num)) {
#'     num <- as.integer(num)
#'     if (is.na(num) || num < 1L || num > nrow(x))
#'       stop("`num` must be a positive integer and at most `nrow(x)`.")
#'   }
#'   if (! is.null(radius)) {
#'     if (is.na(radius) || radius <= 0 || radius == Inf)
#'       stop("`radius` must be a finite non-negative number.")
#'   }
#'
#'   # permute rows of `x` according to `pick_method`
#'   if (pick_method != "first") {
#'     shuffle_idx <- switch (
#'       pick_method,
#'       first = seq(nrow(x)),
#'       last = seq(nrow(x), 1L),
#'       random = sample(nrow(x))
#'     )
#'     x <- x[shuffle_idx, , drop = FALSE]
#'   }
#'   # handle seed selection
#'   if (is.character(seed_index)) {
#'     seed_index <- switch (
#'       match.arg(seed_index, c("random", "minmax")),
#'       random = sample(nrow(x), size = 1L),
#'       minmax = {
#'         mm_idx <- minmax(x,
#'                          dist_method = dist_method)
#'         mm_idx[[1L]]
#'       }
#'     )
#'   } else {
#'     # reset input seed index accordingly
#'     if (pick_method != "first") {
#'       seed_index <- switch (
#'         pick_method,
#'         first = seed_index,
#'         last = nrow(x) + 1L - seed_index,
#'         random = which(shuffle_idx == seed_index)
#'       )
#'     }
#'   }
#'   stopifnot(seed_index >= 1L, seed_index <= nrow(x))
#'
#'   # dispatch to implementations
#'   res <- switch (
#'     engine,
#'     original = landmarks_maxmin_orig(
#'       x = x,
#'       dist_method = dist_method,
#'       num = num, radius = radius,
#'       seed_index = seed_index
#'     ),
#'     `C++` = landmarks_maxmin_cpp(
#'       x = x,
#'       num = if (is.null(num)) 0L else num,
#'       radius = if (is.null(radius)) -1L else radius,
#'       seed_index = seed_index, cover = cover
#'     ),
#'     R = landmarks_maxmin_R(
#'       x = x,
#'       dist_method = dist_method,
#'       num = num, radius = radius,
#'       seed_index = seed_index, cover = cover,
#'       mult_num = extend_num[[1L]], add_num = extend_num[[2L]],
#'       mult_radius = extend_radius[[1L]], add_radius = extend_radius[[2L]]
#'     )
#'   )
#'
#'   # format list as a data frame
#'   stopifnot(is.list(res))
#'   if (length(res) == 1L) {
#'     res <- res[[1]]
#'   } else {
#'     res <- data.frame(landmark = res[[1]], cover_set = I(res[[2]]))
#'   }
#'   # correct for permutation
#'   if (pick_method != "first") {
#'     if (is.list(res)) {
#'       res[[1]] <- shuffle_idx[res[[1]]]
#'       res[[2]] <- lapply(res[[2]], function(set) shuffle_idx[set])
#'     } else {
#'       res <- shuffle_idx[res]
#'     }
#'   }
#'
#'   # print warnings if a parameter was adjusted
#'   ext_num <- num * (1 + extend_num[[1L]]) + extend_num[[2L]]
#'   if (! is.null(num)) {
#'     if (NROW(res) > ext_num) {
#'       warning("Required ", NROW(res),
#'               " (> num = ", ext_num, ") ",
#'               "sets of radius ", radius, ".")
#'     } else if (NROW(res) < ext_num) {
#'       warning("Only ", NROW(res),
#'               " (< num = ", ext_num, ") ",
#'               "distinct landmark points were found.")
#'     }
#'   }
#'
#'   # return landmarks
#'   res
#' }
#'
#' landmarks_maxmin_orig <- function(
#'   x,
#'   dist_method = "euclidean", pick_method = "first",
#'   num = NULL, radius = NULL,
#'   seed_index = 1L
#' ) {
#'   stopifnot(is.matrix(x),
#'             seed_index >= 1L, seed_index <= nrow(x),
#'             # must specify a number of balls or a radius
#'             ! is.null(num) || ! is.null(radius))
#'
#'   shuffle_idx <- switch (
#'     pick_method,
#'     first = NA,
#'     last = rev(seq(nrow(x))),
#'     random = sample(seq(nrow(x)))
#'   )
#'
#'   if (!is.null(num)) {
#'     if (missing(dist_method) || toupper(dist_method) == "EUCLIDEAN") {
#'       lmk_idx <- landmark_maxmin(x, num, seed_index)
#'     } else if (requireNamespace("proxy", quietly = TRUE)) {
#'       stopifnot(toupper(dist_method) %in%
#'                   toupper(proxy::pr_DB$get_entry_names()))
#'       lmk_idx <- vector(mode="integer", num)
#'       lmk_idx[1L] <- seed_index
#'       lmk_min_dist <- rep(Inf, nrow(x))
#'       for (i in 2L:num) {
#'         lmk_dist <- proxy::dist(x, x[lmk_idx[i - 1L], , drop = FALSE],
#'                                      method = dist_method)
#'         lmk_min_dist <- pmin(lmk_dist, lmk_min_dist)
#'         potential_idx <- setdiff(seq(nrow(x)), lmk_idx[c(1:i)])
#'         lmk_idx[i] <- potential_idx[which.max(lmk_min_dist[potential_idx])]
#'       }
#'     } else {
#'       stop(sprintf("Unsupported distance method passed: %s\n", dist_method))
#'     }
#'   } else if(! is.null(radius)) {
#'     stopifnot(toupper(dist_method) %in% toupper(proxy::pr_DB$get_entry_names()))
#'     f_dim <- ncol(x)
#'     f_size <- nrow(x)
#'
#'     # algorithm and variable names as specified in Dlotko paper
#'     C = c() # create a list to store indices of centers/landmarks
#'     next_pt = seed_index # first landmark should be the seed point
#'     while(TRUE){
#'       C = append(C, next_pt) # add new point to list of landmarks
#'
#'       # compute distance between landmark set and each point in the space
#'       dists = proxy::dist(matrix(x[C, ], ncol = f_dim), x, method = dist_method)
#'       sortedDists = matrix(apply(dists, 2, sort), ncol = f_size)
#'
#'       # the next landmark is the point with greatest distance from the current
#'       # landmark set
#'       next_pt = which.max(sortedDists[1, ])
#'       # done when this max distance is < radius, i.e. when all pts are contained
#'       # in an radius-ball
#'       d = sortedDists[1, next_pt]
#'       if(d < radius) { break }
#'     }
#'     lmk_idx = C
#'   }
#'   list(if (is.na(shuffle_idx)) { lmk_idx } else { shuffle_idx[lmk_idx] })
#' }
#'
#' landmarks_maxmin_R <- function(
#'   x,
#'   dist_method = "euclidean", pick_method = "first",
#'   num = NULL, radius = NULL, frac = FALSE,
#'   seed_index = 1L, cover = FALSE,
#'   mult_num = 0, add_num = 0, mult_radius = 0, add_radius = 0
#' ) {
#'
#'   # initialize maxmin, free, and landmark index sets
#'   mm_idx <- seed_index
#'   lmk_idx <- vector(mode = "integer", nrow(x))
#'   free_idx <- seq(nrow(x))
#'   # strike seed and (other) duplicates index from free indices
#'   perm_idx <- c(mm_idx, free_idx[-mm_idx])
#'   free_idx[perm_idx[duplicated(x[perm_idx])]] <- 0L
#'   # initialize distance vector and membership list
#'   lmk_dist <- rep(Inf, times = nrow(x))
#'   # initialize minimum radius and associated number of sets to cover `x`
#'   cover_rad <- Inf
#'   cover_num <- 0L
#'   if (cover) cover_idx <- list()
#'
#'   for (i in seq_along(free_idx)) {
#'
#'     # update vector of landmark points
#'     # -+- assumes data have been permuted according to `pick_method` -+-
#'     lmk_idx[[i]] <- mm_idx[[1L]]
#'
#'     # update vector of free points
#'     if (free_idx[[lmk_idx[[i]]]] == 0L)
#'       stop("A duplicate landmark point was selected, in error.")
#'     free_idx[[lmk_idx[[i]]]] <- 0L
#'
#'     # update landmark distances with distances from new landmark point
#'     lmk_dist <- cbind(
#'       # minimum distances from previous landmark points
#'       lmk_dist,
#'       # distances of all points from newest landmark point
#'       proxy::dist(
#'         x[lmk_idx[[i]], , drop = FALSE],
#'         x,
#'         method = dist_method
#'       )[1, ]
#'     )
#'
#'     # minimum distance from landmarks to `x`
#'     min_dist <- max(pmin(lmk_dist[, 1L], lmk_dist[, 2L]))
#'     # update the minimum radius necessary to cover `x`
#'     if (is.null(radius) && ! is.null(num) && i <= num)
#'       cover_rad <- min_dist
#'     # update membership list
#'     if (cover) {
#'       cover_idx <- if (is.null(radius)) {
#'         # -+- will need to parse later -+-
#'         wh_idx <- which(lmk_dist[, 2L] <=
#'                           cover_rad * (1 + mult_radius) + add_radius)
#'         c(cover_idx,
#'           list(cbind(idx = wh_idx, dist = lmk_dist[wh_idx, 2L])))
#'       } else {
#'         # -+- will not need to parse later -+-
#'         c(cover_idx, list(which(lmk_dist[, 2L] <=
#'                                   radius * (1 + mult_radius) + add_radius)))
#'       }
#'     }
#'
#'     # exhaustion breaks
#'     if (all(free_idx == 0L)) break
#'     # update the minimum number of sets necessary to cover `x`
#'     if (is.null(num) && ! is.null(radius) && min_dist > radius)
#'       cover_num <- i + 1L
#'     # parameter breaks
#'     if ((is.null(num) || i >= num * (1L + mult_num) + add_num) &&
#'         (cover_num == 0L || i >= cover_num * (1L + mult_num) + add_num) &&
#'         (is.null(radius) || min_dist <= radius)) break
#'
#'     # collapse distances to the minimum to each point
#'     lmk_dist <- pmin(lmk_dist[, 1L], lmk_dist[, 2L])
#'     # obtain the maxmin subset
#'     mm_idx <- free_idx[free_idx != 0L]
#'     mm_idx <- mm_idx[lmk_dist[mm_idx] == max(lmk_dist[mm_idx])]
#'
#'   }
#'
#'   # restrict to selected landmarks
#'   lmk_idx <- lmk_idx[seq(i)]
#'   # return data
#'   if (cover) {
#'     if (is.null(radius)) {
#'       # parse extraneous members
#'       cover_idx <- lapply(cover_idx, function(mat) {
#'         unname(mat[mat[, "dist"] <=
#'                      cover_rad * (1 + mult_radius) + add_radius, "idx"])
#'       })
#'     }
#'     # return list of landmark indices and cover membership vectors
#'     list(lmk_idx, cover_idx)
#'   } else {
#'     # return vector of landmark indices
#'     list(lmk_idx)
#'   }
#' }
