set.seed(3)
# small circle sample
X <- tdaunif::sample_circle(n = 12L)
# random seed index
l <- landmarks_maxmin(X, seed_index = "random")
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))
# minmax seed index
l <- landmarks_maxmin(X, seed_index = "minmax")
# plot landmark order at point positions
plot(X, asp = 1, pch = NA)
text(X, labels = order(l))

## Iris data set example
iris_pca <- prcomp(iris[,1:4], retx = TRUE, rank. = 2)
lm <- landmarks_maxmin(iris_pca$x, num = 15, cover = TRUE)

## Helper function
draw_circle <- function(center, radius, ...){
  theta <- seq(0, 2 * pi, length = 200)
  lines(x = radius * cos(theta) + center[1], y = radius * sin(theta) + center[2], ...)
}

## Landmark colors
pt_colors <- sample(rainbow(15))

## Plot the points + landmarks
plot(iris_pca$x, asp = 1)
points(iris_pca$x[lm$landmark,], pch = 20, col = pt_colors)

## Draw colored balls around each landmark
for (i in seq(length(lm$landmark))){
  lm_idx <- lm$landmark[i]
  draw_circle(iris_pca$x[lm_idx,], radius = attr(lm, "cover_radius"), col = pt_colors[i])
}

## Draw a segment moving from each point to its landmark
for (i in seq(length(lm$landmark))){
  cover <- lm$cover_set[[i]]
  for (pt_idx in cover){
    pt <- iris_pca$x[pt_idx,]
    landmark_pt <- iris_pca$x[lm$landmark[i],]
    segments(x0 = pt[1], x1 = landmark_pt[1], y0 = pt[2], y1 = landmark_pt[2], col = adjustcolor("gray", alpha.f = 0.20))
  }
}
points(iris_pca$x[lm$landmark,], pch = 20, col = pt_colors)
