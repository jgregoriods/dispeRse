r <- 0.02
K <- 1
t <- 30

D <- data.frame(dist=c(1,2,3), p=c(0.2, 0.1, 0.05))

grow_pop <- function(N) {
    return((K * N) / ((K - N) * (exp(-r * t)) + N))
}

get_dist <- function(a, b) {
    return (round(sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)))
}

neighborhood <- function(cell, dist) {
    cells <- list()
    k <- 1
    for (i in -dist:dist) {
        for (j in -dist:dist) {
            xi <- cell[1] + i
            yj <- cell[2] + j
            if (get_dist(cell, c(xi, yj)) == dist) {
                cells[[k]] <- c(xi, yj)
                k <- k + 1
            }
        }
    }
    return (cells)
}

disperse <- function() {
    a <- raster(extent(c(0, 100, 0, 100)))
    res(a) <- 1
    values(a) <- 0

    cells <- data.frame(x=20, y=25)
    a[20, 25] <- 0.5

    for (i in 1:2) {
        for (j in 1:nrow(cells)) {
            
        }
    }
}