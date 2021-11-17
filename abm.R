library(raster)
library(viridisLite)

dyn.load("abm.so")

run_model <- function(height, width, base_grid, env_grid, start_x, start_y, num_iter) {
    ret_val <- .C("run_model", height=as.integer(height),
                               width=as.integer(width),
                               base_grid=as.integer(base_grid),
                               env_grid=as.double(env_grid),
                               start_x=as.integer(start_x),
                               start_y=as.integer(start_y),
                               num_iter=as.integer(num_iter))
    gc()
    return(ret_val$base_grid)
}

b <- as.matrix(raster("b.asc"))
env_grid <- matrix(b, nrow=265, byrow=T)
base_grid <- matrix(-1, nrow=340, ncol=265)

res <- run_model(340, 265, base_grid, env_grid, 100, 85, 3000)
r <- raster(matrix(res, nrow=340, ncol=265, byrow=T))
r[values(r) < 0] <- NA
plot(3000 - r, col=plasma(100))
contour(3000 - r, add=T)