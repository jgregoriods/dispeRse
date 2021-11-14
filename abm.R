dyn.load("abm.so")

run_model <- function(height, width, grid, num_iter) {
    ret_val <- .C("run_model", height=as.integer(height),
                               width=as.integer(width),
                               grid=as.integer(grid),
                               num_iter=as.integer(num_iter))
    return(ret_val$grid)
}

grid <- matrix(0, ncol=50, nrow=25)
res <- run_model(50, 25, grid, 200)
r <- raster(matrix(res, nrow=50, ncol=25, byrow=T))
plot(r)