dyn.load("abm.so")

run_model <- function(height, width, base_grid, env_grid, start_x, start_y, num_iter) {
    t0 <- Sys.time()
    ret_val <- .C("run_model", height=as.integer(height),
                               width=as.integer(width),
                               base_grid=as.integer(base_grid),
                               env_grid=as.double(env_grid),
                               start_x=as.integer(start_x),
                               start_y=as.integer(start_y),
                               num_iter=as.integer(num_iter))
    gc()
    print(Sys.time() - t0)
    return(ret_val$base_grid)
}