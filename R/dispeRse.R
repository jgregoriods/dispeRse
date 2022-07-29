#' Simulates first arrival times from one or more origins.
#'
#' The model uses density-dependent growth and emigration. The carrying
#' capacity, growth rates and mobility are allowed to vary with the environment.
#' 
#' The simulation starts with n populated cells at coordinates and start times
#' defined by the parameter coords and runs for a number of time. Each time step
#' corresponds to a generation (defined by parameter t). Growth is applied to
#' every populated cell using a logistic model and emigration to a neighboring
#' cell is calculated from an asymptotic threshold model.
#'
#' The carrying capacity (K) in the density-dependend growth and emigration
#' models is determined by an environment raster, which normally represents a
#' variable assumed to affect population density (e.g. net primary production)
#' scaled to 0-1 range.
#'
#' While carrying capacity depends linearly on the environment, the dependence
#' of the growth rate is allowed to be controlled by a power gamma.
#'
#' For emigration, a threshold phi is considered, expressed as a fraction of
#' carrying capacity. Migrants move to the cell with the highest environmental
#' value in their Moore neighborhood (8 cells).
#'
#' Terrain can be represented by a raster specifying barriers (e.g. mountains),
#' which block movement, and corridors (e.g. rivers), which accelerate movement.
#' For the latter, an acceleration factor can be specified.
#' 
#' @import raster
#' @import sp
#' @param environment A RasterLayer. The local K and r depend on its values,
#' which preferably range from 0 to 1 (fraction of maximum K).
#' @param terrain A RasterLayer. Values must be 0 = no effect, 1 = barrier,
#' 2 = corridor.
#' @param r Numeric. The annual growth rate as a fraction.
#' @param phi Numeric. The fission threshold as a fraction of carrying capacity.
#' @param coords A DataFrame. Must contain columns x, y, and date with the
#' coordinates and starting date (yr BP) of each origin. Coordinates must be in
#' the same system as the environment and terrain layers.
#' @param num_iter Numeric. Number of iterations (steps) of the model.
#' @param t Numeric. The duration, in years, of a generation (model time step).
#' @param dist Numeric. The distance, in km, that migrants move over a
#' generation.
#' @param accel Numeric. The factor by which the usual distance is increased
#' along corridors. E.g. if dist = 50 km and accel = 3, migrants can move up to
#' 150 km along a corridor.
#' @param gamma Numeric. A power that controls the shape of the dependency
#' between r and the environment.
#' @return A RasterLayer with simulated arrival times.
#' @export
#' @useDynLib dispeRse, .registration = TRUE
simulate_dispersal <- function(environment, terrain, coords, num_iter, r=0.025,
                               phi=0.5, t=25, dist=50, accel=3, gamma=1) {
    print("Preparing rasters...")

    ROBINSON <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
    WGS84 <- CRS("+init=epsg:4326")

    coordinates(coords) <- ~x+y
    proj4string(coords) <- proj4string(environment)

    environment <- projectRaster(environment, res=dist*1000, crs=ROBINSON)
    terrain <- projectRaster(terrain, environment, method="ngb")
    coords <- spTransform(coords, ROBINSON)

    NROW <- nrow(environment)
    NCOL <- ncol(environment)

    population <- rep(0, NROW*NCOL)
    arrival <- rep(0, NROW*NCOL)

    environment[is.na(values(environment))] <- -1
    env_values <- values(environment)

    terrain[is.na(values(terrain))] <- -1
    terr_values <- values(terrain)

    grid_coords <- .to_grid(coords, environment)

    x <- grid_coords$x
    y <- grid_coords$y
    start <- grid_coords$date

    print("Running model...")
    ret_val <- .C("run_model", nrow=as.integer(NROW), ncol=as.integer(NCOL),
                population=as.double(population), env=as.double(env_values),
                arrival=as.integer(arrival), r=as.double(r), phi=as.double(phi),
                start=as.integer(start), x=as.integer(x), y=as.integer(y),
                iter=as.integer(num_iter), num_origins=as.integer(length(x)),
                t=as.double(t), terrain=as.integer(terr_values),
                accel=as.integer(accel), gamma=as.double(gamma),
                PACKAGE="dispeRse")
    print("Done.")

    res <- raster(matrix(ret_val$arrival, nrow=NROW, ncol=NCOL, byrow=TRUE))
    res[values(res) == 0] <- NA
    proj4string(res) <- proj4string(environment)
    extent(res) <- extent(environment)
    return(res)
}

#' Lorem ipsum dolor.
#'
#' @import raster
#' @param coords A SpatialPointsDataFrame.
#' @param grid A RasterLayer.
#' @return A DataFrame.
.to_grid <- function(coords, grid) {
    grid_coords <- data.frame(matrix(ncol=3, nrow=0))
    colnames(grid_coords) <- c("x", "y", "date")
    for (i in 1:nrow(coords)) {
        grid_coords[i,1] <- round((coords$x[i] - xmin(grid)) / res(grid)[1])
        grid_coords[i,2] <- round((ymax(grid) - coords$y[i]) / res(grid)[2])
        grid_coords[i,3] <- coords$date[i]
    }
    return(grid_coords)
}
