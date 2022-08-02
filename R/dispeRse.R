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
#' variable or combination of variables assumed to affect population density
#' (e.g. net primary production, elevation) scaled to 0-1 range.
#'
#' The dependence of the carrying capacity and growth rate on the environment
#' is allowed to be controlled by a power gamma. By default, the dependence is
#' linear (gamma = 1).
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
#' @param environment A RasterLayer. Environmental values that affect carrying
#' capacity and growth rate. Typically given as a fraction (0-1) of the max K.
#' @param terrain A RasterLayer. Cells with value 1 are barriers and cells with
#' value 2 are corridors.
#' @param coords A DataFrame. Must contain columns x, y, and date with the
#' coordinates and starting date (yr BP) of each origin. Coordinates must be in
#' the same system as the environment and terrain layers.
#' @param years Numeric. Number of years to run the model for.
#' @param r Numeric. The annual growth rate as a decimal.
#' @param phi Numeric. The fission threshold as a fraction of (local) carrying
#' capacity.
#' @param t Numeric. The duration, in years, of a generation (model time step).
#' @param dist Numeric. The distance, in km, that migrants move over a
#' generation.
#' @param accel Numeric. The factor by which the usual distance is increased
#' along corridors. E.g. if dist = 50 km and accel = 3, migrants can move up to
#' 150 km along a corridor. Must range from 2 to 4.
#' @param gamma Numeric. A power that controls the shape of the dependency
#' between r and the environment.
#' @return A RasterLayer with simulated arrival times.
#' @export
#' @useDynLib dispeRse, .registration = TRUE
simulate_dispersal <- function(environment, terrain, coords, years, r=0.025,
                               phi=0.5, t=25, dist=50, accel=3, gamma=1, updates=NULL) {

    print("Preparing rasters...")

    df <- df[order(-df$date),]

    old_proj <- NULL
    old_res <- NULL
    if (!is.projected(CRS(proj4string(environment)))) {
        old_proj <- proj4string(environment)
        old_res <- res(environment)
    }

    ROBINSON <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

    coordinates(coords) <- ~x+y
    proj4string(coords) <- proj4string(environment)

    environment <- projectRaster(environment, res=dist*1000, crs=ROBINSON)
    terrain <- projectRaster(terrain, environment, method="ngb")
    coords <- spTransform(coords, ROBINSON)

    NROW <- nrow(environment)
    NCOL <- ncol(environment)

    population <- rep(0, NROW*NCOL)
    arrival <- rep(0, NROW*NCOL)

    environment[is.na(as.vector(values(environment)))] <- -1

    env_values <- as.vector(values(environment))

    terrain[is.na(values(terrain))] <- -1
    terr_values <- values(terrain)

    grid_coords <- .to_grid(coords, environment)

    x <- grid_coords$x
    y <- grid_coords$y
    start <- grid_coords$date

    num_iter <- ceiling(years / t)
    update_step = -1
    if (!is.null(updates)) update_step <- c(ceiling((start[1] - updates) / t), -1)

    print("Running model...")
    ret_val <- .C("run_model", nrow=as.integer(NROW), ncol=as.integer(NCOL),
                environment=as.double(env_values), terrain=as.integer(terr_values),
                population=as.double(population), arrival=as.integer(arrival),
                x=as.integer(x), y=as.integer(y), start=as.integer(start),
                num_origins=as.integer(length(x)), num_iter=as.integer(num_iter), 
                r=as.double(r), phi=as.double(phi), t=as.double(t),
                accel=as.integer(accel), gamma=as.double(gamma),
                updates=as.integer(update_step), PACKAGE="dispeRse")
    print("Done.")

    res <- raster(matrix(ret_val$arrival, nrow=NROW, ncol=NCOL, byrow=TRUE))
    res[values(res) == 0] <- NA
    proj4string(res) <- proj4string(environment)
    extent(res) <- extent(environment)

    if (!is.null(old_proj)) {
        res <- projectRaster(res, crs=CRS(old_proj), res=old_res)
    }

    return(res)
}

#' Convert from geographic coordinates in a given projection system to the
#' relative position in rows and columns of a grid.
#'
#' @import raster
#' @param coords A SpatialPointsDataFrame. Must contain a column dates with the
#' start date of the dispersal from each point.
#' @param grid A RasterLayer. The coordinates will be converted to this grid.
#' @return A DataFrame with the converted coordinates.
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
