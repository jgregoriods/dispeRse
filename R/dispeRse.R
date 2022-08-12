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
#' carrying capacity. Migrants are distributed among the available cells of the
#' Moore neighborhood proportionally to the inverse of the distance squared.
#'
#' Terrain can be represented by a raster specifying barriers (e.g. mountains),
#' which block movement, and corridors (e.g. rivers), which accelerate movement.
#' For the latter, an acceleration factor can be specified.
#' 
#' @import raster
#' @import sp
#' @param environment A Raster*. Environmental values that affect carrying
#' capacity and growth rate. Typically given as a fraction (0-1) of the max K.
#' If the environment will be updated during the experiment, this parameter
#' must be a raster stack.
#' @param terrain A Raster*. Cells with value 1 are barriers and cells with
#' value 2 are corridors. If the terrain will be updated during the experiment,
#' this parameter must be a raster stack.
#' @param coords A DataFrame. Must contain columns x, y, and date with the
#' coordinates and starting date (yr BP) of each origin. Coordinates must be in
#' the same system as the environment and terrain layers.
#' @param years Numeric. Number of years to run the model for.
#' @param r Numeric. The annual growth rate as a decimal.
#' @param phi Numeric. The emigration threshold as a fraction of carrying
#' capacity.
#' @param t Numeric. The duration, in years, of a generation (model time step).
#' @param dist Numeric. The distance, in km, that migrants move over a
#' generation.
#' @param accel Numeric. The factor by which the usual distance is increased
#' along corridors. E.g. if dist = 50 km and accel = 3, migrants can move up to
#' 150 km along a corridor. Must range from 2 to 4.
#' @param gamma Numeric. A power that controls the shape of the dependency
#' between r and the environment.
#' @param updates Numeric. Optional vector with the years at which the
#' environment and terrain grids will be updated.
#' @param verbose Boolean. If TRUE, write messages to the console.
#' @return A RasterLayer with simulated arrival times.
#' @export
#' @useDynLib dispeRse, .registration = TRUE
#' @examples
#' terr <- raster::stack(replicate(8, euro_terr))
#' sim <- simulate_dispersal(euro_npp, terr, ppnb, 5500, phi=0.33, updates=seq(10000,4000,-1000))
simulate_dispersal <- function(environment, terrain, coords, years, r=0.025,
                               phi=0.5, t=30, dist=50, accel=3, gamma=1,
                               updates=NULL, verbose=TRUE) {

    if (!accel %in% c(2,3,4)) {
        stop("Acceleration factor must be an integer in {2,3,4}.")  
    }

    if (!is.null(updates) && ((nlayers(environment) < length(updates) + 1) ||
                              (nlayers(terrain) < length(updates) + 1))) {
        stop(paste("The number of layers in the environment and terrain",
                   "raster stacks must be at least equal to the number of",
                   "updates + 1."))
    }

    if (dist <= 0) {
        stop("Migration distance must be > 0.")
    }

    if (phi <= 0 || phi >= 1) {
        stop("Emigration threshold must be in the interval (0,1).")
    }

    if (r <= 0 || r >= 1) {
        stop("Annual growth rate must be in the interval (0,1).")
    }

    if (is.na(proj4string(environment)) || is.na(proj4string(terrain))) {
        stop("No projection associated with input layer(s).")
    }

    if (!"x" %in% colnames(coords) || !"y" %in% colnames(coords) || !"date" %in% colnames(coords)) {
        stop("coords must be a DataFrame with at least columns x, y, date.")
    }

    if (verbose) message("Preparing rasters...")

    coords <- coords[order(-coords$date),]

    if (coords$date[1] - years < 0) {
        warning("Simulation will be stopped at 0 bp")
        years <- coords$date[1]
    }

    old_proj <- NULL
    old_res <- NULL
    if (!is.projected(CRS(proj4string(environment)))) {
        old_proj <- proj4string(environment)
        old_res <- res(environment)
    }

    lon_0 <- round(mean(extent(environment)[1:2]))
    lat_0 <- round(mean(extent(environment)[3:4]))
    lim <- abs(extent(environment)[4] - extent(environment)[3]) / 6
    lat_1 <- round(extent(environment)[3] + lim)
    lat_2 <- round(extent(environment)[4] - lim)
    crs <- CRS(paste("+proj=aea +lat_0=", lat_0, " +lon_0=", lon_0, " +lat_1=",
                     lat_1, " +lat_2=", lat_2,
                     " +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs",
                     sep=""))

    coordinates(coords) <- ~x+y
    proj4string(coords) <- proj4string(environment)

    environment <- projectRaster(environment, res=dist*1000, crs=crs)
    terrain <- projectRaster(terrain, environment, method="ngb")
    coords <- spTransform(coords, crs)

    NROW <- nrow(environment)
    NCOL <- ncol(environment)

    population <- rep(0, NROW*NCOL)
    arrival <- rep(0, NROW*NCOL)

    environment[is.na(as.vector(values(environment)))] <- -1

    env_values <- as.vector(values(environment))

    terrain[is.na(values(terrain))] <- -1
    terr_values <- as.vector(values(terrain))

    grid_coords <- .to_grid(coords, environment)

    x <- grid_coords$x
    y <- grid_coords$y
    start <- grid_coords$date

    num_iter <- floor(years / t)
    update_step = -1
    if (!is.null(updates)) update_step <- c(floor((start[1] - updates) / t), -1)

    if (verbose) message("Running model...")
    ret_val <- .C("run_model", nrow=as.integer(NROW), ncol=as.integer(NCOL),
                environment=as.double(env_values), terrain=as.integer(terr_values),
                population=as.double(population), arrival=as.integer(arrival),
                x=as.integer(x), y=as.integer(y), start=as.integer(start),
                num_origins=as.integer(length(x)), num_iter=as.integer(num_iter), 
                r=as.double(r), phi=as.double(phi), t=as.double(t),
                accel=as.integer(accel), gamma=as.double(gamma),
                updates=as.integer(update_step), PACKAGE="dispeRse")
    if (verbose) message("Done.")

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
#' @param grid A Raster*. The coordinates will be converted to this grid.
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

#' Coordinates and earliest dates (med cal BP) for European Neolithic sites.
#' Dataset adapted from the supplementary information in Pinhasi et al. 2005
#' (https://doi.org/10.1371/journal.pbio.0030410).
#' 
#' @format A SpatialPointsDataFrame object.
"euro_dates"

#' Coordinates and earliest dates (med cal BP) for Late Pre-Pottery Neolithic B
#' sites in the Near East. Dataset adapted from the supplementary information
#' in Pinhasi et al. 2005 (https://doi.org/10.1371/journal.pbio.0030410).
#' 
#' @format A DataFrame object.
"ppnb"

#' Transformed Net Primary Production (NPP) from 11 ka to 4 ka at 1000 yr steps.
#' Calculated with the Miami formula using paleoclimatic data from Beyer et al. 2020.
#' (https://doi.org/10.1038/s41597-020-0552-1). Clipped to max=1350, squared and
#' scaled to [0,1].
#' 
#' @format A RasterStack object.
#' 
"euro_npp"

#' Reclassified terrain layer with elevation > 1750 m as barriers and rivers
#' and coastline as corridors. Terrain reclassified from SRTM. Rivers rasterized
#' from GSHHG (https://www.soest.hawaii.edu/pwessel/gshhg/) and coastline
#' rasterized from rnaturalearth.
#' 
#' @format A RasterLayer object.
"euro_terr"