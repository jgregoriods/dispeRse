#' Simulates first arrival times from one or more origins.
#' The model uses density-dependent growth and emigration.
#' The carrying capacity, growth rates and mobility are
#' allowed to vary with the environment.
#' 
#' The simulation starts with n populated cells at
#' coordinates and start times defined by the parameter
#' coords and runs for a number of time. Each time step
#' corresponds to a generation (defined by parameter t).
#' Growth is applied to every populated cell using a
#' logistic model and emigration to a neighboring cell is
#' calculated from an asymptotic threshold model.
#'
#' The carrying capacity in the density-dependend growth
#' and emigration models is determined by an environment
#' raster, which normally will represent a variable assumed
#' to affect population density (e.g. net primary
#' production) scaled to 0-1 range. Migrants are moved to
#' the cell with the highest environmental value
#' 
#' @import raster
#' @import sp
#' @param environment lorem
#' @param terrain lorem
#' @param r lorem
#' @param phi lorem
#' @param coords lorem
#' @param iter lorem
#' @param t lorem
#' @param dist lorem
#' @param accel lorem
#' @param gamma lorem
#' @return lorem
#' @export
#' @useDynLib dispeRse, .registration = TRUE
run_disp <- function(environment, terrain, r, phi, coords, iter, t, dist, accel, gamma) {

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

    grid_coords <- to_grid(coords, environment)

    x <- grid_coords$x
    y <- grid_coords$y
    start <- grid_coords$date

    ret_val <- .C("run_model", nrow=as.integer(NROW), ncol=as.integer(NCOL),
                population=as.double(population), env=as.double(env_values), arrival=as.integer(arrival),
                r=as.double(r), phi=as.double(phi),
                start=as.integer(start), x=as.integer(x), y=as.integer(y), iter=as.integer(iter),
                num_origins=as.integer(length(x)), t=as.double(t), terrain=as.integer(terr_values),
                accel=as.integer(accel), gamma=as.double(gamma),
                PACKAGE="dispeRse")

    res <- raster(matrix(ret_val$arrival, nrow=NROW, ncol=NCOL, byrow=TRUE))
    res[values(res) == 0] <- NA
    proj4string(res) <- proj4string(environment)
    extent(res) <- extent(environment)
    return(res)
}

#' Lorem ipsum dolor.
#'
#' @import gdistance
#' @import raster
#' @param start_coords lorem
#' @param start_dates lorem
#' @param environment lorem
#' @param base_speed lorem
#' @return lorem
#' @export
sim_dispersal <- function(start_coords, start_dates, environment, base_speed=1) {
    res <- vector(mode="list")
    for (i in 1:nrow(start_coords)) {
        tr <- transition(environment, function(x) 1/mean(x), directions=16)
        tr <- geoCorrection(tr)
        ac <- accCost(tr, as.numeric(start_coords[i,]))
        res[[i]] <- start_dates[i] - ((ac / 1000) * base_speed)
    }
    if (nrow(start_coords) == 1) {
        return(res[[1]])
    } else {
        s <- stack(res)
        max_dates <- max(s, na.rm=T)
        return(max_dates)
    }
}

#' Lorem ipsum dolor.
#'
#' @import raster
#' @param coords lorem
#' @param grid lorem
#' @return lorem
#' @export
to_grid <- function(coords, grid) {
    #coordinates(coords) <- ~x+y
    #proj4string(coords) <- proj4string(grid)
    grid_coords <- data.frame(matrix(ncol=3, nrow=0))
    colnames(grid_coords) <- c("x", "y", "date")
    for (i in 1:nrow(coords)) {
        grid_coords[i,1] <- round((coords$x[i] - xmin(grid)) / res(grid)[1])
        grid_coords[i,2] <- round((ymax(grid) - coords$y[i]) / res(grid)[2])
        grid_coords[i,3] <- coords$date[i]
    }
    return(grid_coords)
}

#' Lorem ipsum dolor.
#'
#' @import raster
#' @param r lorem
#' @return lorem
.scale <- function(r) {
    return ((r - min(values(r), na.rm=T)) / (max(values(r), na.rm=T) - min(values(r), na.rm=T)))
}

#' Lorem ipsum dolor.
#'
#' @import raster
#' @param rain lorem
#' @param temp lorem
#' @return lorem
#' @export
npp <- function(rain, temp) {
    rain_npp <- 3000 * (1 - exp(-0.000664 * rain))
    temp_npp <- 3000 / (1 + exp(1.315-0.119 * temp))
    s <- stack(rain_npp, temp_npp)
    return(min(s))
}
