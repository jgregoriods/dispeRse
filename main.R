library(ADMUR)
library(raster)
library(viridisLite)

source("abm.R")

wgs <- CRS("+init=epsg:4326")
albers <- CRS("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

b <- raster("b.asc")
env_grid <- matrix(as.matrix(b), nrow=265, byrow=T)
base_grid <- matrix(-1, nrow=340, ncol=265)

res <- run_model(340, 265, base_grid, env_grid, 97, 36, 5000)
r <- raster(matrix(res, nrow=340, ncol=265, byrow=T))
r[values(r) < 0] <- NA
#plot(3000 - r, col=plasma(100))
#contour(3000 - r, add=T)

extent(r) <- extent(b)
proj4string(r) <- albers
r <- 5500 - r

sites <- read.csv("./data/data.csv")
coordinates(sites) <- ~Longitude+Latitude
proj4string(sites) <- wgs
sites.m <- spTransform(sites, albers)
#sites.m$max_age <- extract(r, sites.m)
#sites.m$min_age <- 500

if (FALSE) {

Plot <- function(x) {
    plot(x$timeseries$calBP, x$timeseries[["97.5%"]], type="l")
    lines(x$timeseries$calBP, x$timeseries[["2.5%"]])
    lines(x$timeseries$calBP, x$timeseries$SPD, col="red")
}

simTest <- SPDsimulationTest(data=data.frame(age=c(900,600,650,700), sd=c(40,40,40,40), site="Alpha", datingType="14C"),
                             calcurve=shcal20,
                             calrange=c(500, 1100),
                             pars=0.01,
                             type="exp",
                             N=1000)
simTest$pvalue
Plot(simTest)


}

foo <- function(sites) {
    tryCatch (
        {
            invisible(capture.output(simTest <- SPDsimulationTest(sites@data, shcal20, c(max(sites$max_age), min(sites$min_age)), 0.01, "exp", N=50)))
            return(simTest$pvalue)
        },
        error = function(e) {
            return(0)
        }
    )
}

test_dispersal <- function(model_raster, spdf) {
    spdf$max_age <- extract(model_raster, spdf)
    spdf$min_age <- 500
    spdf$datingType <- "14C"

    spdf$t <- NA
    spdf$pvals <- NA

    names <- unique(spdf$site)
    N <- length(names)

    for (i in 1:N) {
        selected <- spdf[spdf$site == names[i],]

        spdf$pvals[which(spdf$site == names[i])] <- foo(selected)

        cat(paste("\r", i, "of", N))
    }

    return(spdf)
}

#x <- test_dispersal(r, sites.m)

to_grid <- function(coords, r) {
    x <- round((coords[1] - extent(r)[1]) / res(r))
    y <- round((extent(r)[4] - coords[2]) / res(r))
    return(c(x, y))
}

exp_growth <- function(x, r) {
    return ((-r * exp(-r * x)) / (exp(-r*max(x)) - exp(-r * min(x))))
}