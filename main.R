library(ADMUR)
library(raster)
library(viridisLite)

source("abm.R")

wgs <- CRS("+init=epsg:4326")
albers <- CRS("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

b <- as.matrix(raster("b.asc"))
env_grid <- matrix(b, nrow=265, byrow=T)
base_grid <- matrix(-1, nrow=340, ncol=265)

res <- run_model(340, 265, base_grid, env_grid, 100, 85, 3000)
r <- raster(matrix(res, nrow=340, ncol=265, byrow=T))
r[values(r) < 0] <- NA
#plot(3000 - r, col=plasma(100))
#contour(3000 - r, add=T)

sites <- read.csv("./data/data.csv")
coordinates(sites) <- ~Longitude+Latitude
proj4string(sites) <- wgs
sites.m <- spTransform(sites, albers)

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