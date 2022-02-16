library(raster)
library(viridisLite)

source("abm.R")

wgs <- CRS("+init=epsg:4326")
albers <- CRS("+proj=aea +lat_0=-32 +lon_0=-60 +lat_1=-5 +lat_2=-42 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

#b <- raster("b.asc")
b <- raster("env2.asc")
b[is.na(values(b))] <- 0
env_grid <- matrix(as.matrix(b), nrow=265, byrow=T)
base_grid <- matrix(-1, nrow=340, ncol=265)

height <- 340
width <- 265
start_x <- 97
start_y <- 36
k <- 1
r <- 0.025
cta <- 0.5
dist <- 2
num_iter <- round(1000 / 30)

res <- run_model(height, width, base_grid, env_grid, start_x, start_y, k, r, cta, dist, num_iter)
res.r <- raster(matrix(res, nrow=height, ncol=width, byrow=T))
res.r[values(res.r) < 0] <- NA

extent(res.r) <- extent(b)
proj4string(res.r) <- albers
res.r <- 5500 - (res.r * 30)

plot(res.r, col=plasma(50))

sites <- read.csv("./data/data.csv")
coordinates(sites) <- ~Longitude+Latitude
proj4string(sites) <- wgs
sites.m <- spTransform(sites, albers)

to_grid <- function(coords, r) {
    x <- round((coords[1] - extent(r)[1]) / res(r))
    y <- round((extent(r)[4] - coords[2]) / res(r))
    return(c(x, y))
}
