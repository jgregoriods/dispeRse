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
#x <- test_dispersal(r, sites.m)

to_grid <- function(coords, r) {
    x <- round((coords[1] - extent(r)[1]) / res(r))
    y <- round((extent(r)[4] - coords[2]) / res(r))
    return(c(x, y))
}
