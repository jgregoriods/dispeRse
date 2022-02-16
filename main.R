source("spatial.R")


data <- read.csv("old/old_dispeRse/data/data.csv")
coordinates(data) <- ~Longitude+Latitude
proj4string(data) <- CRS("+init=epsg:4326")

r <- raster(extent(data))
proj4string(r) <- proj4string(data)
res(r) <- 0.25
values(r) <- 1

simulateDispersal <- function(costRaster, origin, date) {
    tr <- transition(costRaster, function(x) 1 / mean(x), 16)
    tr <- geoCorrection(tr)
    ac <- accCost(tr, origin) / 1000
    ac[values(ac) == Inf] <- NA
    simDates <- date - ac
    return(simDates)
}

s <- simulateDispersal(r, data[323,], 6000)

data$sim_start <- round(extract(s, data))
data$sim_end <- 500

d <- fpc::dbscan(coordinates(data), eps=2.5, MinPts=1)
data$cluster <- d$cluster

# -------------------------------------------------------------


x <- spatial_test(data, model="exp", a=-0.001)