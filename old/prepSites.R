library(raster)
library(rcarbon)

wgs <- CRS("+init=epsg:4326")
robinson <- CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m")

sites <- read.csv("./data/data.csv")
coordinates(sites) <- ~Longitude+Latitude
proj4string(sites) <- wgs

local_spds <- function(spdf, bandwidth) {
	cluster <- zerodist(spdf, zero=bandwidth, unique.ID=T)
	spdf$cluster <- cluster

    N <- length(unique(cluster))
    spds <- vector(mode="list", length=N)

    for (i in 1:N) {
        selected <- spdf[spdf$cluster == unique(cluster)[i],]
        cal <- calibrate(selected$C14Age, selected$C14SD, datenormalised=F)
        minAge <- min(selected$C14Age) - 500
        maxAge <- max(selected$C14Age) + 500
        spd <- spd(cal, timeRange=c(maxAge, minAge))
        spds[[i]] <- spd$grid
    }

    return(spds)
}
