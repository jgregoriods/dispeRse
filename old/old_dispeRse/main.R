source("src.R")

sites <- read.csv("./data/data.csv")
coordinates(sites) <- ~Longitude+Latitude
proj4string(sites) <- wgs

origin <- sites[323,]

model <- dispersal_model(origin, sites, 5000, 500, 1)
res <- test_dispersal(sites, 250, model, 100)
