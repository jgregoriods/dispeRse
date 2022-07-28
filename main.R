library(raster)
library(viridisLite)
library(dispeRse)

b <- shapefile("layers/sasia.shp")
npp <- raster("layers/sasia_npp_50.tif")

env <- npp / 1250
env[values(env) > 1] <- 1

rnd <- (raster("null_50.tif")-1)
env <- env+rnd

terr <- raster("layers/sasia_slo_50.tif")
terr[terr>=0.9] <- 0
terr[terr>0] <- 1
rv <- raster("rivers_.tif")
rv[rv>2] <- 2
terr <- terr+rv
terr[terr>=3] <- 1

#df <- data.frame(x=c(38,72), y=c(36,25), date=c(10500,7700))
df <- data.frame(x=38, y=36, date=10500)

print("running model")
res <- run_disp(env, terr, 0.02, 0.5, df, 250, 20, 50, 2, 0)

plot(res, col=viridis(255))
contour(res, levels=seq(10000, 0, -1000), add=T)

b <- spTransform(b, proj4string(res))
plot(b, add=T)
