library(gdistance)
library(raster)
library(fpc)
library(rcarbon)
library(parallel)

NSIM <- 10
NCORES <- 3



foo <- function(sim) {
    mean <- apply(sim, 1, mean)
    sd <- apply(sim, 1, sd)
    hi <- apply(sim, 1, function(x) {return(quantile(x, .975))})
    lo <- apply(sim, 1, function(x) {return(quantile(x, .025))})
    return(data.frame(mean=mean, sd=sd, hi=hi, lo=lo))
}

sumZscore <- function(sim, stat) {
    z <- lapply(1:length(sim), function(i) {
            if (sim[i] > stat$hi[i] || sim[i] < stat$lo[i]) {
                return(abs((sim[i] - stat$mean[i]) / stat$sd[i]))
            } else {
                return(0)
            }
         })
    z <- filterspd(unlist(z))
    return(sum(z))
    #return(sum(unlist(z)))
}

filterspd <- function(z) {
    thr <- round(length(z) * 0.05)

    # single points
    zpad <- c(0,z,0)
    s1 <- sapply(3:length(zpad)-1, function(i) {
        return(zpad[i] && !zpad[i-1] && !zpad[i+1])
    })

    # double points
    zpad <- c(0,0,z,0,0)
    s2 <- sapply(5:length(zpad)-2, function(i) {
        return(zpad[i] && ((!zpad[i-1] && zpad[i+1] && !zpad[i+2]) || (!zpad[i-2] && zpad[i-1] && !zpad[i+1])))
    })

    s <- s1 | s2

    if (sum(s) > thr) {
        s[sample(which(s), sum(s) - thr)] <- FALSE
    }

    z[s] <- 0

    return(z)
}

tst <- function(real, sim) {
    m <- apply(sim$sim, 2, sumZscore, stat=sim$stat)
    #print(m)
    #z <- sumZscore(real$grid$PrDens, f)
    z <- sumZscore(real, sim$stat)
    #t <- (z - mean(m)) / (sd(m) / (length(m) - 1))
    #pval <- pt(t, df=length(m) - 1, lower.tail=F)
    p.value <- mean(m >= z)
    return(p.value)
}

Plot <- function(real, sim, start, end) {
    plot(start:end, sim$stat$hi, col="grey", type="l")
    lines(start:end, sim$stat$lo, col="grey")
    lines(start:end, real, col="red")
}

plot.SimResult <- function(s) {
    start <- s$timeRange[1]
    end <- s$timeRange[2]
    plot(start:end, s$sim_spd$stat$hi, col="grey", type="l")
    lines(start:end, s$sim_spd$stat$lo, col="grey")
    lines(start:end, s$real_spd, col="red")
}




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

s <- simulateDispersal(r, data[323,], 4000)

data$sim_start <- round(extract(s, data))
data$sim_end <- 500

d <- fpc::dbscan(coordinates(data), eps=1, MinPts=1)
data$cluster <- d$cluster

sampleDates <- function(nsim, start, end) {
    dates <- end:start
    probs <- exp(-0.001 * dates)
    sampled <- sample(dates, nsim, prob=probs, replace=TRUE)
    uncal <- uncalibrate(sampled)
    return(uncal$rCRA)
}


gc()



spatial_test <- function(spdf) {
    t0 <- Sys.time()

    m <- t(mapply(sampleDates, nsim=NSIM, spdf$sim_start, spdf$sim_end))
    spdf$pval <- NA

    clusters <- unique(spdf$cluster)


    resList <- vector("list", length=length(clusters))

    for (i in clusters) {
        sdata <- spdf[spdf$cluster == i,]
        sm <- m[spdf$cluster == i, , drop=FALSE]

        START <- max(max(sdata$age), max(sdata$sim_start))
        END <- min(min(sdata$age), min(sdata$sim_end))
        
        cal <- calibrate(sdata$age, sdata$sd, normalised=F, verbose=F)
        real_spd <- spd(cal, sdata$sd, timeRange=c(START, END), verbose=F)

        cl <- makeCluster(NCORES)
        clusterEvalQ(cl, library(rcarbon))
        clusterExport(cl, c("START", "END", "sm"), envir=environment())

        res <- parLapply(cl, 1:NSIM, function(x) {
            sim_cal <- calibrate(sm[,x], rep(30, nrow(sm)), normalised=F, verbose=F)
            sim_spd <- spd(sim_cal, rep(30, nrow(sm)), timeRange=c(START, END), verbose=FALSE)

            gc()

            return(sim_spd$grid$PrDens)
        })

        stopCluster(cl)

        mat <- matrix(unlist(res), nrow=START-END+1)
        f <- foo(mat)

        sim_data <- list("sim"=mat, "stat"=f)

        #Plot(real_spd$grid$PrDens, sim_data, START, END)
        p <- tst(real_spd$grid$PrDens, sim_data)

        spdf$pval[spdf$cluster == i] <- p

        final <- list("real_spd"=real_spd$grid$PrDens, "sim_spd"=sim_data, "timeRange"=c(START, END), "p"=p)
        class(final) <- "SimResult"
        resList[[i]] <- final
    }

    print(Sys.time() - t0)

    return(resList)
}