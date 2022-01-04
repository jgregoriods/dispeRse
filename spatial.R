library(gdistance)
library(raster)
library(fpc)
library(rcarbon)
library(parallel)
library(rworldmap)


source("calib.R")


NSIM <- 500
NCORES <- 3

BORDERS <- getMap(resolution="low")

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

s <- simulateDispersal(r, data[323,], 5000)

data$sim_start <- round(extract(s, data))
data$sim_end <- 500

d <- fpc::dbscan(coordinates(data), eps=2.5, MinPts=1)
data$cluster <- d$cluster

sampleDates <- function(nsim, start, end) {
    dates <- end:start
    probs <- exp(-0.001 * dates)
    sampled <- sample(dates, nsim, prob=probs, replace=TRUE)
    uncal <- uncalibrate(sampled)
    return(uncal$rCRA)
}


gc()



spatial_test <- function(mySpdf) {
    t0 <- Sys.time()

    myMat <- t(mapply(sampleDates, nsim=NSIM, mySpdf$sim_start, mySpdf$sim_end))
    mySpdf$pval <- NA

    clusters <- unique(mySpdf$cluster)


    #resList <- vector("list", length=length(clusters))

    #pb <- txtProgressBar(0, length(clusters), 0, style=3)

    cl <- makeCluster(NCORES, outfile="")
    clusterEvalQ(cl, library(sp))
    clusterEvalQ(cl, dyn.load('calibc.so'))
    clusterExport(cl, c("Spd", "myMat", "mySpdf", "Calib", "CalibC", "CALCURVE", "NSIM", "foo", "tst", "sumZscore", "filterspd"), envir=environment())


    res <- parLapply(cl, clusters, function(i) {
        sdata <- mySpdf[mySpdf$cluster == i,]
        sm <- myMat[mySpdf$cluster == i, , drop=FALSE]


        START <- max(max(sdata$age), max(sdata$sim_start))
        END <- min(min(sdata$age), min(sdata$sim_end))

        #cal <- calibrate(sdata$age, sdata$sd, normalised=F, verbose=F)
        #real_spd <- spd(cal, sdata$sd, timeRange=c(START, END), verbose=F)

        real_spd <- Spd(sdata$age, sdata$sd)

        
        #clusterEvalQ(cl, library(rcarbon))
        

        spds <- lapply(1:NSIM, function(x) {
            #sim_cal <- calibrate(sm[,x], rep(30, nrow(sm)), normalised=F, verbose=F)
            #sim_spd <- spd(sim_cal, rep(30, nrow(sm)), timeRange=c(START, END), verbose=FALSE)
            sim_spd <- Spd(sm[,x], rep(30, nrow(sm)))

            gc()

            #return(sim_spd$grid$PrDens)
            return(sim_spd)
        })

        #mat <- matrix(unlist(res), nrow=START-END+1)

        mat <- matrix(unlist(spds), nrow=length(real_spd))
        mat <- mat[which(CALCURVE$cal < START & CALCURVE$cal > END),]
        
        f <- foo(mat)

        sim_data <- list("sim"=mat, "stat"=f)

        #p <- tst(real_spd$grid$PrDens, sim_data)
        real_spd <- real_spd[which(CALCURVE$cal < START & CALCURVE$cal > END)]
        p <- tst(real_spd, sim_data)

        #mySpdf$pval[mySpdf$cluster == i] <- p

        #df <- data.frame("calBP"=real_spd$grid$calBP, "PrDens"=real_spd$grid$PrDens, "lo"=f$lo, "hi"=f$hi)
        df <- data.frame("calBP"=CALCURVE[which(CALCURVE$cal < START & CALCURVE$cal > END),1], "PrDens"=real_spd, "lo"=f$lo, "hi"=f$hi)

        spdmodeltest <- list("result"=df, "pval"=p, "siteCluster"=i)
        class(spdmodeltest) <- c(class(spdmodeltest), "SpdModelTest")

        #resList[[i]] <- spdmodeltest
        return(spdmodeltest)
    })

    stopCluster(cl)

    mySpdf$pval <- NA
    for (i in res) {
        mySpdf$pval[mySpdf$cluster == i$siteCluster] <- i$pval
    }

    #new_spdf <- SpatialPointsDataFrame(coordinates(mySpdf), data.frame("pval"=mySpdf$pval, "cluster"=mySpdf$cluster))
    #proj4string(new_spdf) <- proj4string(mySpdf)

    print(Sys.time() - t0)

    #return(list(new_spdf, resList))
    return(list(mySpdf, res))
}

Plot <- function(x) {
    xlim <- c(extent(x)[1]-1, extent(x)[2]+1)
    ylim <- c(extent(x)[3]-1, extent(x)[4]+1)
    plot(BORDERS, col="lightgrey", border="white", xlim=xlim, ylim=ylim)
    cls <- s[[1]][which(!duplicated(s[[1]]$cluster)),]
    plot(x, pch=20, add=T)
    sig <- x[x$pval < 0.05,]
    plot(sig, pch=20, col="red", add=T)
    text(coordinates(cls)[,1]+.5, coordinates(cls)[,2]+.5, cls$cluster, cex=.8)
}


Plt <- function(x) {
    par(mfrow=c(ceiling(length(x) / 2),2), mar=rep(2, 4))
    for (i in 1:length(x)) {
        plot(x[[i]], main=i)
    }
}

