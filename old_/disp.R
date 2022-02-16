library(rcarbon)
library(parallel)

START <- 2500
END <- 500

sampleDates <- function(n) {
    dates <- END:START
    probs <- exp(-0.005 * dates)
    sampled <- sample(dates, n, prob=probs, replace=TRUE)
    #sampled <- sample(dates, n, replace=TRUE)
    uncal <- uncalibrate(sampled)
    return(data.frame(C14Age=uncal$rCRA, C14SD=30))
}

simSPD <- function(ndates, nsim) {
    #m <- matrix(nrow=START-END+1, ncol=nsim)
    
    cl <- makeCluster(4)
    clusterEvalQ(cl, library(rcarbon))
    clusterExport(cl, c("sampleDates", "START", "END"))

    res <- parLapply(cl, 1:nsim, function(x) {
        s <- sampleDates(ndates)
        cal <- calibrate(s$C14Age, s$C14SD, verbose=FALSE, normalised=F)
        spd <- spd(cal, s$C14SD, timeRange=c(START, END), verbose=FALSE)
        
        gc()
        
        return(spd$grid$PrDens)
    })

    stopCluster(cl)

    m <- matrix(unlist(res), nrow=START-END+1)

    #for (i in 1:nsim) {
    #    s <- sampleDates(ndates)
    #    cal <- calibrate(s$C14Age, s$C14SD, verbose=FALSE, normalised=F)
    #    spd <- spd(cal, s$C14SD, timeRange=c(START, END), verbose=FALSE)
    #    m[,i] <- spd$grid$PrDens
    #}
    
    f <- foo(m)
    return(list("sim" = m, "stat" = f))
}

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

Plot <- function(real, sim) {
    plot(START:END, sim$stat$hi, col="grey", type="l")
    lines(START:END, sim$stat$lo, col="grey")
    lines(START:END, real, col="red")
}



t0 <- Sys.time()

cal <- calibrate(round(runif(10, 500, 1100)), rep(30, 10), normalised=F)
real <- spd(cal, rep(30, 10), timeRange=c(START, END))
sim <- simSPD(10, 1000)

p <- tst(real$grid$PrDens, sim)

print(Sys.time() - t0)

print(p)

Plot(real$grid$PrDens, sim)

