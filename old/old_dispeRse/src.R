library(parallel)
library(raster)
library(rcarbon)
library(rgdal)
library(ADMUR)

wgs <- CRS("+init=epsg:4326")

dispersal_model <- function(origin, sites, start_date, end_date, speed) {
    dists <- spDistsN1(sites, origin, longlat=T)
    max_ages <- round(start_date - (dists / speed))
    min_ages <- end_date
    sites$max_age <- max_ages
    sites$min_age <- min_ages
    return(sites)
}

testDispersal <- function(spdf, nsim=1000) {
    spdf$datingType <- "14C"
    spdf$pvals <- NA
    names <- unique(spdf$site)
    N <- length(names)
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    for (i in 1:N) {
        selected <- spdf[spdf$site == names[i],]

        tryCatch(
            expr = {
                invisible(capture.output(simTest <- SPDsimulationTest(selected@data, shcal20, c(selected$max_age[1], selected$min_age[1]), 0.01, "exp", N=nsim)))
                spdf$pvals[which(spdf$site == names[i])] <- simTest$pvalue
            },
            error = function(e) {
                spdf$pvals[which(spdf$site == names[i])] <- 0
            }
        )
        setTxtProgressBar(pb, i)
    }
    close(pb)
    return(spdf)
}

PlotModel <- function(spdf) {
    sig_05 = spdf[spdf$pvals < 0.05,]
    sig_01 = spdf[spdf$pvals < 0.01,]
    plot(spdf)
    plot(sig_05, col="red", add=T)
    plot(sig_01, col="brown", add=T)
}

test_dispersal <- function(spdf, bandwidth, model, nsim=1000) {
	cluster <- zerodist(spdf, zero=bandwidth, unique.ID=T)
	spdf$cluster <- cluster

    spdf$max_age <- extract(model$max_age, spdf)
    spdf$min_age <- extract(model$min_age, spdf)

    spdf$datingType <- "14C"

    spdf$t <- NA
    spdf$pvals <- NA

    N <- length(unique(cluster))
    pb <- txtProgressBar(min = 0, max = N, style = 3)
    for (i in 1:N) {
        selected <- spdf[spdf$cluster == unique(cluster)[i],]
        invisible(capture.output(simTest <- SPDsimulationTest(selected@data, shcal20, c(max(selected$max_age), min(selected$min_age)), 0.01, "exp", N=nsim)))
        spdf$pvals[which(spdf$cluster == unique(cluster)[i])] <- simTest$pvalue
        setTxtProgressBar(pb, i)
    }
    close(pb)
    test_result <- list("sites"=spdf, "model"=model)
    class(test_result) <- "dispersal_test"
    return(test_result)
}

Plot <- function(x) {
    plot(x$timeseries$calBP, x$timeseries[["97.5%"]], type="l")
    lines(x$timeseries$calBP, x$timeseries[["2.5%"]])
    lines(x$timeseries$calBP, x$timeseries$SPD, col="blue")
}

sample_dates <- function(n, time_range, model='exp') {
    model <- convertPars(0.001, time_range[2]:time_range[1], 'exp')
    cal_dates <- simulateCalendarDates(model, n)
    c14 <- uncalibrate(cal_dates, rep(30, length(cal_dates)))
    return (c14)
}

sum_dates <- function(x, time_range) {
    cal_dates <- calibrate(x$rCRA, x$rError, datenormalised=F, verbose=F)
    summed <- spd(cal_dates, timeRange=c(time_range[1], time_range[2]), verbose=F)
    return (summed$grid)
}

z_score <- function(x, mean, sd, upper, lower) {
    z <- abs((x - mean) / sd)
    selected <- z[x > upper | x < lower]
    selected <- selected[!is.infinite(selected)]
    return(sum(na.omit(selected)))
}

simulate_spd <- function(n, time_range, model='exp', nsim=999, ncores=1) {
    #t1 <- Sys.time()
    m <- matrix(nrow=time_range[1]-time_range[2]+1, ncol=nsim+1)
    m[,1] <- time_range[1]:time_range[2]
    #cat("Sampling, calibrating and summing dates...\n")
    cl <- makeCluster(ncores)
    clusterEvalQ(cl, library(rcarbon))
    clusterExport(cl=cl, varlist=c("sample_dates", "sum_dates"))
    l <- parLapply(cl, 1:nsim, function(x) {
        dates <- sample_dates(n, time_range)
        summed <- sum_dates(dates, time_range)
        return(summed$PrDens)
    })
    stopCluster(cl)
    gc()
    for (i in 1:nsim) {
        m[,i+1] <- l[[i]]
    }
    #for (i in 1:nsim) {
    #    dates <- sample_dates(n, time_range)
    #    summed <- sum_dates(dates, time_range)
    #    m[,i+1] <- summed$PrDens
    #    cat("\r", floor(i / nsim * 100), "%")
    #    gc()
    #}
    #cat("\n")
    df <- data.frame(m)
    df <- transform(df, upper=apply(df[,2:nsim+1], 1, function(x) quantile(x, 0.975)),
                        lower=apply(df[,2:nsim+1], 1, function(x) quantile(x, 0.025)),
                        mean=apply(df[,2:nsim+1], 1, mean),
                        sd=apply(df[,2:nsim+1], 1, sd))
    colnames(df)[1] <- "calBP"
    #print(Sys.time() - t1)
    class(df) <- c("simulated_spd", "data.frame")
    return(df)
}

compare_spds <- function(sim_spd, real_spd, nsim) {
    merged <- merge(sim_spd, real_spd, by="calBP", all=T)
    sim_scores <- c()
    for (i in 1:nsim) {
        sim_scores[i] <- z_score(merged[,i+1], merged$mean, merged$sd, merged$upper, merged$lower)
    }
    real_score <- z_score(merged$PrDens, merged$mean, merged$sd, merged$upper, merged$lower)
    t <- (real_score - mean(sim_scores)) / (sd(sim_scores) / sqrt(nsim))
    pval <- pt(t, df=nsim-1, lower.tail=F)
    return(c(t, pval))
}

plot.dispersal_test <- function(dispersal_test) {
    sig_05 = dispersal_test$sites[dispersal_test$sites$pvals < 0.05,]
    sig_01 = dispersal_test$sites[dispersal_test$sites$pvals < 0.01,]
    plot(dispersal_test$sites)
    plot(sig_05, col="red", add=T)
    plot(sig_01, col="brown", add=T)
}

plot.simulated_spd <- function(sspd, real_spd=NULL) {
    plot(sspd$calBP, sspd$upper, type="l", lwd=0)
    lines(sspd$calBP, sspd$lower, lwd=0)
    polygon(c(sspd$calBP, rev(sspd$calBP)), c(sspd$lower, rev(sspd$upper)), col="ivory2", border="ivory2")
    if (!is.null(real_spd)) {
        lines(real_spd$calBP, real_spd$PrDens)
    }
}