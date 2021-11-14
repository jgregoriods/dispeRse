library(rcarbon)

sampleDates <- function(n, start, end, model='exp') {
    if (model == 'exp') {
        caldates <- round(start + (rexp(n) * (end - start) / 5))
    }
    c14 <- uncalibrate(caldates, rep(30, length(caldates)))
    return (c14)
}

calibDates <- function(x, timeRange) {
    xcal <- calibrate(x$rCRA, x$rError, datenormalised=F, verbose=F)
    xsum <- spd(xcal, timeRange=timeRange, verbose=F)
    return (xsum$grid)
}

sim <- function(n, start, end, model='exp', nsim=10) {
    m <- matrix(nrow=end-start+1, ncol=nsim+1)
    m[,1] <- end:start
    cat("Sampling, calibrating and summing dates...\n")
    for (i in 1:nsim) {
        dates <- sampleDates(n, start, end)
        summed <- calibDates(dates, c(end, start))
        m[,i+1] <- summed$PrDens
        cat("\r", floor(i / nsim * 100), "%")
        gc()
    }
    cat("\n")
    df <- data.frame(m)
    df <- transform(df, upper=apply(df[,2:nsim+1], 1, function(x) quantile(x, 0.975)), lower=apply(df[,2:nsim+1], 1, function(x) quantile(x, 0.025)))
    df <- data.frame(calBP=df$X1, upper=df$upper, lower=df$lower)
    return(df)
}

x <- sim(5, 1200, 4000, nsim=99)
plot(x$calBP, x$upper, type="l")
lines(x$calBP, x$lower)