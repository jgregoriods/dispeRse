library(rcarbon)

START <- 2500
END <- 500

sampleDates <- function(n) {
    dates <- END:START
    probs <- exp(-0.005 * dates)
    sampled <- sample(dates, n, prob=probs, replace=TRUE)
    uncal <- uncalibrate(sampled)
    return(data.frame(C14Age=uncal$rCRA, C14SD=30))
}

simSPD <- function(ndates, nsim) {
    m <- matrix(nrow=START-END+1, ncol=nsim)
    for (i in 1:nsim) {
        s <- sampleDates(ndates)
        cal <- calibrate(s$C14Age, s$C14SD, verbose=FALSE)
        spd <- spd(cal, s$C14SD, timeRange=c(START, END), verbose=FALSE)
        m[,i] <- spd$grid$PrDens
    }
    hi <- apply(m, 1, function(x) {return(quantile(x, .975))})
    lo <- apply(m, 1, function(x) {return(quantile(x, .025))})
    m <- cbind(m, hi, lo)
    return(m)
}

compareSPD <- function(real, sim) {
    sim <- cbind(sim, "PrDens"=real$grid$PrDens)
    a <- apply(sim, 1, function(x) {return(x["PrDens"] >= x["lo"] && x["PrDens"] <= x["hi"])})
    plot(real$grid$calBP, sim[,"hi"], type="l", col="gray")
    lines(real$grid$calBP, sim[,"lo"], col="gray")
    lines(real$grid$calBP, real$grid$PrDens, col="red")
    return(a)
}

a <- simSPD(5, 100)
