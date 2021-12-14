library(rcarbon)

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
    m <- matrix(nrow=START-END+1, ncol=nsim)
    for (i in 1:nsim) {
        s <- sampleDates(ndates)
        cal <- calibrate(s$C14Age, s$C14SD, verbose=FALSE, normalised=F)
        spd <- spd(cal, s$C14SD, timeRange=c(START, END), verbose=FALSE)
        m[,i] <- spd$grid$PrDens
    }
    return(m)
}

compareSPD <- function(real, sim) {
    sim <- cbind(sim, "PrDens"=real$grid$PrDens)
    a <- apply(sim, 1, function(x) {
            if (x["PrDens"] >= x["lo"] && x["PrDens"] <= x["hi"]) {
                return(x["PrDens"])
            } else {
                return(0)
            }
         })
    plot(real$grid$calBP, sim[,"hi"], type="l", col="gray")
    lines(real$grid$calBP, sim[,"lo"], col="gray")
    lines(real$grid$calBP, real$grid$PrDens, col="red")
    return(a)
}

foo <- function(sim) {
    mean <- apply(sim, 1, mean)
    sd <- apply(sim, 1, sd)
    hi <- apply(sim, 1, function(x) {return(quantile(x, .975))})
    lo <- apply(sim, 1, function(x) {return(quantile(x, .025))})
    return(data.frame(mean=mean, sd=sd, hi=hi, lo=lo))
}

sumZscore <- function(data, stats) {
    z <- lapply(1:length(data), function(i) {
            if (data[i] > stats$hi[i] || data[i] < stats$lo[i]) {
                return(abs((data[i] - stats$mean[i]) / stats$sd[i]))
            } else {
                return(0)
            }
         })
    return(sum(unlist(z)))
}

filterspd <- function(z) {
    res <- sapply(2:2000, function(i) {z[i] && !z[i-1] && !z[i+1]})
}

tst <- function(real, sim) {
    f <- foo(sim)
    m <- apply(sim, 2, sumZscore, stats=f)
    #z <- sumZscore(real$grid$PrDens, f)
    z <- sumZscore(real, f)
    t <- (z - mean(m)) / (sd(m) / (length(m) - 1))
    pval <- pt(t, df=length(m) - 1, lower.tail=F)
    plot(START:END, f$hi, col="grey", type="l")
    lines(START:END, f$lo, col="grey")
    lines(START:END, real, col="red")
    return(pval)
}

cal <- calibrate(c(1000, 1100, 1400, 1000, 1050), rep(30, 5), normalised=F)
real <- spd(cal, rep(30, 5), timeRange=c(START, END))
sim <- simSPD(5, 100)

tst(real$grid$PrDens, sim)