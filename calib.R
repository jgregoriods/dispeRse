dyn.load('calibc.so')

CALCURVE <- read.csv("intcal20.csv")
calcurve_c14 <- as.integer(CALCURVE$C14)
calcurve_error <- as.integer(CALCURVE$error)
len <- as.integer(nrow(CALCURVE))

Calib <- function(mu, sigma) {
    aux <- sigma^2 + CALCURVE$error^2
    return ( exp( -  ((mu - CALCURVE$C14)^2  / (2 * ( aux )) ) )  /  sqrt(2*pi*aux)  )
}

CalibC <- function(mu, sigma) {
    res <- vector(mode='numeric', length=nrow(CALCURVE))
    x <- .C('calibc', as.integer(CALCURVE$C14), as.integer(CALCURVE$error), as.integer(mu), as.integer(sigma), as.double(res), as.integer(nrow(CALCURVE)))
    return(x[[5]])
}

Spd <- function(ages, errors) {
    #m <- mapply(Calib, ages, errors, SIMPLIFY=F)
    m <- mapply(CalibC, ages, errors, SIMPLIFY=F)
    return(Reduce("+", m))
}

ASD <- function(ages, errors) {
    df <- cbind(ages, errors)
    s <- sapply(1:nrow(df), function(x) Calib(df[x,1], df[x,2]))
    return(rowSums(s))
}