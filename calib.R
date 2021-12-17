CALCURVE <- read.csv("intcal20.csv")

Calib <- function(mu, sigma) {
    aux <- sigma^2 + CALCURVE$error^2
    return ( exp( -  ((mu - CALCURVE$C14)^2  / (2 * ( aux )) ) )  /  sqrt(2*pi*aux)  )
}

Spd <- function(ages, errors) {
    m <- mapply(Calib, ages, errors, SIMPLIFY=F)
    return(Reduce("+", m))
}

ASD <- function(ages, errors) {
    df <- cbind(ages, errors)
    s <- sapply(1:nrow(df), function(x) Calib(df[x,1], df[x,2]))
    return(rowSums(s))
}