ncomb <- function(n, y) {
    out <- factorial(n)/
        (factorial(n-y) * factorial(y))
    return(out)
}

beta_sum <- function(abar, bbar) {
    bmode <- (abar - 1)/(abar + bbar -2)
    bmean <- abar/(abar + bbar)
    bvar <- abar*bbar/((abar + bbar)^2 * (abar + bbar + 1))
    out <- c(bmode=bmode, bmean=bmean, bvar=bvar)
    return(out)
}

dbb <- function(m, y, alpha, beta) {
    ## Density of the beta-binomial distribution.
    ## See p.155 in textbook.
    out <- ncomb(m, y) * beta(y + alpha, m-y+beta)/
        beta(alpha, beta)
    return(out)
}

bb_getCI <- function(m, alpha, beta, tail_prob=0.025) {
    ## Equal tail credible interval for the Beta-Bin PPD.
    ds <- sapply(0:m, function(y) {
        dbb(m, y, alpha, beta)
    })
    cum_ds <- cumsum(ds)
    i_begin <- which(cum_ds == cum_ds[cum_ds > tail_prob][1])
    valid_end <- cum_ds[cum_ds < 1 - tail_prob]
    i_end <- which(cum_ds == valid_end[length(valid_end)])
    out <- c(i_begin, i_end)
    return(out)
}

beta_hpd <- function(alpha, beta) {
    f <- function(theta, alpha, beta) {
        ## Difference between the upper and lower
        ## limit densities. The cumulative densities
        ## f(b) - f(theta) are always 0.95
        b <- qbeta(pbeta(theta, alpha, beta) + 0.95, alpha, beta)
        out <- (dbeta(theta, alpha, beta) - dbeta(b, alpha, beta))^2
        return(out)
    }
    ## Finding theta so that f(b) - f(theta) = 0
    hpdmin <- optimize(f, lower=0, upper=qbeta(0.05, alpha, beta),
                       alpha=alpha, beta=beta)$minimum
    hpdmax <- qbeta(p=pbeta(hpdmin,alpha,beta)+0.95,
                    alpha, beta)
    return(c(lo=hpdmin, hi=hpdmax))
}

pb <- function(t0, alpha, beta) {
    f <- function(to, t0, alpha, beta) {
        d1 <- dbeta(t0, alpha, beta)
        d2 <- dbeta(to, alpha, beta)
        out <- (d1 - d2)^2
        return(out)
    }
    bool <- pbeta(t0, alpha, beta) > 0.5
    if(bool) {
        to <- optimize(f, lower=0, upper=t0, t0=t0,
                       alpha=alpha, beta=beta)$minimum
        pb_compl <- pbeta(t0, alpha, beta) -
            pbeta(to, alpha, beta)
    }
    else {
        to <- optimize(f, lower=t0, upper=1, t0=t0,
                       alpha=alpha, beta=beta)$minimum
        pb_compl <- pbeta(to, alpha, beta) -
            pbeta(t0, alpha, beta)
    }
    return(c(t0=t0, to=to, pb=1-pb_compl, pb_compl=pb_compl))
}






