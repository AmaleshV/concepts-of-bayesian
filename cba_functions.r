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
    mode_b <- (alpha - 1)/(alpha + beta - 2)
    bool <- t0 > mode_b
    ## Check for extreme cases 0 and 1
    if(mode_b == 1) {
        to <- 1
        obj <- NA
        pb_compl <- 1 - pbeta(t0, alpha, beta)
        return(c(t0=t0, to=to, pb=1-pb_compl,
                 pb_compl=pb_compl, obj=obj))
    }
    if(mode_b == 0) {
        to <- 0
        obj <- NA
        pb_compl <- pbeta(t0, alpha, beta)
        return(c(t0=t0, to=to, pb=1-pb_compl,
                 pb_compl=pb_compl, obj=obj))
    }

    ## Compute HPD if no extreme case
    if(bool) {
        to <- optimize(f, lower=0, upper=mode_b, t0=t0,
                       alpha=alpha, beta=beta)$minimum
        obj <- optimize(f, lower=0, upper=mode_b, t0=t0,
                       alpha=alpha, beta=beta)$objective
        pb_compl <- pbeta(t0, alpha, beta) -
            pbeta(to, alpha, beta)
    }
    else {
        to <- optimize(f, lower=mode_b, upper=1, t0=t0,
                       alpha=alpha, beta=beta)$minimum
        obj <- optimize(f, lower=mode_b, upper=1, t0=t0,
                        alpha=alpha, beta=beta)$objective
        pb_compl <- pbeta(to, alpha, beta) -
            pbeta(t0, alpha, beta)
    }
    return(c(t0=t0, to=to, pb=1-pb_compl, pb_compl=pb_compl, obj=obj))
}

pb2 <- function(t0=0.25, alpha, beta)  {
    f <- function(theta) {
        d <- dbeta(0.25, alpha, beta)
        out <- 1/beta(alpha, beta) * theta^(alpha-1) * (1-theta)^(beta-1) - d
        return(out)
    }

    d <- function(theta) {
        dbeta(theta, shape1=alpha, shape2=beta)
    }

    q5 <- qbeta(0.5, alpha, beta)
    bool <- q5 > t0
    if(bool) {
        to <- uniroot(f, c(q5, 1))$root
        pb_compl <- integrate(d, t0, to)$value
    }
    else {
        to <- uniroot(f, c(0, q5))$root
        pb_compl <- integrate(d, to, t0)$value
    }
    
    return(c(t0=t0,
             to=to,
             pb=1-pb_compl,
             pb_compl=pb_compl))
}

contour_plot <- function(t0=0.25, alpha, beta, lo=0, hi=1) {
    dfun <- function(theta) {
        dbeta(theta, shape1=alpha, shape2=beta)
    }
    input <- pb(t0, alpha, beta)
    pb <- input[3]
    pb <- round(pb, 3)
    to <- input[2]
    bool <- to > t0
    if(bool) {
        x <- seq(t0, to, length.out=1e3)
        y <- dbeta(x, alpha, beta)
        x <- c(t0, x, to)
        y <- c(0, y , 0)
    }
    else {
        x <- seq(to, t0, length.out=1e3)
        y <- dbeta(x, alpha, beta)
        x <- c(to, x, t0)
        y <- c(0, y , 0)
    }

    curve(dfun, lo, hi, xlab='Theta', ylab='Beta density',
          main=paste0('Alpha=', alpha, ', Beta=', beta,
                      ', pb=', pb))
    abline(h=dfun(t0), lty=3)
    polygon(x, y, col='grey40')
}



