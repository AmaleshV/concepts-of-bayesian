rm(list=ls())

setwd('/home/marc/School/Concepts_Bayesian_Inference/cba_assigment/')

source("./cba_functions.r")
library(rjags)


dyme.dat <- cbind(dose=c(0, 62.5, 125, 250, 500),
                  nfoetuses=c(282, 225, 290, 261, 141),
                  nmalform=c(67, 43, 193, 250, 141))
dyme.dat <- as.data.frame(dyme.dat)

## Questions 1-3, conjugate analyses and analytical results.
## Question 4-10 based on MCMC. 

## Question 1
## The malformations are treated as successes
prev <- mapply(function(x, y) y/x, dyme.dat$nfoetuse,
               dyme.dat$nmalform)
names(prev) <- dyme.dat$dose

y <- 30
n <- 224

abar_vec <- dyme.dat$nmalform + y
bbar_vec <- dyme.dat$nfoetuses - dyme.dat$nmalform + n - y

## Posteriors using each dose data to obtain
## a prior distribution
post_sum <- mapply(beta_sum, abar_vec, bbar_vec)
colnames(post_sum) <- dyme.dat$dose

## Question 2
## Naive approach
## Using the modes (equivalent of mle) for each dose group
## 25 animals per dose group. 
n <- 25
post_dat <- as.data.frame(t(post_sum))
## Probably best to round these answers
sapply(post_dat$bmode, function(x) x * n)
## PPD approach: beta binomial, 95% interval with equal tails
## Uses bar_vec's from question 1
bb_int <- mapply(function(alpha, beta) {
    bb_getCI(25, alpha, beta)
}, abar_vec, bbar_vec)

colnames(bb_int) <- dyme.dat$dose
rownames(bb_int) <- c('lo', 'hi')

## Question 3
p <- 0.25

contour.dat <- mapply(function(alpha, beta) {
    pb(p, alpha, beta)
}, abar_vec, bbar_vec)
colnames(contour.dat) <- dyme.dat$dose
contour.dat <- round(contour.dat, 3)

## Question 4
logit.dat <- list(y=dyme.dat$nmalform,
                  n=dyme.dat$nfoetuses,
                  sdose=c(scale(dyme.dat$dose)),
                  ## dose=dyme.dat$dose, 
                  g=nrow(dyme.dat))

params=c("alpha", "beta", "logitp", "prevp")

logit.inits <- list(
    list(alpha=-5, beta=-5),
    list(alpha=10, beta=10))

logit.mod1 <- jags.model('./cba_logit.jag', logit.dat,
                         n.chains=2)

## Burnin
update(logit.mod1, 2e3)

logit.sim1 <- coda.samples(logit.mod1, n.iter=1e4,
                           variable.names=params)

logit.csim1 <- as.mcmc(do.call('rbind', logit.sim1))

## plot(logit.sim1, ask=TRUE)
## gelman.diag(logit.sim1)
## gelman.plot(logit.sim1)

summary(logit.sim1)

## Question 5
logit.mod2 <- jags.model('./cba_logit_t.jag', logit.dat,
                         n.chains=2, inits=logit.inits)

## Burnin
update(logit.mod2, 2e3)

logit.sim2 <- coda.samples(logit.mod2, variable.names=params,
                           n.iter=1e4)

logit.csim2 <- as.mcmc(do.call('rbind', logit.sim2))

## plot(logit.sim2)
## gelman.diag(logit.sim2)
## gelman.plot(logit.sim2)

summary(logit.sim2)

## Fitting a frequentist model for comparison
dyme.dat$sdose <- c(scale(dyme.dat$dose))

Y <- cbind(dyme.dat$nmalform, dyme.dat$nfoetuses)

fit1 <- glm(Y ~ sdose, data=dyme.dat, family='binomial')
fit2 <- glm(Y ~ dose, data=dyme.dat, family='binomial')

pi <- sum(dyme.dat$nmalform)/sum(dyme.dat$nfoetuses)
predict(fit1, type='response')


## Question 6
## Graphical and summary statistics of results

## Question 7
keep <- grepl('prevp', colnames(logit.csim1))
prevp.dat <- logit.csim1[, keep]

prevp.dat <- apply(prevp.dat, 2, function(x) {
    quantile(x, probs=c(0.025, .5, 0.975))
})
prevp.dat <- t(prevp.dat)

plot(dyme.dat$dose, prev, type='n')
polygon(c(rev(dyme.dat$dose), dyme.dat$dose),
        c(rev(prevp.dat[, 1]), prevp.dat[, 3]), col='grey80')
lines(dyme.dat$dose, prevp.dat[, 2], type='l', lty=1)
lines(dyme.dat$dose, prevp.dat[, 1], type='l', lty=2)
lines(dyme.dat$dose, prevp.dat[, 3], type='l', lty=2)
points(dyme.dat$dose, prev)

## Question 8
adv_ratio <- sapply(prev[-1], function(x) x/prev[1])

## Question 9
sapply(prev[-1], function(x) {
    bool <- rbinom(1e5, 1e3, x) > rbinom(1e5, 1e3, prev[1])
    sum(bool)/1e5
})

## Question 10
param.dat <- logit.csim1[, c(1, 2)]
qs <- 0.01 * (1 - prev[1]) + prev[1]

bmds <- (log(qs/(1 - qs)) - param.dat[, 1])/param.dat[, 2]
bmd <- bmds * sd(dyme.dat$dose) + mean(dyme.dat$dose)
summary(bmd)




