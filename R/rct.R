#' Pretest-posttest RCT for quantitative observations with
#'
#' Show both the head and tail of an R object XXXX
#'
#' In a typical pretest-posttest RCT, subjects are randomized to two treatments, and response is measured at baseline,
#' prior to intervention with the randomized treatment (pretest), and at prespecified follow-up time (posttest).
#' Interest focuses on the effect of treatments on the change between mean baseline and follow-up response.
#' Missing posttest response for some subjects is routine, and disregarding missing cases can lead to invalid inference.
prepost.test <- function(baseline, post, treatment) {

    ## Check factor
    if ("factor" %in% class(treatment)) {
        treat <- !(treatment==levels(treatment)[1])
    }
    else {
        treat <- !(treatment==0)
    }

    ## Handle missing baseline/treatment

    ##

    DF <- data.frame(baseline, post, treatment, treat)

    ## Complete case analysis
    ccanalysis <- lm(post ~ treatment + baseline, data=DF)

    miss <- is.na(DF$post)
    R <- !miss           # Observed post indicator

#    Z <- is.na(post)
#    R <- 1-Z
    N <- length(baseline)
    N1 <- sum(treatment)
    N0 <- N - N1

    m1 <- lm(post ~ baseline, subset=(treat), data=DF)
    m0 <- lm(post ~ baseline, subset=(!treat), data=DF)

#    m1 <- lm(post ~ treatment*baseline)

#    eq1 <- predict(m1)
#    eq0 <- predict(m0)

    eq1 <- predict(m1, newdata=data.frame(baseline, treat=rep(TRUE, N)))
    eq0 <- predict(m1, newdata=data.frame(baseline, treat=rep(FALSE, N)))
#    eq1 <- 1
#    eq0 <- 1

    delta <- N1 / N


    pi1 <- 1
    pi0 <- 1
    if ( sum(miss[DF$treat])>0 ) {
        missmodel <- glm( R  ~ baseline, family="binomial", subset=(treat==TRUE), data=DF)
        pi1 <- predict(missmodel, newdata=data.frame(DF$baseline, treat=rep(TRUE, nrow(DF))), type="response")
    }
    if ( sum(miss[!DF$treat])>0 ) {
        missmodel <- glm( R  ~ baseline, family="binomial", subset=(treat==FALSE), data=DF)
        pi0 <- predict(missmodel, newdata=data.frame(DF$baseline, treat=rep(FALSE, nrow(DF))), type="response")
    }

    post[is.na(post)] <- 0

    mu1 <- 1/N1*( sum((R * treat * post ) / pi1) -
                  sum((treat - delta)*eq1) -
                  sum((R - pi1)*treat*eq1/pi1)
                 )


    mu2 <- 1/N0*(sum(R*(1-treat)*post/pi0) +
                     sum((treat - delta)*eq0) -
                     sum((R-pi0)*(1-treat)*eq0/pi0)
                 )
    betahat <- mu1 - mu2

#    betahat <- coef(m1)[2] - coef(m2)[2]

#    rval <- list(statistic = tstat, parameter = df, p.value = pval,
#                 conf.int = cint, estimate = betahat, null.value = mu,
#                 alternative = alternative, method = method, data.name = dname)
#    class(rval) <- "htest"
#    return(rval)
    betahat
}


ppt <- function(x, y, treat) {
    summary(lm(y ~ treat+x))
}
