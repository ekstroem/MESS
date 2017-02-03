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

    print(DF)
    ## Complete case analysis
    ccanalysis <- lm(post ~ treatment + baseline, data=DF)

    miss <- is.na(DF$post)
    R <- !miss           # Observed post indicator

#    Z <- is.na(post)
#    R <- 1-Z
    N <- length(baseline)
    N1 <- sum(treatment)
    N0 <- N - N1

    m1 <- lm(post ~ baseline, subset=(treatment==1))
    m0 <- lm(post ~ baseline, subset=(treatment==0))

#    m1 <- lm(post ~ treatment*baseline)

    eq1 <- predict(m1)
    eq0 <- predict(m0)

    eq1 <- predict(m1, newdata=data.frame(baseline, treatment=rep(1,N)))
    eq0 <- predict(m1, newdata=data.frame(baseline, treatment=rep(0,N)))

    delta <- N1 / N


    pi1 <- 1
    pi0 <- 1
    print(sum(DF$miss[DF$treat]))
    if ( sum(DF$miss[DF$treat])>0 ) {
        print(DF$R)
        missmodel <- glm( R  ~ baseline*treat, family="binomial", subset=(treat==TRUE), data=DF)
        pi1 <- predict(missmodel, newdata=data.frame(x, treat=rep(1, 100)), type="response")
    }
    if ( sum(DF$miss[!DF$treat])>0 ) {
        missmodel <- glm( R  ~ baseline*treat, family="binomial", subset=(treat==FALSE), data=DF)
        pi0 <- predict(missmodel, newdata=data.frame(x, treat=rep(1, 100)), type="response")
    }

    print(pi0)
    print(pi1)

    print("asdasd")

#    pi0 <- predict(missmodel, newdata=data.frame(x, treat=rep(0, 100)), type="response")
#    pi0
#    pi1 <- predict(missmodel, newdata=data.frame(x, treat=rep(1, 100)), type="response")
#    pi1
#    } else if sum(Z) {

#    }


#    cat("kjhkjh")
#    print(pi1)


    post[is.na(post)] <- 0

    mu1 <- 1/N1*( sum((R * treatment * post ) / pi1) -
                  sum((treatment - delta)*eq1) -
                  sum((R - pi1)*treatment*eq1/pi1)
                 )

#    mu1 <- 1/N1*( sum((treatment[R] * post[R] ) / pi1[R]) -
#                  sum((treatment - delta)*eq1) -
#                  sum((R - pi1)*treatment*eq1/pi1)
#                 )

#    print(length(eq1))
#    print(eq1)
#    print(length(treatment-delta))
#    print((treatment-delta)*eq1)

    mu2 <- 1/N0*(sum(R*(1-treatment)*post/pi0) +
                     sum((treatment - delta)*eq0) -
                     sum((R-pi0)*(1-treatment)*eq0/pi0)
                 )
    betahat <- mu1 - mu2

#    betahat <- coef(m1)[2] - coef(m2)[2]

    betahat
}


ppt <- function(x, y, treat) {
    summary(lm(y ~ treat+x))
}
