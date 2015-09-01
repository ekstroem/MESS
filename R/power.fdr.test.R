#' Power calculations for two sample t tests with unequal sample size when controlled by the false discovery rate
#'
#' Compute power of a single test, or determine parameters to obtain target
#' power for equal and unequal sample sizes when the significance is determined by the false discovery rate.
#'
#' @param n Number of observations (per group)
#' @param delta True difference in means
#' @param sd Standard deviation
#' @param fdr False discovery rate
#' @param power Power of test (1 minus Type II error probability)
#' @param pi0 The proportion of tests with with no difference  (true nulls)
#' @param ratio The ratio n2/n1 between the larger group and the smaller group. Should be a value equal to or greater than 1 since n2 is the larger group. Defaults to 1 (equal group sizes)
#' @param type Type of t test. Currently only two-sample tests are implemented
#' @param alternative One- or two-sided test. Currently only two-sided tests are implemented.
#' @param df.method Method for calculating the degrees of default. Possibilities are welch (the default) or classical.
#' @param strict Use strict interpretation in two-sided case
#'
#' @return  Object of class \code{power.htest}, a list of the arguments (including the computed one)
#' augmented with \code{method} and \code{note} elements.
#' @details Exactly one of the parameters \code{n}, \code{delta}, \code{power}, \code{sd}, \code{fdr}, \code{ratio}
#' must be passed as NULL,
#' and that parameter is determined from the others. Notice that the last two have non-NULL defaults
#' so NULL must be explicitly passed if you want to compute them.
#'
#' If \code{strict = TRUE} is used, the power will include the probability
#' of rejection in the opposite direction of the true effect, in the
#' two-sided case. Without this the power will be half the
#' significance level if the true difference is zero.
#' @note \code{uniroot} is used to solve power equation for unknowns, so you may
#' see errors from it, notably about inability to bracket the root
#' when invalid arguments are given.
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{\link{power.t.test}} \code{\link{power.prop.test}}
#' @keywords htest
#' @examples
#' power.fdr.test(n=NULL, delta=2, sd=1, fdr=0.05, power=.8, pi0=0.95)
#' @export
power.fdr.test <-
  function (n = NULL, delta = NULL, sd = 1, fdr = 0.05, power = NULL, pi0 = NULL, sd.ratio=1,
            ratio = 1,
            type = c("two.sample"),
            alternative = c("two.sided"),
            df.method = c("welch", "classical"),
            strict = FALSE)
{
  type <- match.arg(type)
  if (type == "two.sample") {

    if (sum(sapply(list(n, delta, fdr, sd, power, ratio), is.null)) != 1)
      stop("exactly one of n, delta, sd, power, fdr, and ratio must be NULL")

    if (!is.null(ratio) && ratio <= 0)
      stop("ratio between group sizes must be positive")
  }
  else {
      ratio <- 1
      sd.ratio <- 1
      if (sum(sapply(list(n, delta, sd, power, sig.level), is.null)) != 1)
          stop("exactly one of n, delta, sd, power, and fdr must be NULL")
  }

  alternative <- match.arg(alternative)
  df.method <- match.arg(df.method)
  tsample <- switch(type, one.sample = 1, two.sample = 2, paired = 1)
  tside <- switch(alternative, one.sided = 1, two.sided = 2)
  if (tside == 2 && !is.null(delta))
      delta <- abs(delta)
  p.body <- quote({
      ## Compute the df
      nu <- switch(tsample, n-1, switch(df.method,
                                        welch=(sd^2/n + (sd*sd.ratio)^2/(n*ratio))^2/((sd^2/n)^2/(n-1) + ((sd*sd.ratio)^2/(ratio*n))^2/(n*ratio-1)),
                                        classical=(1+ratio)*n-2))
      ## Compute nu (the non-centrality parameter) - theta
      ncp <- switch(tsample, sqrt(n/tsample), sqrt(n/(1+sd.ratio/ratio))) * delta/sd  # theta
      Lambda <- fdr/(1-fdr)*(1-pi0)/pi0

      print(c("NU", nu))
      print(ncp)

      cg <- uniroot(function(x) {
                        print(x)
                        print(2*pt(-x, df=nu)/(pt(x, df=nu, ncp=ncp, lower.tail=FALSE) + pt(-x, df=nu, ncp=ncp))-Lambda)

                        2*pt(-x, df=nu)/(pt(x, df=nu, ncp=ncp, lower.tail=FALSE) + pt(-x, df=nu, ncp=ncp))-Lambda }, c(0, 1e+02))$root
      print(c("cg", cg))
      pt(cg, df=nu, ncp=ncp, lower.tail=FALSE)
      ##    pt(qt(sig.level/tside, nu, lower = FALSE), nu, ncp = switch(tsample, sqrt(n/tsample), sqrt(n/(1+sd.ratio/ratio))) * delta/sd, lower = FALSE)
  })
  if (strict & tside == 2)
      p.body <- quote({
          ## Compute the df
          nu <- switch(tsample, n-1, switch(df.method,
                                            welch=(sd^2/n + (sd*sd.ratio)^2/(n*ratio))^2/((sd^2/n)^2/(n-1) + ((sd*sd.ratio)^2/(ratio*n))^2/(n*ratio-1)),
                                            classical=(1+ratio)*n-2))
          ## Compute nu (the non-centrality parameter) - theta
          ncp <- switch(tsample, sqrt(n/tsample), sqrt(n/(1+sd.ratio/ratio))) * delta/sd  # theta
          Lambda <- fdr/(1-fdr)*(1-pi0)/pi0

          cg <- uniroot(function(x) {2*pt(-x, df=nu)/(1 - pt(x, df=nu, ncp=ncp) + pt(-x, df=nu, ncp=ncp))-Lambda }, c(0, 1e+01))$root
          pt(cg, df=nu, ncp=ncp, lower.tail=FALSE) + pt(-cg, df=nu, ncp=ncp)
      })

#  print(p.body)
  fff <- function(n) {eval(p.body) - power}
  fff(6)
  fff(10)
  print("HEJ")

  if (is.null(power))
    power <- eval(p.body)
  else if (is.null(n))
    n <- uniroot(function(n) eval(p.body) - power, c(4, 1e+07))$root
  else if (is.null(sd))
    sd <- uniroot(function(sd) eval(p.body) - power, delta * c(1e-07, 1e+07))$root
  else if (is.null(delta))
    delta <- uniroot(function(delta) eval(p.body) - power, sd * c(1e-07, 1e+07))$root
  else if (is.null(fdr))
    sig.level <- uniroot(function(fdr) eval(p.body) - power, c(1e-10, 1 - 1e-10))$root
  else if (is.null(ratio))
    ratio <- uniroot(function(ratio) eval(p.body) - power, c(2/n, 1e+07))$root
  else if (is.null(sd.ratio))
    sd.ratio <- uniroot(function(sd.ratio) eval(p.body) - power, c(1e-07, 1e+07))$root
  else stop("internal error")
  NOTE <- switch(type,
                 paired = "n is number of *pairs*, sd is std.dev. of *differences* within pairs",
                 two.sample = ifelse(ratio==1, "n is number in *each* group", "n is vector of number in each group"),
                 NULL)
#  n <- switch(type, paired=n, two.sample=c(n, ifelse(ratio==1, NULL, n*ratio)), one.sample=n)
#  sd <- switch(type, paired=sd, two.sample=c(sd, ifelse(ratio==1, NULL, sd*sd.ratio)), one.sample=sd)

  if (type=="two.sample" & ratio!=1) {
      n <- c(n, n*ratio)
      sd <- c(sd*sd.ratio)
  }


  METHOD <- paste(switch(type, one.sample = "One-sample", two.sample = ifelse(ratio==1, "Two-sample", "Two-sample with unequal sizes"),
                         paired = "Paired"), "t test power calculation when controlling false discovery rate")
  structure(list(n = n, delta = delta, sd = sd, fdr = fdr, pi0 = pi0,
                 power = power, alternative = alternative, note = NOTE,
                 method = METHOD), class = "power.htest")
}


##
##
##

findcrit <- function(fdr=0.05, pi0=NULL, delta=2, sigma=1, n=10) {

    p.body <- quote({
    theta <- abs(delta)/(sigma*sqrt(2/n))
    Lambda <- fdr/(1-fdr)*(1-pi0)/pi0

    2*pt(-x, df=2*n-2)/(1 - pt(x, df=2*n-2, ncp=theta) + pt(-x, df=2*n-2, ncp=theta))-Lambda
    })

    cg <- uniroot(function(x) eval(p.body) , c(0, 5))$root

    pt(cg, df=2*n-2, ncp=2/sqrt(2/n), lower.tail=FALSE) + pt(-cg, df=2*n-2, ncp=2/sqrt(2/n))
}

ppp <- function(x, n, pi0, gamma, lambda=.5, delta) {

    p.body <- quote({
        pi0*x-gamma*(pi0*x + (1-pi0)*(pnorm(sqrt(n*lambda*(1-lambda))*abs(delta) + qnorm(x/2)) + pnorm(-sqrt(n*lambda*(1-lambda))*abs(delta) + qnorm(x/2)))   )
    })

    alphamax <- eval(p.body)

    alphamax <- x
    (1-gamma)*pi0*alphamax/gamma/(1-pi0)
}
