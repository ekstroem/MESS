#' Compute Hajnal-Curie-Cohen model
#'
#' XXXSplit a matrix into block diagonal sub matrices according to clusters and
#' XXXcombine the lower triangular parts into a vector
#'
#'
#' @param limit the number of AID children for each donor
#' @param S the number of donors
#' @param malefreq the frequency of male births (defaults to 0.51). Set to NULL to force equal male/female distribution among donor offspring
#' @param NaturalChildren A vector of numbers representing the probabilities of that number of natural children (first element corresponds to zero children, second to 1 etc.).
#' @param diag logical. Should the diagonal be included?
#' @return Returns a numeric vector containing the elements of the lower
#' triangular sub matrices.
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @seealso \code{\link{lower.tri}}
#' @keywords manip
#' @examples
#'
#' m <- matrix(1:64, ncol=8)
#' cluster <- c(1, 1, 1, 1, 2, 2, 3, 3)
#' lower.tri.vector(m, cluster)
#'
#' @export hcc
hcc <- function(limit = 12,
                S = 1,
                l = 0.86,
                A = 60000,
                Q = 6/7,
                CC = 4,
                agediff = setNames(rep(1, 41), -20:20),
                FUN.nat.age = hcc.nat.age,
                FUN.aid.age = hcc.aid.age,
                FUN.aid.fam.age = hcc.aid.family.age,
                NaturalChildren = c(.2, .3, .3, .2),
                malefreq = .51,
                familyfreq = 0,
                criteria = c("offspring", "families")) {

    ## Age
    ## Pop age difference
    ## CC


    criteria <- match.arg(criteria)

    ## Sanity check. The geography vectors should have the same length or be scalar

    ## NaturalChildren
    if (sum(abs(NaturalChildren) < .Machine$double.eps^2)) {
        stop("NaturalChildren must be a vector of probabilities")
    }
    NaturalChildren <- abs(NaturalChildren)/sum(abs(NaturalChildren))
    n.nat <- sample(seq(0,length(NaturalChildren)-1), 1, prob=NaturalChildren)

    ## AID children
    n.aid <- limit
    aid.family <- rbinom(n.aid, size=1, prob=familyfreq)
    n.aid.family <- sum(aid.family)

    age.donor <- FUN.aid.age(n.aid)

    ## Glue everything together
    n <- n.aid + n.nat + n.aid.family
    donor <- data.frame(id=seq(n),
                        type=c(rep(1,n.aid), rep(2,n.nat), rep(1,n.aid.family)),
                        sex=rbinom(n, size=1, prob=malefreq),  # 1 male, 0 female
                        age = c(age.donor, FUN.nat.age(n.nat), age.donor[aid.family==1] + FUN.aid.fam.age(n.aid.family)),
                        momid=c(seq(n.aid+n.nat), seq(n.aid)[aid.family==1])
                        )

    ## Find all possible combos
    ## Keep only relevant gender combinations
    keep <- abs(outer(donor$sex, donor$sex, FUN="-"))

    ## Remove nat vs nat
    keep[donor$type==2,donor$type==2] <-  0

    ## Remove donor vs donor from same mother
    keep[outer(donor$momid, donor$momid, FUN="==")] <-  0

    ## Remove the duplicates
    keep[lower.tri(keep)] <- 0

    RR <- outer(donor$age, donor$age, "-")[keep]

    ## This could be sped up by being computed outside the function
    agespan <- seq(floor(min(as.numeric(names(agediff)))), ceiling(max(as.numeric(names(agediff)))), 1)
    dr <- hist(agediff, breaks=agespan, plot=FALSE)$density

    dbar <- max(dr)

    ## Check span her
    res <- 2* sum(dr[as.vector(RR) - agespan[1] +1]) * 2 * l * CC * Q * S/ A

    list(probability=res, df=donor, criteria=criteria)
}

hcc.nat.age <- function(n) {
    ## Draw from the marginal distributions - then order
    ## Not ideal
    rep(33.5, n)
}

hcc.aid.age <- function(n) {
    ## Draw first age from the donor information
    probs <- c(0.00615603614623044, 0.0184558127425089, 0.0396434133634452,
               0.0649260847393735, 0.0863569855166136, 0.100869544200651, 0.10595658647087,
               0.099507320379802, 0.0865748210929907, 0.0707606214357199, 0.0546002111552864,
               0.0430105148606771, 0.0335987463285334, 0.0256644391133457, 0.0212556002279848,
               0.0195926191617051, 0.0200303311379039, 0.0204637575051394, 0.018234457992895,
               0.0137984877765597, 0.00976582770276006, 0.00782930985731277,
               0.00778332018555664, 0.00830227622716194, 0.00732420653449786,
               0.00440308177934223, 0.00174704329404207, 0.000792347200327918,
               0.0010166226945827, 0.00100755023965697, 0.000466022238166756,
               9.69403339694676e-05, 9.06036438747496e-06)

    startage <- sample(seq(18, 55, 1), 1, prob=probs)

    ## Then draw the remaining age difference from the following
    ## distribution (based on the donor data)
    c(startage,
      startage + cumsum(rpois(n-1, lambda=exp(-0.65 + rnorm(n-1, sd=0.34)))))
}

hcc.aid.family.age <- function(n) {
    rnbinom(n, size=10, mu=2)+1
}
