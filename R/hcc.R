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
                agediff = setNames(c(1.6e-06, 1.6e-06, 1.6e-06, 7.8e-06, 1.57e-05, 1.72e-05, 3.92e-05,   # Computed from 2001- onwards
                    7.83e-05, 0.0001332, 0.0001864, 0.0002773, 0.0004982, 0.0007113,
                    0.0011828, 0.0018972, 0.0027981, 0.0046013, 0.0070798, 0.0114273,
                    0.0178694, 0.0279525, 0.0451279, 0.0751312, 0.1178276, 0.1257393,
                    0.1151502, 0.0992046, 0.0802933, 0.0629127, 0.0480795, 0.0367776,
                    0.0277645, 0.021161, 0.0161132, 0.0119725, 0.0090146, 0.0069952,
                    0.0052374, 0.0041564, 0.0031286, 0.002386, 0.0019317, 0.0014868,
                    0.0011703, 0.0009807, 0.0006737, 0.000622, 0.0005013, 0.0003885,
                    0.0002632, 0.0002225, 0.0001661, 0.0001379, 0.0001018, 9.24e-05,
                    7.05e-05, 6.89e-05, 3.45e-05, 2.35e-05, 3.92e-05, 1.41e-05, 1.41e-05,
                    1.1e-05, 4.7e-06, 3.1e-06, 9.4e-06, 3.1e-06, 4.7e-06, 3.1e-06,
                  3.1e-06, 3.1e-06), -23:47),
                FUN.nat.age = hcc.nat.age,
                FUN.aid.age = hcc.aid.age,
                FUN.aid.fam.age = hcc.aid.family.age,
                NaturalChildren = c(.2, .3, .3, .2),
                malefreq = .51,
                familyfreq = 0.14,
                criteria = c("offspring", "families")) {

    ## Age
    ## Pop age difference
    ## CC

    ##
    ## Sanity checks.
    ##
    criteria <- match.arg(criteria)
    ## The geography vectors should have the same length or be scalar
    if (length(S)==1 & max(length(A), length(Q))> 1)
        warning("S is donors for *each* region")
    inputl <- sapply(list(A, Q, S), length)
    if(length(inputl[inputl>1])>1)
        stop("S, A, and Q should either be scalars or have the same length")

    ## NaturalChildren
    if (sum(abs(NaturalChildren) < .Machine$double.eps^2)) {
        stop("NaturalChildren must be a vector of probabilities")
    }
    NaturalChildren <- abs(NaturalChildren)/sum(abs(NaturalChildren))
    n.nat <- sample(seq(0,length(NaturalChildren)-1), 1, prob=NaturalChildren) ## Number of natural children

    ## AID children
    n.aid <- limit
    aid.family <- rbinom(n.aid, size=1, prob=familyfreq)     # Find out which aid children have a sibling. Can only have 1
    n.aid.family <- sum(aid.family)
    if (criteria != "families") {   # Do not make donor families unless criteria is families
        n.aid.family <- 0
        aid.family <- rep(n.aid, 0)
    }

    age.donor <- FUN.aid.age(n.aid)     # Donor age when getting AID children


    ## Glue everything together
    n <- n.aid + n.nat + n.aid.family

    donor <- data.frame(id=seq(n),                             # Individual id within family
                        type=c(rep(1,n.aid), rep(2,n.nat), rep(1,n.aid.family)),   # 1 is AID, 2 is natural
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

    # print(outer(donor$age, donor$age, "-"))

    ## This could be sped up by being computed outside the function
    agespan <- seq(floor(min(as.numeric(names(agediff)))), ceiling(max(as.numeric(names(agediff)))), 1)
#    print(agespan)
                                        #    dr <- hist(agediff, breaks=agespan, plot=FALSE)$density

    dr <- agediff

    dbar <- max(dr)

    ## Check span her
    res <- 2* sum(dr[as.vector(RR) - agespan[1] +1]) * 2 * l * CC * sum(Q * S / A)

    list(probability=res, df=donor, criteria=criteria, regions=max(sapply(list(S, Q, A), length)), DonorsPerRegion=S)
}

hcc.nat.age <- function(n) {
    ## Draw from the marginal distributions - then order
    ## Not ideal
    res <- integer(0)

    if (n>0) {
        res <- rep(31.5, n) + (1:n)*2
    }
    res
}

hcc.aid.age <- function(n) {

    if (n==0)
        stop("Must provide donor age")
    ## Draw first age from the donor information

    probs <- c(0.0125, 0.0273, 0.0496, 0.0570, 0.0718, 0.1014, 0.0903, 0.0792, 0.0755,
               0.0755, 0.0533, 0.0236, 0.0570, 0.0162, 0.0310, 0.0199, 0.0162, 0.0125,
               0.0236, 0.0199, 0.0162, 0.0125, 0.0088, 0.0051, 0.0051, 0.0051, 0.0125,
               0.0088, 0.0019, 0.0014, 0.0014, 0.0014, 0.0051, 0.0014)

    startage <- sample(seq(19, 52, 1), 1, prob=probs)

    ## Then draw the remaining age difference from the following
    ## distribution (based on the donor data)
    c(startage,
      startage + cumsum(rpois(n-1, lambda=exp(-0.65 + rnorm(n-1, sd=0.34)))))
}

## Return the age difference between first donor child and 2nd donor child in a family
hcc.aid.family.age <- function(n) {
    rnbinom(n, size=10, mu=2)+1
}
