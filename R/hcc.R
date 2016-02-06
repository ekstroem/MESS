#' Compute Hajnal-Curie-Cohen model
#'
#' XXXSplit a matrix into block diagonal sub matrices according to clusters and
#' XXXcombine the lower triangular parts into a vector
#'
#'
#' @param limit a square matrix
#' @param malefreq the frequency of male births (defaults to 0.51). Set to NULL to force equal male/female distribution among donor offspring
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
#' @export lower.tri.vector

hcc <- function(limit = 12,
                S = 1,
                l = 0.86,
                A = NULL,
                Q = 6/7,
                C = 4,
                PopAgeDifference = 1,
                DonorChildAge,
                DonorAge,
                NaturalChildren = c(.2, .3, .3, .2),
                malefreq = .51,
                familyfreq = 0,
                criteria = c("offspring", "families")) {

    ## Sanity check. The geography vectors should have the same length or be scalar

    ## NaturalChildren
    if (sum(abs(NaturalChildren) < .Machine.eps^2)) {
        stop("NaturalChildren must be a vector of probabilities")
    }
    NaturalChildren <- abs(NaturalChildren)/sum(abs(NaturalChildren))
    n.nat <- sample(seq(0,length(NaturalChildren)), 1, prob=NaturalChildren)
    gender.nat <- rbinom(nnat, size=1, prob=malefreq)

    ## AID children
    n.aid <- limit
    gender.aid <- rbinom(naid, size=1, prob=malefreq)

    ## Number of combinations


    ## Fathers age

    ## Sample AID as from our data
    ##  res <- tmp[Donor==sample(Donor, 1)]
    ##  n <- res$n



}
