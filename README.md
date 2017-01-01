# MESS

Development version of the R package **MESS** (Miscellaneous Esoteric Statistical Scripts)

To install the development version of **MESS** run the following
command from within R (this requires that the devtools package is
already installed on the system.

```r
devtools::install_github('ekstroem/MESS')
```

![Download counter](http://cranlogs.r-pkg.org/badges/grand-total/MESS)


# Package overview

The package contains a collection of various functions that have
repeatedly. Several of them are shortcuts to combinations of standard
R function that I keep using (and forgetting) repeatedly.

The list below is far from complete.


## Statistical functions

* cmd - Correlation matrix distnce
* drop1.geeglm - 
* gkgamma - compute Goodman-Kruskal's gamma statistic for a
two-dimensional table of ordered categories
* mfastLmCpp - fast computation of simple regression slopes for each
predictor represented by a column in a matrix

Power calculations are - for better or worse - part of my job, and
while I generally recommend researchers to simulate their design and
test procedure, there are some benefits to be gained from asymptotic
approximations in standard designs.

In time these will be moved to the Austin package, but for now they
still reside in MESS.

* power_t_test extends the standard power.t.test function to
accommodate different group sizes and/or variances.
* power_mcnemar_test
* power_binomial_test


## Graphical functions

* wallyplot
* col.shade
* col.tint
* col.alpha

## Utility functions

* age - compute the age of a person from two date (birth date and
  current date)
* auc - computes the area under the curve for two vectors (x-values
and y-values)
* expand - Expand table or matrix to data frame where each observation
in the table becomes a single observation in the data frame with
corresponding information for each for each combination of the table
dimensions.
* fac2num - convert a factor to numerical
* filldown - Fill down missing values with the latest non-missing
value
* lower.tri.vector - 
* write.xml - save a data frame as an xml file


## Computational functions

* qdiag - fast extraction of matrix diagonal
* quadform - fast computation of a quadratic form $t(X) %*% X$
* repmat - fast replication of a matrix
* sinv - invert a symmetric positive-definite matrix
* tracemp - fast computation of trace of matrix product trace(t(A) %*% B)


