# MESS

Development version of the R package **MESS** (Miscellaneous Esoteric
Statistical Scripts). This package contains a collection of various
semi-useful functions that I have written and over the years. 

There is no real overall theme to the functions in the package but I
have tried to group the contents in the overview below.

To install the development version of **MESS** run the following
command from within R (this requires that the devtools package is
already installed on the system.)

```r
devtools::install_github('ekstroem/MESS')
```

[![Travis-CI Build Status](https://travis-ci.org/ekstroem/MESS.svg?branch=master)](https://travis-ci.org/ekstroem/MESS) ![Download counter](http://cranlogs.r-pkg.org/badges/grand-total/MESS)


# Package overview

The list below is far from complete.

## Statistical functions

* `cmd` - Correlation matrix distance. A measure of the similarity of two matrices of equal dimensions.
* `drop1.geeglm` - A `drop1` extension for `geeglm` objects
* `gkgamma` - compute Goodman-Kruskal's gamma statistic for a
two-dimensional table of ordered categories
* `monte_carlo_chisq_test` - Monte Carlo tests of 2x2 tables with fixed margin(s). Also works for r x c tables.
* `mfastLmCpp` - fast computation of simple regression slopes for each
predictor represented by a column in a matrix
* `ks_cumtest` - One-sample Kolmogorov-Smirnov discrete cumulative comparison.
* qic - 
* `rud` - Randomized treatments for an RCT based on an urn model
* screenr -

Power calculations are - for better or worse - part of my job, and
while I generally recommend researchers to simulate their design and
test procedure, there are some benefits to be gained from asymptotic
approximations in standard designs.

In time these will be moved to the `Austin` package, but for now they
still reside in `MESS`.

* `power_prop_test` extends the standard power.prop.test function to
accommodate different group sizes.
* `power_t_test` extends the standard power.t.test function to
accommodate different group sizes and/or variances.
* `power_mcnemar_test` - power calculations for exact and asymptotic McNemar test in a 2 by 2 table 
* `power_binom_test` - power calculations for exact test of a simple null hypothesis in a Bernoulli experiment

## Graphical functions

* col.shade
* col.tint
* col.alpha
* `rootonorm` - plot Tukey's hanging root-o-gram for comparison of a histogram to a normal distribution.
* `wallyplot` - plot a Wally plot for evaluation of a residual plot.


## Utility functions

* `age` - computes the age in years of a person from two date (birth date and current date) 
* `auc` - computes the area under the curve for two vectors (x-values and y-values). Can handle ranges and missing observations
* `categorize` - produce tables using a data argument 
* `expand_table` - Expand table or matrix to data frame where each observation
in the table becomes a single observation in the data frame with
corresponding information for each for each combination of the table
dimensions.
* `fac2num` - convert a factor to numerical. I keep forgetting the simple code for this so ended up writing a function that could do it.
* `filldown` - Fill down missing values in a vector with the latest non-missing
value. A last measurement carried forward function. And fast.
* lower.tri.vector - 
* `write.xml` - save a data frame as an xml file


## Computational functions

* `conditional_rowMeans` - a function similar to `rowMeans` but it only returns the mean if a prespecified number of observations is available.
* `pairwise_Schur_product` - compute Schur products (element-wise) of all pairwise combinations of columns in matrix
* `qdiag` - fast extraction of matrix diagonal
* `quadform` - fast computation of a quadratic form `t(X) %*% M %*% X`
* `repmat` - fast replication of a matrix. Can be replicated both row-wise and columns-wise.
* `sinv` - invert a symmetric positive-definite matrix. Some speedup can be gained if we know the matrix is symmetric
* `tracemp` - fast computation of trace of matrix product `trace(t(A) %*% B)`


## Datasets

* `bdstat` - Monthly live births and deaths in Denmark from January 1901 to March 2013
* `earthquakes` - Information on earthquakes worldwide in 2015 with a magnitude greater than 3 on the Richter scale.
* `happiness` - Data on happiness, taxation rates, countris and continents
* `kwdata` - example data to show that Krushkal-Wallis' test examines more than differences in medians.
* `soccer` - List of players on the Danish national soccer team
