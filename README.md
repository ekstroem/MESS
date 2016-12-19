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


Power calculations are - for better or worse - part of my job, and
while I generally recommend researchers to simulate their design and
test procedure, there are some benefits to be gained from asymptotic
approximations in standard designs.

In time these will be moved to the Austin package, but for now they
still reside in MESS.


## Utility functions

* fac2num - convert a factor to numerical
* write.xml - save a data frame as an xml file


## Computational functions

* quadform - fast computation of a quadratic form $t(X) %*% X$
* repmat - fast replication of a matrix
* sinv - invert a symmetric positive-definite matrix
* tracemp - fast computation of trace of matrix product


