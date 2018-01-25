# MESS 0.5.1

## Bug fixes

*  Fixed a case in `screen_variables` where R currently only emits a
   warning when if/while statement is used with a condition of length
   greater than one.


# MESS 0.5.0

## Changes and bug fixes

*  Updated the limits of the `residualplot` function to ensure that the blurred area are always within the plot region
*  Updated the `residualplot` function to accommodate glm objects with binomial and poisson families
*  Fixed bug in `auc` when spline interpolation was used. Thanks to james-jenkins!
*  Fixed the exact test for `power_mcnemar_test` including improved documentation.


## New functions

*  Added the `cumsumbinning` to create binned variables based on cumulative sums with threshold.
*  Added the `pairwise_Schur_product` to compute the products (element-wise) of all pairwise combinations of columns in matrix.
*  Added the `ks_cumtest` to compute a Kolmoggorov-Schmirnov test for ordinal cumulate proportions.
*  Added the `rud` function to simulate a statistical design from the random urn design model.
*  Added the `power_prop_test` to allow for different sample sizes.

