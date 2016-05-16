# MESS

Development version of the R package MESS

To install the development version of MESS run the following commands
from within R
```{r}
library(devtools)
install_github('ekstroem/MESS')
```

The stable version of MESS can always be obtained from CRAN

```{r}
install.packages("MESS")
```

## Showcase

### Utility functions

* age
* auc
* fac2num
* repmat


## Planned future functions

* A function to extract the estimated working correlation matrix from
geeglm objects
* A function to extract the naive standard errors from geepack
* A fast age calculator. The `age` function already implements this
  
