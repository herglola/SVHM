# SVHM
Implements a support vector machine algorithm for cencored data to predict event times in survival analysis. The algorithm is described in detail in [[1]](#1).

Here is an example for the dataset bmt from the R package KMsurv

```sh
library(KMsurv)
library(SVHM)

##############
# Parameters #
##############

gamma_squared <- .005
k <- 3
cross_validation_val <- 4
test_size=.3
cost_grid <- 2^c(-12:12)

covariates <- c('z1', 'z2', 'z3', 'z4', 'z5', 'z6', 'z7', 'z8', 'z9', 'z10')

######################
#  Model prediction  #
######################

data(bmt)

model <- SVHM:::create_svhm(bmt, 
                            covariates, 
                            cross_validation_val, 
                            cost_grid,
                            gamma_squared, 
                            k, 
                            test_size, 
                            varName_cencored="d3",
                            varName_futime = "t2", 
                            opt='mosek')
```

## References
<a id="1">[1]</a> 
Yuanjia Wang and Tianle Chen and Donglin Zeng (2016)
Support Vector Hazards Machine: A Counting Process Framework for Learning Risk Scores for Censored Outcomes
Journal of Machine Learning Research, 17(167), 1-37.
http://jmlr.org/papers/v17/16-007.html
