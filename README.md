# SVHM
Implements a support vector machine algorithm for cencored data to predict event times in survival analysis. The algorithm is described in detail in [[1]](#1).

# Installation
The package can be installed via the following command:
```sh
devtools::install_github("herglola/SVHM")
```
The functions 
```sh
create_svhm()
create_time_svhm()
train_svhm()
train_time_svhm()
```
can be executed directly. All other functions are meant for internal usage and can be accessed via the 
```sh
SVHM:::func_name()
```
command.

# Minimal examples
Here is an example for the dataset bmt from the R package KMsurv

```sh
library(KMsurv)
library(SVHM)

##############
# Parameters #
##############

opt = 'osqp'
gamma_squared <- 100
k <- 3
cross_validation_val <- 4
test_size=.3
cost_grid <- 2^c(-6:6)

covariates <- c('z1', 'z3', 'z7')

######################
#  Model prediction  #
######################

data(bmt, package='KMsurv')

model <- create_svhm(bmt,
                     covariates,
                     cross_validation_val=cross_validation_val,
                     cost_grid=cost_grid,
                     varName_cencored="d3",
                     varName_futime="t2",
                     k=k,
                     test_size=test_size,
                     opt=opt,
                     gamma_squared=gamma_squared,
                     choose = 'c')
```
The following example is for the time dependent model of the SVHM:
```sh
library(timereg)
library(SVHM)
##############
# Parameters #
##############

opt = 'osqp'
gamma_squared <- 200
test_size=.3
cost <- .5

covariates <- c('prot')

######################
#  Model prediction  #
######################

data(csl)

time_model <- create_time_svhm(csl, 
                                c("prot", "prot.base"), 
                                .5, 
                                varName_cencored='dc',
                                varName_futime='eventT', 
                                start_interval='lt', 
                                end_interval='rt',
                                test_size=.3,
                                opt='osqp,
                                gamma_squared=gamma_squared)

```
# Notes and Remarks
Optimization is done either via the Rmosek package or the osqp package in R. And while the osqp package is open to use for everybody the Rmosek package requires a license which can be aquired for example via an academic institution.

## References
<a id="1">[1]</a> 
Yuanjia Wang and Tianle Chen and Donglin Zeng (2016)
Support Vector Hazards Machine: A Counting Process Framework for Learning Risk Scores for Censored Outcomes
Journal of Machine Learning Research, 17(167), 1-37.
http://jmlr.org/papers/v17/16-007.html
