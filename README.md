# lcphR
lcphR: Latent Class Proportional Hazards Regression

# Installation Guide
```{r}
# install.packages("devtools")
devtools::install_github("tengfei-emory/lcphR")
library(lcphR)
```
Currently `lcphR` supports R version >= 4.1.1.

# Example: analyze a simulated dataset

## Data simulation

By default, the function `simulation(n)` generates a dataset with n observations under the scenario 1 described by Fei, Hanfelt and Peng (submitted).
```{r}
# generate a dataset with 500 individuals
dat <- simulation(500)
```
Specifically, it returns a data frame of 2 latent classes with 2 baseline covariates (`Xcov1` and `Xcov2`), time-to-event (`tildet`), censoring indicator (`delta`) and true latent class labels (`latent`).

## Model fitting

The analysis for the dataset `dat` can be conducted by running `LSCA` function:

  ```{r}
lcfit <- lcphR(dat,num_class=2,covx=c('Xcov1','Xcov2'),tolEM=1e-3,varest=T,traceplot=T,initial='kmeans')
```

The output list `lcfit` contains the following information:

  `alpha`: estimated unknown parameters under the class membership probability submodel.

`zeta`: estimated unknown parameters under the class-specific proportional hazards submodel.

`tevent`: uncensored event times.

`chaz`: estimated baseline hazard function for class 1.

`ASE`: standard error estimates for alpha and zeta, by the close-form estimator derived from observed-data log likelihood.

`nASESr`: standard error estimates for alpha and zeta, by the numerical estimator derived from observed-data profile log likelihood.

`chaztgt` and `chaztgtASE`: point estimates and corresponding standard error estimates for the baseline hazard function at specified times (under development).

`tau`: posterior membership probabilities

`p`: baseline membership probabilities

`loglik`, `obsloglik`: complete-data log likelihood and observed-data log likelihood after the last iteration of the EM algorithm.

`AIC`, `BIC`, `CEBIC`: Akaike information criterion, Bayesian information criterion, and Classification entropy incorporated BIC.

`timediff`: Total time used in point estimation and variance estimation.

`diffEM`: L-2 norm of the difference of posterior membership probability vector between the last and one prior to the last iteration.

`diffPAR`: L-2 norm of the difference of point estimates (alpha, beta, hazard increments) between the last and one prior to the last iteration.

`difflA`: convergence criterion (Aitken acceleration).

`numiter`: total number of iterations.

`entropy`: standardized entropy R-square ranging from 0 to 1. A value close one indicates better separation of classes.

`censor`: censoring rate of the observations.

## Reference

Fei, T, Hanfelt, J, Peng, L. [Latent Class Proportional Hazards Regression with Heterogeneous Survival Data](https://arxiv.org/pdf/2202.00775.pdf). Statistics and its Interface, accepted.
