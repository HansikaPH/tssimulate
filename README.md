
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tssimulate

<!-- badges: start -->

<!-- badges: end -->

tssimulate is a library for simulating time series from known Data
Generating Processes.

## Installation

<!-- You can install the released version of tssimulate from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->

<!-- install.packages("tssimulate") -->

<!-- ``` -->

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("HansikaPH/tssimulate")
```

## Example

This is a basic example which shows you how to simulate 5 time series of
length 10, using an AR(3) DGP:

``` r
library(tssimulate)
## basic example code

sim <- sim_ar(ar_order = 3, length = 10, no_series = 5)
sim$parameters
#>        ar1        ar2        ar3 
#> 0.30785840 0.34958885 0.05190761

sim$series
#>              1          2          3          4            5          6
#> ts1  1.5568155  0.2707576  2.1837464  0.1933141  2.368602772  2.2988208
#> ts2  2.4328491  1.4857462  1.8310515  2.8889276  1.841875181  2.2023139
#> ts3 -0.5567058 -0.4812871 -1.3628601  0.4366138 -1.103635228 -1.0264660
#> ts4  0.4577363  0.3026447  0.3504894  2.1333015 -0.003492689  0.8400392
#> ts5 -0.1047165 -0.1143659 -0.6000691 -1.5226410 -1.877783146  0.4052015
#>              7          8          9         10
#> ts1  0.8681080  0.8315542  0.1689837 -0.8593711
#> ts2  1.9758661  2.4245447  1.2681592  0.6638493
#> ts3 -1.2334030  1.6019393 -0.8785626  0.1883955
#> ts4  1.2166242 -1.4299060 -0.3416499 -0.6351076
#> ts5 -0.9690013 -1.0751606 -3.4468447 -1.5273845
```
