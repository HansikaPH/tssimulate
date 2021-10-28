source("./ts_simulate.R")

plot_series <- function (series){
  plot(series, type="l")
}

plot_pacf <- function (series){
  ggPacf(series, lag.max = 50)
}

plot_acf <- function (series){
  ggAcf(series, lag.max = 50)
}

plot_series(result$series[1,])
plot_pacf(result$series[1,])
plot_acf(result$series[1,])

# test ar series generation
result <- sim_ar(ar_order = 5, length = 200)
result <- sim_ar(ar_order = 5, length = 200, no_series = 5)
result <- sim_ar(ar_order = 5, length = 200, no_series = 10, burn_in = 200)
result <- sim_ar(ar_order = 5, length = 200, no_series = 10, max_root = 0.5)
result <- sim_ar(ar_order = 5, length = 200, no_series = 10, max_root = 10)
result <- sim_ar(ar_order = 5, length = 200, no_series = 1, max_root = 10)
result <- sim_ar(ar_order = 5, params=0.5, length = 200, no_series = 10)
result <- sim_ar(ar_order = 2, params=c(0.5, 0.2), length = 200, no_series = 10)
result <- sim_ar(ar_order = 2, params=c(0.5, 0.2), length = 200, no_series = 1)
result <- sim_ar(ar_order = 2, params=c(0.5, 0.2), sample_series = USAccDeaths, length = 200, no_series = 10)
result <- sim_ar(ar_order = 5, sample_series = USAccDeaths, length = 200, no_series = 10)
result <- sim_ar(ar_order = 5, sample_series = rep(3,100), length = 200, no_series = 10)
result <- sim_ar(ar_order = 5, sample_series = rnorm(100), length = 200, no_series = 10)
result <- sim_ar(ar_order = 5, sample_series = rnorm(100), length = 200, no_series = 1)


# test sar series generation
result <- sim_sar(sar_order = 4, length = 200, no_series = 10)
result <- sim_sar(sar_order = 4, length = 200, no_series = 10, frequency = 12)
result <- sim_sar(sar_order = 4, length = 200, no_series = 1, frequency = 12)
result <- sim_sar(sar_order = 4, length = 200, no_series = 10, frequency = 7, burn_in = 200)
result <- sim_sar(sar_order = 4, params=c(0.3, 0.2, 0.1, 0.2), length = 200, no_series = 10, frequency = 7)
result <- sim_sar(sar_order = 3, params=c(0.3, 0.2, 0.1, 0.2), length = 200, no_series = 10, frequency = 7)
result <- sim_sar(sar_order = 4, params=c(0.3, 0.2, 0.1, 0.2), length = 200, no_series = 1, frequency = 7)
result <- sim_sar(sar_order = 4, params=c(0.3, 0.2, 0.1, 0.2), sample_series = USAccDeaths, length = 200, no_series = 1, frequency = 12)
result <- sim_sar(sar_order = 2, sample_series = USAccDeaths, length = 200, no_series = 1, frequency = 12)
result <- sim_sar(sar_order = 2, sample_series = USAccDeaths, length = 200, no_series = 10, frequency = 12)
result <- sim_sar(sar_order = 2, sample_series = rnorm(100), length = 200, no_series = 10, frequency = 12)

# test chaotic logistic map series generation
result <- sim_chaotic_logistic_map(length = 200)
result <- sim_chaotic_logistic_map(length = 200, no_series = 10)
result <- sim_chaotic_logistic_map(length = 200, no_series = 10, burn_in = 300)
result <- sim_chaotic_logistic_map(length = 200, params = c(0.5, 0.6))
result <- sim_chaotic_logistic_map(length = 200, params = 0.6)
result <- sim_chaotic_logistic_map(length = 200, params = c(3.5))
result <- sim_chaotic_logistic_map(length = 200, params = 3.5)

# test mackey glass series generation
result <- sim_mackey_glass_equations(length = 200)
result <- sim_mackey_glass_equations(length = 200, tau=5)
result <- sim_mackey_glass_equations(length = 200, tau=50)
result <- sim_mackey_glass_equations(length = 200, sample=0.2)
result <- sim_mackey_glass_equations(length = 200, tau=50, sample=10)
result <- sim_mackey_glass_equations(length = 200, no_series = 10)
result <- sim_mackey_glass_equations(length = 200, no_series = 10, burn_in = 300)

# test fourier terms series generation
result <- sim_fourier_terms(length = 200, seasonal_periods = 7)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7))
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7), no_series = 10)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7, 24*15), params=c(0.5, 0.3), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = 24, params=c(0.5, 0.2), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = 24, params=0.5, seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7, 24*15), params=list(c(0.5, 0.3), 0.3), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7, 24*15), params=list(c(0.5, 0.3), c(0.5, 0.3, 0.3, 0.4), c(0.5, 0.3)), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7, 24*15), no_fourier_terms = 4, seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7, 24*15), no_fourier_terms = c(4, 3, 1), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = 24, no_fourier_terms = 4, seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = 24, no_fourier_terms = 1, params=c(0.5, 0.3, 0.2), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = 24, no_fourier_terms = 1, params=c(0.5, 0.3), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7), no_fourier_terms = 1, params=list(c(0.5, 0.3), c(0.5,0.3)), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7), no_fourier_terms = c(1, 1), params=list(c(0.5, 0.3), c(0.5,0.3)), seas_sd = 5)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7), no_fourier_terms = c(1, 1), params=list(c(0.5, 0.3), c(0.5,0.3)), seas_sd = 5, no_series = 10)
result <- sim_fourier_terms(length = 200, seasonal_periods = c(24, 24*7), no_fourier_terms = c(1, 1), params=list(c(0.5, 0.3), c(0.5,0.3)), seas_sd = 5, no_series = 10, burn_in = 100)
