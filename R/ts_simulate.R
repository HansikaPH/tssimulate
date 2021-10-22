library(forecast)
library(smooth)

set.seed(1)

#' Create the format of the dataframe returned from all the functions
#' @return The `series_set` which is the formatted dataframe
#' @noRd
format_df <- function (series_set){
  series_set <- as.matrix(series_set)
  rownames(series_set) <- paste0("ts", 1:nrow(series_set))
  colnames(series_set) <- c(1:ncol(series_set))
  series_set
}

#' Simulating from AR model
#'
#' @description This function simulates a defined number of series from a defined AR model.
#'
#' @param ar_order The order of the AR model to be used for simulation
#' @param length The length of the simulated series
#'
#' @return A list with named arguments `parameters` and the `series` generated
#' @example sim_ar(ar_order=3, length=100)
#'
#' @export
sim_ar <- function(ar_order, length, params=NULL, sample_series=NULL, no_series=1, max_root=5, burn_in=100){
  series_set <- NULL
  if(!is.null(params) & !is.null(sample_series)){
    stop("please specify only one out of the arguments 'params' or 'sample_series'")
  }
  else if (!is.null(sample_series)){
    if (class(sample_series) != "ts"){
      sample_series <- ts(sample_series, frequency = 1)
    }
    tryCatch(
    {
      arima_mod <- Arima(sample_series, order=c(ar_order,0,0))
    },
    error=function(cond) {
      message(paste0("Error when fitting AR model of order ", ar_order, " to the provided series."))
      message(cond)
    }
    )
    params <- arima_mod$coef
    for (i in 1:no_series){
      series <- simulate(arima_mod, nsim=length + burn_in, seed=i)
      series <- series[(burn_in + 1): length(series)]
      series_set <- rbind(series_set, series)
    }
  } else{
    if (is.null(params)){ # generate coefficients automatically if not specified by user
      if(max_root <= 1.1)
        stop("max_root has to be bigger than 1.1")

      l <- ar_order
      s <- sign(runif(l,-1,1))

      # the AR process is stationary if the absolute value of all the roots of the characteristic polynomial are greater than 1
      polyRoots <- s*runif(l,1.1,max_root)

      #calculate coefficients from the roots of the characteristic polynomial
      coeff <- 1
      for(root in polyRoots)
        coeff <- c(0,coeff) - c(root*coeff,0)

      nCoeff <- coeff / coeff[1]
      params <- -nCoeff[2:length(nCoeff)]
    }else if(ar_order != length(params)){
      stop("the ar_order needs to be equal to the length of the params vector")
    }
    for (i in 1:no_series){
      ts <- arima.sim(model=list(ar=params), n=length, n.start = burn_in)
      series_set <- rbind(series_set, ts)
    }
    # set the names for parameters
    names(params) <- paste0("ar", 1:ar_order)
  }
  series_set <- format_df(series_set)
  list("parameters" = params, "series"=series_set)
}

# simulating from SAR model
sim_sar <- function(sar_order, length, frequency=12, params=NULL, sample_series=NULL, no_series=1, burn_in=100){
  series_set <- NULL
  if(!is.null(params) & !is.null(sample_series)){
    stop("please specify only one out of the arguments 'params' or 'sample_series'")
  }else if (!is.null(sample_series)){
    if (class(sample_series) != "ts"){
      sample_series <- ts(sample_series, frequency = frequency)
    }
    tryCatch(
    {
      seasonal_arima_mod <- Arima(sample_series, seasonal=c(sar_order,0,0))
    },
    error=function(cond) {
      message(paste0("Error when fitting SAR model of order ", sar_order, " to the provided series."))
      message(cond)
    }
    )
    params <- seasonal_arima_mod$coef
    for (i in 1:no_series){
      series <- simulate(seasonal_arima_mod, nsim=length + burn_in, seed=i)
      series <- series[(burn_in + 1): length]
      series_set <- rbind(series_set, series)
    }
  }else{
    if (is.null(params)){
      params <- as.numeric(sim.ssarima(ar.orders=c(0,sar_order), i.orders=c(0,0), ma.orders=c(0,0), lags=c(0,frequency),
                          frequency = frequency, bounds = "admissible", obs=(length+burn_in))$AR)
    }else if(sar_order != length(params)){
      stop("the sar_order needs to be equal to the length of the params vector")
    }
    sim <- sim.ssarima(ar.orders=c(0,sar_order), i.orders=c(0,0), ma.orders=c(0,0), AR=params, lags=c(0,frequency),
                            frequency = frequency, bounds = "admissible", obs=(length+burn_in), nsim=no_series)
    if(no_series == 1){
      series_set <-  t(sim$data[(burn_in+1):(length+burn_in)])
    }else{
     series_set <- as.matrix(t(sim$data)[,(burn_in+1):(length+burn_in)])
    }
    names(params) <- paste0("sar", 1:sar_order)
  }
  series_set <- format_df(series_set)
  list("parameters" = params, "series"=series_set)
}

# simulating from Chaotic Logistic Map DGP
sim_chaotic_logistic_map <- function(length, params=NULL, no_series=1, burn_in=100){
  if (!is.null(params)){
    if (length(params) != 1){
      stop("please specify only one param for the Chaotic Logistic Map DGP")
    }else if(params <= 1 | params >= 4){
      message("param value not in the range (1,4). please note that the generated series may not show the intended chaotic behaviour.")
    }
  }else{
    params <- runif(1, min=1, max=4)
  }
  series <- numeric(length + burn_in + 1)
  initial_value <- runif(1, min=0, max=1)
  series[1] <- initial_value
  series_set <- NULL

  for (i in 1:no_series){
    for (j in 2:(length + burn_in + 1)){
      series[j] <- params * series[j-1] * (1 - series[j-1])
    }
    noise <- rnorm(length + burn_in + 1)/10
    series <- series + noise
    series <- pmax(series, 0)
    series <- series[(burn_in + 1):length(series)]
    series_set <- rbind(series_set, series)
  }
  names(params) <- "r"
  series_set <- format_df(series_set)
  list("parameters" = params, "series"=series_set)
}

# simulating from Mackey-Glass Equations DGP
sim_mackey_glass_equations <- function(length, tau=NULL, sample=NULL, no_series=1, burn_in=100){
  # default values of coefficients to ensure chaoticity in series
  beta <- 0.2
  gamma <- 0.1
  n <- 10

  no_discrete_steps <- 1000

  # generate tau value if null
  if (is.null(tau)){
    tau <- runif(1, min=17, max=30)
  }else if (tau < 17){
    message("tau value less than 17. please note that the generated series may be perfectly periodic.")
  }
  if (is.null(sample)){
    sample <- runif(1, tau/100, tau/10) # calculate sample from tau and no_discrete_steps
  }
  sample <- floor(no_discrete_steps * sample / tau)
  # burn-in period + actual length of the series sampled at a particular rate
  grids <- (no_discrete_steps * burn_in) + (sample * length)

  A <- (2 * no_discrete_steps - gamma * tau) / (2 * no_discrete_steps + gamma * tau)
  B <- beta * tau / (2 * no_discrete_steps + gamma * tau)

  series_set <- NULL
  for (i in 1:no_series){
    series <- numeric(grids)
    series[1:no_discrete_steps] <- 0.5 + 0.05 * (-1 + 2 * runif(no_discrete_steps, min=0, max=1))
    series[no_discrete_steps+1] <- A * series[no_discrete_steps] + B * (series[1] / (1 + series[1] ^ n))
    for (j in (no_discrete_steps+1):(grids-1)){
      series[(j + 1)] <- A * series[j] + B * (series[(j - no_discrete_steps)] / (1 + series[(j - no_discrete_steps)] ^ n) +
                                 series[(j - no_discrete_steps + 1)] / (1 + series[j - no_discrete_steps + 1] ^ n))
    }
    series <- series[(no_discrete_steps * burn_in + 1):length(series)] # discard the burn-in
    sampled_indices <- seq(from=1, to=length(series), by=sample)
    series <- series[sampled_indices]
    series_set <- rbind(series_set, series)
  }
  params <- c(tau, beta, gamma, n, sample)
  names(params) <- c("tau", "beta", "gamma", "n", "sample")
  series_set <- format_df(series_set)
  list("parameters" = params, "series"=series_set)
}

# simulating from Fourier Terms DGP
sim_fourier_terms <- function(length, seasonal_periods, no_fourier_terms=NULL, params=NULL, no_series=1, burn_in=100, seas_sd=0.3){
  # input sanity checking
  if (!is.null(no_fourier_terms)){
    # check no_fourier_terms against params
    if(!is.null(params)){
      if(class(params) != "list"){
        params <- list(params)
      }
      if (length(no_fourier_terms) != length(params)){
        stop("length of the no_fourier_terms vector and params list are nost consistent.
        please specify the params as a list where each list element corresponds to the fourier term coefficients for a single period")
      }
      for (i in seq_along(params)){
        if (length(params[[i]]) != 2*no_fourier_terms[i]){
          stop("number of fourier coefficients for each seasonality period should be twice the no_fourier_terms for that seasonality")
        }
      }
    }else{ # check no_fourier_terms against seasonal_periods
      if(length(seasonal_periods) != length(no_fourier_terms))
        stop("length of the no_fourier_terms vector needs to be equal to the length of the seasonal_periods vector")
    }
  }
  if (!is.null(params)){
    if(class(params) != "list"){
      params <- list(params)
    }
    # check params against seasonal_periods
    if (length(seasonal_periods) != length(params))
      stop("if specifying params, please specify them correctly for every seasonality period")

    # check if the number of coefficients for each period is a multiple of 2
    for (i in seq_along(params)){
      if(length(params[[i]]) %% 2 != 0){
        stop("the number of coefficients for each period should be a multiple of 2. please specifiy the same number of coefficients for both sine and cosine terms.")
      }
    }
  }

  # set parameters after input sanity checking
  if(is.null(no_fourier_terms)){
    if (is.null(params))
      no_fourier_terms <- rep(5, length(seasonal_periods))
    else{
      no_fourier_terms<-c()
      for (i in seq_along(params)){
        no_fourier_terms <- c(no_fourier_terms, length(params[[i]])/2)
      }
    }
  }
  if(is.null(params)){
    params <- list()
    for (i in seq_along(no_fourier_terms)){
      params[[i]] <- rnorm(2*no_fourier_terms[i])
    }
  }

  # generate the series
  series_set <- NULL
  for (i in 1:no_series){
    series <- numeric(length + burn_in)
    for(j in seq_along(seasonal_periods)){
      period <- seasonal_periods[j]
      n_seasons <- (length + burn_in) %/% period + 1
      coeff <- params[[j]]
      n_coeff <- no_fourier_terms[j]

      sinx <- cosx <- matrix(0, ncol=n_coeff, nrow=period)

      for(k in seq(n_coeff)) {
        sinx[,k] <- sin(2*pi*k*seq(period)/period)
        cosx[,k] <- cos(2*pi*k*seq(period)/period)
      }
      terms <- cbind(sinx,cosx)
      adjusted_coeff <- coeff + rnorm(length(coeff), sd=seas_sd)

      seasonality <- as.numeric(terms %*% adjusted_coeff)
      seasonality <- rep(seasonality, n_seasons)
      seasonality <- seasonality %>% head(length + burn_in) %>% scale() %>% as.numeric()
      series <- series + seasonality

      # set the names of the param properly
      if(i==1){
        names(params[[j]]) <- c(paste0("sin", 1:no_fourier_terms[[j]]), paste0("cos", 1:no_fourier_terms[[j]]))
      }
    }
    # add remainder noise to the seasonality
    remainder <- rnorm((length + burn_in), sd=0.5)
    series <- series + remainder
    series <- series[(burn_in+1):length(series)]
    series_set <- rbind(series_set, series)
  }
  series_set <- format_df(series_set)
  params <- list("fourier_term_coefficients"=params, "seas_sd"=seas_sd)
  list("parameters" = params, "series"=series_set)
}
