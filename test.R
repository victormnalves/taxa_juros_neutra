library(tis)
library(tidyverse)
library(mFilter)
library(nloptr)

shiftQuarter <- function(original.start,shift){
  #################################################################
  # This function takes in a (year,quarter) date in time series format
  # and a shift number, and returns the (year,quarter) date corresponding
  # to the shift. Positive values of shift produce leads and negative values
  # of shift produce lags.
  # For example, entering 2014q1 with a shift of -1 would return 2013q4.
  # Entering 2014q1 with a shift of 1 would return 2014q2.
  # In each case, the first argument of the function must be entered as
  # a two-element vector, where the first element corresponds to the year
  # and the second element corresponds to the quarter.
  # For example, Q12014 must be entered as "c(2014,1)".
  ################################################################

  # Leads (positive values of shift)
  if (shift > 0) {
    new.start = c(0,0)
    sum = original.start[2] + shift

    # Get the year value
    if (sum <= 4) {
      new.start[1] = original.start[1]
    }
    else {
      new.start[1] = original.start[1] + ceiling(sum/4) - 1
    }

    # Get the quarter value
    if (sum %% 4 > 0) {
      new.start[2] = sum %% 4
    }
    else {
      new.start[2] = sum %% 4 + 4
    }
  }

  # Lags (negative values of shift)
  else {
    new.start = c(0,0)
    diff = original.start[2] - abs(shift)

    # Get the year value
    if (diff > 0) {
      new.start[1] = original.start[1]
    }
    else {
      new.start[1] = original.start[1] - (1 + floor(abs(diff)/4))
    }

    # Get the quarter value
    if (diff %% 4 > 0) {
      new.start[2] = diff %% 4
    }
    else {
      new.start[2] = diff %% 4 + 4
    }
  }

  return(new.start)}

calculate.covariance <-
  function(initial.parameters,theta.lb,theta.ub,
           y.data,x.data,stage,lambda.g=NA,lambda.z=NA,xi.00){

    n.state.vars <- length(xi.00)

    # Set covariance matrix equal to 0.2 times the identity matrix
    P.00 <- diag(0.2, n.state.vars, n.state.vars)

    # Get parameter estimates via maximum likelihood
    f <- function(theta) {return(-log.likelihood.wrapper(theta,
                                                         y.data, x.data, stage,
                                                         lambda.g, lambda.z,
                                                         xi.00, P.00)$ll.cum)}

    nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                         lb=theta.lb,ub=theta.ub,opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))

    theta <- nloptr.out$solution

    # Run Kalman filter with above covariance matrix and corresponding parameter estimates
    states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)

    # Save initial covariance matrix
    P.00 <- states$filtered$P.ttm1[1:n.state.vars,]

    return(P.00)
  }

kalman.log.likelihood <- function(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x) {
  T <- dim(y)[1]
  n <- dim(y)[2]
  ll.vec <- matrix(NA,T,1)
  ll.cum <- 0
  xi.tt <- xi.tm1tm1
  P.tt  <- P.tm1tm1
  for (t in 1:T){

    xi.ttm1 <- F_matrix %*% xi.tt
    P.ttm1 <- F_matrix %*% P.tt %*% t(F_matrix) + Q
    prediction.error <- (as.vector(y[t,]) - as.vector(t(A) %*% as.vector(x[t,])) - as.vector(t(H) %*% xi.ttm1))
    HPHR <- t(H) %*% P.ttm1 %*% H + R
    ll.vec[t] <- drop(-(n / 2) * log(2 * atan(1) * 4) - 0.5 * log(det(HPHR))
                      -0.5 * prediction.error %*% solve(HPHR, prediction.error, tol = 0))
    ll.cum <- ll.cum + ll.vec[t]

    xi.tt <- xi.ttm1 + P.ttm1 %*% H %*% solve(HPHR, prediction.error, tol = 0)
    P.tt  <- P.ttm1 - P.ttm1 %*% H %*% solve(HPHR, t(H) %*% P.ttm1, tol = 0)
  }
  return(list("ll.vec"=ll.vec,"ll.cum"=ll.cum))
}

shiftMonth <- function(original.start,shift){
  #################################################################
  # This function takes in a (year,month) date in time series format
  # and a shift number, and returns the (year,month) date corresponding
  # to the shift. Positive values of shift produce leads and negative values
  # of shift produce lags.
  # For example, entering 2014m1 with a shift of -1 would return 2013m12.
  # Entering 2014m1 with a shift of 1 would return 2014m2.
  # In each case, the first argument of the function must be entered as
  # a two-element vector, where the first element corresponds to the year
  # and the second element corresponds to the month.
  # This function is analogous to shiftQuarter().
  ################################################################

  # Leads (positive values of shift)
  if (shift > 0) {
    new.start = c(0,0)
    sum = original.start[2] + shift

    # Get the year value
    if (sum <= 12) {
      new.start[1] = original.start[1]
    }
    else {
      new.start[1] = original.start[1] + ceiling(sum/12) - 1
    }

    # Get the month value
    if (sum %% 12 > 0) {
      new.start[2] = sum %% 12
    }
    else {
      new.start[2] = sum %% 12 + 12
    }
  }

  # Lags (negative values of shift)
  else {
    new.start = c(0,0)
    diff = original.start[2] - abs(shift)

    # Get the year value
    if (diff > 0) {
      new.start[1] = original.start[1]
    }
    else {
      new.start[1] = original.start[1] - (1 + floor(abs(diff)/12))
    }

    # Get the month value
    if (diff %% 12 > 0) {
      new.start[2] = diff %% 12
    }
    else {
      new.start[2] = diff %% 12 + 12
    }
  }

  return(new.start)}


getFRED <- function(url, freq = "Quarterly") {
  ##########################################################################################
  # This function downloads data from FRED. It returns quarterly data.
  # User must provide the FRED url.
  ###########################################################################################

  FREDraw <- readLines(url)
  # download.file(url, destfile = paste0('FREDtemp.txt'))
  #
  # txt.file.name <- paste0(folder,substr(url, regexpr('[a-zA-z0-9]*.txt',url),1000))
  # if (!file.exists(txt.file.name)){
  #     # Download the data from FRED
  #     download.file(url, destfile = paste0(folder,'FREDtemp.txt'))
  #     system(paste0('wget --no-check-certificate "', url, '"'))
  #     system(paste('mv',substr(url, regexpr('[a-zA-z0-9]*.txt',url),1000),txt.file.name))
  # }
  # FREDraw <- readLines(txt.file.name)

  # Frequency
  freq.FRED <- gsub(' ', '',substr(FREDraw[which(regexpr('Frequency', FREDraw)==1)],
                                   (nchar('Frequency')+2),100))

  # Where does the data start
  datastart = which(gsub(' ', '',FREDraw)=='DATEVALUE') - 2

  #data <- read.table('FREDtemp.txt', skip = datastart, header = TRUE)
  data <- FREDraw
  data <- read.table(textConnection(data), skip = datastart, header = TRUE)

  # get starting date
  first.year  <- as.numeric(format(as.Date(data$DATE[1]),'%Y'))
  first.month <- as.numeric(format(as.Date(data$DATE[1]),'%m'))

  # Adjust frequency
  if (freq.FRED == 'Quarterly'){
    first.q  <- (first.month-1)/3 + 1
    data.tis <- tis(data$VALUE, start = c(first.year, first.q), tif = 'quarterly')
  } else if (freq.FRED == 'Monthly') {
    data.tis <- tis(data$VALUE, start = c(first.year, first.month), tif = 'monthly')
  }

  # Convert frequency
  if (freq.FRED == 'Monthly' & freq == 'Quarterly') {
    data.tis <- convert(data.tis, tif = 'quarterly', method = 'constant', observed. = 'averaged')
  }

  return(data.tis)
}


splice <- function(s1, s2, splice.date, freq) {
  ##########################################################################################
  # This function splices two series, with the series s2 beginning at splice.date
  # and extended back using the growth rate at the splice.date times series s1
  # The freq argument accepts two values - 'quarterly' and 'monthly' -
  # but it could be modified to take more.
  ##########################################################################################
  t <- splice.date #renaming for convenience
  if (freq == "quarterly" | freq == "Quarterly") {
    t.minus.1 <- shiftQuarter(t,-1)
  }
  else if (freq == "monthly" | freq == "Monthly") {
    t.minus.1 <- shiftMonth(t,-1)
  }
  else { stop("You must enter 'quarterly' or 'monthly' for freq.") }
  ratio <- as.numeric(window(s2,start = t, end = t)/
                        window(s1,start = t, end = t))

  return(mergeSeries(ratio*window(s1,end = t.minus.1),window(s2, start = t)))
}


gradient <- function(f, x, delta = x * 0 + 1.0e-5) {
  ##########################################################################################
  # This function computes the gradient of a function f given a vector input x.
  ##########################################################################################
  g <- x * 0
  for (i in 1:length(x)) {
    x1 <- x
    x1[i] <- x1[i] + delta[i]
    f1 <- f(x1)
    x2 <- x
    x2[i] <- x2[i] - delta[i]
    f2 <- f(x2)
    g[i] <- (f1 - f2) / delta[i] / 2
  }
  return(g)
}

folder <- "data/input/"
file <- "rstar.data.us.csv"
pcGDP <- FALSE

get.US.data <- function(folder = "data/input/" , file ="rstar.data.us.csv", pcGDP = FALSE) {
  # Import data using the function getFRED() in utilities.R
  # If the connection does not work, try the old URL:
  # "https://research.stlouisfed.org/fred2/data/NAME.txt"

  data.start <- c(1960,1)
  data.start.year <- 1960
  data.start.quarter <- 1

  # create alternate output measures
  # choice between GDP and GDP per capita
  if (pcGDP == FALSE) {
    gdp <- getFRED('https://fred.stlouisfed.org/data/GDPC1.txt')
  }
  else {
    gdp <- getFRED('https://fred.stlouisfed.org/data/A939RX0Q048SBEA.txt')
  }


  price.index     <- getFRED('https://fred.stlouisfed.org/data/PCEPILFE.txt')
  #ny.discount     <- getFRED('https://fred.stlouisfed.org/data/INTDSRUSM193N.txt')
  fed.funds       <- getFRED('https://fred.stlouisfed.org/data/FEDFUNDS.txt')

  #------------------------------------------------------------------------------#
  # Prepare Data
  #------------------------------------------------------------------------------#

  # Take log of real GDP
  gdp.log <- log(gdp)

  # Create an annualized inflation series using the price index
  inflation <- 400*log(price.index/Lag(price.index, k=1))

  # Inflation expectations measure: 4-quarter moving average of past inflation
  inflation.expectations <- (inflation + Lag(inflation, k=1) + Lag(inflation, k=2) + Lag(inflation, k=3))/4

  # Express interest rate data on a 365-day basis
  #ny.discount.eff <- 100*((1+ny.discount/36000)^365 -1)
  fed.funds.eff   <- 100*((1+fed.funds/36000)^365 -1)

  # NY Fed discount rate is used prior to 1965; thereafter, use the effective federal funds rate
  #interest.rate <- mergeSeries(window(ny.discount.eff, end = c(1964,4)),window(fed.funds.eff, start = c(1965,1)))
  interest.rate <- fed.funds.eff
  #------------------------------------------------------------------------------#
  # Output Data
  #------------------------------------------------------------------------------#
  #get the end of the data
  data_t <- window(cbind(gdp.log, inflation, inflation.expectations, interest.rate))
  quarters <- nrow(data_t)
  data.end.year <- data.start.year + quarters %/% 4
  data.end.quarter <- (data.start.quarter + quarters %% 4) - 1
  data.end <- c(data.end.year, data.end.quarter)

  message("Dateset spans from ", data.start.year, " Q", data.start.quarter,
          " to ", data.end.year, " Q", data.end.quarter)

  # create date column sequence
  date <- seq(from = (as.Date(ti(shiftQuarter(data.start,-1),'quarterly'))+1),
              to = (as.Date(ti(shiftQuarter(data.end,-1),tif='quarterly'))+1),
              by = 'quarter')
  # merge time series into a dataframe
  data.out <- (window(cbind(gdp.log, inflation, inflation.expectations, interest.rate), start = data.start, end = data.end))
  data.out <- cbind(date, as.data.frame(data.out))
  write.table(data.out, file = paste0(folder, file), sep = ',',
              col.names = TRUE, quote = FALSE, na = '.', row.names = FALSE)
}

# unpack.parameters.stage1
unpack.parameters.stage1 <- function(parameters, y.data, x.data, xi.00, P.00) {
  A         <- matrix(0, 4, 2)
  A[1:2, 1] <- parameters[1:2]  # a_y,1, a_y,2
  A[1, 2]   <- parameters[4]    # b_y
  A[3, 2]   <- parameters[3]    # b_pi
  A[4, 2]   <- 1-parameters[3]  # 1 - b_pi

  H         <- matrix(0, 3, 2)
  H[1, 1]   <- 1
  H[2:3, 1] <- -parameters[1:2] # -a_y,1, -a_y,2
  H[2, 2]   <- -parameters[4]   # -b_y

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) # sigma_y~, sigma_pi
  Q         <- matrix(0, 3, 3)
  Q[1, 1]   <- parameters[8]^2  # sigma_y*

  F_matrix <- matrix(0, 3, 3)
  F_matrix[1, 1] <- F_matrix[2, 1] <- F_matrix[3, 2] <- 1

  # Make the data stationary
  y.data[, 1] <- y.data[, 1] - 1:dim(y.data)[1] * parameters[5] # g
  x.data[, 1] <- x.data[, 1] - 0:(dim(x.data)[1]-1) * parameters[5]
  x.data[, 2] <- x.data[, 2] - -1:(dim(x.data)[1]-2) * parameters[5]

  return(list("xi.00"=xi.00, "P.00"=P.00, "F_matrix"=F_matrix, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

# unpack.parameters.stage2
unpack.parameters.stage2 <- function(parameters, y.data, x.data, lambda.g, xi.00=NA, P.00=NA) {
  A         <- matrix(0, 2, 7)
  A[1, 1:2] <- parameters[1:2]   # a_y,1, a_y,2
  A[1, 3:4] <- parameters[3]/2   # a_r/2
  A[1, 7]   <- parameters[4]     # a_0
  A[2, 1]   <- parameters[7]     # b_y
  A[2, 5]   <- parameters[6]     # b_pi
  A[2, 6]   <- 1 - parameters[6] # 1 - b_pi
  A         <- t(A)

  H         <- matrix(0, 2, 4)
  H[1, 1  ] <- 1
  H[1, 2:3] <- -parameters[1:2] # -a_y,1, -a_y,2
  H[1, 4  ] <- parameters[5]    # a_g
  H[2, 2]   <- -parameters[7]   # -b_y
  H         <- t(H)

  R         <- diag(c(parameters[8]^2, parameters[9]^2)) # sigma_y~, sigma_pi
  Q         <- matrix(0, 4, 4)
  Q[1, 1]   <- parameters[10]^2              # sigma_y*
  Q[4, 4]   <- (lambda.g * parameters[10])^2 # sigma_y*

  F_matrix <- matrix(0, 4, 4)
  F_matrix[1, 1] <- F_matrix[1, 4] <- F_matrix[2, 1] <- F_matrix[3, 2] <- F_matrix[4,4] <- 1

  return(list("xi.00"=xi.00, "P.00"=P.00, "F_matrix"=F_matrix, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

# unpack.parameters.stage3
unpack.parameters.stage3 <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00, P.00) {
  A         <- matrix(0, 2, 6)
  A[1, 1:2] <- parameters[1:2]   # a_y,1, a_y,2
  A[1, 3:4] <- parameters[3]/2   # a_r/2
  A[2, 1]   <- parameters[5]     # b_y
  A[2, 5]   <- parameters[4]     # b_pi
  A[2, 6]   <- 1 - parameters[4] # 1 - b_pi
  A         <- t(A)


  H         <- matrix(0, 2, 7)
  H[1, 1]   <- 1
  H[1, 2:3] <- -parameters[1:2]       # a_y,1, a_y,2
  H[1, 4:5] <- -parameters[3] * 2     # -a_r/2 (annualized)
  H[1, 6:7] <- -parameters[3]/2       # -a_r/2
  H[2, 2]   <- -parameters[5]         # -b_y
  H         <- t(H)

  R         <- diag(c(parameters[6]^2, parameters[7]^2)) # sigma_y~, sigma_pi

  Q         <- matrix(0, 7, 7)
  Q[1, 1]   <- (1+lambda.g^2)*parameters[8]^2                  # sigma_y*
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4]<- (lambda.g*parameters[8])^2 # sigma_y*
  Q[6, 6]   <- (lambda.z*parameters[6]/parameters[3])^2        # sigma_y~/a_r

  F_matrix <- matrix(0, 7, 7)
  F_matrix[1, 1] <- F_matrix[1, 4] <- F_matrix[2, 1] <- F_matrix[3, 2] <- F_matrix[4,4] <- F_matrix[5,4]<- F_matrix[6,6] <- F_matrix[7,6] <- 1

  return(list("xi.00"=xi.00, "P.00"=P.00, "F_matrix"=F_matrix, "Q"=Q, "A"=A, "H"=H, "R"=R, "x.data"=x.data, "y.data"=y.data))
}

kalman.states <- function(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x) {
  filtered <- kalman.states.filtered(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x)
  smoothed <- kalman.states.smoothed(filtered$xi.ttm1, filtered$P.ttm1, filtered$xi.tt, filtered$P.tt,
                                     F_matrix, Q, A, H, R, y, x)
  return(list("filtered"=filtered, "smoothed"=smoothed))
}

# kalman.states.filtered
kalman.states.filtered <- function(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x, t=1) {
  xi.ttm1 <- as.vector(F_matrix %*% xi.tm1tm1)
  P.ttm1 <- F_matrix %*% P.tm1tm1 %*% t(F_matrix) + Q
  prediction.error <- (as.vector(y[t,]) - as.vector(t(A) %*% as.vector(x[t,])) - as.vector(t(H) %*% xi.ttm1))
  HPHR <- t(H) %*% P.ttm1 %*% H + R
  xi.tt <- xi.ttm1 + as.vector(P.ttm1 %*% H %*% solve(HPHR, prediction.error, tol = 0))
  P.tt <- P.ttm1 - P.ttm1 %*% H %*% solve(HPHR, t(H) %*% P.ttm1, tol = 0)
  if (t == dim(y)[1]) {
    return(list("xi.ttm1"=xi.ttm1, "P.ttm1"=P.ttm1, "xi.tt"=xi.tt, "P.tt"=P.tt))
  } else {
    tmp <- kalman.states.filtered(xi.tt, P.tt, F_matrix, Q, A, H, R, y, x, t+1)
    return(list("xi.ttm1"=rbind(xi.ttm1, tmp$xi.ttm1),
                "P.ttm1"=rbind(P.ttm1, tmp$P.ttm1),
                "xi.tt"=rbind(xi.tt, tmp$xi.tt),
                "P.tt"=rbind(P.tt, tmp$P.tt)))
  }
}

# kalman.states.smoothed
kalman.states.smoothed <- function(xi.ttm1.array, P.ttm1.array, xi.tt.array, P.tt.array,
                                   F_matrix, Q, A, H, R, y, x, t=dim(y)[1], xi.tp1T=NA, P.tp1T=NA) {
  n <- dim(xi.ttm1.array)[2]
  if (t == dim(y)[1]) {
    xi.tT <- xi.tt.array[t,]
    P.tT <- P.tt.array[((t-1)*n+1):(t*n),]
    tmp <- kalman.states.smoothed(xi.ttm1.array, P.ttm1.array, xi.tt.array, P.tt.array,
                                  F_matrix, Q, A, H, R, y, x, t-1, xi.tT, P.tT)
    return(list("xi.tT"=rbind(tmp$xi.tT, xi.tT),
                "P.tT" =rbind(tmp$P.tT, P.tT)))
  } else {
    P.tt <- P.tt.array[((t-1)*n+1):(t*n),]
    P.tp1t <- P.ttm1.array[(t*n+1):((t+1)*n),]

    #added for debug
    # P.tp1t matrix is singular. And the source of my pain
    # #-------------------------------------------------------------------------
    if (!matrixcalc::is.singular.matrix(P.tp1t)) {
      P.tp1t <- Matrix::nearPD(P.tp1t)$mat
      print("P.tp1t matrix singular. Found closest PD")
    }

    J.t <- P.tt %*% t(F_matrix) %*% solve(P.tp1t, tol = 0)
    #-------------------------------------------------------------------------
    xi.tt <- xi.tt.array[t,]
    xi.tp1t <- xi.ttm1.array[t+1,]
    xi.tT <- xi.tt + as.vector(J.t %*% (xi.tp1T - xi.tp1t))
    P.tT <- P.tt + J.t %*% (P.tp1T - P.tp1t) %*% t(J.t)
    if (t > 1) {
      tmp <- kalman.states.smoothed(xi.ttm1.array, P.ttm1.array, xi.tt.array, P.tt.array,
                                    F_matrix, Q, A, H, R, y, x, t-1, xi.tT, P.tT)
      return(list("xi.tT"=rbind(tmp$xi.tT, xi.tT),
                  "P.tT" =rbind(tmp$P.tT, P.tT)))
    } else {
      return(list("xi.tT"=xi.tT, "P.tT"=P.tT))
    }
  }
}

# kalman.states
kalman.states <- function(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x) {
  filtered <- kalman.states.filtered(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x)
  smoothed <- kalman.states.smoothed(filtered$xi.ttm1, filtered$P.ttm1, filtered$xi.tt, filtered$P.tt,
                                     F_matrix, Q, A, H, R, y, x)
  return(list("filtered"=filtered, "smoothed"=smoothed))
}

# kalman.states.filtered
kalman.states.filtered <- function(xi.tm1tm1, P.tm1tm1, F_matrix, Q, A, H, R, y, x, t=1) {
  xi.ttm1 <- as.vector(F_matrix %*% xi.tm1tm1)
  P.ttm1 <- F_matrix %*% P.tm1tm1 %*% t(F_matrix) + Q
  prediction.error <- (as.vector(y[t,]) - as.vector(t(A) %*% as.vector(x[t,])) - as.vector(t(H) %*% xi.ttm1))
  HPHR <- t(H) %*% P.ttm1 %*% H + R
  xi.tt <- xi.ttm1 + as.vector(P.ttm1 %*% H %*% solve(HPHR, prediction.error, tol = 0))
  P.tt <- P.ttm1 - P.ttm1 %*% H %*% solve(HPHR, t(H) %*% P.ttm1, tol = 0)
  if (t == dim(y)[1]) {
    return(list("xi.ttm1"=xi.ttm1, "P.ttm1"=P.ttm1, "xi.tt"=xi.tt, "P.tt"=P.tt))
  } else {
    tmp <- kalman.states.filtered(xi.tt, P.tt, F_matrix, Q, A, H, R, y, x, t+1)
    return(list("xi.ttm1"=rbind(xi.ttm1, tmp$xi.ttm1),
                "P.ttm1"=rbind(P.ttm1, tmp$P.ttm1),
                "xi.tt"=rbind(xi.tt, tmp$xi.tt),
                "P.tt"=rbind(P.tt, tmp$P.tt)))
  }
}

kalman.standard.errors <- function(T, states, theta, y.data, x.data, stage,
                                   lambda.g, lambda.z, xi.00, P.00, niter = niter,
                                   a3.constraint=NA, b2.constraint=NA) {
  message('Computing Standard Errors')

  # Set a3.constraint to -0.0025 if a constraint is not specified in stage 3
  if (is.na(a3.constraint)) {
    a3.constraint <- -0.0025
  }
  # Set b2.constraint to 0.025 if a constraint is not specified in stage 3
  if (is.na(b2.constraint)) {
    b2.constraint <- 0.025
  }

  message("Standard Error Procedure: a3.constraint")
  message(a3.constraint)

  message("Standard Error Procedure: b2.constraint")
  message(b2.constraint)

  n.params <- length(theta)
  n.state.vars <- length(xi.00)

  # Return vector of log likelihood values at each time t
  log.likelihood.estimated.vector <- log.likelihood.wrapper(theta, y.data, x.data, stage = 3,lambda.g=lambda.g, lambda.z=lambda.z, xi.00=xi.00, P.00=P.00)$ll.vec
  stin <- states$smoothed$xi.tT[1,] # First smoothed state vector
  pp1  <- states$filtered$P.ttm1[1:n.state.vars,] # First covariance matrix
  eigenstuff.pp1   <- eigen(pp1)
  eigenvectors.pp1 <- eigenstuff.pp1$vectors # Eigenvectors of first covariance matrix
  # Eigenvectors without a positive first entry are multiplied by -1 to ensure
  # consistency across different versions of R, which choose the sign differently
  for (l in 1:n.state.vars) {
    if (eigenvectors.pp1[1,l] < 0 ) { eigenvectors.pp1[,l] <- -eigenvectors.pp1[,l] }
  }
  eigenvalues.pp1  <- eigenstuff.pp1$value   # Eigenvalues of first covariance matrix
  dg   <- diag(x = eigenvalues.pp1)
  hh2  <- eigenvectors.pp1 %*% sqrt(dg)

  # Compute information matrix from difference in gradients of the likelihood function
  # from varying theta (parameter vector) values
  likelihood.gradient <- matrix(NA,T,n.params)
  for (i in 1:n.params){
    delta   <- max(theta[i]*1e-6, 1e-6)
    d.theta <- theta
    d.theta[i] <- theta[i] + delta
    likelihood.gradient[,i] <-  (log.likelihood.wrapper(d.theta, y.data, x.data, stage = 3,lambda.g=lambda.g, lambda.z=lambda.z, xi.00=xi.00, P.00=P.00)$ll.vec -
                                   log.likelihood.estimated.vector)/delta
  }
  info <- solve(t(likelihood.gradient) %*% likelihood.gradient, tol = 0) # Information matrix
  bse <- sqrt(diag(info))
  t.stats <- abs(theta) / bse

  # Smoothed estimates
  g      <- 4 * states$smoothed$xi.tT[,4]
  ypot   <- states$smoothed$xi.tT[,1]
  z      <- states$smoothed$xi.tT[,6]
  rstar  <- g + z

  # cum1 cumulates terms for parameter uncertainty;
  # cum2 cumulates terms for filter uncertainty
  cum1 <- matrix(0,T,3)
  cum2 <- matrix(0,T,3)
  eigenstuff.info   <- eigen(info)
  eigenvectors.info <- eigenstuff.info$vectors # Eigenvectors of information matrix
  # Eigenvectors without a positive first entry are multiplied by -1 to ensure
  # consistency across different versions of R, which choose the sign differently
  for (l in 1:n.params) {
    if (eigenvectors.info[1,l] < 0 ) { eigenvectors.info[,l] <- -eigenvectors.info[,l] }
  }
  eigenvalues.info  <- eigenstuff.info$value # Eigenvalues of information matrix
  dg <- diag(x = eigenvalues.info)
  hh <- eigenvectors.info %*% sqrt(dg)

  set.seed(50)

  # Store the number of draws excluded for violating constraints
  good.draws                 <- 0
  excluded.draw.counter      <- 0
  excluded.draw.counter.a3   <- 0
  excluded.draw.counter.b2   <- 0
  excluded.draw.counter.a1a2 <- 0

  # See HLW footnote 7 for description of procedure
  # niter is the number of iterations; we discard draws that violate constraints
  while (good.draws < niter) {
    theta.i <- (hh %*% rnorm(n.params) + theta)[,1]
    if ( (theta.i[3] <= a3.constraint) & (theta.i[5] >= b2.constraint) & (theta.i[1] + theta.i[2] < 1) ) {
      xi.00.i  <- c(t(hh2 %*% rnorm(n.state.vars) + stin))
      states.i <- kalman.states.wrapper(theta.i, y.data, x.data, stage, lambda.g, lambda.z, xi.00.i, pp1)

      g.i    <- 4 * states.i$smoothed$xi.tT[,4]
      ypot.i <- states.i$smoothed$xi.tT[,1]
      z.i    <- states.i$smoothed$xi.tT[,6]
      r.i    <- g.i + z.i

      cum1[,1] <- cum1[,1]+(ypot.i-ypot)^2
      cum1[,2] <- cum1[,2]+(r.i-rstar)^2
      cum1[,3] <- cum1[,3]+(g.i-g)^2

      P.ttm1.i   <- states.i$smoothed$P.tT
      P.ttm1.i.f <- states.i$filtered$P.tt
      for (j in 1:(T-1)){
        cum2[j,1]  <- cum2[j,1] + P.ttm1.i[(j * n.state.vars +1),1]
        cum2[j,2]  <- cum2[j,2] + 16 * P.ttm1.i[(j*n.state.vars+4),4] + P.ttm1.i[(j*n.state.vars+6),6]
        cum2[j,3]  <- cum2[j,3] + P.ttm1.i[(j*n.state.vars+4),4]
      }
      cum2[T,1] <- cum2[T,1] + P.ttm1.i.f[((T-1)*n.state.vars+1),1]
      cum2[T,2] <- cum2[T,2] + (16 * P.ttm1.i.f[((T-1)*n.state.vars+4),4] + P.ttm1.i.f[((T-1)*n.state.vars+6),6])
      cum2[T,3] <- cum2[T,3] + P.ttm1.i.f[((T-1)*n.state.vars+4),4]
      good.draws <- good.draws + 1
    } else {
      excluded.draw.counter <- excluded.draw.counter + 1
      if (theta.i[3] > a3.constraint) {
        excluded.draw.counter.a3 <- excluded.draw.counter.a3 + 1
      }
      if (theta.i[5] < b2.constraint) {
        excluded.draw.counter.b2 <- excluded.draw.counter.b2 + 1
      }
      if ((theta.i[1] + theta.i[2]) >= 1) {
        excluded.draw.counter.a1a2 <- excluded.draw.counter.a1a2 + 1
      }
    }
  }
  cum1 <- cum1/niter # Measure of parameter uncertainty
  cum2 <- cum2/niter # Measure of filter uncertainty
  cum2[,3] <- 16*cum2[,3] # Variance for growth at an annualized rate

  # Standard errors for estimates of the states
  # Order: y*, r*, g
  se <- sqrt(cum1 + cum2)

  #rm(.Random.seed)
  rm(.Random.seed, envir=.GlobalEnv)
  return(list("se.mean"=colMeans(se),
              "se"=se,"t.stats"=t.stats,"bse"=bse,
              "number.excluded"=excluded.draw.counter,
              "number.excluded.a3"=excluded.draw.counter.a3,
              "number.excluded.b2"=excluded.draw.counter.b2,
              "number.excluded.a1a2"=excluded.draw.counter.a1a2))
}

log.likelihood.wrapper <- function(parameters, y.data, x.data, stage = NA,
                                   lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){

  if (stage == 1) {
    out <- unpack.parameters.stage1(parameters, y.data, x.data,
                                    xi.00, P.00)
  } else if (stage == 2) {
    out <- unpack.parameters.stage2(parameters, y.data, x.data,
                                    lambda.g, xi.00, P.00)
  } else if (stage == 3) {
    out <- unpack.parameters.stage3(parameters, y.data, x.data,
                                    lambda.g, lambda.z, xi.00, P.00)
  } else {
    stop('You need to enter a stage number in log.likelihood.wrapper.')
  }


  for (n in names(out)) {
    eval(parse(text=paste0(n, "<-out$", n)))
  }
  return(kalman.log.likelihood(xi.00, P.00, F_matrix, Q, A, H, R, y.data, x.data))
}


kalman.states.wrapper <- function(parameters, y.data, x.data, stage = NA,
                                  lambda.g=NA, lambda.z=NA, xi.00=NA, P.00=NA){
  if (stage == 1) {
    out <- unpack.parameters.stage1(parameters, y.data, x.data,
                                    xi.00, P.00)
  } else if (stage == 2) {
    out <- unpack.parameters.stage2(parameters, y.data, x.data,
                                    lambda.g, xi.00, P.00)
  } else if (stage == 3) {
    out <- unpack.parameters.stage3(parameters, y.data, x.data,
                                    lambda.g, lambda.z, xi.00, P.00)
  } else {
    stop('You need to enter a stage number in kalman.states.wrapper.')
  }

  for (n in names(out)) {
    eval(parse(text=paste0(n, "<-out$", n)))
  }
  T <- dim(y.data)[1]
  states <- kalman.states(xi.00, P.00, F_matrix, Q, A, H, R, y.data, x.data)
  if (stage == 1) {
    states$filtered$xi.tt <- states$filtered$xi.tt + cbind(1:T,0:(T-1),-1:(T-2)) * parameters[5]
    states$smoothed$xi.tT <- states$smoothed$xi.tT + cbind(1:T,0:(T-1),-1:(T-2)) * parameters[5]
  }
  return(states)
}

rstar.stage1 <- function(log.output,
                         inflation,
                         b2.constraint=NA,
                         g.pot.start.index = g.pot.start.index) {
  message("Running rstar.stage1...")
  stage <- 1

  # Data must start 4 quarters before the estimation period
  T <- length(log.output) - 4

  # Original output gap estimate
  x.og <- cbind(rep(1,T+4), 1:(T+4))
  y.og <- log.output
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og, tol = 0)) * 100

  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  xi.00 <- c(100*g.pot[3:1])

  # IS curve
  y.is <- output.gap[5:(T+4)]
  x.is <- cbind(output.gap[4:(T+3)], output.gap[3:(T+2)])
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is, tol = 0)
  r.is <- y.is - x.is %*% b.is
  s.is <- sqrt(sum(r.is^2) / (length(r.is)-(dim(x.is)[2])))

  # Phillips curve
  y.ph <- inflation[5:(T+4)]
  x.ph <- cbind(inflation[4:(T+3)],
                (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                output.gap[4:(T+3)])
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph, tol = 0)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph)-(dim(x.ph)[2])))

  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)

  # Starting values for the parameter vector
  initial.parameters <- c(b.is, b.ph[1], b.ph[3], 0.85, s.is, s.ph, 0.5)

  # Set an upper and lower bound on the parameter vectors:
  # The vector is unbounded unless values are otherwise specified
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))

  # Set a lower bound for the Phillips curve slope (b_2) of b2.constraint, if not NA
  # In HLW, b2.constraint = 0.025
  if (!is.na(b2.constraint)) {
    if (initial.parameters[4] < b2.constraint) {
      initial.parameters[4] <- b2.constraint
    }
    theta.lb[4] <- b2.constraint
  }

  # Set the initial covariance matrix (see footnote 6)
  P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, NA, NA, xi.00)

  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage, NA, NA, xi.00, P.00)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))
  theta <- nloptr.out$solution

  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, NA, NA, xi.00, P.00)$ll.cum

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- kalman.states.wrapper(theta, y.data, x.data, stage, NA, NA, xi.00, P.00)

  # One-sided (filtered) estimates
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)

  # Two-sided (smoothed) estimates
  potential.smoothed  <- as.vector(states$smoothed$xi.tT[,1])/100
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100)

  # Save variables to return
  return.list                <- list()
  return.list$theta          <- theta
  return.list$log.likelihood <- log.likelihood
  return.list$states         <- states
  return.list$xi.00          <- xi.00
  return.list$P.00           <- P.00
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return(return.list)
}

median.unbiased.estimator.stage1 <- function(series) {
  message("Running median.unbiased.estimator.stage1...")
  T <- length(series)
  y <- 400 * diff(series)

  stat <- rep(T-2*4)
  for (i in 4:(T-5)) {
    xr <- cbind(rep(1, T-1), c(rep(0,i),rep(1,T-i-1)))
    xi <- solve(t(xr) %*% xr, tol = 0)
    b  <- solve(t(xr) %*% xr, t(xr) %*% y, tol = 0)
    s3 <- sum((y-xr%*%b)^2)/(T-2-1)
    stat[i+1-4] = b[2]/sqrt(s3*xi[2,2])
  }

  ew <- 0
  for (i in 1:length(stat)) {
    ew <- ew+exp(stat[i]^2/2)
  }
  ew  <- log(ew/length(stat))
  mw  <- sum(stat^2) / length(stat)
  qlr <- max(stat^2)

  # Values are from Table 3 in Stock and Watson (1998)
  # Test Statistic: Exponential Wald (EW)
  valew <- c(0.426, 0.476, 0.516, 0.661, 0.826, 1.111,
             1.419, 1.762, 2.355, 2.91,  3.413, 3.868, 4.925,
             5.684, 6.670, 7.690, 8.477, 9.191, 10.693, 12.024,
             13.089, 14.440, 16.191, 17.332, 18.699, 20.464,
             21.667, 23.851, 25.538, 26.762, 27.874)
  # Test Statistic: Mean Wald (MW)
  valmw <- c(0.689, 0.757, 0.806, 1.015, 1.234, 1.632,
             2.018, 2.390, 3.081, 3.699, 4.222, 4.776, 5.767,
             6.586, 7.703, 8.683, 9.467, 10.101, 11.639, 13.039,
             13.900, 15.214, 16.806, 18.330, 19.020, 20.562,
             21.837, 24.350, 26.248, 27.089, 27.758)
  # Test Statistic: QLR
  valql <- c(3.198, 3.416, 3.594, 4.106, 4.848, 5.689,
             6.682, 7.626, 9.16,  10.66, 11.841, 13.098, 15.451,
             17.094, 19.423, 21.682, 23.342, 24.920, 28.174, 30.736,
             33.313, 36.109, 39.673, 41.955, 45.056, 48.647, 50.983,
             55.514, 59.278, 61.311, 64.016)

  lame <- NA
  lamm <- NA
  lamq <- NA

  # Median-unbiased estimator of lambda_g for given values of the test
  # statistics are obtained using the procedure described in the
  # footnote to Stock and Watson (1998) Table 3.
  if (ew <= valew[1]) {
    lame <- 0
  } else {
    for (i in 1:(length(valew)-1)) {
      if ((ew > valew[i]) & (ew <= valew[i+1])) {
        lame <- i-1+(ew-valew[i])/(valew[i+1]-valew[i])
      }
    }
  }
  if (mw <= valmw[1]) {
    lamm <- 0
  } else {
    for (i in 1:(length(valmw)-1)) {
      if ((mw > valmw[i]) & (mw <= valmw[i+1])) {
        lamm <- i-1+(mw-valmw[i])/(valmw[i+1]-valmw[i])
      }
    }
  }
  if (qlr <= valql[1]) {
    lamq <- 0
  } else {
    for (i in 1:(length(valql)-1)) {
      if ((qlr > valql[i]) & (qlr <= valql[i+1])) {
        lamq <- i-1+(qlr-valql[i])/(valql[i+1]-valql[i])
      }
    }
  }
  if (is.na(lame) | is.na(lamm) | is.na(lamq)) {
    message("At least one statistic has an NA value. Check to see if your EW, MW, and/or QLR value is outside of Table 3.")
  }

  stats <- c(ew, mw, qlr)
  lams  <- c(lame, lamm, lamq)
  return(lame/(T-1))
}

rstar.stage2 <- function(log.output,
                         inflation,
                         real.interest.rate,
                         lambda.g,
                         a3.constraint=NA,
                         b2.constraint=NA,
                         g.pot.start.index = g.pot.start.index) {
  message("Running rstar.stage2...")
  stage <- 2

  # Data must start 4 quarters before the estimation period
  T <- length(log.output) - 4

  # Original output gap estimate
  x.og <- cbind(rep(1,T+4), 1:(T+4))
  y.og <- log.output
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og, tol = 0)) * 100

  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)
  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2])

  # IS curve
  y.is <- output.gap[5:(T+4)]
  x.is <- cbind(output.gap[4:(T+3)], output.gap[3:(T+2)],
                (real.interest.rate[4:(T+3)] + real.interest.rate[3:(T+2)])/2,
                rep(1,T))
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is, tol = 0)
  r.is <- as.vector(y.is - x.is %*% b.is)
  s.is <- sqrt(sum(r.is^2) / (length(r.is)-(dim(x.is)[2])))

  # Phillips curve
  y.ph <- inflation[5:(T+4)]
  x.ph <- cbind(inflation[4:(T+3)],
                (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                output.gap[4:(T+3)])
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph, tol = 0)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph)-(dim(x.ph)[2])))

  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                  rep(1,T))

  # Starting values for the parameter vector
  initial.parameters <- c(b.is, -b.is[3], b.ph[1], b.ph[3], s.is, s.ph, 0.5)

  # Set an upper and lower bound on the parameter vectors:
  # The vector is unbounded unless values are otherwise specified
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))

  # Set a lower bound for the Phillips curve slope (b_2) of b2.constraint, if not NA
  # In HLW, b2.constraint = 0.025
  if (!is.na(b2.constraint)) {
    if (initial.parameters[7] < b2.constraint) {
      initial.parameters[7] <- b2.constraint
    }
    theta.lb[7] <- b2.constraint
  }

  # Set an upper bound for the IS curve slope (a_3) of a3.constraint, if not NA
  # In HLW, a3.constraint = -0.0025
  if (!is.na(a3.constraint)) {
    if (initial.parameters[3] > a3.constraint) {
      initial.parameters[3] <- a3.constraint
    }
    theta.ub[3] <- a3.constraint
  }

  # Set the initial covariance matrix (see footnote 6)
  P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, NA, xi.00)

  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))
  theta <- nloptr.out$solution

  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)$ll.cum

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, NA, xi.00, P.00)

  # Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tt[,4] * 4
  potential.smoothed  <- c(states$smoothed$xi.tT[1, 3:2], states$smoothed$xi.tT[,1])
  output.gap.smoothed <- 100 * log.output[3:(T+4)] - potential.smoothed

  # Inputs for median.unbiased.estimator.stage2.R
  y <- output.gap.smoothed[3:length(output.gap.smoothed)]
  x <- cbind(output.gap.smoothed[2:(length(output.gap.smoothed)-1)],
             output.gap.smoothed[1:(length(output.gap.smoothed)-2)],
             (x.data[,3]+x.data[,4])/2,
             states$smoothed$xi.tT[,4],
             rep(1,T))

  # One-sided (filtered) estimates
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)

  # Save variables to return
  return.list <- list()
  return.list$y              <- y
  return.list$x              <- x
  return.list$theta          <- theta
  return.list$log.likelihood <- log.likelihood
  return.list$states         <- states
  return.list$xi.00          <- xi.00
  return.list$P.00           <- P.00
  return.list$trend.filtered      <- trend.filtered
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$trend.smoothed      <- trend.smoothed
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return(return.list)
}

median.unbiased.estimator.stage2 <- function(y, x) {
  message("Running median.unbiased.estimator.stage2...")
  T <- dim(x)[1]
  stat <- rep(0, T-2*4+1)
  for (i in 4:(T-4)) {
    xr <- cbind(x, c(rep(0, i), rep(1, T-i)))
    xi <- solve(t(xr)%*%xr, tol = 0)
    b <- solve(t(xr)%*%xr,t(xr)%*%y, tol = 0)
    s3 <- sum((y-xr%*%b)^2)/(T-dim(xr)[2])
    stat[i+1-4] <- b[dim(xr)[2]]/sqrt(s3*xi[dim(xr)[2],dim(xr)[2]])
  }
  ew <- 0
  for (i in 1:length(stat)) {
    ew <- ew+exp((stat[i]^2)/2)
  }
  ew <- log(ew/length(stat))
  mw <- mean(stat^2)
  qlr <- max(stat^2)

  # Values are from Table 3 in Stock and Watson (1998)
  # Test Statistic: Exponential Wald (EW)
  valew <- c(0.426, 0.476, 0.516, 0.661, 0.826, 1.111,
             1.419, 1.762, 2.355, 2.91,  3.413, 3.868, 4.925,
             5.684, 6.670, 7.690, 8.477, 9.191, 10.693, 12.024,
             13.089, 14.440, 16.191, 17.332, 18.699, 20.464,
             21.667, 23.851, 25.538, 26.762, 27.874)
  # Test Statistic: Mean Wald (MW)
  valmw <- c(0.689, 0.757, 0.806, 1.015, 1.234, 1.632,
             2.018, 2.390, 3.081, 3.699, 4.222, 4.776, 5.767,
             6.586, 7.703, 8.683, 9.467, 10.101, 11.639, 13.039,
             13.900, 15.214, 16.806, 18.330, 19.020, 20.562,
             21.837, 24.350, 26.248, 27.089, 27.758)
  # Test Statistic: QLR
  valql <- c(3.198, 3.416, 3.594, 4.106, 4.848, 5.689,
             6.682, 7.626, 9.16,  10.66, 11.841, 13.098, 15.451,
             17.094, 19.423, 21.682, 23.342, 24.920, 28.174, 30.736,
             33.313, 36.109, 39.673, 41.955, 45.056, 48.647, 50.983,
             55.514, 59.278, 61.311, 64.016)

  lame <- NA
  lamm <- NA
  lamq <- NA

  # Median-unbiased estimator of lambda_g for given values of the test
  # statistics are obtained using the procedure described in the
  # footnote to Stock and Watson (1998) Table 3.
  if (ew <= valew[1]) {
    lame <- 0
  } else {
    for (i in 1:(length(valew)-1)) {
      if ((ew > valew[i]) & (ew <= valew[i+1])) {
        lame <- i-1+(ew-valew[i])/(valew[i+1]-valew[i])
      }
    }
  }

  if (mw <= valmw[1]) {
    lamm <- 0
  } else {
    for (i in 1:(length(valmw)-1)) {
      if ((mw > valmw[i]) & (mw <= valmw[i+1])) {
        lamm <- i-1+(mw-valmw[i])/(valmw[i+1]-valmw[i])
      }
    }
  }

  if (qlr <= valql[1]) {
    lamq <- 0
  } else {
    for (i in 1:(length(valql)-1)) {
      if ((qlr > valql[i]) & (qlr <= valql[i+1])) {
        lamq <- i-1+(qlr-valql[i])/(valql[i+1]-valql[i])
      }
    }
  }
  if (is.na(lame) | is.na(lamm) | is.na(lamq)) {
    message("At least one statistic has an NA value. Check to see if your EW, MW, and/or QLR value is outside of Table 3.")
  }

  stats <- c(ew,mw,qlr)
  lams  <- c(lame,lamm,lamq)
  return(lame/T)
}

rstar.stage3 <- function(log.output,
                         inflation,
                         real.interest.rate,
                         nominal.interest.rate,
                         lambda.g,
                         lambda.z,
                         a3.constraint=NA,
                         b2.constraint=NA,
                         run.se = TRUE,
                         g.pot.start.index = g.pot.start.index,
                         niter = niter) {
  message("Running rstar.stage3...")
  stage <- 3

  # Data must start 4 quarters before the estimation period
  T <- length(log.output) - 4

  # Original output gap estimate
  x.og <- cbind(rep(1,T+4), 1:(T+4))
  y.og <- log.output
  output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og, tol = 0)) * 100

  # Initialization of state vector for Kalman filter using HP trend of log output
  log.output.hp.trend <- hpfilter(log.output,freq=36000,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)
  xi.00 <- c(100*g.pot[3:1],100*g.pot.diff[2:1],0,0)

  # IS curve
  y.is <- output.gap[5:(T+4)]
  x.is <- cbind(output.gap[4:(T+3)], output.gap[3:(T+2)],
                (real.interest.rate[4:(T+3)] + real.interest.rate[3:(T+2)])/2,
                rep(1,T))
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is, tol = 0)
  r.is <- as.vector(y.is - x.is %*% b.is)
  s.is <- sqrt(sum(r.is^2) / (length(r.is)-(dim(x.is)[2])))

  # Phillips curve
  y.ph <- inflation[5:(T+4)]
  x.ph <- cbind(inflation[4:(T+3)],
                (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3,
                output.gap[4:(T+3)])
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph, tol = 0)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (length(r.ph)-(dim(x.ph)[2])))

  y.data <- cbind(100 * log.output[5:(T+4)],
                  inflation[5:(T+4)])
  x.data <- cbind(100 * log.output[4:(T+3)],
                  100 * log.output[3:(T+2)],
                  real.interest.rate[4:(T+3)],
                  real.interest.rate[3:(T+2)],
                  inflation[4:(T+3)],
                  (inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/3)

  # Starting values for the parameter vector
  initial.parameters <- c(b.is[1:3], b.ph[1], b.ph[3], s.is, s.ph, 0.7)

  # Set an upper and lower bound on the parameter vectors:
  # The vector is unbounded unless values are otherwise specified
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))

  # Set a lower bound for the Phillips curve slope (b_2) of b2.constraint, if not NA
  # In HLW, b2.constraint = 0.025
  if (!is.na(b2.constraint)) {
    message("Setting a lower bound on b_2 of ", as.character(b2.constraint))
    if (initial.parameters[5] < b2.constraint) {
      initial.parameters[5] <- b2.constraint
    }
    theta.lb[5] <- b2.constraint
  }

  # Set an upper bound for the IS curve slope (a_3) of a3.constraint, if not NA
  # In HLW, a3.constraint = -0.0025
  if (!is.na(a3.constraint)) {
    message("Setting an upper bound on a_3 of ", as.character(a3.constraint))
    if (initial.parameters[3] > a3.constraint) {
      initial.parameters[3] <- a3.constraint
    }
    theta.ub[3] <- a3.constraint
  }

  # Set the initial covariance matrix (see footnote 6)
  P.00 <- calculate.covariance(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, lambda.z, xi.00)

  # Get parameter estimates via maximum likelihood
  f <- function(theta) {return(-log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_LBFGS","xtol_rel"=1.0e-8))
  theta <- nloptr.out$solution

  log.likelihood <- log.likelihood.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum

  # Get state vectors (xi.tt, xi.ttm1, xi.tT, P.tt, P.ttm1, P.tT) via Kalman filter
  states <- kalman.states.wrapper(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)

  # If run.se = TRUE, compute standard errors for estimates of the states (see footnote 7) and report run time
  if (run.se) {
    ptm <- proc.time()
    se <- kalman.standard.errors(T, states, theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00, niter, a3.constraint, b2.constraint)
    message("Standard error procedure run time")
    message(proc.time() - ptm)
  }

  # One-sided (filtered) estimates
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  z.filtered          <- states$filtered$xi.tt[,6]
  rstar.filtered      <- trend.filtered + z.filtered
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)

  # Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tt[,4] * 4
  z.smoothed          <- states$smoothed$xi.tt[,6]
  rstar.smoothed      <- trend.smoothed + z.smoothed
  potential.smoothed  <- states$smoothed$xi.tt[,1]/100
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100)

  # Save variables to return
  return.list <- list()
  return.list$rstar.filtered      <- rstar.filtered
  return.list$trend.filtered      <- trend.filtered
  return.list$z.filtered          <- z.filtered
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$rstar.smoothed      <- rstar.smoothed
  return.list$trend.smoothed      <- trend.smoothed
  return.list$z.smoothed          <- z.smoothed
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return.list$theta               <- theta
  return.list$log.likelihood      <- log.likelihood
  return.list$states              <- states
  return.list$xi.00               <- xi.00
  return.list$P.00                <- P.00
  return.list$y.data              <- y.data
  return.list$initial.parameters  <- initial.parameters
  if (run.se) { return.list$se    <- se }
  return(return.list)
}

run.hlw.estimation <- function(log.output, inflation, real.interest.rate, nominal.interest.rate,
                               a3.constraint = NA, b2.constraint = NA, run.se = TRUE,
                               g.pot.start.index = g.pot.start.index, niter = niter) {
  message("Running HLW ESTIMATION...")
  # Running the stage 1 model
  out.stage1 <- rstar.stage1(log.output,
                             inflation,
                             b2.constraint,
                             g.pot.start.index)

  # Median unbiased estimate of lambda_g
  lambda.g <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)

  # Running the stage 2 model
  out.stage2 <- rstar.stage2(log.output,
                             inflation,
                             real.interest.rate,
                             lambda.g,
                             a3.constraint,
                             b2.constraint,
                             g.pot.start.index)

  # Median unbiased estimate of lambda_z
  lambda.z <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x)

  # Running the stage 3 model
  out.stage3 <- rstar.stage3(log.output,
                             inflation,
                             real.interest.rate,
                             nominal.interest.rate,
                             lambda.g,
                             lambda.z,
                             a3.constraint,
                             b2.constraint,
                             run.se,
                             g.pot.start.index,
                             niter = niter)

  return(list(out.stage1=out.stage1,out.stage2=out.stage2,out.stage3=out.stage3,
              lambda.g=lambda.g,lambda.z=lambda.z))
}

hlw_nri <- function(input_file,
                    country_name,
                    output_folder,
                    niter = 5000,
                    run.se = TRUE) {


  # Read in input file for a country
  data <- read.table(input_file,
                     sep = ',',
                     na.strings = ".",
                     header=TRUE,
                     stringsAsFactors=FALSE) %>%
    mutate(date = as.Date(parse_date_time(date, c("Ymd", "mdY"))))

  # perform data quality checks
  #-------------------------------------------------------------------------
  # minimum quarters is 98 otherwise singular matrices occur
  message("Checking that data is at least 98 quarters...")
  if (nrow(data) < 98) {
    message("Your data is of length ", nrow(data),
            " please make sure it spans at least 98 quarters for the estimation to run.")
    stop("You should get a longer sample period")
  } else {
    message("Data is of length ", nrow(data), " continuing to estimation precedure")
  }

  # zeros or NAs break the code
  message("Please check yourself that the data doesn't contain non-sensible zeros...")
  message("Checking that data does not contain illegal values...")
  for (data_column in c(2,3,4,5)) {
    if (sum(is.na(data[data_column]))>5) {
      stop(paste0(colnames(data[data_column])), " contains ", sum(data[data_column]), " NAs")
    }
    if (sum(data[data_column] == 0)>0) {
      stop(paste0(colnames(data[data_column])), " contains ", sum(data[data_column]), " zeros")
    }
  }

  # provide additional info
  if(run.se == TRUE) {message("Running with standard error calculations")}
  if(run.se == FALSE) {message("Running without standard error calculations")}
  message("Iterations for Monte Carlo simulation are set to: ", niter, ". Default is 5000.")

  # Automatically set the start and end dates of the estimation sample
  #------------------------------------------------------------------------------
  # (format is c(year,quarter))
  data.start.quarter <- as.numeric(quarter(data$date[1]))
  data.start.year <- as.numeric(year(data$date[1]))
  data.start <- c(data.start.year, data.start.quarter)

  # The estimation process uses data beginning 4 quarters prior to the sample start
  sample.start    <- shiftQuarter(data.start,+4)

  data.end.quarter <- as.numeric(quarter(data$date[length(data$date)]))
  data.end.year <- as.numeric(year(data$date[length(data$date)]))
  sample.end   <- c(data.end.year, data.end.quarter)

  message("Data starts in ", data.start.year, " quarter ", data.start.quarter)
  message("Sample starts 4 quarters after data start")
  message("Data ends in ", data.end.year, " quarter ", data.end.quarter )

  #------------------------------------------------------------------------------#
  # Define variables
  #------------------------------------------------------------------------------#
  # Upper bound on a_3 parameter (slope of the IS curve)
  a3.constraint <- -0.0025

  # Lower bound on b_2 parameter (slope of the Phillips curve)
  b2.constraint <- 0.025


  # Set start index for y
  g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(data.start,'quarterly')

  # Set column names for CSV output
  output.col.names <- c("Date","rstar","g","z","output gap","","All results are output from the Stage 3 model.",rep("",8),"Standard Errors","Date","y*","r*","g","","rrgap")

  # Set number of iterations for Monte Carlo standard error procedure
  # original is 5000
  niter <- niter

  # Because the MC standard error procedure is time consuming, we include a run switch
  # Set run.se to TRUE to run the procedure
  run.se <- run.se

  # create variable arrays from data
  #--------------------------------------------------------------------------
  #the source of my pain. Using an index. the matrix is singular when logged.
  # but estimates are wonky when output is logged,


  log.output            <- as.numeric(data$gdp.log)
  inflation              <- as.numeric(data$inflation)
  # Inflation expectations measure: 4-quarter moving average of past inflation
  inflation.expectations <- as.numeric(data$inflation.expectations)
  nominal.interest.rate  <- as.numeric(data$interest.rate)
  real.interest.rate     <- as.numeric(nominal.interest.rate - inflation.expectations)


  # Run HLW estimation
  estimation <- run.hlw.estimation(log.output, inflation, real.interest.rate, nominal.interest.rate,
                                   a3.constraint = a3.constraint, b2.constraint = b2.constraint,
                                   run.se = run.se, g.pot.start.index = g.pot.start.index,
                                   niter = niter)

  # One-sided (filtered) estimates
  one.sided.est <- cbind(estimation$out.stage3$rstar.filtered,
                         estimation$out.stage3$trend.filtered,
                         estimation$out.stage3$z.filtered,
                         estimation$out.stage3$output.gap.filtered)

  # Two-sided (smoothed) estimates
  two.sided.est <- cbind(estimation$out.stage3$rstar.smoothed,
                         estimation$out.stage3$trend.smoothed,
                         estimation$out.stage3$z.smoothed,
                         estimation$out.stage3$output.gap.smoothed)

  # Save one-sided estimates to CSV
  #------------------------ change url to dynamically change
  write.table(one.sided.est, paste0(output_folder, 'one.sided.est.', country_name, '.csv'), row.names = FALSE, col.names = c("rstar","g","z","output gap"), quote = FALSE, sep = ',', na = ".")

  # Save two-sided estimates to CSV
  #------------------------ change url to dynamically change
  write.table(two.sided.est, paste0(output_folder, 'two.sided.est.', country_name, '.csv'), row.names = FALSE, col.names = c("rstar","g","z","output gap"), quote = FALSE, sep = ',', na = ".")

  # Save output to CSV
  #------------------------ change url to dynamically change
  output <- format.output(estimation, one.sided.est, real.interest.rate, sample.start, sample.end, run.se = run.se, niter = niter)
  write.table(output, paste0(output_folder, 'output.', country_name, '.csv'), col.names = output.col.names, quote=FALSE, row.names=FALSE, sep = ',', na = '')
  message("Done. You can find the output in: ", output_folder)
}

format.output <- function(country.estimation, one.sided.est.country,
                          real.rate.country, start, end, run.se = TRUE, niter=niter) {
    output.country <- data.frame(matrix(NA,dim(one.sided.est.country)[1],22))

    output.country[,1]   <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1),
                                to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1),
                                by = 'quarter')
    output.country[,2:5] <- one.sided.est.country

    output.country[1,7]    <- "Parameter Point Estimates"
    output.country[2,7:15] <- c("a_1","a_2","a_3","b_1","b_2",
                                "sigma_1","sigma_2","sigma_4","a_1 + a_2")
    output.country[3,7:14] <- country.estimation$out.stage3$theta
    output.country[3,15]   <- country.estimation$out.stage3$theta[1]
                            + country.estimation$out.stage3$theta[2]
    # Include standard errors in output only if run.se switch is TRUE
    if (run.se) {
        output.country[4,7]    <- "T Statistics"
        output.country[5,7:14] <- country.estimation$out.stage3$se$t.stats

        output.country[8,7]    <- "Average Standard Errors"
        output.country[9,7:9]  <- c("y*","r*","g")
        output.country[10,7:9] <- country.estimation$out.stage3$se$se.mean

        output.country[12,7]   <- "Restrictions on MC draws: a_3 < -0.0025; b_2 > 0.025; a_1 + a_2 < 1"
        output.country[13,7]   <- "Draws excluded:"; output.country[13,9] <-
          country.estimation$out.stage3$se$number.excluded
        output.country[13,10]  <- "Total:"; output.country[13,11] <- niter
        output.country[14,7]   <- "Percent excluded:"; output.country[14,9] <-
          as.numeric(output.country[13,9]) / (as.numeric(output.country[13,9]) + as.numeric(output.country[13,11]))
        output.country[15,7]   <- "Draws excluded because a_3 > -0.0025:"; output.country[15,11] <- country.estimation$out.stage3$se$number.excluded.a3
        output.country[16,7]   <- "Draws excluded because b_2 <  0.025:"; output.country[16,11] <- country.estimation$out.stage3$se$number.excluded.b2
        output.country[17,7]   <- "Draws excluded because a_1 + a_2 < 1:"; output.country[17,11] <- country.estimation$out.stage3$se$number.excluded.a1a2
    }

    output.country[19,7] <- "Signal-to-noise Ratios"
    output.country[20,7] <- "lambda_g"; output.country[20,8] <- country.estimation$lambda.g
    output.country[21,7] <- "lambda_z"; output.country[21,8] <- country.estimation$lambda.z
    output.country[19,11] <- "Log Likelihood"; output.country[20,11] <- country.estimation$out.stage3$log.likelihood

    output.country[24,7] <- "State vector: [y_{t}* y_{t-1}* y_{t-2}* g_{t-1} g_{t-2} z_{t-1} z_{t-2}]"
    output.country[25,7] <- "Initial State Vector"
    output.country[26,7:13] <- country.estimation$out.stage3$xi.00
    output.country[28,7] <- "Initial Covariance Matrix"
    output.country[29:35,7:13] <- country.estimation$out.stage3$P.00

    if (run.se) {
        output.country[,17]   <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
        output.country[,18:20] <- country.estimation$out.stage3$se$se
    }

    output.country[,22] <- real.rate.country[5:length(real.rate.country)] - country.estimation$out.stage3$rstar.filtered

    return(output.country)
}

input_file <-  "data/input/rstar.data.us.csv"  # location of input file
country_name <-  "US.test"                     #country name for output file
output_folder <-  "data/output/"                  # output folder for estimation
run.se <-  FALSE                                    # turn on standard errors
niter <-5000       #iterations for monte carlo simulation

hlw_nri(input_file, country_name, output_folder, run.se, niter)

