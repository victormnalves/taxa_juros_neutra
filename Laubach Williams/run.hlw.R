#------------------------------------------------------------------------------#
# File:        run.hlw.R
#
# Description: This the main file for HLW, which does the following:
#              (1) Prepares data to be used in estimation
#              (2) Runs the three-stage HLW estimation for each economy
#             (3) Saves output.
#
# COVID-19 Update: This code was updated in May 2020 to adjust the HLW model
#              and estimation procedure to allow for estimation during and
#              after the COVID-19 pandemic. Please see the note on modifying
#              the HLW and LW models during the COVID-19 pandemic at
# https://www.newyorkfed.org/medialibrary/media/research/policy/rstar/LW_HLW_COVID_note
#------------------------------------------------------------------------------#
rm(list=ls())

#------------------------------------------------------------------------------#
# Prepare data to be used in estimation.
#
# Output will be saved in the inputData folder. Manually store the provided
# COVID-19 index in the inputData folder as well.
#
# IMPORTANT: Set the data start and end dates manually in each prepare.rstar.data file
# See the included code guide for instructions on downloading and storing raw data.
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
# Load required packages and source all programs to be used in HLW estimation.
#------------------------------------------------------------------------------#
if (!require("tis")) {install.packages("tis"); library("tis")} # Time series package
if (!require("mFilter")) {install.packages("mFilter"); library("mFilter")} # HP filter
if (!require("nloptr")) {install.packages("nloptr"); library("nloptr")} # Optimization

# Source all R programs; see code guide for details of each
source("calculate.covariance.R")
source("format.output.R")
source("kalman.log.likelihood.R")
source("kalman.standard.errors.R")
source("kalman.states.R")
source("kalman.states.wrapper.R")
source("log.likelihood.wrapper.R")
source("median.unbiased.estimator.stage1.R")
source("median.unbiased.estimator.stage2.R")
source("rstar.stage1.R")
source("rstar.stage2.R")
source("rstar.stage3.R")
source("run.hlw.estimation.R")
source("unpack.parameters.stage1.R")
source("unpack.parameters.stage2.R")
source("unpack.parameters.stage3.R")
source("utilities.R")

#------------------------------------------------------------------------------#
# Define variables
#------------------------------------------------------------------------------#

# Upper bound on a_3 parameter (slope of the IS curve)
a3.constraint <- -0.0025

# Lower bound on b_2 parameter (slope of the Phillips curve)
b2.constraint <- 0.25

# Set the start and end dates of the estimation sample (format is c(year,quarter))
sample.start <- c(2002,1)
sample.end   <- c(2023,4)

# Calculate number of quarters to omit in calculating intial parameters through 19q4
# See COVID-adjustment note for explanation
T.og0.omit.yq <-  sample.end - c(2019,4)
T.og0.omit    <-  max(T.og0.omit.yq[1]*4 + T.og0.omit.yq[2],0)

# The estimation process uses data beginning 4 quarters prior to the sample start
data.start    <- shiftQuarter(sample.start,-4)

# Set start index for y
g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(data.start,'quarterly')

# Set column names for CSV output
output.col.names = c("Date","rstar","g","z","output gap","","All results are output from the Stage 3 model.",rep("",9),"Standard Errors","Date","y*","r*","g","","rrgap","adj. output gap")

# Set number of iterations for Monte Carlo standard error procedure
niter <- 5000

# Because the MC standard error procedure is time consuming, we include a run switch
# Set run.se to TRUE to run the procedure
run.se <- TRUE

#------------------------------------------------------------------------------#
# United States: Read in data, run estimation, and save output
#------------------------------------------------------------------------------#
# Read in output of prepare.rstar.data.us.R
us.data <- readxl::read_xlsx('dados.xlsx')

us.log.output             <- us.data$ln_gdp
us.inflation              <- us.data$inflation
us.inflation.expectations <- us.data$expectations
us.nominal.interest.rate  <- us.data$rates
us.real.interest.rate     <- us.nominal.interest.rate - us.inflation.expectations

us.covid.dummy            <- us.data$covid

# Run HLW estimation for the US
us.estimation <- run.hlw.estimation(us.log.output, us.inflation, us.real.interest.rate, us.nominal.interest.rate,
                                    us.covid.dummy, T.og0.omit = T.og0.omit,
                                    a3.constraint = a3.constraint, b2.constraint = b2.constraint, run.se = run.se,
                                    sample.end)

# One-sided (filtered) estimates
one.sided.est.us <- cbind(us.estimation$out.stage3$rstar.filtered,
                          us.estimation$out.stage3$trend.filtered,
                          us.estimation$out.stage3$z.filtered,
                          us.estimation$out.stage3$output.gap.filtered)



# Save output to CSV
output.us <- format.output(us.estimation, one.sided.est.us, us.real.interest.rate, sample.start, sample.end, run.se = run.se)
write.table(one.sided.est.us, 'output/output.us.csv', quote=FALSE, row.names=FALSE, sep = ',', na = '')
save.image('output/outdata.us.RData')

# Save one-sided estimates to CSV
#write.table(one.sided.est.us, 'output/one.sided.est.us.csv', row.names = FALSE, col.names = c("rstar","g","z","output gap"), quote = FALSE, sep = ',', na = ".")

