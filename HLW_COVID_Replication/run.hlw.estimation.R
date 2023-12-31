#------------------------------------------------------------------------------
# File:        run.hlw.estimation.R
#
# Description: Runs the three stages of the HLW estimation using country-specific
#              inputs and returns output from each stage.
#------------------------------------------------------------------------------#
run.hlw.estimation <- function(log.output, inflation, real.interest.rate, nominal.interest.rate,
                               covid.dummy, T.og0.omit = 0,
                               a3.constraint = NA, b2.constraint = NA, run.se = TRUE, sample.end) {


    # Running the stage 1 model
    out.stage1 <- rstar.stage1(log.output,
                                inflation,
                                covid.dummy,
                                T.og0.omit,
                                b2.constraint)

    # Median unbiased estimate of lambda_g
    lambda.g <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)

    # Running the stage 2 model
    out.stage2 <- rstar.stage2(log.output,
                               inflation,
                               real.interest.rate,
                               covid.dummy,
                               lambda.g,
                               T.og0.omit,
                               a3.constraint,
                               b2.constraint)

    # Median unbiased estimate of lambda_z
    lambda.z <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x)

    # Running the stage 3 model
    out.stage3 <- rstar.stage3(log.output,
                               inflation,
                               real.interest.rate,
                               nominal.interest.rate,
                               covid.dummy,
                               lambda.g,
                               lambda.z,
                               T.og0.omit,
                               a3.constraint,
                               b2.constraint,
                               run.se,
                               sample.end)

    return(list(out.stage1=out.stage1,out.stage2=out.stage2,out.stage3=out.stage3,
                lambda.g=lambda.g,lambda.z=lambda.z))
}
