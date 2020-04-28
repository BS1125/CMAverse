get_coef <- function(formulas, regressions, model, yreg, mreg, data) {

  if (model == "rb") {

    mediator_regression <- regressions$mediator_regression[[1]]

    outcome_regression <- regressions$outcome_regression

    ## Store coefficients from regressions
    betas  <- as.vector(t(coefficients(mediator_regression)))

    thetas <- coefficients(outcome_regression)

    ## Store covariances from regressions

    vcov_betas <- vcov(mediator_regression)

    vcov_thetas <- vcov(outcome_regression)

    if (mreg == "linear") { variance <- sigma(mediator_regression)^2
    } else variance <- NULL

    if (yreg == "aft_weibull")
      vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

    ## Build block diagonal matrix
    vcov_block <- Matrix::bdiag(vcov_thetas, vcov_betas)

    coef <- list(betas = betas, thetas = thetas,
                 vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
                 variance = variance)

  } else if (model == "iorw") {

    tot_regression = regressions$tot_regression
    dir_regression = regressions$dir_regression

    if (is.factor(data[, exposure])) {

      dir_coef <- unname(coef(dir_regression)[2+(0:(length(levels(data[, exposure]))-2))])

      tot_coef <- unname(coef(tot_regression)[2+(0:(length(levels(data[, exposure]))-2))])

    } else {

      dir_coef <- unname(coef(dir_regression)[2])

      tot_coef <- unname(coef(tot_regression)[2])

    }

    coef <- c(tot_coef = tot_coef, dir_coef = dir_coef)

  }

  return(coef)

}
