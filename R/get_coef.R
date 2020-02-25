get_coef <- function(formulas = NULL, regressions = NULL, mreg, yreg, model) {

if (model == "rb") {

  assign("mediator_formula", formulas$mediator_formula[[1]], envir = globalenv())

  mediator_regression <- regressions$mediator_regression

  outcome_regression <- regressions$outcome_regression

  ## Store coefficients from regressions
  betas  <- as.vector(t(coefficients(mediator_regression[[1]])))

  thetas <- coefficients(outcome_regression)

  ## Store covariances from regressions

  vcov_betas <- vcov(mediator_regression[[1]])

  vcov_thetas <- vcov(outcome_regression)

  if (mreg == "linear") { variance <- sigma(mediator_regression[[1]])^2
    } else variance <- NULL

  if (yreg == "aft_weibull")
    vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

  ## Build block diagonal matrix
  vcov_block <- Matrix::bdiag(vcov_thetas, vcov_betas)

  coef <- list(betas = betas, thetas = thetas,
               vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
               variance = variance)

  rm(mediator_formula, pos = ".GlobalEnv")

} else if (model == "iorw") {

  tot_outcome_regression = regressions$tot_outcome_regression
  dir_outcome_regression = regressions$dir_outcome_regression

    dir_coef <- unname(coef(dir_outcome_regression)[2])

    tot_coef <- unname(coef(tot_outcome_regression)[2])

    coef <- c(tot_coef = tot_coef, dir_coef = dir_coef)

  }

  return(coef)

}
