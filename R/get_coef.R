get_coef <- function(regressions = NULL, mreg, yreg) {

  mediator_regression <- regressions$mediator_regression

  outcome_regression <- regressions$outcome_regression

  ## Store coefficients from regressions
  betas  <- lapply(1:length(mediator_regression), FUN = function(x) coefficients(mediator_regression[[x]]))

  thetas <- coefficients(outcome_regression)

  ## Store covariances from regressions
  vcov_betas <- lapply(1:length(mediator_regression), FUN = function(x) vcov(mediator_regression[[x]]))

  vcov_thetas <- vcov(outcome_regression)

  variance <- lapply(1:length(mediator_regression), FUN = function(x) {
    if (mreg[x] == "linear") { sigma(mediator_regression[[x]])^2
    } else if (mreg[x] == "logistic") NULL })

  if (yreg == "aft_weibull")
    vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

  ## Build block diagonal matrix
  vcov_block <- do.call(Matrix::bdiag, append(list(vcov_thetas), vcov_betas))

  coef <- list(betas = betas, thetas = thetas,
               vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
               variance = variance)

  return(coef)
}
