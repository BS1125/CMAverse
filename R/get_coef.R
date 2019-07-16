get_coef <- function(regressions = NULL, outcome, treatment, mediator, covariates,
                     interaction, event, mreg, yreg) {

  mediator_regression <- regressions$mediator_regression

  outcome_regression <- regressions$outcome_regression

  ## Store coefficients from regressions
  betas  <- coefficients(mediator_regression)

  thetas <- coefficients(outcome_regression)

  ## Store covariances from regressions
  vcov_betas <- vcov(mediator_regression)

  vcov_thetas <- vcov(outcome_regression)

  if (mreg == "linear")
    variance = summary(mediator_regression)$sigma^2
  else variance = NULL

  if (yreg == "aft_weibull")
    vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

  ## Build block diagonal matrix
  vcov_block <- Matrix::bdiag(vcov_thetas, vcov_betas)

  coef <- list(outcome = outcome, treatment = treatment, mediator = mediator,
               covariates = covariates, interaction = interaction,
               event = event, mediator_reg = mreg, outcome_reg = yreg,betas = betas, thetas=thetas,
               vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
               variance = variance)

  return(coef)
}

