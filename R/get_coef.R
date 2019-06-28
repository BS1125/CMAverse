get_coef <- function(regressions = NULL) {

  outcome <- regressions$outcome

  treatment <- regressions$treatment

  mediator <- regressions$mediator

  covariates <- regressions$covariates

  interaction <- regressions$interaction

  event <- regressions$event

  yreg <- regressions$outcome_reg

  mreg <- regressions$mediator_reg

  mediator_regression <- regressions$mediator_regression

  outcome_regression <- regressions$outcome_regression

  ## Store coefficients from regressions
  betas  <- coefficients(mediator_regression)

  thetas <- coefficients(outcome_regression)

  ## Store covariances from regressions
  vcov_betas <- vcov(mediator_regression)

  vcov_thetas <- vcov(outcome_regression)

  if (yreg == "aft_weibull")
    vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

  ## Build block diagonal matrix
  vcov_block <- Matrix::bdiag(vcov_thetas, vcov_betas)

  variance <- (summary(mediator_regression)$sigma)^2

  coef <- list(outcome = outcome, treatment = treatment, mediator = mediator,
               covariates = covariates, interaction = interaction, event = event,
               outcome_reg = yreg, mediator_reg = mreg,
               betas = betas, thetas=thetas,
               vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
               variance = variance)

  return(coef)
}

#coef = get_coef(regressions)
