run_regressions <- function(formulas = NULL, data = NULL, vecc = NULL) {

  outcome <- formulas$outcome

  treatment <- formulas$treatment

  mediator <- formulas$mediator

  covariates <- formulas$covariates

  interaction <- formulas$interaction

  event <- formulas$event

  yreg <- formulas$outcome_reg

  outcome_formula <- formulas$outcome_formula

  mreg <- formulas$mediator_reg

  mediator_formula<- formulas$mediator_formula

  if (is.null(covariates) & !is.null(vecc)) {
    warning("Incompatible arguments")
  } else if (!is.null(covariates) & is.null(vecc)) {
    vecc <- colMeans(as.data.frame(data[, covariates]))
  }

  if (yreg == "linear") {
    outcome_regression  <- lm(outcome_formula, data = data)
  }

  if (yreg == "logistic") {
    outcome_regression  <- glm(outcome_formula, family = binomial(), data = data)
  }

  if (yreg == "loglinear") {
    outcome_regression  <- glm(outcome_formula, family = binomial("log"), data = data)
  }

  if (yreg == "poisson") {
    outcome_regression  <- glm(outcome_formula, family = poisson(), data = data)
  }

  if (yreg == "quasipoisson") {
    outcome_regression  <- glm(outcome_formula, family = quasipoisson(), data = data)
  }

  if (yreg == "negbin") {
    outcome_regression  <- glm.nb(outcome_formula, data = data)
  }

  if (yreg == "coxph") {
    outcome_regression <- coxph(as.formula(outcome_formula), data = data)
  }

  if (yreg == "aft_exp") {
    outcome_regression <- survreg(as.formula(outcome_formula), dist = "exponential", data = data)
  }

  if (yreg == "aft_weibull") {
    outcome_regression <- survreg(as.formula(outcome_formula), dist = "weibull", data = data)
  }

  if (mreg == "linear") {
    mediator_regression <- lm(mediator_formula, data = data)
  }

  if (mreg == "logistic") {
    mediator_regression <- glm(mediator_formula, family = binomial(), data = data)
  }

  regressions <- list(outcome = outcome, treatment = treatment, mediator = mediator,
              covariates = covariates, interaction = interaction, event = event,
              outcome_reg = yreg,
              mediator_reg = mreg,
              outcome_formula = outcome_formula,
              mediator_formula = mediator_formula,
              mediator_regression = mediator_regression,
              outcome_regression = outcome_regression)

  return(regressions)
}

#regressions = run_regressions(formulas, data)
