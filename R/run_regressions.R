run_regressions <- function(formulas = NULL, outcome, treatment, mediator, covariates,
                            interaction, event, yreg, mreg, data = NULL) {

  outcome_formula <- formulas$outcome_formula

  mediator_formula<- formulas$mediator_formula

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

  regressions <- list(mediator_regression = mediator_regression,
              outcome_regression = outcome_regression)

  return(regressions)
}
