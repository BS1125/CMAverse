mediation <- function(regressions = NULL, robustSE = TRUE, sims = 100, ...) {

  treatment <- regressions$treatment
  mediator <- regressions$mediator
  mediator_regression <- regressions$mediator_regression
  outcome_regression <- regressions$outcome_regression


  result <- mediation::mediate(model.m = mediator_regression,
                               model.y = outcome_regression,
                               treat = treatment,
                               mediator = mediator,
                               robustSE = robustSE,
                               sims = sims, ...)
  return(summary(result))
}

#mediation(regressions)
