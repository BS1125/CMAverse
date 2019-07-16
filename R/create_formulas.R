create_formulas <- function(outcome = NULL, treatment = NULL, mediator = NULL, covariates = c(),
                           interaction = NULL, event = NULL,
                           mreg = c("linear", "logistic"),
                           yreg = c("linear", "logistic", "loglinear", "poisson","quasipoisson",
                                    "negbin", "coxph", "aft_exp", "aft_weibull")) {


  mediator_formula_basic <- paste(mediator, treatment, sep = " ~ ")

  outcome_formula_basic  <- paste(paste(outcome, treatment, sep = " ~ "),
                                  mediator,
                                  sep = " + ")

  if (interaction) {

    outcome_formula_basic <- paste(outcome_formula_basic,
                                   paste(treatment, mediator, sep = " * "),
                                   sep = " + ")
  }

  if (length(covariates) == 0) {

    mediator_formula <- mediator_formula_basic

    outcome_formula  <- outcome_formula_basic

  } else {

    mediator_formula <- paste(mediator_formula_basic,
                              paste(covariates, collapse = " + "),
                              sep = " + ")

    outcome_formula  <- paste(outcome_formula_basic,
                              paste(covariates, collapse = " + "),
                              sep = " + ")
  }

  if (yreg %in% c("coxph","aft_exp","aft_weibull")) {

    l <- strsplit(outcome_formula, split = " ~ ")

    l[[1]][1] <- paste0("Surv(", outcome, ", ", event, ")")

    outcome_formula <- paste(l[[1]][1], l[[1]][2], sep = " ~ ")
  }

  formulas <- list(outcome_formula = outcome_formula,
                  mediator_formula = mediator_formula)

  return(formulas)
}
