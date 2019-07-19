causal_mediation <- function(data, outcome, treatment, mediator, covariates, vecc = NULL,
                             interaction = TRUE, type = c("delta", "bootstrap"), nboot = 100,
                             mreg = c("linear", "logistic"),
                             yreg = c("linear", "logistic", "loglinear", "poisson",
                                      "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull"),
                             event = NULL, m = 0, a_star = 1, a = 0) {

  if (is.null(covariates) & !is.null(vecc)) {
    warning("Incompatible arguments")
  } else if (!is.null(covariates) & is.null(vecc)) {
    vecc <- colMeans(as.data.frame(data[, covariates]))
  }

  formulas <- create_formulas(outcome = outcome, treatment = treatment, mediator = mediator,
                              covariates = covariates, interaction = interaction, event = event,
                              mreg = mreg, yreg = yreg)

  regressions <- run_regressions(formulas = formulas, outcome = outcome, treatment = treatment,
                                 mediator = mediator, covariates = covariates, interaction = interaction,
                                 event = event, mreg = mreg, yreg = yreg, data = data)

  coef <- get_coef(regressions = regressions, outcome = outcome, treatment = treatment,
                   mediator = mediator, covariates = covariates, interaction = interaction,
                   event = event, mreg = mreg, yreg = yreg)

  n=nrow(data)

  if (type == "delta") {

    effect_estimate <- format_df_delta(delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a),
                                       conf = conf, n = n, yreg  = yreg)

    class(effect_estimate)=c("delta_out","data.frame")

  }else if (type == "bootstrap") {

    effect_estimate <- format_df_boot(bootstrap(formulas = formulas, outcome = outcome, treatment = treatment,
                                 mediator = mediator, covariates = covariates, interaction = interaction,
                                 event = event, mreg = mreg, yreg = yreg, data = data,
                                 nboot = nboot, vecc = vecc,
                                 m = m, a_star = a_star, a = a),
                                 conf = conf, n = n, yreg  = yreg)

    class(effect_estimate)=c("boot_out","data.frame")

    }



  return(effect_estimate)
}
