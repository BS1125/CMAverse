bootstrap_step <- function(formulas = NULL, outcome, treatment, mediator, covariates,
                          interaction, event, yreg, mreg, data = NULL, indices = NULL, vecc = NULL,
                           m = NULL, a_star = NULL, a = NULL) {

  regressions <- run_regressions(formulas = formulas, outcome = outcome, treatment = treatment,
                                 mediator = mediator, covariates = covariates, interaction = interaction,
                                 event = event, mreg = mreg, yreg = yreg, data = data)

  coef <- get_coef(regressions = regressions, outcome = outcome, treatment = treatment,
                   mediator = mediator, covariates = covariates, interaction = interaction,
                   event = event, mreg = mreg, yreg = yreg)

  cde_boot <- CDE_boot(coef = coef, m = m, a_star = a_star, a = a)

  nde_boot <- NDE_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  nie_boot <- NIE_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  te_boot <- total_effect_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  pm_boot <- proportion_mediated_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  return(as.numeric(c(cde = cde_boot,
                      pnde = nde_boot$pnde, tnde = nde_boot$tnde,
                      pnie = nie_boot$pnie, tnie = nie_boot$tnie,
                      te = te_boot,
                      pm = pm_boot)))
}

bootstrap <- function(formulas = NULL, outcome, treatment, mediator, covariates,
                      interaction, event = NULL, yreg, mreg, data = NULL, nboot = NULL, vecc = c(),
                      m = NULL, a_star = NULL, a = NULL) {

  estimate <- bootstrap_step(formulas = formulas, outcome = outcome, treatment = treatment,
                             mediator = mediator, covariates = covariates,
                             interaction = interaction, event = event,
                             mreg = mreg, yreg = yreg, data = data,
                             vecc=vecc, m = m, a_star = a_star, a = a)
  boot_results <- NULL

  for (i in 1:nboot) {

    indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
    data_boot <- data[indices,]

    boot_results <- rbind(boot_results,
                          bootstrap_step(formulas = formulas, outcome = outcome, treatment = treatment,
                                         mediator = mediator, covariates = covariates,
                                         interaction = interaction, event = event,
                                         mreg = mreg, yreg = yreg, data = data_boot,
                                         vecc=vecc, m = m, a_star = a_star, a = a))
  }

  out <- list(estimate = estimate, boot_results = boot_results)

  class(out) <- "boot_out"

  return(out)
}


bootstrap_decom_4way <- function(data, outcome, treatment, mediator, covariates, vecc = NULL,
                           interaction = TRUE, nboot = 100,conf=0.95,
                           mreg = c("linear", "logistic"),
                           yreg = c("linear", "logistic", "loglinear", "poisson",
                                    "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull"),
                           event = NULL, m = 0, a_star = 1, a = 0) {

  estimate <- decomp_4way_step(data = data, outcome = outcome, treatment = treatment,
                             mediator = mediator, covariates = covariates,
                             interaction = interaction, event = event,
                             mreg = mreg, yreg = yreg,
                             vecc=vecc, m = m, a_star = a_star, a = a)
  boot_results <- NULL

  for (i in 1:nboot) {

    indices <- sample(1:nrow(data), nrow(data), replace = TRUE)
    data_boot <- data[indices,]

    boot_results <- rbind(boot_results,
                          decomp_4way_step(data = data_boot, outcome = outcome, treatment = treatment,
                                           mediator = mediator, covariates = covariates,
                                           interaction = interaction, event = event,
                                           mreg = mreg, yreg = yreg,
                                           vecc=vecc, m = m, a_star = a_star, a = a))
  }

  CI_lower <- apply(boot_results, 2, function(x) quantile(x, (1-conf)/2))

  CI_upper <- apply(boot_results, 2, function(x) quantile(x, conf+(1-conf)/2))

  bias <- apply(boot_results, 2, mean) - estimate

  std.error <- apply(boot_results, 2, sd)

  out <- cbind(estimate, bias, std.error, CI_lower, CI_upper)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(out) <- c("estimate", "bias", "std.error", label_CI[1], label_CI[2])

  rownames(out) <- c(names(estimate))

  class(out) <- "boot_out"

  return(out)
}

print.boot_out <- function(boot_out, digits = 2) {

  printCoefmat(boot_out, digits = digits)

}
