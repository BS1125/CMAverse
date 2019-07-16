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

print.boot_out <- function(boot_out, digits = 2, conf = 0.95) {

  printCoefmat(format_df_boot(boot_out, conf = conf), digits = digits)

}
