boostrap_step <- function(formulas = NULL, data = NULL, indices = NULL, vecc = NULL,
                         m = NULL, a_star = NULL, a = NULL) {
  data_boot <- data[indices, ]

  coef <- get_coef(formulas = formulas, data = data_boot)

  outcome <- coef$outcome
  treatment <- coef$treatment
  mediator <- coef$mediator
  covariates <- coef$covariates
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

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

bootstrap <- function(formulas = NULL, data = NULL, nboot = NULL, vecc = c(),
                     m = NULL, a_star = NULL, a = NULL) {

  out <- boot::boot(data = data, statistic = boostrap_step, R = nboot,
                    formulas = formulas, vecc = vecc, m = m, a_star = a_star, a = a)

  class(out) <- "boot_out"

  return(out)
}

print.boot_out <- function(boot_out, digits = 2, conf = 0.95) {

  printCoefmat(format_df_boot(boot_out, conf = conf, n = nrow(data)), digits = digits)

}
