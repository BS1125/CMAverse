delta <- function(coef = NULL, vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  cde_boot <- CDE_boot(coef = coef, m = m, a_star = a_star, a = a)

  nde_boot <- NDE_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  nie_boot <- NIE_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  te_boot <- total_effect_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  nde_delta <- NDE_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  nie_delta <- NIE_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  te_delta <- total_effect_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  cde_delta <- CDE_delta(coef = coef, m = m, a_star = a_star, a = a)

  pm_boot <- proportion_mediated_boot(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  pm_delta <- proportion_mediated_delta(coef = coef,vecc = vecc, m = m, a_star = a_star, a = a)

  out <- list( cded = cde_boot, se_cded = cde_delta,
               pnde = nde_boot$pnde, se_pnde = nde_delta$se_pnde_delta,
               tnde = nde_boot$tnde, se_tnde = nie_delta$se_pnie_delta,
               pnie = nie_boot$pnie, se_pnie = nie_delta$se_pnie_delta,
               tnie = nie_boot$tnie, se_tnie = nie_delta$se_tnie_delta,
               te = te_boot, se_te = te_delta,
               pm = pm_boot, se_pm = pm_delta)

  class(out) <- "delta_out"

  return(out)
}

#try=delta(coef = coef, m = 1, a_star = 2, a = 3)

print.delta_out <- function(delta_out, digits = 2, conf = 0.95) {

  printCoefmat(format_df_delta(delta_out, conf = conf, n = nrow(data)),
               digits = digits, has.Pvalue = TRUE)

}
