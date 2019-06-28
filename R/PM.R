proportion_mediated_boot_function <- function(pnde, tnie, te, ycont = FALSE) {
  res <- 0
  if (ycont) {
    res <- tnie / (pnde + te)
  } else {
    res <-  (pnde * (tnie - 1)) / (pnde * tnie - 1)
  }
  res
}

proportion_mediated_delta_function  <- function(ycont = FALSE) {
  res <- 0
  if (ycont) {
    res <- "~x2 / (x1 + x3)"
  } else {
    res <-  "~(x1 * (x2 - 1)) / (x1 * x2 - 1)"
  }
  as.formula(res)
}


proportion_mediated_boot = function(coef = list(), vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment = coef$treatment
  mediator = coef$mediator
  covariates = coef$covariates
  interaction = coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  nde_boot <- NDE_boot(coef, vecc, m, a_star, a)

  nie_boot <- NIE_boot(coef, vecc, m, a_star, a)

  te_boot <- total_effect_boot(coef, vecc, m, a_star, a)

  pm_boot <- proportion_mediated_boot_function(nde_boot$pnde, nie_boot$tnie, te_boot,
                                               ycont = (yreg == "linear"))

  pm_boot
}

#proportion_mediated_boot(coef, m = 1, a_star = 2, a = 3)

proportion_mediated_delta = function(coef, vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment = coef$treatment
  mediator = coef$mediator
  covariates = coef$covariates
  interaction = coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  nde_boot <- NDE_boot(coef, vecc, m, a_star, a)

  nie_boot <- NIE_boot(coef, vecc, m, a_star, a)

  te_boot <- total_effect_boot(coef, vecc, m, a_star, a)

  se_pnde_delta <- NDE_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)$se_pnde_delta

  se_tnie_delta <- NIE_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)$se_tnie_delta

  se_te_delta <- total_effect_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)

  pm_delta <- proportion_mediated_delta_function(ycont = (yreg == "linear"))

  se_pm_delta <- msm::deltamethod(pm_delta, c(nde_boot$pnde, nie_boot$tnie, te_boot),
                                  Matrix::bdiag(se_pnde_delta^2, se_tnie_delta^2, se_te_delta^2))

  se_pm_delta
}

#proportion_mediated_delta(coef = coef, m = 1, a_star = 2, a = 3)
