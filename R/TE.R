total_effect_boot_function <- function(pnde, tnie, ycont = FALSE) {
  res <- 0
  if (ycont) {
    res <- tnie + pnde
  } else {
    res <- tnie * pnde
  }
  res
}


total_effect_delta_function <- function(ycont = FALSE) {
  res <- 0
  if (ycont) {
    res <- "~x1 + x2"
  } else {
    res <- "~x1 * x2"
  }
  as.formula(res)
}


total_effect_boot = function(coef = list(), vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment
  mediator <- coef$mediator
  covariates <- coef$covariates
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg


  nde_boot <- NDE_boot(coef, vecc, m, a_star, a)

  nie_boot <- NIE_boot(coef, vecc, m, a_star, a)

  te_boot <- total_effect_boot_function(nde_boot$pnde, nie_boot$tnie, ycont = (yreg == "linear"))

  return(te_boot)
}

#total_effect_boot(coef, m = 1, a_star = 2, a = 3)

total_effect_delta = function(coef = list(), vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment
  mediator <- coef$mediator
  covariates <- coef$covariates
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  nde_boot <- NDE_boot(coef, vecc, m, a_star, a)

  nie_boot <- NIE_boot(coef, vecc, m, a_star, a)

  se_pnde_delta <- NDE_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)$se_pnde_delta

  se_tnie_delta <- NIE_delta(coef = coef, vecc = vecc, m = m, a_star = a_star, a = a)$se_tnie_delta

  te_delta <- total_effect_delta_function(ycont = (yreg == "linear"))

  se_te_delta <- msm::deltamethod(te_delta, c(nde_boot$pnde, nie_boot$tnie),
                                  Matrix::bdiag(se_pnde_delta^2, se_tnie_delta^2))

  return(se_te_delta)
}

#total_effect_delta(coef, m = 1, a_star = 2, a = 3)
