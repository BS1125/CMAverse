te_function <- function(pnde, tnie, yreg) {

  if (yreg == "linear") {

    te <- tnie + pnde

    return(c(te = te))

  } else {

    ORte <- tnie * pnde

    return(c(ORte = ORte))
  }
}


te_se_delta <- function(betas, thetas, variance, vcov_block, treatment, mediator,
                        covariates, vecc, interaction, m,  a_star, a, mreg, yreg) {

  if (yreg == "linear") {

    te_formula <- as.formula("~x1 + x2")

    pnde <- nde_function(betas, thetas, variance, vcov_block, treatment, mediator,
                         covariates, vecc, interaction, m,  a_star, a, mreg, yreg)["pnde"]

    tnie <- nie_function(betas, thetas, treatment, mediator, covariates, vecc,
                         m, interaction, a_star, a, mreg, yreg)["tnie"]

  } else {

    te_formula <- as.formula("~x1 * x2")

    pnde <- nde_function(betas, thetas, variance, vcov_block, treatment, mediator,
                         covariates, vecc, interaction, m,  a_star, a, mreg, yreg)["ORpnde"]

    tnie <- nie_function(betas, thetas, treatment, mediator, covariates, vecc = vecc,
                         m, interaction, a_star, a, mreg, yreg)["ORtnie"]

  }

  se_pnde_delta <- nde_se_delta(thetas, betas, vcov_block, treatment, mediator, interaction,
                                m, vecc, a_star, a, variance, mreg, yreg)["pnde_se_delta"]

  se_tnie_delta <- nie_se_delta(thetas, betas, vcov_block, treatment, mediator, m, vecc, interaction,
                                a_star, a, mreg, yreg)["tnie_se_delta"]

  te_se_delta <- msm::deltamethod(te_formula, c(pnde, tnie),
                                  Matrix::bdiag(se_pnde_delta^2, se_tnie_delta^2))

  return(c(te_se_delta = te_se_delta))
}


te_se_delta(thetas=coef1$thetas, betas=coef1$betas,variance=coef1$variance,vcov_block=coef1$vcov_block,
            treatment, mediator,
            covariates, vecc, interaction, m,  a_star, a, mreg, yreg)
