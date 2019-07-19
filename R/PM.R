pm_function <- function(betas, thetas, variance, vcov_block, treatment, mediator,
                        covariates, vecc, interaction, a_star, a, mreg, yreg) {


  pnde <- nde_function(betas, thetas, variance, vcov_block, treatment, mediator,
                         covariates, vecc, interaction, a_star, a, mreg, yreg)$pnde

  tnie <- nie_function(betas, thetas, treatment, mediator, covariates, vecc,
                         interaction, a_star, a, mreg, yreg)$tnie

  te <- te_function(betas, thetas, variance, vcov_block, treatment, mediator,
                      covariates, vecc, interaction, a_star, a, mreg, yreg)

  if (yreg == "linear") { pm <- tnie / (pnde + te)
  } else { pm <-  (pnde * (tnie - 1)) / (pnde * tnie - 1) }

  return(pm)

}

pm_se_delta  <- function(betas, thetas, variance, vcov_block, treatment, mediator,
                         covariates, vecc, interaction, a_star, a, mreg, yreg) {


    pnde <- nde_function(betas, thetas, variance, vcov_block, treatment, mediator,
                         covariates, vecc, interaction, a_star, a, mreg, yreg)$pnde

    tnie <- nie_function(betas, thetas, treatment, mediator, covariates, vecc,
                         interaction, a_star, a, mreg, yreg)$tnie

    te <- te_function(betas, thetas, variance, vcov_block, treatment, mediator,
                      covariates, vecc, interaction, a_star, a, mreg, yreg)

    pnde_se_delta <- nde_se_delta(thetas, betas, vcov_block, treatment, mediator, interaction,
                                  vecc, a_star, a, variance, mreg, yreg)$pnde_se_delta

    tnie_se_delta <- nie_se_delta(thetas, betas, vcov_block, treatment, mediator, vecc, interaction,
                                  a_star, a, mreg, yreg)$tnie_se_delta

    te_se_delta <- te_se_delta(betas, thetas, variance, vcov_block, treatment, mediator,
                               covariates, vecc, interaction, a_star, a, mreg, yreg)

    if (yreg == "linear") { pm_formula <- as.formula("~x2 / (x1 + x3)")
    } else { pm_formula <-  as.formula("~(x1 * (x2 - 1)) / (x1 * x2 - 1)") }

    pm_se_delta <- msm::deltamethod(pm_formula, c(pnde, tnie, te),
                                  Matrix::bdiag(pnde_se_delta^2, tnie_se_delta^2, te_se_delta^2))

    return(pm_se_delta)

}

pm_se_delta(thetas=coef1$thetas, betas=coef1$betas,variance=coef1$variance,vcov_block=coef1$vcov_block,
            treatment, mediator, covariates, vecc, interaction, a_star, a, mreg, yreg)
