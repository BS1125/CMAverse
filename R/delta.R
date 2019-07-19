delta <- function(coef = NULL, vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  betas <- coef$betas

  thetas <- coef$thetas

  variance <- coef$variance

  vcov_block <- coef$vcov_block

  treatment <- coef$treatment

  mediator <- coef$mediator

  covariates <- coef$covariates

  interaction <- coef$interaction

  mreg <- coef$mediator_reg

  yreg <- coef$outcome_reg

  cde <- cde_function(thetas, treatment, mediator, interaction, m, a_star, a, yreg)

  nde <- nde_function(betas, thetas, variance, vcov_block, treatment, mediator,
                      covariates, vecc, interaction, m,  a_star, a, mreg, yreg)

  nie <- nie_function(betas, thetas, treatment, mediator, covariates, vecc = vecc,
                      m, interaction, a_star, a, mreg, yreg)

  te <- te_function(betas, thetas, variance, vcov_block, treatment, mediator,
                    covariates, vecc, interaction, m,  a_star, a, mreg, yreg)

  pm <- pm_function(betas, thetas, variance, vcov_block, treatment, mediator,
                    covariates, vecc, interaction, m,  a_star, a, mreg, yreg)

  se_cde <- cde_se_delta(thetas, vcov_thetas, treatment, mediator, m, a_star, a,
                            interaction, yreg)

  se_nde <- nde_se_delta(thetas, betas, vcov_block, treatment, mediator, interaction,
                            m, vecc, a_star, a, variance, mreg, yreg)

  se_nie <- nie_se_delta(thetas, betas, vcov_block, treatment, mediator, m, vecc, interaction,
                         a_star, a, mreg, yreg)

  se_te <- te_se_delta(betas, thetas, variance, vcov_block, treatment, mediator,
                                 covariates, vecc, interaction, m,  a_star, a, mreg, yreg)

  se_pm <- pm_se_delta(betas, thetas, variance, vcov_block, treatment, mediator,
                             covariates, vecc, interaction, m,  a_star, a, mreg, yreg)

  out <- list(cde = cde, se_cde = se_cde,
           pnde = nde$pnde, se_pnde = se_nde$pnde_se_delta,
           tnde = nde$tnde, se_tnde = se_nde$tnde_se_delta,
           pnie = nie$pnie, se_pnie = se_nie$pnie_se_delta,
           tnie = nie$tnie, se_tnie = se_nie$tnie_se_delta,
           te = te, se_te = se_te,
           pm = pm, se_pm = se_pm)

  return(out)

}

#try=delta(coef = coef1, vecc=c(1,1,1),m = 1, a_star = 0, a = 1)

