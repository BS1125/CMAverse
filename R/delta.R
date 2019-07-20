delta <- function(coef = NULL, vecc = c(), m_star = NULL, a_star = NULL, a = NULL) {

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

  cde <- cde_function(thetas, treatment, mediator, interaction, m_star, a_star, a, yreg)

  nde <- nde_function(betas, thetas, variance, vcov_block, treatment, mediator,
                      covariates, vecc, interaction, a_star, a, mreg, yreg)

  nie <- nie_function(betas, thetas, treatment, mediator, covariates, vecc = vecc,
                      interaction, a_star, a, mreg, yreg)

  te <- te_function(betas, thetas, variance, vcov_block, treatment, mediator,
                    covariates, vecc, interaction, a_star, a, mreg, yreg)

  pm <- pm_function(betas, thetas, variance, vcov_block, treatment, mediator,
                    covariates, vecc, interaction, a_star, a, mreg, yreg)

  se_cde <- cde_se_delta(thetas, vcov_thetas, treatment, mediator, m_star, a_star, a,
                            interaction, yreg)

  se_nde <- nde_se_delta(thetas, betas, vcov_block, treatment, mediator, interaction,
                         vecc, a_star, a, variance, mreg, yreg)

  se_nie <- nie_se_delta(thetas, betas, vcov_block, treatment, mediator, vecc, interaction,
                         a_star, a, mreg, yreg)

  se_te <- te_se_delta(betas, thetas, variance, vcov_block, treatment, mediator,
                                 covariates, vecc, interaction, a_star, a, mreg, yreg)

  se_pm <- pm_se_delta(betas, thetas, variance, vcov_block, treatment, mediator,
                             covariates, vecc, interaction, a_star, a, mreg, yreg)

  out <- list(cde = cde, se_cde = se_cde,
           pnde = nde$pnde, se_pnde = se_nde$pnde_se_delta,
           tnde = nde$tnde, se_tnde = se_nde$tnde_se_delta,
           pnie = nie$pnie, se_pnie = se_nie$pnie_se_delta,
           tnie = nie$tnie, se_tnie = se_nie$tnie_se_delta,
           te = te, se_te = se_te,
           pm = pm, se_pm = se_pm)

  return(out)

}


delta_4way <- function(coef = NULL, vecc = c(), m_star = NULL, a_star = NULL, a = NULL) {

  betas <- coef$betas

  thetas <- coef$thetas

  variance <- coef$variance

  vcov_thetas <- coef$vcov_thetas

  vcov_betas <- coef$vcov_betas

  vcov_block <- coef$vcov_block

  outcome <- coef$outcome

  treatment <- coef$treatment

  mediator <- coef$mediator

  covariates <- coef$covariates

  interaction <- coef$interaction

  event <- coef$event

  mreg <- coef$mediator_reg

  yreg <- coef$outcome_reg

  effect_estimates <- bootstrap_4way_step(data = data, indices = c(1:nrow(data)), outcome = outcome,
                                          treatment = treatment, mediator = mediator, covariates = covariates,
                                          vecc = vecc, interaction = interaction, event = event,
                                          m_star = m_star, a_star = a_star, a = a,
                                          mreg = mreg, yreg = yreg)

  theta0 <- "x1"

  theta1 <- "x2"

  theta2 <- "x3"

  theta3 <- ifelse(interaction, paste0("x", length(thetas)), "0")

  beta0 <- paste0("x", length(thetas)+1)

  beta1 <- paste0("x", length(thetas)+2)

  XC <- ifelse(length(vecc)  > 0,
               paste0("x", k + 2 + 1:length(vecc), "*", "vecc_", 1:length(vecc), collapse = "+"),
               "0")

  if (yreg == "linear") {

    if (mreg == "linear") {

      cde_formula <- as.formula(stringr::str_replace_all(
        paste0("~(", theta1, "+", theta3, "*m_star)*(a-a_star)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      intref_formula <- as.formula(stringr::str_replace_all(
        paste0("~", theta3, "*(", beta0, "+", beta1, "*a_star", "+", XC,
               "-m_star)*(a-a_star)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      intmed_formula <- as.formula(stringr::str_replace_all(
        paste0("~", theta3, "*", beta1, "*(a-a_star)^2"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      pie_formula <- as.formula(stringr::str_replace_all(
        paste0("~(", theta2, "*", beta1, "+", theta3, "*", beta1, "*a_star)*(a-a_star)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

    } else if (mreg=="logistic") {

      cde_formula <- as.formula(stringr::str_replace_all(
        paste0("~(", theta1, "+", theta3, "*m_star)*(a-a_star)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      intref_formula <- as.formula(stringr::str_replace_all(
        paste0("~", theta3, "*(a-a_star)*(exp(", beta0, "+", beta1, "*a_star",
               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, "))-m_star)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))


      intmed_formula <- as.formula(stringr::str_replace_all(
        paste0("~", theta3, "*(a-a_star)*(exp(", beta0, "+", beta1, "*a",
               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a", "+", XC, "))-exp(",
               beta0, "+", beta1, "*a_star",
               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, ")))"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      pie_formula <- as.formula(stringr::str_replace_all(
        paste0("~(", theta2, "+", theta3, "*a_star)*(exp(", beta0, "+", beta1, "*a",
               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a", "+", XC, "))-exp(",
               beta0, "+", beta1, "*a_star",
               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, ")))"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

    }

    cde_se_delta <- msm::deltamethod(cde_formula, thetas, vcov_thetas)

    intref_se_delta <- msm::deltamethod(intref_formula, c(thetas, betas), vcov_block)

    intmed_se_delta <- msm::deltamethod(intmed_formula, c(thetas, betas), vcov_block)

    pie_se_delta <- msm::deltamethod(pie_formula, c(thetas, betas), vcov_block)

    te_se_delta <- msm::deltamethod(as.formula("~x1+x2+x3+x4"),
                                    c(effect_estimates["cde"], effect_estimates["intref"],
                                      effect_estimates["intmed"], effect_estimates["pie"]),
                                    Matrix::bdiag(cde_se_delta^2, intref_se_delta^2,
                                                  intmed_se_delta^2, pie_se_delta^2))

    cde_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                          c(effect_estimates["cde"], effect_estimates["te"]),
                                          Matrix::bdiag(cde_se_delta^2, te_se_delta^2))

    intref_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                             c(effect_estimates["intref"], effect_estimates["te"]),
                                             Matrix::bdiag(intref_se_delta^2, te_se_delta^2))

    intmed_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                             c(effect_estimates["intmed"], effect_estimates["te"]),
                                             Matrix::bdiag(intmed_se_delta^2, te_se_delta^2))

    pie_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                          c(effect_estimates["pie"], effect_estimates["te"]),
                                          Matrix::bdiag(pie_se_delta^2, te_se_delta^2))

    overall_pm_se_delta <- msm::deltamethod(as.formula("~(x1+x2)/x3"),
                                            c(effect_estimates["pie"], effect_estimates["intmed"],
                                              effect_estimates["te"]),
                                            Matrix::bdiag(pie_se_delta^2, intmed_se_delta^2,
                                                          te_se_delta^2))

    overall_int_se_delta <- msm::deltamethod(as.formula("~(x1+x2)/x3"),
                                             c(effect_estimates["intref"], effect_estimates["intmed"],
                                               effect_estimates["te"]),
                                             Matrix::bdiag(intref_se_delta^2, intmed_se_delta^2,
                                                           te_se_delta^2))

    overall_pe_se_delta <- msm::deltamethod(as.formula("~(x1+x2+x3)/x4"),
                                            c(effect_estimates["intref"], effect_estimates["intmed"],
                                              effect_estimates["pie"],effect_estimates["te"]),
                                            Matrix::bdiag(intref_se_delta^2, intmed_se_delta^2,
                                                          pie_se_delta^2, te_se_delta^2))

    return(list(cde = effect_estimates["cde"], cde_se_delta = cde_se_delta,
                intref = effect_estimates["intref"], intref_se_delta = intref_se_delta,
                intmed = effect_estimates["intmed"], intmed_se_delta = intmed_se_delta,
                pie = effect_estimates["pie"], pie_se_delta = pie_se_delta,
                te = effect_estimates["te"], te_se_delta = te_se_delta,
                cde_prop = effect_estimates["cde_prop"], cde_prop_se_delta = cde_prop_se_delta,
                intref_prop = effect_estimates["intref_prop"], intref_prop_se_delta = intref_prop_se_delta,
                intmed_prop = effect_estimates["intmed_prop"], intmed_prop_se_delta = intmed_prop_se_delta,
                pie_prop = effect_estimates["pie_prop"], pie_prop_se_delta = pie_prop_se_delta,
                overall_pm = effect_estimates["overall_pm"], overall_pm_se_delta = overall_pm_se_delta,
                overall_int = effect_estimates["overall_int"], overall_int_se_delta = overall_int_se_delta,
                overall_pe = effect_estimates["overall_pe"], overall_pe_se_delta = overall_pe_se_delta))

  } else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                          "negbin", "coxph", "aft_exp", "aft_weibull")) {

    if (mreg=="linear") {
      cde_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp(", theta1, "*(a-a_star)+", theta2, "*m_star+", theta3,
               "*a*m_star-(", theta2, "+", theta3, "*a_star)(", beta0, "+",
               beta1,"*a_star+", XC, ")-0.5*(", theta2, "+", theta3,
               "*a_star)^2*variance)-exp(", theta2, "*m_star+", theta3,
               "*a_star*m_star-(", theta2, "+", theta3, "*a_star)(", beta0,
               "+", beta1, "*a_star", "+", XC, ")-0.5*(",theta2, "+", theta3,
               "*a_star)^2*variance)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star),
                    "\\bvariance\\b" = as.character(variance))))

      intref_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+",
               beta1,"*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
               theta3, "^2*variance*(a^2-a_star^2))-1-exp(", theta1, "*(a-a_star)+",
               theta2, "*m_star+", theta3, "*a*m_star-(", theta2, "+", theta3,
               "*a_star)(", beta0, "+", beta1, "*a_star+", XC, ")-0.5*(", theta2,
               "+", theta3, "*a_star)^2*variance)+exp(", theta2, "*m_star+",
               theta3, "*a_star*m_star-(", theta2, "+", theta3, "*a_star)(", beta0,
               "+", beta1,"*a_star+", XC, ")-0.5*(", theta2, "+", theta3,
               "*a_star)^2*variance)"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star),
                    "\\bvariance\\b" = as.character(variance))))

      intmed_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp((", theta1, "+", theta2, "*", beta1, "+", theta3,
               "*(", beta0, "+", beta1,"*a_star+", "+", beta1,"*a+", XC,
               "+", theta2, "*variance))*(a-a_star)+0.5*", theta3,
               "^2*variance*(a^2-a_star^2))-exp((", theta2, "*", beta1,
               "+", theta3, "*", beta1, "*a_star)*(a-a_star))-exp((",
               theta1, "+", theta3, "*(", beta0, "+", beta1, "*a_star+", XC,
               "+", theta2, "*variance))*(a-a_star)+0.5*", theta3,
               "^2*variance*(a^2-a_star^2))+1"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star),
                    "\\bvariance\\b" = as.character(variance))))

      pie_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp((", theta2, "*", beta1, "+", theta3, "*", beta1,
               "*a_star)*(a-a_star))-1"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      total_rr_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
               "*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
               theta3, "^2*variance*(a^2-a_star^2))*exp((", theta2, "*",
               beta1, "+", theta3, "*", beta1, "*a)*(a-a_star))"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star),
                    "\\bvariance\\b" = as.character(variance))))

    } else if (mreg=="logistic") {

      cde_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp(", theta1, "*(a-a_star)+", theta2, "*m_star+", theta3,
               "*a*m_star)*(1+exp(", beta0, "+", beta1, "*a_star+", XC,
               "))/(1+exp(",  beta0, "+", beta1,  "*a_star+", XC, "+",
               theta2, "+", theta3, "*a_star))-exp(", theta2, "*m_star+",
               theta3, "*a_star*m_star)*(1+exp(", beta0, "+", beta1,
               "*a_star+", XC, "))/(1+exp(", beta0, "+", beta1, "*a_star+",
               XC, "+", theta2, "+", theta3, "*a_star))"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      intref_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp(", theta1, "*(a-a_star))*(1+exp(", beta0, "+",
               beta1, "*a_star+", XC, "+", theta2, "+", theta3,
               "*a))/(1+exp(", beta0, "+", beta1, "*a_star+", XC, "+",
               theta2, "+", theta3, "*a_star))-1-exp(", theta1,
               "*(a-a_star)+", theta2, "*m_star+", theta3, "*a*m_star)*(1+exp(",
               beta0, "+", beta1, "*a_star+", XC, "))*exp((", theta1, "+",
               theta3, "*m_star)*(a-a_star))/(1+exp(", beta0, "+", beta1,
               "*a_star+", XC, "+", theta2, "+", theta3, "*a_star))+exp(",
               theta2, "*m_star+", theta3, "*a_star*m_star)*(1+exp(",
               beta0, "+", beta1, "*a_star+", XC, "))/(1+exp(", beta0, "+",
               beta1, "*a_star+", XC, "+", theta2, "+", theta3, "*a_star))"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      intmed_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~(exp(", theta1, "*(a-a_star))*(1+exp(", beta0, "+",
               beta1, "*a+", XC, "+", theta2, "+", theta3, "*a))*(1+exp(",
               beta0, "+", beta1, "*a_star+", XC, ")))/((1+exp(",
               beta0, "+", beta1, "*a_star+", XC, "+", theta2, "+", theta3,
               "*a_star))*(1+exp(", beta0, "+", beta1, "*a+", XC,
               ")))-((1+exp(", beta0, "+", beta1, "*a+", XC, "+", theta2,
               "+", theta3, "*a_star))*(1+exp(", beta0, "+", beta1,
               "*a_star+", XC, ")))/((1+exp(", beta0, "+", beta1, "*a_star+",
               XC, "+", theta2, "+", theta3, "*a_star))*(1+exp(", beta0, "+",
               beta1, "*a+", XC, ")))-(exp(", theta1, "*(a-a_star))*(1+exp(",
               beta0, "+", beta1, "*a_star+", XC, "+", theta2, "+", theta3,
               "*a)))/(1+exp(", beta0, "+", beta1, "*a_star+", XC, "+",
               theta2, "+", theta3, "*a_star))+1"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      pie_comp_formula <- as.formula(stringr::str_replace_all(
        paste0("~((1+exp(", beta0, "+", beta1, "*a_star+", XC, "))*(1+exp(",
               beta0, "+", beta1, "*a+", XC, "+", theta2, "+", theta3,
               "*a_star)))/((1+exp(", beta0, "+", beta1, "*a+", XC,
               "))*(1+exp(", beta0, "+", beta1, "*a_star+", XC, "+", theta2,
               "+", theta3, "*a_star)))-1"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))

      total_rr_formula <- as.formula(stringr::str_replace_all(
        paste0("~exp(", theta1, "*a)*(1+exp(", beta0, "+", beta1, "*a_star",
               "+", XC, "))*(1+exp(", beta0, "+", beta1, "*a", "+", XC, "+",
               theta2, "+", theta3, "*a))/(exp(", theta1, "*a_star)*(1+exp(",
               beta0, "+", beta1, "*a+", XC, "))*(1+exp(", beta0, "+",
               beta1, "*a_star+", XC, "+", theta2, "+", theta3, "*a_star)))"),
        pattern = c("\\ba_star\\b" = as.character(a_star),
                    "\\ba\\b" = as.character(a),
                    "\\bm_star\\b" = as.character(m_star))))
    }

    cde_comp_se_delta <- msm::deltamethod(cde_comp_formula, c(thetas, betas), vcov_block)

    intref_comp_se_delta <- msm::deltamethod(intref_comp_formula, c(thetas, betas), vcov_block)

    intmed_comp_se_delta <- msm::deltamethod(intmed_comp_formula, c(thetas, betas), vcov_block)

    pie_comp_se_delta <- msm::deltamethod(pie_comp_formula, c(thetas, betas), vcov_block)

    total_rr_se_delta <- msm::deltamethod(total_rr_formula, c(thetas, betas), vcov_block)

    total_err_se_delta <- total_rr_se_delta

    tcomp_se_delta <- msm::deltamethod(as.formula("~x1+x2+x3+x4"),
                                       c(effect_estimates["cde_comp"], effect_estimates["intref_comp"],
                                         effect_estimates["intmed_comp"], effect_estimates["pie_comp"]),
                                       Matrix::bdiag(cde_comp_se_delta^2, intref_comp_se_delta^2,
                                                     intmed_comp_se_delta^2, pie_comp_se_delta^2))

    cde_err_se_delta <- msm::deltamethod(as.formula("~x1*(x2-1)/x3"),
                                         c(effect_estimates["cde_comp"], effect_estimates["total_rr"],
                                           effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                             effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                         Matrix::bdiag(cde_comp_se_delta^2, total_rr_se_delta^2,
                                                       tcomp_se_delta^2))

    intmed_err_se_delta <- msm::deltamethod(as.formula("~x1*(x2-1)/x3"),
                                            c(effect_estimates["intmed_comp"], effect_estimates["total_rr"],
                                              effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                            Matrix::bdiag(intmed_comp_se_delta^2, total_rr_se_delta^2,
                                                          tcomp_se_delta^2))

    intref_err_se_delta <- msm::deltamethod(as.formula("~x1*(x2-1)/x3"),
                                            c(effect_estimates["intref_comp"], effect_estimates["total_rr"],
                                              effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                            Matrix::bdiag(intref_comp_se_delta^2, total_rr_se_delta^2,
                                                          tcomp_se_delta^2))

    pie_err_se_delta <- msm::deltamethod(as.formula("~x1*(x2-1)/x3"),
                                         c(effect_estimates["pie_comp"], effect_estimates["total_rr"],
                                           effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                             effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                         Matrix::bdiag(pie_comp_se_delta^2, total_rr_se_delta^2,
                                                       tcomp_se_delta^2))

    cde_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                          c(effect_estimates["cde_comp"],
                                            effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                              effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                          Matrix::bdiag(cde_comp_se_delta^2, tcomp_se_delta^2))

    intmed_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                             c(effect_estimates["intmed_comp"],
                                               effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                 effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                             Matrix::bdiag(intmed_comp_se_delta^2, tcomp_se_delta^2))

    intref_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                             c(effect_estimates["intref_comp"],
                                               effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                 effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                             Matrix::bdiag(intref_comp_se_delta^2, tcomp_se_delta^2))

    pie_prop_se_delta <- msm::deltamethod(as.formula("~x1/x2"),
                                          c(effect_estimates["pie_comp"],
                                            effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                              effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                          Matrix::bdiag(pie_comp_se_delta^2, tcomp_se_delta^2))

    overall_pm_se_delta <- msm::deltamethod(as.formula("~(x1+x2)/x3"),
                                            c(effect_estimates["pie_comp"], effect_estimates["intmed_comp"],
                                              effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                            Matrix::bdiag(pie_comp_se_delta^2, intmed_comp_se_delta^2,
                                                          tcomp_se_delta^2))

    overall_int_se_delta <- msm::deltamethod(as.formula("~(x1+x2)/x3"),
                                             c(effect_estimates["intref_comp"], effect_estimates["intmed_comp"],
                                               effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                 effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                             Matrix::bdiag(intref_comp_se_delta^2, intmed_comp_se_delta^2,
                                                           tcomp_se_delta^2))

    overall_pe_se_delta <- msm::deltamethod(as.formula("~(x1+x2+x3)/x4"),
                                            c(effect_estimates["intref_comp"], effect_estimates["intmed_comp"],
                                              effect_estimates["pie_comp"],
                                              effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                            Matrix::bdiag(intref_comp_se_delta^2, intmed_comp_se_delta^2,
                                                          pie_comp_se_delta^2, tcomp_se_delta^2))

    return(list(cde_comp = unname(effect_estimates["cde_comp"]), cde_comp_se_delta = cde_comp_se_delta,
                intref_comp = unname(effect_estimates["intref_comp"]), intref_comp_se_delta = intref_comp_se_delta,
                intmed_comp = unname(effect_estimates["intmed_comp"]), intmed_comp_se_delta = intmed_comp_se_delta,
                pie_comp = unname(effect_estimates["pie_comp"]), pie_comp_se_delta = pie_comp_se_delta,
                total_rr = unname(effect_estimates["total_rr"]), total_rr_se_delta = total_rr_se_delta,
                cde_err = unname(effect_estimates["cde_err"]), cde_err_se_delta = cde_err_se_delta,
                intref_err = unname(effect_estimates["intref_err"]), intref_err_se_delta = intref_err_se_delta,
                intmed_err = unname(effect_estimates["intmed_err"]), intmed_err_se_delta = intmed_err_se_delta,
                pie_err = unname(effect_estimates["pie_err"]), pie_err_se_delta = pie_err_se_delta,
                cde_prop = unname(effect_estimates["cde_prop"]), cde_prop_se_delta = cde_prop_se_delta,
                intref_prop = unname(effect_estimates["intref_prop"]), intref_prop_se_delta = intref_prop_se_delta,
                intmed_prop = unname(effect_estimates["intmed_prop"]), intmed_prop_se_delta = intmed_prop_se_delta,
                pie_prop = unname(effect_estimates["pie_prop"]), pie_prop_se_delta = pie_prop_se_delta,
                overall_pm = unname(effect_estimates["overall_pm"]), overall_pm_se_delta = overall_pm_se_delta,
                overall_int = unname(effect_estimates["overall_int"]), overall_int_se_delta = overall_int_se_delta,
                overall_pe = unname(effect_estimates["overall_pe"]), overall_pe_se_delta = overall_pe_se_delta))
  }
}


