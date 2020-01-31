delta_method <- function(mediator, thetas, betas, variance, vcov_block,
                         EMint, mreg, yreg, vecc, a, a_star, m_star) {

  theta0 <- "x1"

  theta1 <- "x2"

  theta2 <- "x3"

  theta3 <- ifelse(EMint, paste0("x", length(thetas)), "0")

  beta0 <- paste0("x", length(thetas)+1)

  beta1 <- paste0("x", length(thetas)+2)

  XC <- ifelse(length(vecc)  > 0,
               paste0("x", length(thetas) + 2 + 1:length(vecc), "*", vecc,
                      collapse = "+"),
               "0")

  if (yreg == "linear") {

    cde_formula <- tnde_formula <- pnde_formula <- paste0(theta1, "*(a-a_star)")

    tnie_formula <- pnie_formula <- ifelse(mreg == "linear",
                                           paste(beta1, theta2, "(a-a_star)", sep = "*"), paste0("exp(", beta0,
                                                                                                 "+", beta1, "*a+", XC, ")*", theta2, "/(1+exp(", beta0, "+",
                                                                                                 beta1, "*a+", XC, "))-exp(", beta0, "+", beta1,
                                                                                                 "*a_star+", XC, ")*", theta2, "/(1+exp(", beta0, "+", beta1,
                                                                                                 "*a_star+", XC, "))"))

    if (EMint == TRUE) {

      cde_formula <- paste(cde_formula, paste(theta3, as.character(m_star), "(a-a_star)",
                                              sep = "*",collapse = "+"), sep = "+")


      pnde_formula <- paste(pnde_formula, ifelse(mreg == "linear",
                                                 paste0(theta3, "*(a-a_star)*(", beta0, "+", beta1,
                                                        "*a_star+", XC, ")", collapse = "+"), paste0(theta3, "*(a-a_star)*", "exp(", beta0,
                                                                                                     "+", beta1, "*a_star+", XC, ")", "/(1+exp(",
                                                                                                     beta0, "+", beta1, "*a_star+", XC, "))")), sep = "+")

      tnde_formula <- paste(tnde_formula, ifelse(mreg == "linear",
                                                 paste0(theta3, "*(a-a_star)*(", beta0, "+", beta1,
                                                        "*a+", XC, ")", collapse = "+"), paste0(theta3, "*(a-a_star)*", "exp(", beta0,
                                                                                                "+", beta1, "*a+", XC, ")", "/(1+exp(",
                                                                                                beta0, "+", beta1, "*a+", XC, "))")), sep = "+")

      pnie_formula <- paste(pnie_formula, ifelse(mreg == "linear",
                                                 paste0(theta3, "*", beta1, "*a_star*(a-a_star)", collapse = "+"),
                                                 paste0("exp(", beta0, "+", beta1, "*a+", XC, ")*", theta3,
                                                        "*a_star/(1+exp(", beta0, "+",
                                                        beta1, "*a+", XC, "))-exp(", beta0, "+",
                                                        beta1, "*a_star+", XC, ")*", theta3,
                                                        "*a_star/(1+exp(", beta0, "+", beta1,
                                                        "*a_star+", XC, "))")), sep = "+")

      tnie_formula <- paste(tnie_formula, ifelse(mreg == "linear",
                                                 paste0(theta3, "*", beta1, "*a*(a-a_star)", collapse = "+"),
                                                 paste0("exp(", beta0, "+", beta1, "*a+", XC, ")*", theta3,
                                                        "*a/(1+exp(", beta0, "+",
                                                        beta1, "*a+", XC, "))-exp(", beta0, "+",
                                                        beta1, "*a_star+", XC, ")*", theta3,
                                                        "*a/(1+exp(", beta0, "+", beta1,
                                                        "*a_star+", XC, "))")), sep = "+")
    }

    te_formula <- paste0("(", tnie_formula, ")+(", pnde_formula, ")")

    pm_formula <- paste0("(", tnie_formula, ")/((", pnde_formula, ")+(", te_formula, "))")

    if (mreg == "linear") {

      intref_formula <- paste0(theta3, "*(", beta0, "+", beta1, "*a_star", "+", XC,
                               "-m_star)*(a-a_star)")

      intmed_formula <- paste0(theta3, "*", beta1, "*(a-a_star)^2")

    } else if (mreg=="logistic") {

      intref_formula <- paste0(theta3, "*(a-a_star)*(exp(", beta0, "+", beta1, "*a_star",
                               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, "))-m_star)")


      intmed_formula <- paste0(theta3, "*(a-a_star)*(exp(", beta0, "+", beta1, "*a",
                               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a", "+", XC, "))-exp(",
                               beta0, "+", beta1, "*a_star",
                               "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, ")))")

    }

    pie_formula <- pnie_formula

    cde_prop_formula <- paste0("(", cde_formula, ")/(", te_formula, ")")

    intref_prop_formula <- paste0("(", intref_formula, ")/(", te_formula, ")")

    intmed_prop_formula <- paste0("(", intmed_formula, ")/(", te_formula, ")")

    pie_prop_formula <- paste0("(", pie_formula, ")/(", te_formula, ")")

    overall_pm_formula <- paste0("((", pnie_formula, ")+(", intmed_formula, "))/(", te_formula, ")")

    overall_int_formula <- paste0("((", intref_formula, ")+(", intmed_formula, "))/(", te_formula, ")")

    overall_pe_formula <- paste0("((", intref_formula, ")+(", intmed_formula, ")+(", pnie_formula,
                                 "))/(", te_formula, ")")

    delta_formula <- list(cde_formula = cde_formula, pnde_formula = pnde_formula,
                          tnde_formula = tnde_formula, pnie_formula = pnie_formula,
                          tnie_formula = tnie_formula, te_formula = te_formula, pm_formula = pm_formula,
                          intref_formula = intref_formula, intmed_formula = intmed_formula, pie_formula = pie_formula,
                          cde_prop_formula = cde_prop_formula, intref_prop_formula = intref_prop_formula,
                          intmed_prop_formula = intmed_prop_formula, pie_prop_formula = pie_prop_formula,
                          overall_pm_formula = overall_pm_formula, overall_int_formula = overall_int_formula,
                          overall_pe_formula = overall_pe_formula)

  } else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                          "negbin", "coxph", "aft_exp", "aft_weibull")) {

    if (mreg == "linear") {

      cde_rr_formula <- paste0("exp((", theta1, "+", theta3, "*m_star)*(a-a_star))")

      pnde_rr_formula <- paste0("exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                                "*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                                theta3, "^2*variance*(a^2-a_star^2))")

      tnde_rr_formula <- paste0("exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                                "*a+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                                theta3, "^2*variance*(a^2-a_star^2))")

      pnie_rr_formula <- paste0("exp((", theta2, "*", beta1, "+", theta2, "*",
                                beta1, "*a_star)*(a-a_star))")

      tnie_rr_formula <- paste0("exp((", theta2, "*", beta1, "+", theta3, "*",
                                beta1, "*a)*(a-a_star))")

      cde_err_formula <- paste0("exp(", theta1, "*(a-a_star)+", theta2, "*m_star+", theta3,
                                "*a*m_star-(", theta2, "+", theta3, "*a_star)*(", beta0, "+",
                                beta1,"*a_star+", XC, ")-0.5*(", theta2, "+", theta3,
                                "*a_star)^2*variance)-exp(", theta2, "*m_star+", theta3,
                                "*a_star*m_star-(", theta2, "+", theta3, "*a_star)*(", beta0,
                                "+", beta1, "*a_star", "+", XC, ")-0.5*(",theta2, "+", theta3,
                                "*a_star)^2*variance)")

      intref_err_formula <- paste0("exp((", theta1, "+", theta3, "*(", beta0, "+",
                                   beta1,"*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                                   theta3, "^2*variance*(a^2-a_star^2))-1-exp(", theta1, "*(a-a_star)+",
                                   theta2, "*m_star+", theta3, "*a*m_star-(", theta2, "+", theta3,
                                   "*a_star)*(", beta0, "+", beta1, "*a_star+", XC, ")-0.5*(", theta2,
                                   "+", theta3, "*a_star)^2*variance)+exp(", theta2, "*m_star+",
                                   theta3, "*a_star*m_star-(", theta2, "+", theta3, "*a_star)*(", beta0,
                                   "+", beta1,"*a_star+", XC, ")-0.5*(", theta2, "+", theta3,
                                   "*a_star)^2*variance)")

      intmed_err_formula <- paste0("exp((", theta1, "+", theta2, "*", beta1, "+", theta3,
                                   "*(", beta0, "+", beta1,"*a_star+", "+", beta1,"*a+", XC,
                                   "+", theta2, "*variance))*(a-a_star)+0.5*", theta3,
                                   "^2*variance*(a^2-a_star^2))-exp((", theta2, "*", beta1,
                                   "+", theta3, "*", beta1, "*a_star)*(a-a_star))-exp((",
                                   theta1, "+", theta3, "*(", beta0, "+", beta1, "*a_star+", XC,
                                   "+", theta2, "*variance))*(a-a_star)+0.5*", theta3,
                                   "^2*variance*(a^2-a_star^2))+1")

    } else if (mreg=="logistic") {

      cde_rr_formula <- paste0("exp((", theta1, "+", theta3, "*m_star)*(a-a_star))")

      pnde_rr_formula <- paste0("exp(", theta1, "*(a-a_star))*(1+exp(", theta2, "+",
                                theta3, "*a+", beta0, "+", beta1, "*a_star+", XC,
                                "))/(1+exp(", theta2, "+", theta3, "*a_star+", beta0,
                                "+", beta1, "*a_star+", XC, "))")

      tnde_rr_formula <- paste0("exp(", theta1, "*(a-a_star))*(1+exp(", theta2, "+",
                                theta3, "*a+", beta0, "+", beta1, "*a+", XC,
                                "))/(1+exp(", theta2, "+", theta3, "*a_star+", beta0,
                                "+", beta1, "*a+", XC, "))")

      pnie_rr_formula <- paste0("(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                                "))*(1+exp(", theta2, "+", theta3, "*a_star+",
                                beta0, "+", beta1, "*a+", XC, "))/((1+exp(",
                                beta0, "+", beta1, "*a+", XC, "))*(1+exp(",
                                theta2, "+", theta3, "*a_star+", beta0, "+", beta1,
                                "*a_star+", XC, ")))")

      tnie_rr_formula <- paste0("(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                                "))*(1+exp(", theta2, "+", theta3, "*a+",
                                beta0, "+", beta1, "*a+", XC, "))/((1+exp(",
                                beta0, "+", beta1, "*a+", XC, "))*(1+exp(",
                                theta2, "+", theta3, "*a+", beta0, "+", beta1,
                                "*a_star+", XC, ")))")

      cde_err_formula <- paste0("exp(", theta1, "*(a-a_star)+", theta2, "*m_star+", theta3,
                                "*a*m_star)*(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                                "))/(1+exp(",  beta0, "+", beta1,  "*a_star+", XC, "+",
                                theta2, "+", theta3, "*a_star))-exp(", theta2, "*m_star+",
                                theta3, "*a_star*m_star)*(1+exp(", beta0, "+", beta1,
                                "*a_star+", XC, "))/(1+exp(", beta0, "+", beta1, "*a_star+",
                                XC, "+", theta2, "+", theta3, "*a_star))")


      intref_err_formula <- paste0("exp(", theta1, "*(a-a_star))*(1+exp(", beta0, "+",
                                   beta1, "*a_star+", XC, "+", theta2, "+", theta3,
                                   "*a))/(1+exp(", beta0, "+", beta1, "*a_star+", XC, "+",
                                   theta2, "+", theta3, "*a_star))-1-exp(", theta1,
                                   "*(a-a_star)+", theta2, "*m_star+", theta3, "*a*m_star)*(1+exp(",
                                   beta0, "+", beta1, "*a_star+", XC, "))*exp((", theta1, "+",
                                   theta3, "*m_star)*(a-a_star))/(1+exp(", beta0, "+", beta1,
                                   "*a_star+", XC, "+", theta2, "+", theta3, "*a_star))+exp(",
                                   theta2, "*m_star+", theta3, "*a_star*m_star)*(1+exp(",
                                   beta0, "+", beta1, "*a_star+", XC, "))/(1+exp(", beta0, "+",
                                   beta1, "*a_star+", XC, "+", theta2, "+", theta3, "*a_star))")

      intmed_err_formula <- paste0("(exp(", theta1, "*(a-a_star))*(1+exp(", beta0, "+",
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
                                   theta2, "+", theta3, "*a_star))+1")

    }

    te_rr_formula <- paste0("(", tnie_rr_formula, ")*(", pnde_rr_formula, ")")

    pm_formula <- paste0("((", pnde_rr_formula, ")*(", tnie_rr_formula, "-1))/(", te_rr_formula, "-1)")

    pie_err_formula <- paste0(pnie_rr_formula, "-1")

    cde_err_prop_formula <- paste0("(", cde_err_formula, ")/(", te_rr_formula, "-1)")

    intmed_err_prop_formula <- paste0("(", intmed_err_formula, ")/(", te_rr_formula, "-1)")

    intref_err_prop_formula <-paste0("(", intref_err_formula, ")/(", te_rr_formula, "-1)")

    pie_err_prop_formula <- paste0("(", pie_err_formula, ")/(", te_rr_formula, "-1)")

    overall_pm_formula <- paste0("((", pnie_rr_formula, "-1)+(", intmed_err_formula, "))/(",
                                 te_rr_formula, "-1)")

    overall_int_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, "))/(",
                                  te_rr_formula, "-1)")

    overall_pe_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, ")+(",
                                 ORpnie_formula,
                                 "-1))/(", te_rr_formula, "-1)")

    delta_formula <- list(cde_rr_formula = cde_rr_formula, pnde_rr_formula = pnde_rr_formula,
                          tnde_rr_formula = tnde_rr_formula, pnie_rr_formula = pnie_rr_formula,
                          tnie_rr_formula = tnie_rr_formula, te_rr_formula = te_rr_formula, pm_formula = pm_formula,
                          intref_err_formula = intref_err_formula, intmed_err_formula = intmed_err_formula,
                          pie_err_formula = pie_err_formula, cde_err_prop_formula = cde_err_prop_formula,
                          intref_err_prop_formula = intref_err_prop_formula,
                          intmed_err_prop_formula = intmed_err_prop_formula, pie_err_prop_formula = pie_err_prop_formula,
                          overall_pm_formula = overall_pm_formula, overall_int_formula = overall_int_formula,
                          overall_pe_formula = overall_pe_formula)

  }

  effect_se <- c()

  for (formula in names(delta_formula)) {

    delta_formula[[formula]] <- as.formula(paste0("~", stringr::str_replace_all(
      delta_formula[[formula]],
      pattern = c("\\ba_star\\b" = as.character(a_star),
                  "\\ba\\b" = as.character(a),
                  "\\bm_star\\b" = as.character(m_star),
                  "\\bvariance\\b" = as.character(variance)))))

    effect_se <- c(effect_se, msm::deltamethod(delta_formula[[formula]], c(thetas, betas), vcov_block))

  }

  names(effect_se) = stringr::str_replace(names(delta_formula), "_formula", "_se")

  return(effect_se)

}
