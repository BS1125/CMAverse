causal_mediation <- function(data, method = c("delta", "bootstrap", "simulation"),
                             nboot = 100, conf = 0.95, nsims = 1000,
                             outcome, treatment, mediator, covariates,
                             vecc = NULL, interaction = TRUE, event = NULL,
                             m_star = 0, a_star = 1, a = 0,
                             mreg = c("linear", "logistic"),
                             yreg = c("linear", "logistic", "loglinear", "poisson",
                                      "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) {

  if (is.null(covariates) & !is.null(vecc)) {
    warning("Incompatible arguments")
  } else if (!is.null(covariates) & is.null(vecc)) {
    vecc <- colMeans(as.data.frame(data[, covariates]))
  }

  formulas <- create_formulas(outcome = outcome, treatment = treatment, mediator = mediator,
                              covariates = covariates, interaction = interaction, event = event,
                              mreg = mreg, yreg = yreg)

  regressions <- run_regressions(formulas = formulas, outcome = outcome, treatment = treatment,
                                 mediator = mediator, covariates = covariates, interaction = interaction,
                                 event = event, mreg = mreg, yreg = yreg, data = data)

  coef <- get_coef(regressions = regressions, outcome = outcome, treatment = treatment,
                   mediator = mediator, covariates = covariates, interaction = interaction,
                   event = event, mreg = mreg, yreg = yreg)

  n <- nrow(data)

  betas <- coef$betas

  thetas <- coef$thetas

  variance <- coef$variance

  vcov_thetas <- coef$vcov_thetas

  vcov_betas <- coef$vcov_betas

  vcov_block <- coef$vcov_block

  effect_estimates <- bootstrap_step(data = data, indices = c(1:nrow(data)), outcome = outcome,
                                          treatment = treatment, mediator = mediator, covariates = covariates,
                                          vecc = vecc, interaction = interaction, event = event,
                                          m_star = m_star, a_star = a_star, a = a,
                                          mreg = mreg, yreg = yreg)

###############################Delta method: 3 way and 4 way decomposition################################

  if (method == "delta") {

    theta0 <- "x1"

    theta1 <- "x2"

    theta2 <- "x3"

    theta3 <- ifelse(interaction, paste0("x", length(thetas)), "0")

    beta0 <- paste0("x", length(thetas)+1)

    beta1 <- paste0("x", length(thetas)+2)

    XC <- ifelse(length(vecc)  > 0,
                 paste0("x", length(thetas) + 2 + 1:length(vecc), "*", "vecc_", 1:length(vecc),
                        collapse = "+"),
                 "0")

    if (yreg == "linear") {

      if (mreg == "linear") {

        cde_formula <- stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*m_star)*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnde_formula <- stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a_star+", XC, "))*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnde_formula <- stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a+", XC, "))*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnie_formula <- stringr::str_replace_all(
          paste0("~(", theta2, "*", beta1, "+", theta3, "*", beta1,
                 "*a_star)*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnie_formula <- stringr::str_replace_all(
          paste0("~(", theta2, "*", beta1, "+", theta3, "*", beta1,
                 "*a)*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

      } else if (mreg=="logistic") {

        cde_formula <- stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*m_star)*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnde_formula <- stringr::str_replace_all(
          paste0("~", theta1, "*(a-a_star)+(", theta3, "*(a-a_star))*exp(",
                 beta0, "+", beta1, "*a_star+", XC, ")/(1+exp(",
                 beta0, "+", beta1, "*a_star+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnde_formula <- stringr::str_replace_all(
          paste0("~", theta1, "*(a-a_star)+(", theta3, "*(a-a_star))*exp(",
                 beta0, "+", beta1, "*a+", XC, ")/(1+exp(",
                 beta0, "+", beta1, "*a+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnie_formula <- stringr::str_replace_all(
          paste0("~(", theta2, "+", theta3, "*a_star)*(exp(", beta0,
                 "+", beta1, "*a+", XC, ")/(1+exp(", beta0, "+",
                 beta1, "*a+", XC, "))-exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")/(1+exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnie_formula <- stringr::str_replace_all(
          paste0("~(", theta2, "+", theta3, "*a)*(exp(", beta0,
                 "+", beta1, "*a+", XC, ")/(1+exp(", beta0, "+",
                 beta1, "*a+", XC, "))-exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")/(1+exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

      }

      j <- length(vecc)

      if (j > 0) {

        for (formula in c("cde_formula", "pnde_formula", "tnde_formula",
                          "pnie_formula", "tnie_formula")) {

          mid <- get(formula)

          for (i in 1:j) {

            mid <-   stringr::str_replace_all(mid,paste("vecc", i, sep = "_"),
                                              as.character(vecc[i]))

          }

          assign(formula, as.formula(mid))

        }
      }

      cde_se <- msm::deltamethod(cde_formula, thetas, vcov_thetas)

      pnde_se <- msm::deltamethod(pnde_formula, c(thetas, betas), vcov_block)

      tnde_se <- msm::deltamethod(tnde_formula, c(thetas, betas), vcov_block)

      pnie_se <- msm::deltamethod(pnie_formula, c(thetas, betas), vcov_block)

      tnie_se <- msm::deltamethod(tnie_formula, c(thetas, betas), vcov_block)

      te_formula <- as.formula("~x1+x2")

      te_se <- msm::deltamethod(te_formula, c(effect_estimates["pnde"], effect_estimates["tnie"]),
                                Matrix::bdiag(pnde_se^2, tnie_se^2))

      pm_formula <- as.formula("~x2/(x1+x3)")

      pm_se <- msm::deltamethod(pm_formula, c(effect_estimates["pnde"], effect_estimates["tnie"],
                                              effect_estimates["te"]),
                                Matrix::bdiag(pnde_se^2, tnie_se^2, te_se^2))

      intref_formula <- intmed_formula <- as.formula("~x1-x2")

      intref_se <- msm::deltamethod(intref_formula, c(effect_estimates["pnde"], effect_estimates["cde"]),
                                Matrix::bdiag(pnde_se^2, cde_se^2))

      intmed_se <- msm::deltamethod(intmed_formula, c(effect_estimates["tnie"], effect_estimates["pnie"]),
                                    Matrix::bdiag(tnie_se^2, pnie_se^2))

      cde_prop_formula <- intref_prop_formula <- intmed_prop_formula <- pie_prop_formula <- as.formula("~x1/x2")

      cde_prop_se <- msm::deltamethod(cde_prop_formula, c(effect_estimates["cde"], effect_estimates["te"]),
                                    Matrix::bdiag(cde_se^2, te_se^2))

      intref_prop_se <- msm::deltamethod(intref_prop_formula, c(effect_estimates["intref"], effect_estimates["te"]),
                                      Matrix::bdiag(intref_se^2, te_se^2))

      intmed_prop_se <- msm::deltamethod(intmed_prop_formula, c(effect_estimates["intmed"], effect_estimates["te"]),
                                         Matrix::bdiag(intmed_se^2, te_se^2))

      pie_prop_se <- msm::deltamethod(pie_prop_formula, c(effect_estimates["pie"], effect_estimates["te"]),
                                         Matrix::bdiag(pie_se^2, te_se^2))

      overall_pm_formula <- overall_int_formula <- as.formula("~(x1+x2)/x3")

      overall_pm_se <- msm::deltamethod(overall_pm_formula, c(effect_estimates["pie"], effect_estimates["intmed"],
                                          effect_estimates["te"]),
                                      Matrix::bdiag(pie_se^2, intmed_se^2, te_se^2))

      overall_int_se <- msm::deltamethod(overall_int_formula, c(effect_estimates["intref"], effect_estimates["intmed"],
                                                              effect_estimates["te"]),
                                        Matrix::bdiag(intref_se^2, intmed_se^2, te_se^2))

      overall_pe_formula <- as.formula("~(x1+x2+x3)/x4")

      overall_pe_se <- msm::deltamethod(overall_pe_formula, c(effect_estimates["intref"], effect_estimates["intmed"],
                                                              effect_estimates["pie"], effect_estimates["te"]),
                                         Matrix::bdiag(intref_se^2, intmed_se^2, pie_se^2, te_se^2))

      # 3 way decomposition effects and se's
      decomp3way <-  c(cde = unname(effect_estimates["cde"]), cde_se = cde_se,
                       pnde = unname(effect_estimates["pnde"]), pnde_se = pnde_se,
                       tnde = unname(effect_estimates["tnde"]), tnde_se = tnde_se,
                       pnie = unname(effect_estimates["pnie"]), pnie_se = pnie_se,
                       tnie = unname(effect_estimates["tnie"]), tnie_se = tnie_se,
                       te = unname(effect_estimates["te"]), te_se = te_se,
                       pm = unname(effect_estimates["pm"]), pm_se = pm_se)

      # 4 way decomposition effects and se's
      decomp4way <-  c(cde = unname(effect_estimates["cde"]), cde_se = cde_se,
                       intref = unname(effect_estimates["intref"]), intref_se = intref_se,
                       intmed = unname(effect_estimates["intmed"]), intmed_se = intmed_se,
                       pie = unname(effect_estimates["pie"]), pie_se = pie_se,
                       te = unname(effect_estimates["te"]), te_se = te_se,
                       cde_prop = unname(effect_estimates["cde_prop"]), cde_prop_se = cde_prop_se,
                       intref_prop = unname(effect_estimates["intref_prop"]), intref_prop_se = intref_prop_se,
                       intmed_prop = unname(effect_estimates["intmed_prop"]), intmed_prop_se = intmed_prop_se,
                       pie_prop = unname(effect_estimates["pie_prop"]), pie_prop_se = pie_prop_se,
                       overall_pm = unname(effect_estimates["overall_pm"]), overall_pm_se = overall_pm_se,
                       overall_int = unname(effect_estimates["overall_int"]), overall_int_se = overall_int_se,
                       overall_pe = unname(effect_estimates["overall_pe"]), overall_pe_se = overall_pe_se)

    } else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                            "negbin", "coxph", "aft_exp", "aft_weibull")) {

      if (mreg=="linear") {

        cde_formula <- stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*m_star)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnde_formula <- stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                 theta3, "^2*variance*(a^2-a_star^2))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star),
                      "\\bvariance\\b" = as.character(variance)))

        tnde_formula <- stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                 theta3, "^2*variance*(a^2-a_star^2))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star),
                      "\\bvariance\\b" = as.character(variance)))

        pnie_formula <- stringr::str_replace_all(
          paste0("~exp((", theta2, "*", beta1, "+", theta2, "*",
                 beta1, "*a_star)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnie_formula <- stringr::str_replace_all(
          paste0("~exp((", theta2, "*", beta1, "+", theta2, "*",
                 beta1, "*a)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        cde_err_formula <- stringr::str_replace_all(
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
                      "\\bvariance\\b" = as.character(variance)))

      } else if (mreg=="logistic") {

        cde_formula <- stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*m_star)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnde_formula <- stringr::str_replace_all(
          paste0("~exp(", theta1, "*(a-a_star))*(1+exp(", theta2, "+",
                 theta3, "*a+", beta0, "+", beta1, "*a_star+", XC,
                 "))/(1+exp(", theta2, "+", theta3, "*a_star+", beta0,
                 "+", beta1, "*a_star+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnde_formula <- stringr::str_replace_all(
          paste0("~exp(", theta1, "*(a-a_star))*(1+exp(", theta2, "+",
                 theta3, "*a+", beta0, "+", beta1, "*a+", XC,
                 "))/(1+exp(", theta2, "+", theta3, "*a_star+", beta0,
                 "+", beta1, "*a+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        pnie_formula <- stringr::str_replace_all(
          paste0("~(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                 "))*(1+exp(", theta2, "+", theta3, "*a_star+",
                 beta0, "+", beta1, "*a+", XC, "))/((1+exp(",
                 beta0, "+", beta1, "*a+", XC, "))*(1+exp(",
                 theta2, "+", theta3, "*a_star+", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        tnie_formula <- stringr::str_replace_all(
          paste0("~(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                 "))*(1+exp(", theta2, "+", theta3, "*a+",
                 beta0, "+", beta1, "*a+", XC, "))/((1+exp(",
                 beta0, "+", beta1, "*a+", XC, "))*(1+exp(",
                 theta2, "+", theta3, "*a+", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

        cde_err_formula <- stringr::str_replace_all(
          paste0("~exp(", theta1, "*(a-a_star)+", theta2, "*m_star+", theta3,
                 "*a*m_star)*(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                 "))/(1+exp(",  beta0, "+", beta1,  "*a_star+", XC, "+",
                 theta2, "+", theta3, "*a_star))-exp(", theta2, "*m_star+",
                 theta3, "*a_star*m_star)*(1+exp(", beta0, "+", beta1,
                 "*a_star+", XC, "))/(1+exp(", beta0, "+", beta1, "*a_star+",
                 XC, "+", theta2, "+", theta3, "*a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star)))

      }

      j <- length(vecc)

      if (j > 0) {

        for (formula in c("cde_formula", "pnde_formula", "tnde_formula",
                          "pnie_formula", "tnie_formula", "cde_err_formula")) {

          mid <- get(formula)

          for (i in 1:j) {

            mid <-   stringr::str_replace_all(mid,paste("vecc", i, sep = "_"),
                                              as.character(vecc[i]))

          }

          assign(formula, as.formula(mid))

        }
      }

      cde_rr_se <- msm::deltamethod(cde_formula, thetas, vcov_thetas)

      pnde_rr_se <- msm::deltamethod(pnde_formula, c(thetas, betas), vcov_block)

      tnde_rr_se <- msm::deltamethod(tnde_formula, c(thetas, betas), vcov_block)

      pnie_rr_se <- pie_err_se <- msm::deltamethod(pnie_formula, c(thetas, betas), vcov_block)

      tnie_rr_se <- msm::deltamethod(tnie_formula, c(thetas, betas), vcov_block)

      te_formula <- as.formula("~x1*x2")

      te_rr_se <- total_err_se <- msm::deltamethod(te_formula, c(effect_estimates["pnde_rr"],
                                                              effect_estimates["tnie_rr"]),
                                Matrix::bdiag(pnde_rr_se^2, tnie_rr_se^2))

      pm_formula <- as.formula("~x1*(x2-1)/(x1*x2-1)")

      pm_se <- msm::deltamethod(pm_formula, c(effect_estimates["pnde_rr"], effect_estimates["tnie_rr"],
                                              effect_estimates["te_rr"]),
                                Matrix::bdiag(pnde_rr_se^2, tnie_rr_se^2, te_rr_se^2))

      cde_err_se <-  msm::deltamethod(cde_err_formula, c(thetas, betas), vcov_block)

      intref_err_formula <- as.formula("~x1-x2-1")

      intref_err_se <- msm::deltamethod(intref_err_formula, c(effect_estimates["pnde_rr"],
                                                              effect_estimates["cde_rr"]),
                       Matrix::bdiag(pnde_rr_se^2, cde_rr_se^2))

      intmed_err_formula <- as.formula("~x1*x2-x2-x3+1")

      intmed_err_se <- msm::deltamethod(intmed_err_formula, c(effect_estimates["tnie_rr"],
                                                              effect_estimates["pnde_rr"],
                                                              effect_estimates["pnie_rr"]),
                                    Matrix::bdiag(tnie_rr_se^2, pnde_rr_se^2, pnie_rr_se^2))

      cde_err_prop_formula <- intmed_err_prop_formula <- intref_err_prop_formula <-
        pie_err_prop_formula <- as.formula("~x1/x2")

      cde_err_prop_se <- msm::deltamethod(cde_err_prop_formula, c(effect_estimates["cde_err"],
                                                          effect_estimates["total_err"]),
                                      Matrix::bdiag(cde_err_se^2, total_err_se^2))

      intmed_err_prop_se <- msm::deltamethod(intmed_err_prop_formula, c(effect_estimates["intmed_err"],
                                                          effect_estimates["total_err"]),
                                      Matrix::bdiag(intmed_err_se^2, total_err_se^2))

      intref_err_prop_se <- msm::deltamethod(intref_err_prop_formula, c(effect_estimates["intref_err"],
                                                          effect_estimates["total_err"]),
                                      Matrix::bdiag(intref_err_se^2, total_err_se^2))

      pie_err_prop_se <- msm::deltamethod(pie_err_prop_formula, c(effect_estimates["pie_err"],
                                                          effect_estimates["total_err"]),
                                      Matrix::bdiag(pie_err_se^2, total_err_se^2))

      overall_pm_formula <- overall_int_formula <- as.formula("~(x1+x2)/x3")

      overall_pm_se <- msm::deltamethod(overall_pm_formula, c(effect_estimates["pie_err"],
                                                              effect_estimates["intmed_err"],
                                                              effect_estimates["total_err"]),
                                      Matrix::bdiag(pie_err_se^2, intmed_err_se^2, total_err_se^2))

      overall_int_se <- msm::deltamethod(overall_int_formula, c(effect_estimates["intref_err"],
                                                                effect_estimates["intmed_err"],
                                                                effect_estimates["total_err"]),
                                        Matrix::bdiag(intref_err_se^2, intmed_err_se^2, total_err_se^2))

      overall_pe_formula <- as.formula("~(x1+x2+x3)/x4")

      overall_pe_se <- msm::deltamethod(overall_pe_formula, c(effect_estimates["intref_err"],
                                                              effect_estimates["intmed_err"],
                                                              effect_estimates["pie_err"],
                                                              effect_estimates["total_err"]),
                                         Matrix::bdiag(intref_err_se^2, intmed_err_se^2,
                                                       pie_err_se^2, total_err_se^2))

      # 3 way decomposition effects and se's
      decomp3way <-  c(cde_rr = unname(effect_estimates["cde_rr"]), cde_rr_se = cde_rr_se,
                       pnde_rr = unname(effect_estimates["pnde_rr"]), pnde_rr_se = pnde_rr_se,
                       tnde_rr = unname(effect_estimates["tnde_rr"]), tnde_rr_se = tnde_rr_se,
                       pnie_rr = unname(effect_estimates["pnie_rr"]), pnie_rr_se = pnie_rr_se,
                       tnie_rr = unname(effect_estimates["tnie_rr"]), tnie_rr_se = tnie_rr_se,
                       te_rr = unname(effect_estimates["te_rr"]), te_rr_se = te_rr_se,
                       pm = unname(effect_estimates["pm"]), pm_se = pm_se)

      # 4 way decomposition effects and se's
      decomp4way <-  c(cde_err = unname(effect_estimates["cde_err"]), cde_err_se = cde_err_se,
                       intref_err = unname(effect_estimates["intref_err"]), intref_err_se = intref_err_se,
                       intmed_err = unname(effect_estimates["intmed_err"]), intmed_err_se = intmed_err_se,
                       pie_err = unname(effect_estimates["pie_err"]), pie_err_se = pie_err_se,
                       total_err = unname(effect_estimates["total_err"]), total_err_se = total_err_se,
                       cde_err_prop = unname(effect_estimates["cde_err_prop"]), cde_err_prop_se = cde_err_prop_se,
                       intref_err_prop = unname(effect_estimates["intref_err_prop"]),
                       intref_err_prop_se = intref_err_prop_se,
                       intmed_err_prop = unname(effect_estimates["intmed_err_prop"]),
                       intmed_err_prop_se = intmed_err_prop_se,
                       pie_err_prop = unname(effect_estimates["pie_err_prop"]), pie_err_prop_se = pie_err_prop_se,
                       overall_pm = unname(effect_estimates["overall_pm"]), overall_pm_se = overall_pm_se,
                       overall_int = unname(effect_estimates["overall_int"]), overall_int_se = overall_int_se,
                       overall_pe = unname(effect_estimates["overall_pe"]), overall_pe_se = overall_pe_se)

    }

    out <- list(decomp3way = decomp3way, decomp4way = decomp4way)

  }


#####################################Bootstrapping: 3 way and 4 way decomposition####################################

  if (method == "bootstrap") {

     out <- boot::boot(data = data, statistic = bootstrap_step, R = nboot,
                        outcome = outcome, treatment = treatment, mediator = mediator,
                        covariates = covariates, vecc = vecc, interaction = interaction,
                        event = event, mreg = mreg, yreg = yreg,
                        m_star = m_star, a_star = a_star, a = a)

     se <- apply(out$t, 2, sd)

     out <- effect_estimates

     for (i in 1:length(effect_estimates))
       out[paste0(names(effect_estimates)[i], "_se")] <- se[i]

     if (yreg == "linear") {

       # 3 way decomposition effects and se's
       decomp3way <-  c(cde = unname(out["cde"]), cde_se =  unname(out["cde_se"]),
                      pnde = unname(out["pnde"]), pnde_se =  unname(out["pnde_se"]),
                      tnde = unname(out["tnde"]), tnde_se =  unname(out["tnde_se"]),
                      pnie = unname(out["pnie"]), pnie_se =  unname(out["pnie_se"]),
                      tnie = unname(out["tnie"]), tnie_se =  unname(out["tnie_se"]),
                      te = unname(out["te"]), te_se =  unname(out["te_se"]),
                      pm = unname(out["pm"]), pm_se =  unname(out["pm_se"]))

       # 4 way decomposition effects and se's
       decomp4way <-  c(cde = unname(out["cde"]), cde_se = unname(out["cde_se"]),
                      intref = unname(out["intref"]), intref_se = unname(out["intref_se"]),
                      intmed = unname(out["intmed"]), intmed_se = unname(out["intmed_se"]),
                      pie = unname(out["pie"]), pie_se = unname(out["pie_se"]),
                      te = unname(out["te"]), te_se = unname(out["te_se"]),
                      cde_prop = unname(out["cde_prop"]), cde_prop_se = unname(out["cde_prop_se"]),
                      intref_prop = unname(out["intref_prop"]), intref_prop_se = unname(out["intref_prop_se"]),
                      intmed_prop = unname(out["intmed_prop"]), intmed_prop_se = unname(out["intmed_prop_se"]),
                      pie_prop = unname(out["pie_prop"]), pie_prop_se = unname(out["pie_prop_se"]),
                      overall_pm = unname(out["overall_pm"]), overall_pm_se = unname(out["overall_pm_se"]),
                      overall_int = unname(out["overall_int"]), overall_int_se = unname(out["overall_int_se"]),
                      overall_pe = unname(out["overall_pe"]), overall_pe_se = unname(out["overall_pe_se"]))

     } else if (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                     "negbin", "coxph", "aft_exp", "aft_weibull")) {

       # 3 way decomposition effects and se's
       decomp3way <-  c(cde_rr = unname(out["cde_rr"]), cde_rr_se = unname(out["cde_rr_se"]),
                        pnde_rr = unname(out["pnde_rr"]), pnde_rr_se = unname(out["pnde_rr_se"]),
                        tnde_rr = unname(out["tnde_rr"]), tnde_rr_se = unname(out["tnde_rr_se"]),
                        pnie_rr = unname(out["pnie_rr"]), pnie_rr_se = unname(out["pnie_rr_se"]),
                        tnie_rr = unname(out["tnie_rr"]), tnie_rr_se = unname(out["tnie_rr_se"]),
                        te_rr = unname(out["te_rr"]), te_rr_se = unname(out["te_rr_se"]),
                        pm = unname(out["pm"]), pm_se = unname(out["pm_se"]))

       # 4 way decomposition effects and se's
       decomp4way <-  c(cde_err = unname(out["cde_err"]), cde_err_se = unname(out["cde_err_se"]),
                        intref_err = unname(out["intref_err"]), intref_err_se = unname(out["intref_err_se"]),
                        intmed_err = unname(out["intmed_err"]), intmed_err_se = unname(out["intmed_err_se"]),
                        pie_err = unname(out["pie_err"]), pie_err_se = unname(out["pie_err_se"]),
                        total_err = unname(out["total_err"]), total_err_se = unname(out["total_err_se"]),
                        cde_err_prop = unname(out["cde_err_prop"]), cde_err_prop_se = unname(out["cde_err_prop_se"]),
                        intref_err_prop = unname(out["intref_err_prop"]),
                        intref_err_prop_se = unname(out["intref_err_prop_se"]),
                        intmed_err_prop = unname(out["intmed_err_prop"]),
                        intmed_err_prop_se = unname(out["intmed_err_prop_se"]),
                        pie_err_prop = unname(out["pie_err_prop"]), pie_err_prop_se = unname(out["pie_err_prop_se"]),
                        overall_pm = unname(out["overall_pm"]), overall_pm_se = unname(out["overall_pm_se"]),
                        overall_int = unname(out["overall_int"]), overall_int_se = unname(out["overall_int_se"]),
                        overall_pe = unname(out["overall_pe"]), overall_pe_se = unname(out["overall_pe_se"]))

     }

     out <- list(decomp3way = decomp3way, decomp4way = decomp4way)

  }


###############################Simulation-based approach: 3 way and 4 way decomposition##################################

  if (method == "simulation") {

    EY0m_sim <- EY1m_sim <- EY00_sim <- EY01_sim <- EY10_sim <- EY11_sim <-c()

    for (j in 1:nsims) {

      data_sim <- data[sample(1:nrow(data), replace = TRUE),]

      regressions_sim <- run_regressions(formulas = formulas, outcome = outcome, treatment = treatment,
                                         mediator = mediator, covariates = covariates, interaction = interaction,
                                         event = event, mreg = mreg, yreg = yreg, data = data_sim)

      # counterfactual design matrices for the mediator model
      mdesign_a <- cbind(data.frame(treatment = c(rep(a,n))),data[,covariates])

      mdesign_a_star <- cbind(data.frame(treatment = c(rep(a_star,n))),data[,covariates])

      colnames(mdesign_a) <- colnames(mdesign_a_star) <- c(treatment, covariates)

      # simulate counterfactual mediators
      m_a <- predict(regressions_sim$mediator_regression, newdata = mdesign_a, type = "response")

      m_a_star <- predict(regressions_sim$mediator_regression, newdata = mdesign_a_star, type = "response")

      # counterfactual design matrices for the outcome model
      ydesign0m <- cbind(data.frame(treatment = c(rep(a_star,n)), mediator = rep(m, n)), data[,covariates])

      ydesign1m <- cbind(data.frame(treatment = c(rep(a,n)), mediator = rep(m, n)), data[,covariates])

      ydesign00 <- cbind(data.frame(treatment = c(rep(a_star,n)), mediator = m_a_star), data[,covariates])

      ydesign01 <- cbind(data.frame(treatment = c(rep(a_star,n)), mediator = m_a), data[,covariates])

      ydesign10 <- cbind(data.frame(treatment = c(rep(a,n)), mediator = m_a_star), data[,covariates])

      ydesign11 <- cbind(data.frame(treatment = c(rep(a,n)), mediator = m_a), data[,covariates])

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(treatment, mediator, covariates)

      # simulate counterfactual mediators
      y0m <- predict(regressions_sim$outcome_regression, newdata =  ydesign0m, type = "response")

      y1m <- predict(regressions_sim$outcome_regression, newdata =  ydesign1m, type = "response")

      y00 <- predict(regressions_sim$outcome_regression, newdata =  ydesign00, type = "response")

      y01 <- predict(regressions_sim$outcome_regression, newdata =  ydesign01, type = "response")

      y10 <- predict(regressions_sim$outcome_regression, newdata =  ydesign10, type = "response")

      y11 <- predict(regressions_sim$outcome_regression, newdata =  ydesign11, type = "response")

      # calculate causal effect components
      EY0m_sim <- c(EY0m_sim, sum(y0m)/n) #E(Ya0m*)

      EY1m_sim <- c(EY1m_sim, sum(y1m)/n) #E(Ya1m*)

      EY00_sim <- c(EY00_sim, sum(y00)/(n*nsims)) #E(Ya0Ma0)

      EY10_sim <- c(EY10_sim, sum(y10)/(n*nsims)) #E(Ya1Ma0)

      EY01_sim <- c(EY01_sim, sum(y01)/(n*nsims)) #E(Ya0Ma1)

      EY11_sim <- c(EY11_sim, sum(y11)/(n*nsims)) #E(Ya1Ma1)

      }

    # linear outcome
    if (yreg == "linear") {

         cde_sim <- EY1m_sim - EY0m_sim # cde in the 3 way or 4 way decomposition

         pnde_sim <- EY10_sim - EY00_sim # pnde in the 3 way decomposition

         tnie_sim <- EY11_sim - EY10_sim # tnie in the 3 way decomposition

         tnde_sim <- EY11_sim - EY01_sim # tnde in the 3 way decomposition

         pnie_sim <- EY01_sim - EY00_sim # pnie in the 3 way decomposition

         te_sim <- tnie_sim + pnde_sim # te in the 3 way or 4 way decomposition

         pm_sim <- tnie_sim / (pnde_sim + te_sim) # pm in the 3 way decomposition

         intref_sim <- pnde_sim - cde_sim # intref in the 4 way decomposition

         intmed_sim <- tnie_sim - pnie_sim # intmed in the 4 way decomposition

         pie_sim <- pnie_sim # pie in the 4 way decomposition

         cde_prop_sim <- cde_sim/te_sim # proportion cde in the 4 way decomposition

         intref_prop_sim <- intref_sim/te_sim # proportion intref in the 4 way decomposition

         intmed_prop_sim <- intmed_sim/te_sim # proportion intmed in the 4 way decomposition

         pie_prop_sim <- pie_sim/te_sim # proportion pie in the 4 way decomposition

         overall_pm_sim <- (pie_sim + intmed_sim)/te_sim # overall pm in the 4 way decomposition

         overall_int_sim <- (intref_sim + intmed_sim)/te_sim # overall int in the 4 way decomposition

         overall_pe_sim <- (intref_sim + intmed_sim + pie_sim)/te_sim # overall pe in the 4 way decomposition

         # calculate the mean of each simulated effects
         for (effect in c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm", "intref", "intmed",
                          "pie", "cde_prop", "intref_prop", "intmed_prop", "pie_prop",
                          "overall_pm", "overall_int", "overall_pe")) {

           mid <- mean(get(paste0(effect, "_sim")))

           assign(effect, mid)

         }

         # calculate the se of each simulated effects
         for (effect_se in c("cde_se", "pnde_se", "tnde_se", "pnie_se", "tnie_se", "te_se", "pm_se",
                             "intref_se", "intmed_se", "pie_se", "cde_prop_se", "intref_prop_se",
                             "intmed_prop_se", "pie_prop_se", "overall_pm_se",
                             "overall_int_se", "overall_pe_se")) {

           mid <- sd(get(stringr::str_replace(effect_se, "_se", "_sim")))

           assign(effect_se, mid)

         }

         # 3 way decomposition effects and se's
         decomp3way <-  c(cde = cde, cde_se = cde_se,
                          pnde = pnde, pnde_se = pnde_se,
                          tnde = tnde, tnde_se = tnde_se,
                          pnie = pnie, pnie_se = pnie_se,
                          tnie = tnie, tnie_se = tnie_se,
                          te = te, te_se = te_se,
                          pm = pm, pm_se = pm_se)

         # 4 way decomposition effects and se's
         decomp4way <-  c(cde = cde, cde_se = cde_se, intref = intref, intref_se = intref_se,
                          intmed = intmed, intmed_se = intmed_se, pie = pie, pie_se = pie_se,
                          te = te, te_se = te_se, cde_prop = cde_prop, cde_prop_se = cde_prop_se,
                          intref_prop = intref_prop, intref_prop_se = intref_prop_se,
                          intmed_prop = intmed_prop, intmed_prop_se = intmed_prop_se,
                          pie_prop = pie_prop, pie_prop_se = pie_prop_se,
                          overall_pm = overall_pm, overall_pm_se = overall_pm_se,
                          overall_int = overall_int, overall_int_se = overall_int_se,
                          overall_pe = overall_pe, overall_pe_se = overall_pe_se)

         out <- list(decomp3way = decomp3way, decomp4way = decomp4way,
                     sims = data.frame(cde = cde_sim, pnde = pnde_sim, tnde = tnde_sim,
                                       pnie = pnie_sim, tnie = tnie_sim, intref = intref_sim,
                                       intmed = intmed_sim, pie_sim = pnie_sim))

       }

     # nonlinear outcome
     if (yreg != "linear") {

         cde_rr_sim <- EY1m_sim/EY0m_sim # cde in the 3 way decomposition

         pnde_rr_sim <- EY10_sim/EY00_sim # pnde in the 3 way decomposition

         tnie_rr_sim <- EY11_sim/EY10_sim # tnie in the 3 way decomposition

         tnde_rr_sim <- EY11_sim/EY01_sim # tnde in the 3 way decomposition

         pnie_rr_sim <- EY01_sim/EY00_sim # pnie in the 3 way decomposition

         te_rr_sim <- tnie_rr_sim * pnde_rr_sim # te in the 3 way or 4 way decomposition

         pm_sim <- (pnde_rr_sim * (tnie_rr_sim - 1)) / (pnde_rr_sim * tnie_rr_sim - 1) # pm in the 3 way decomposition

         cde_err_sim <- (EY1m_sim-EY0m_sim)/EY00_sim # excess RR due to cde in the 4 way decomposition

         intref_err_sim <- pnde_rr_sim - 1 - cde_err_sim # excess RR due to intref in the 4 way decomposition

         intmed_err_sim <- tnie_rr_sim * pnde_rr_sim - pnde_rr_sim - pnie_rr_sim + 1 # excess RR due to intmed in the 4 way decomposition

         pie_err_sim <- pnie_rr_sim - 1 # excess RR due to pie in the 4 way decomposition

         total_err_sim <- EY11_sim / EY00_sim - 1 # total excess RR in the 4 way decomposition

         cde_prop_sim <- cde_err_sim/total_err_sim # proportion cde

         intmed_prop_sim <- intmed_err_sim/total_err_sim # proportion intmed

         intref_prop_sim <- intref_err_sim/total_err_sim # proportion intref

         pie_prop_sim <- pie_err_sim/total_err_sim # proportion pie

         overall_pm_sim <- (pie_err_sim+intmed_err_sim)/total_err_sim # overall pm in the 4 way decomposition

         overall_int_sim <- (intref_err_sim+intmed_err_sim)/total_err_sim # overall int in the 4 way decomposition

         overall_pe_sim <- (intref_err_sim+intmed_err_sim+pie_err_sim)/total_err_sim # overall pe in the 4 way decomposition

         # calculate the mean of each simulated effects
         for (effect in c("cde_rr", "pnde_rr", "tnde_rr", "pnie_rr", "tnie_rr", "te_rr", "pm",
                          "cde_err", "intref_err", "intmed_err", "pie_err", "total_err",
                          "cde_prop", "intref_prop", "intmed_prop", "pie_prop",
                          "overall_pm", "overall_int", "overall_pe")) {

           mid <- mean(get(paste0(effect, "_sim")))

           assign(effect, mid)

         }

         # calculate the se of each simulated effects
         for (effect_se in c("cde_rr_se", "pnde_rr_se", "tnde_rr_se", "pnie_rr_se", "tnie_rr_se",
                             "te_rr_se", "pm_se", "cde_err_se", "intref_err_se", "intmed_err_se",
                             "pie_err_se", "total_err_se", "cde_prop_se", "intref_prop_se",
                             "intmed_prop_se", "pie_prop_se",
                             "overall_pm_se", "overall_int_se", "overall_pe_se")) {

           mid <- sd(get(stringr::str_replace(effect_se, "_se", "_sim")))

           assign(effect_se, mid)

         }

         # 3 way decomposition effects and se's
         decomp3way <-  c(cde_rr = cde_rr, cde_rr_se = cde_rr_se,
                          pnde_rr = pnde_rr, pnde_rr_se = pnde_rr_se,
                          tnde_rr = tnde_rr, tnde_rr_se = tnde_rr_se,
                          pnie_rr = pnie_rr, pnie_rr_se = pnie_rr_se,
                          tnie_rr = tnie_rr, tnie_rr_se = tnie_rr_se,
                          te_rr = te_rr, te_rr_se = te_rr_se,
                          pm = pm, pm_se = pm_se)

         # 4 way decomposition effects and se's
         decomp4way <-  c(cde_err = cde_err, cde_err_se = cde_err_se, intref_err = intref_err,
                          intref_err_se = intref_err_se, intmed_err = intmed_err,
                          intmed_err_se = intmed_err_se, pie_err = pie_err, pie_err_se = pie_err_se,
                          total_err = total_err, total_err_se = total_err_se,
                          cde_prop = cde_prop, cde_prop_se = cde_prop_se,
                          intref_prop = intref_prop, intref_prop_se = intref_prop_se,
                          intmed_prop = intmed_prop, intmed_prop_se = intmed_prop_se,
                          pie_prop = pie_prop, pie_prop_se = pie_prop_se,
                          overall_pm = overall_pm, overall_pm_se = overall_pm_se,
                          overall_int = overall_int, overall_int_se = overall_int_se,
                          overall_pe = overall_pe, overall_pe_se = overall_pe_se)

         out <- list(decomp3way = decomp3way, decomp4way = decomp4way,
                     sims = data.frame(cde_rr = cde_rr_sim, pnde_rr = pnde_rr_sim, tnde_rr = tnde_rr_sim,
                                       pnie_rr = pnie_rr_sim, tnie_rr = tnie_rr_sim,
                                       cde_err_sim = cde_err_sim, intref_err_sim = intref_err_sim,
                                       intmed_err_sim = intmed_err_sim, pie_err_sim = pie_err_sim))

         }
       }

   return(out)

  }
