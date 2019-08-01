causal_mediation <- function(data, nway = c(3,4), method = c("delta", "bootstrap", "simulation"),
                             nboot = 100, conf = 0.95, outcome, treatment, mediator, covariates,
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
                                 event = event, mreg = mreg, yreg = yreg, data = data_boot)

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

######################################Three way decomposition####################################

  if (nway == 3) {

###############################Three way decomposition: Delta method#############################

    if (method == "delta") {

      effect_estimates <- bootstrap_3way_step(data = data, indices = c(1:nrow(data)), outcome = outcome,
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

        te_formula <- as.formula("~x1+x2")

        pm_formula <- as.formula("~x2/(x1+x3)")

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
        }

        te_formula <- as.formula("~x1*x2")

        pm_formula <- as.formula("~x1*(x2-1)/(x1*x2-1)")

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

      cde_se_delta <- msm::deltamethod(cde_formula, thetas, vcov_thetas)

      pnde_se_delta <- msm::deltamethod(pnde_formula, c(thetas, betas), vcov_block)

      tnde_se_delta <- msm::deltamethod(tnde_formula, c(thetas, betas), vcov_block)

      pnie_se_delta <- msm::deltamethod(pnie_formula, c(thetas, betas), vcov_block)

      tnie_se_delta <- msm::deltamethod(tnie_formula, c(thetas, betas), vcov_block)

      te_se_delta <- msm::deltamethod(te_formula, c(effect_estimates["pnde"], effect_estimates["tnie"]),
                                        Matrix::bdiag(pnde_se_delta^2, tnie_se_delta^2))

      pm_se_delta <- msm::deltamethod(pm_formula, c(effect_estimates["pnde"], effect_estimates["tnie"],
                                                      effect_estimates["te"]),
                                        Matrix::bdiag(pnde_se_delta^2, tnie_se_delta^2, te_se_delta^2))

      out <- c(cde = unname(effect_estimates["cde"]), cde_se = cde_se_delta,
                 pnde = unname(effect_estimates["pnde"]), pnde_se = pnde_se_delta,
                 tnde = unname(effect_estimates["tnde"]), tnde_se = tnde_se_delta,
                 pnie = unname(effect_estimates["pnie"]), pnie_se = pnie_se_delta,
                 tnie = unname(effect_estimates["tnie"]), tnie_se = tnie_se_delta,
                 te = unname(effect_estimates["te"]), te_se = te_se_delta,
                 pm = unname(effect_estimates["pm"]), pm_se = pm_se_delta)

      class(out) <- c("delta_out", "data.frame")

    }

###############################Three way decomposition: Bootstrapping#############################

    if (method == "bootstrap") {

      out <- boot::boot(data = data, statistic = bootstrap_3way_step, R = nboot,
                                      outcome = outcome, treatment = treatment, mediator = mediator,
                                      covariates = covariates, vecc = vecc, interaction = interaction,
                                      event = event, mreg = mreg, yreg = yreg,
                                      m_star = m_star, a_star = a_star, a = a)

      class(out) <- c("boot_out", "data.frame")

    }
}

#######################################Four way decomposition######################################

  if (nway == 4) {

#################################Four way decomposition: delta method##############################

    if (method == "delta") {

      effect_estimates <- bootstrap_4way_step(data = data, indices = c(1:nrow(data)), outcome = outcome,
                                            treatment = treatment, mediator = mediator,
                                            covariates = covariates, vecc = vecc,
                                            interaction = interaction, event = event,
                                            m_star = m_star, a_star = a_star, a = a,
                                            mreg = mreg, yreg = yreg)

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

          intref_formula <- stringr::str_replace_all(
            paste0("~", theta3, "*(", beta0, "+", beta1, "*a_star", "+", XC,
                   "-m_star)*(a-a_star)"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

          intmed_formula <- stringr::str_replace_all(
            paste0("~", theta3, "*", beta1, "*(a-a_star)^2"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

          pie_formula <- stringr::str_replace_all(
            paste0("~(", theta2, "*", beta1, "+", theta3, "*", beta1, "*a_star)*(a-a_star)"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

        } else if (mreg=="logistic") {

          cde_formula <- stringr::str_replace_all(
            paste0("~(", theta1, "+", theta3, "*m_star)*(a-a_star)"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

          intref_formula <- stringr::str_replace_all(
            paste0("~", theta3, "*(a-a_star)*(exp(", beta0, "+", beta1, "*a_star",
                   "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, "))-m_star)"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))


          intmed_formula <- stringr::str_replace_all(
            paste0("~", theta3, "*(a-a_star)*(exp(", beta0, "+", beta1, "*a",
                   "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a", "+", XC, "))-exp(",
                   beta0, "+", beta1, "*a_star",
                   "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, ")))"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

          pie_formula <- stringr::str_replace_all(
            paste0("~(", theta2, "+", theta3, "*a_star)*(exp(", beta0, "+", beta1, "*a",
                   "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a", "+", XC, "))-exp(",
                   beta0, "+", beta1, "*a_star",
                   "+", XC, ")/(1+exp(", beta0, "+", beta1, "*a_star", "+", XC, ")))"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

        }

        j <- length(vecc)

        if (j > 0) {

          for (formula in c("cde_formula", "intref_formula", "intmed_formula", "pie_formula")) {

            mid <- get(formula)

            for (i in 1:j) {

              mid <-   stringr::str_replace_all(mid,paste("vecc", i, sep = "_"),
                                                as.character(vecc[i]))

            }

            assign(formula, as.formula(mid))

          }
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

        out <- c(cde = unname(effect_estimates["cde"]), cde_se_delta = cde_se_delta,
                 intref = unname(effect_estimates["intref"]), intref_se_delta = intref_se_delta,
                 intmed = unname(effect_estimates["intmed"]), intmed_se_delta = intmed_se_delta,
                 pie = unname(effect_estimates["pie"]), pie_se_delta = pie_se_delta,
                 te = unname(effect_estimates["te"]), te_se_delta = te_se_delta,
                 cde_prop = unname(effect_estimates["cde_prop"]), cde_prop_se_delta = cde_prop_se_delta,
                 intref_prop = unname(effect_estimates["intref_prop"]), intref_prop_se_delta = intref_prop_se_delta,
                 intmed_prop = unname(effect_estimates["intmed_prop"]), intmed_prop_se_delta = intmed_prop_se_delta,
                 pie_prop = unname(effect_estimates["pie_prop"]), pie_prop_se_delta = pie_prop_se_delta,
                 overall_pm = unname(effect_estimates["overall_pm"]), overall_pm_se_delta = overall_pm_se_delta,
                 overall_int = unname(effect_estimates["overall_int"]), overall_int_se_delta = overall_int_se_delta,
                 overall_pe = unname(effect_estimates["overall_pe"]), overall_pe_se_delta = overall_pe_se_delta)

        class(out) <- c("delta_out", "data.frame")

      } else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                              "negbin", "coxph", "aft_exp", "aft_weibull")) {

        if (mreg=="linear") {
          cde_comp_formula <- stringr::str_replace_all(
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

          intref_comp_formula <- stringr::str_replace_all(
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
                        "\\bvariance\\b" = as.character(variance)))

          intmed_comp_formula <- stringr::str_replace_all(
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
                        "\\bvariance\\b" = as.character(variance)))

          pie_comp_formula <- stringr::str_replace_all(
            paste0("~exp((", theta2, "*", beta1, "+", theta3, "*", beta1,
                   "*a_star)*(a-a_star))-1"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))

          total_rr_formula <- stringr::str_replace_all(
            paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                   "*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                   theta3, "^2*variance*(a^2-a_star^2))*exp((", theta2, "*",
                   beta1, "+", theta3, "*", beta1, "*a)*(a-a_star))"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star),
                        "\\bvariance\\b" = as.character(variance)))

        } else if (mreg=="logistic") {

          cde_comp_formula <- stringr::str_replace_all(
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

          intref_comp_formula <- stringr::str_replace_all(
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
                        "\\bm_star\\b" = as.character(m_star)))

          intmed_comp_formula <- stringr::str_replace_all(
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
                        "\\bm_star\\b" = as.character(m_star)))

          pie_comp_formula <- stringr::str_replace_all(
            paste0("~((1+exp(", beta0, "+", beta1, "*a_star+", XC, "))*(1+exp(",
                   beta0, "+", beta1, "*a+", XC, "+", theta2, "+", theta3,
                   "*a_star)))/((1+exp(", beta0, "+", beta1, "*a+", XC,
                   "))*(1+exp(", beta0, "+", beta1, "*a_star+", XC, "+", theta2,
                   "+", theta3, "*a_star)))-1"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a)))

          total_rr_formula <- stringr::str_replace_all(
            paste0("~exp(", theta1, "*a)*(1+exp(", beta0, "+", beta1, "*a_star",
                   "+", XC, "))*(1+exp(", beta0, "+", beta1, "*a", "+", XC, "+",
                   theta2, "+", theta3, "*a))/(exp(", theta1, "*a_star)*(1+exp(",
                   beta0, "+", beta1, "*a+", XC, "))*(1+exp(", beta0, "+",
                   beta1, "*a_star+", XC, "+", theta2, "+", theta3, "*a_star)))"),
            pattern = c("\\ba_star\\b" = as.character(a_star),
                        "\\ba\\b" = as.character(a),
                        "\\bm_star\\b" = as.character(m_star)))
        }

        j <- length(vecc)

        if (j > 0) {

          for (formula in c("cde_comp_formula", "intref_comp_formula", "intmed_comp_formula",
                            "pie_comp_formula", "total_rr_formula")) {

            mid <- get(formula)

            for (i in 1:j) {

              mid <-   stringr::str_replace_all(mid,paste("vecc", i, sep = "_"),
                                                as.character(vecc[i]))

            }

            assign(formula, as.formula(mid))

          }
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

        cde_err_se_delta <- msm::deltamethod(as.formula("~x1*x2/x3"),
                                             c(effect_estimates["cde_comp"], effect_estimates["total_err"],
                                               (effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                  effect_estimates["intmed_comp"]+effect_estimates["pie_comp"])),
                                             Matrix::bdiag(cde_comp_se_delta^2, total_err_se_delta^2,
                                                           tcomp_se_delta^2))

        intmed_err_se_delta <- msm::deltamethod(as.formula("~x1*x2/x3"),
                                                c(effect_estimates["intmed_comp"], effect_estimates["total_err"],
                                                  effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                    effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                                Matrix::bdiag(intmed_comp_se_delta^2, total_err_se_delta^2,
                                                              tcomp_se_delta^2))

        intref_err_se_delta <- msm::deltamethod(as.formula("~x1*x2/x3"),
                                                c(effect_estimates["intref_comp"], effect_estimates["total_err"],
                                                  effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                    effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                                Matrix::bdiag(intref_comp_se_delta^2, total_err_se_delta^2,
                                                              tcomp_se_delta^2))

        pie_err_se_delta <- msm::deltamethod(as.formula("~x1*x2/x3"),
                                             c(effect_estimates["pie_comp"], effect_estimates["total_err"],
                                               effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                                 effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
                                             Matrix::bdiag(pie_comp_se_delta^2, total_err_se_delta^2,
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

        out <- c(cde_comp = unname(effect_estimates["cde_comp"]), cde_comp_se_delta = cde_comp_se_delta,
                 intref_comp = unname(effect_estimates["intref_comp"]), intref_comp_se_delta = intref_comp_se_delta,
                 intmed_comp = unname(effect_estimates["intmed_comp"]), intmed_comp_se_delta = intmed_comp_se_delta,
                 pie_comp = unname(effect_estimates["pie_comp"]), pie_comp_se_delta = pie_comp_se_delta,
                 total_err = unname(effect_estimates["total_err"]), total_err_se_delta = total_err_se_delta,
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
                 overall_pe = unname(effect_estimates["overall_pe"]), overall_pe_se_delta = overall_pe_se_delta)

        class(out) <- c("delta_out", "data.frame")

        }
      }

#####################################Four way decomposition: bootstrapping##############################

    if (method == "bootstrap") {

      out <- boot::boot(data = data, statistic = bootstrap_4way_step, R = nboot,
                        outcome = outcome, treatment = treatment, mediator = mediator,
                        covariates = covariates, vecc = vecc, interaction = interaction,
                        event = event, mreg = mreg, yreg = yreg,
                        m_star = m_star, a_star = a_star, a = a)


      class(out) <- c("boot_out", "data.frame")

    }
  }

##########################Simulation-based decomposition: 3 way and 4 way############################

  if (method == "simulation") {

     mdesign_a <- as.matrix(cbind(data.frame(intercept = rep(1,n),treatment = c(rep(a,n))),
                                    data[,covariates]))

     mdesign_a_star <- as.matrix(cbind(data.frame(intercept = rep(1,n),treatment = c(rep(a_star,n))),
                                         data[,covariates]))

     thetas_sim <- as.matrix(mvtnorm::rmvnorm(nsims, mean = thetas, sigma = vcov_thetas))

       betas_sim <- as.matrix(mvtnorm::rmvnorm(nsims, mean = betas, sigma = vcov_betas))

       cde_sim <- c()

       pnde_sim <- c()

       tnie_sim <- c()

       tnde_sim <- c()

       pnie_sim <- c()

       intref_sim <- c()

       intmed_sim <- c()

       for (j in 1:nsims) {

         m_mean_a <- mdesign_a%*%betas_sim[j,]

         m_mean_a_star <- mdesign_a_star%*%betas_sim[j,]

         covariatesTerm_sim <- rowSums(data[,covariates]*thetas_sim[j,covariates])

         #simulate potential values of the mediator
         if (mreg == "linear") {

           m_sims_a <- matrix(rnorm(n*nsims, mean = m_mean_a, sd = sqrt(variance)),
                              ncol = nsims) #n by nsims

           m_sims_a_star <- matrix(rnorm(n*nsims, mean = m_mean_a_star, sd = sqrt(variance)),
                                   ncol = nsims) #n by nsims

         } else if (mreg == "logistic")  {

           m_sims_a <- matrix(rbinom(n*nsims, size = 1,
                                     prob = boot::inv.logit(m_mean_a)), ncol = nsims) #n by nsims

           m_sims_a_star <- matrix(rbinom(n*nsims, size = 1,
                                          prob = boot::inv.logit(m_mean_a_star)), ncol = nsims) #n by nsims

         }

         theta_interaction <- ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")], 0)

         #calculate expected values of the outcome
         y0m <- thetas_sim[j,"(Intercept)"]+a_star*thetas_sim[j,treatment]+
           m_star*thetas_sim[j,mediator]+covariatesTerm_sim+
           a_star*m_star*theta_interaction

         y1m <- thetas_sim[j,"(Intercept)"]+a*thetas_sim[j,treatment]+
           m_star*thetas_sim[j,mediator]+covariatesTerm_sim+
           a*m_star*theta_interaction

         y00 <- sapply(1:nsims,FUN=function(k)
           thetas_sim[j,"(Intercept)"]+a_star*thetas_sim[j,treatment]+
             m_sims_a_star[,k]*thetas_sim[j,mediator]+covariatesTerm_sim+
             a_star*m_sims_a_star[,k]*theta_interaction)

         y01 <- sapply(1:nsims,FUN =function(k)
           thetas_sim[j,"(Intercept)"]+a_star*thetas_sim[j,treatment]+
             m_sims_a[,k]*thetas_sim[j,mediator]+covariatesTerm_sim+
             a_star*m_sims_a[,k]*theta_interaction)

         y10 <- sapply(1:nsims,FUN =function(k)
           thetas_sim[j,"(Intercept)"]+a*thetas_sim[j,treatment]+
             m_sims_a_star[,k]*thetas_sim[j,mediator]+covariatesTerm_sim+
             a*m_sims_a_star[,k]*theta_interaction)

         y11 <- sapply(1:nsims,FUN =function(k)
           thetas_sim[j,"(Intercept)"]+a*thetas_sim[j,treatment]+
             m_sims_a[,k]*thetas_sim[j,mediator]+covariatesTerm_sim+
             a*m_sims_a[,k]*theta_interaction)

         if (yreg == "logistic") {

           y0m <- boot::inv.logit(y0m)

           y1m <- boot::inv.logit(y1m)

           y00 <- boot::inv.logit(y00)

           y01 <- boot::inv.logit(y01)

           y10 <- boot::inv.logit(y10)

           y11 <- boot::inv.logit(y11)

         }

         if (yreg %in% c("loglinear", "poisson", "quasipoisson", "negbin",
                         "coxph", "aft_exp", "aft_weibull")) {

           y0m <- exp(y0m)

           y1m <- exp(y1m)

           y00 <- exp(y00)

           y01 <- exp(y00)

           y10 <- exp(y10)

           y11 <- exp(y11)

         }

         EY0m_sim <- c(EY0m_sim, sum(y0m)/n)

         EY1m_sim <- c(EY1m_sim, sum(y1m)/n)

         EY00_sim <- c(EY00_sim, sum(y00)/(n*nsims))

         EY10_sim <- c(EY10_sim, sum(y10)/(n*nsims))

         EY01_sim <- c(EY01_sim, sum(y01)/(n*nsims))

         EY11_sim <- c(EY11_sim, sum(y11)/(n*nsims))

       }

       if (yreg == "linear") {

         cde_sim <- EY1m_sim - EY0m_sim

         pnde_sim <- EY10_sim - EY00_sim

         tnie_sim <- EY11_sim - EY10_sim

         tnde_sim <- EY11_sim - EY01_sim

         pnie_sim <- EY01_sim - EY00_sim

         te_sim <- tnie_sim + pnde_sim

         pm_sim <- tnie_sim / (pnde_sim + te_sim)

         intref_sim <- pnde_sim - cde_sim

         intmed_sim <- tnie_sim - pnie_sim

         cde_prop_sim <- cde_sim/te_sim

         intref_prop_sim <- intref_sim/te_sim

         intmed_prop_sim <- intmed_sim/te_sim

         pie_prop_sim <- pnie_sim/te_sim

         overall_pm_sim <- (pie_sim + intmed_sim)/te_sim

         overall_int_sim <- (intref_sim + intmed_sim)/te_sim

         overall_pe_sim <- (intref_sim + intmed_sim + pie_sim)/te_sim

         for (effect in c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm", "intref", "intmed",
                          "cde_prop", "intref_prop", "intmed_prop", "pie_prop",
                          "overall_pm", "overall_int", "overall_pe")) {

           mid <- mean(get(paste0(effect, "_sim")))

           assign(effect, mid)

         }

         for (effect_se in c("cde_se", "pnde_se", "tnde_se", "pnie_se", "tnie_se", "te_se", "pm_se",
                             "intref_se", "intmed_se", "cde_prop_se", "intref_prop_se",
                             "intmed_prop_se", "pie_prop_se", "overall_pm_se",
                             "overall_int_se", "overall_pe_se")) {

           mid <- sd(get(stringr::str_replace(effect_se, "_se", "_sim")))

           assign(effect_se, mid)

         }

         decomp3way <-  c(cde = cde, cde_se = cde_se,
                          pnde = pnde, pnde_se = pnde_se,
                          tnde = tnde, tnde_se = tnde_se,
                          pnie = pnie, pnie_se = pnie_se,
                          tnie = tnie, tnie_se = tnie_se,
                          te = te, te_se = te_se,
                          pm = pm, pm_se = pm_se)

         decomp4way <-  c(cde = cde, cde_se = cde_se, intref = intref, intref_se = intref_se,
                          intmed = intmed, intmed_se = intmed_se, pie = pnie, pie_se = pnie_se,
                          te = te, te_se = te_se, cde_prop = cde_prop, cde_prop_se = cde_prop_se,
                          intref_prop = intref_prop, intref_prop_se = intref_prop_se,
                          intmed_prop = intmed_prop, intmed_prop_se = intmed_prop_se,
                          pie_prop = pie_prop, pie_prop_se = pie_prop_se,
                          overall_pm = overall_pm, overall_pm_se = overall_pm_se,
                          overall_int = overall_int, overall_int_se = overall_int_se,
                          overall_pe = overall_pe, overall_pe_se = overall_pe_se)

         class(decomp3way) <- c("simulation_out", "data.frame")

         class(decomp4way) <- c("simulation_out", "data.frame")

         out <- list(decomp3way = decomp3way, decomp4way = decomp4way,
                     sims = data.frame(cde = cde_sim, pnde = pnde_sim, tnde = tnde_sim,
                                       pnie = pnie_sim, tnie = tnie_sim, intref = intref_sim,
                                       intmed = intmed_sim, pie_sim = pnie_sim))

       }

       if (yreg != "linear") {

         cde_sim <- EY1m_sim/EY0m_sim

         pnde_sim <- EY10_sim/EY00_sim

         tnie_sim <- EY11_sim/EY10_sim

         tnde_sim <- EY11_sim/EY01_sim

         pnie_sim <- EY01_sim/EY00_sim

         te_sim <- tnie_sim * pnde_sim

         pm_sim <- (pnde_sim * (tnie_sim - 1)) / (pnde_sim * tnie_sim - 1)

         cde_err_sim <- (EY1m_sim-EY0m_sim)/EY00_sim

         intref_err_sim <- pnde_sim - 1 - cde_comp_sim

         intmed_err_sim <- tnie_sim * pnde_sim - pnde_sim - pnie_sim + 1

         pie_err_sim <- pnie_sim - 1

         total_err_sim <- EY11_sim / EY00_sim - 1

         cde_prop_sim <- cde_err_sim/total_err_sim

         intmed_prop_sim <- intmed_err_sim/total_err_sim

         intref_prop_sim <- intref_err_sim/total_err_sim

         pie_prop_sim <- pie_err_sim/total_err_sim

         overall_pm_sim <- (pie_err_sim+intmed_err_sim)/total_err_sim

         overall_int_sim <- (intref_err_sim+intmed_err_sim)/total_err_sim

         overall_pe_sim <- (intref_err_sim+intmed_err_sim+pie_err_sim)/total_err_sim

         for (effect in c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm", "cde_err",
                          "intref_err", "intmed_err", "pie_err", "total_err",
                          "cde_prop", "intref_prop", "intmed_prop", "pie_prop",
                          "overall_pm", "overall_int", "overall_pe")) {

           mid <- mean(get(paste0(effect, "_sim")))

           assign(effect, mid)

         }

         for (effect_se in c("cde_se", "pnde_se", "tnde_se", "pnie_se", "tnie_se", "te_se",
                             "pm_se", "cde_err_se", "intref_err_se", "intmed_err_se", "pie_err_se",
                             "total_err_se", "cde_prop_se", "intref_prop_se", "intmed_prop_se",
                             "pie_prop_se", "overall_pm_se", "overall_int_se", "overall_pe_se")) {

           mid <- sd(get(stringr::str_replace(effect_se, "_se", "_sim")))

           assign(effect_se, mid)

         }

         decomp3way <-  c(cde = cde, cde_se = cde_se,
                          pnde = pnde, pnde_se = pnde_se,
                          tnde = tnde, tnde_se = tnde_se,
                          pnie = pnie, pnie_se = pnie_se,
                          tnie = tnie, tnie_se = tnie_se,
                          te = te, te_se = te_se,
                          pm = pm, pm_se = pm_se)

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

         class(decomp3way) <- c("simulation_out", "data.frame")

         class(decomp4way) <- c("simulation_out", "data.frame")

         out <- list(decomp3way = decomp3way, decomp4way = decomp4way,
                     sims = data.frame(cde = cde_sim, pnde = pnde_sim, tnde = tnde_sim,
                                       pnie = pnie_sim, tnie = tnie_sim,
                                       cde_err_sim = cde_err_sim, intref_err_sim = intref_err_sim,
                                       intmed_err_sim = intmed_err_sim, pie_err_sim = pie_err_sim))

         }
       }

   return(out)

  }






