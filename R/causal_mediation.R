causal_mediation <- function(data, nway = c(3,4), type = c("delta", "bootstrap"), nboot = 100, conf = 0.95,
                             outcome, treatment, mediator, covariates, vecc = NULL,
                             interaction = TRUE, event = NULL, m_star = 0, a_star = 1, a = 0,
                             mreg = c("linear", "logistic"),
                             yreg = c("linear", "logistic", "loglinear", "poisson",
                                      "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) {

  delta_3way <- function(coef = NULL, vecc = c(), m_star = NULL, a_star = NULL, a = NULL) {

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

    effect_estimates <- bootstrap_step(data = data, indices = c(1:nrow(data)), outcome = outcome,
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

        pnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a_star+", XC, "))*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a+", XC, "))*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        pnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta2, "*", beta1, "+", theta3, "*", beta1,
                 "*a_star*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta2, "*", beta1, "+", theta3, "*", beta1,
                 "*a*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

      } else if (mreg=="logistic") {

        cde_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta1, "+", theta3, "*m_star)*(a-a_star)"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        pnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~", theta1, "*(a-a_star)+(", theta3, "*(a-a_star))*exp(",
                 beta0, "+", beta1, "*a_star+", XC, ")/(1+exp(",
                 beta0, "+", beta1, "*a_star+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~", theta1, "*(a-a_star)+(", theta3, "*(a-a_star))*exp(",
                 beta0, "+", beta1, "*a+", XC, ")/(1+exp(",
                 beta0, "+", beta1, "*a+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        pnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta2, "+", theta3, "*a_star)*(exp(", beta0,
                 "+", beta1, "*a+", XC, ")/(1+exp(", beta0, "+",
                 beta1, "*a+", XC, "))-exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")/(1+exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~(", theta2, "+", theta3, "*a)*(exp(", beta0,
                 "+", beta1, "*a+", XC, ")/(1+exp(", beta0, "+",
                 beta1, "*a+", XC, "))-exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")/(1+exp(", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

      }

      te_formula <- as.formula("~x1+x2")

      pm_formula <- as.formula("~x2/(x1+x3)")

    } else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                            "negbin", "coxph", "aft_exp", "aft_weibull")) {

      if (mreg=="linear") {

        cde_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*m_star)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        pnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a_star+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                 theta3, "^2*variance*(a^2-a_star^2))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star),
                      "\\bvariance\\b" = as.character(variance))))

        tnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*(", beta0, "+", beta1,
                 "*a+", XC, "+", theta2, "*variance))*(a-a_star)+0.5*",
                 theta3, "^2*variance*(a^2-a_star^2))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star),
                      "\\bvariance\\b" = as.character(variance))))

        pnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp((", theta2, "*", beta1, "+", theta2, "*",
                 beta1, "*a_star)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp((", theta2, "*", beta1, "+", theta2, "*",
                 beta1, "*a)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

      } else if (mreg=="logistic") {

        cde_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp((", theta1, "+", theta3, "*m_star)*(a-a_star))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        pnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp(", theta1, "*(a-a_star))*(1+exp(", theta2, "+",
                 theta3, "*a+", beta0, "+", beta1, "*a_star+", XC,
                 "))/(1+exp(", theta2, "+", theta3, "*a_star+", beta0,
                 "+", beta1, "*a_star+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnde_formula <- as.formula(stringr::str_replace_all(
          paste0("~exp(", theta1, "*(a-a_star))*(1+exp(", theta2, "+",
                 theta3, "*a+", beta0, "+", beta1, "*a+", XC,
                 "))/(1+exp(", theta2, "+", theta3, "*a+", beta0,
                 "+", beta1, "*a_star+", XC, "))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        pnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                 "))*(1+exp(", theta2, "+", theta3, "*a_star+",
                 beta0, "+", beta1, "*a+", XC, "))/((1+exp(",
                 beta0, "+", beta1, "*a+", XC, "))*(1+exp(",
                 theta2, "+", theta3, "*a_star+", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))

        tnie_formula <- as.formula(stringr::str_replace_all(
          paste0("~(1+exp(", beta0, "+", beta1, "*a_star+", XC,
                 "))*(1+exp(", theta2, "+", theta3, "*a+",
                 beta0, "+", beta1, "*a+", XC, "))/((1+exp(",
                 beta0, "+", beta1, "*a+", XC, "))*(1+exp(",
                 theta2, "+", theta3, "*a+", beta0, "+", beta1,
                 "*a_star+", XC, ")))"),
          pattern = c("\\ba_star\\b" = as.character(a_star),
                      "\\ba\\b" = as.character(a),
                      "\\bm_star\\b" = as.character(m_star))))
      }

      te_formula <- as.formula("~x1*x2")

      pm_formula <- as.formula("~x1*(x2-1)/(x1*x2-1)")

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

    return(c(cde = effect_estimates["cde"], cde_se_delta = cde_se_delta,
             pnde = effect_estimates["pnde"], pnde_se_delta = pnde_se_delta,
             tnde = effect_estimates["tnde"], tnde_se_delta = tnde_se_delta,
             pnie = effect_estimates["pnie"], pnie_se_delta = pnie_se_delta,
             tnie = effect_estimates["tnie"], tnie_se_delta = tnie_se_delta,
             te = effect_estimates["te"], te_se_delta = te_se_delta,
             pm = effect_estimates["pm"], pm_se_delta = pm_se_delta))

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

      return(c(cde = effect_estimates["cde"], cde_se_delta = cde_se_delta,
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
                      "\\ba\\b" = as.character(a))))

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

      cde_err_se_delta <- msm::deltamethod(as.formula("~x1*x2/x3"),
                                           c(effect_estimates["cde_comp"], effect_estimates["total_err"],
                                             effect_estimates["cde_comp"]+effect_estimates["intref_comp"]+
                                               effect_estimates["intmed_comp"]+effect_estimates["pie_comp"]),
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

      return(c(cde_comp = unname(effect_estimates["cde_comp"]), cde_comp_se_delta = cde_comp_se_delta,
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
               overall_pe = unname(effect_estimates["overall_pe"]), overall_pe_se_delta = overall_pe_se_delta))
    }
  }


  bootstrap_3way_step <- function(data, indices, outcome, treatment, mediator, covariates, vecc,
                                  interaction, event, mreg, yreg, m_star, a_star, a) {

    data_boot <- data[indices, ]

    formulas <- create_formulas(outcome = outcome, treatment = treatment,
                                mediator = mediator, covariates = covariates, interaction = interaction,
                                event = event, mreg = mreg, yreg = yreg)

    regressions <- run_regressions(formulas = formulas, outcome = outcome, treatment = treatment,
                                   mediator = mediator, covariates = covariates, interaction = interaction,
                                   event = event, mreg = mreg, yreg = yreg, data = data_boot)

    coef <- get_coef(regressions = regressions, outcome = outcome, treatment = treatment,
                     mediator = mediator, covariates = covariates, interaction = interaction,
                     event = event, mreg = mreg, yreg = yreg)

    thetas <- coef$thetas

    betas <- coef$betas

    variance <- coef$variance

    covariatesTerm <- ifelse(!is.null(vecc), sum(betas[covariates]*t(vecc)), 0)

    interactionTerm <- ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")], 0)

    if (mreg != "linear" & yreg != "linear") {

      cde <- unname(exp(thetas[treatment] + interactionTerm*m_star * (a - a_star)))

      pnde <- unname((exp(thetas[treatment] * (a - a_star)) * (1 + exp(thetas[mediator] +
              interactionTerm * a + betas[1] + betas[treatment] * a_star +
              covariatesTerm))) /(1 + exp(thetas[mediator] + interactionTerm * a_star +
              betas[1] + betas[treatment] * a_star + covariatesTerm)))

      tnde <- unname((exp(thetas[treatment] * (a - a_star)) * (1 + exp(thetas[mediator] +
              interactionTerm * a + betas[1] + betas[treatment] * a + covariatesTerm))) /
              (1 + exp(thetas[mediator] + interactionTerm * a_star +  betas[1] +
              betas[treatment] * a + covariatesTerm)))

      pnie <- unname(((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
              (1 + exp(thetas[mediator] + interactionTerm * a_star + betas[1] +
              betas[treatment] * a + covariatesTerm))) /
              ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) *
              (1+ exp(thetas[mediator] + interactionTerm * a_star + betas[1] +
              betas[treatment] * a_star + covariatesTerm))))

      tnie <- unname(((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
              (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a +
              covariatesTerm))) / ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) *
              (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a_star +
              covariatesTerm))))

      te <- tnie * pnde

      pm <-  (pnde * (tnie - 1)) / (pnde * tnie - 1)

    } else if (mreg != "linear" & yreg == "linear") {

      cde <- unname((thetas[treatment] + interactionTerm * m_star * (a - a_star)))

      pnde <- unname(thetas[treatment] * (a - a_star) + interactionTerm*(a - a_star) *
              (exp(betas[1] + betas[treatment] * a_star + covariatesTerm) /
              (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))

      tnde <- unname(thetas[treatment] * (a - a_star) + interactionTerm*(a - a_star) *
              (exp(betas[1] + betas[treatment] * a + covariatesTerm) /
              (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm))))

      pnie <- unname((thetas[mediator]+interactionTerm*a_start) * (exp(betas[1] +
              betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] +
              betas[treatment] * a + covariatesTerm)) - exp(betas[1] + betas[treatment] * a_star +
              covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))

      tnie <- unname((thetas[mediator]+interactionTerm*a) * (exp(betas[1] +
              betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] +
              betas[treatment] * a + covariatesTerm)) - exp(betas[1] + betas[treatment] * a_star +
              covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))

      te <- tnie + pnde

      pm <- tnie / (pnde + te)

    } else if (mreg == "linear" & yreg != "linear") {

      cde <- unname(exp(thetas[treatment] + interactionTerm * m_star * (a - a_star)))

      pnde <- unname(exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star +
              covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
              0.5 * interactionTerm ^ 2 * variance * (a^2 - a_star ^ 2)))

      tnde <- unname(exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a +
              covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
              0.5 * interactionTerm ^ 2 * variance * (a^2 - a_star ^ 2)))

      pnie <- unname(exp((thetas[mediator] * betas[treatment] +
              interactionTerm * betas[treatment] * a_star) * (a - a_star)))

      tnie <- unname(exp((thetas[mediator] * betas[treatment] +
              interactionTerm * betas[treatment] * a) * (a - a_star)))

      te <- tnie * pnde

      pm <-  (pnde * (tnie - 1)) / (pnde * tnie - 1)

    } else if (mreg == "linear" & yreg == "linear") {

      cde <- unname((thetas[treatment] + ifelse(interaction,
             thetas[paste(treatment, mediator, sep = ":")] * m_star, 0)) * (a - a_star))

      pnde <- unname((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star +
              covariatesTerm))*(a - a_star))

      tnde <- unname((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a +
              covariatesTerm))*(a - a_star))

      pnie <- unname((thetas[mediator] * betas[treatment] +
              interactionTerm * betas[treatment] * a_star) * (a - a_star))

      tnie <- unname((thetas[mediator] * betas[treatment] +
              interactionTerm * betas[treatment] * a) * (a - a_star))

      te <- tnie + pnde

      pm <- tnie / (pnde + te)

    }

    out <- c(cde = cde,
             pnde = pnde, tnde = tnde,
             pnie = pnie, tnie = tnie,
             te = te,
             pm = pm)

    return(out)

  }


  bootstrap_3way <- function(data, outcome, treatment, mediator, covariates, vecc,
                             interaction, event, mreg, yreg, m_star, a_star, a, nboot) {

    bootstrap_results <- boot::boot(data = data, statistic = bootstrap_3way_step, R = nboot,
                                    outcome = outcome, treatment = treatment, mediator = mediator,
                                    covariates = covariates, vecc = vecc, interaction = interaction,
                                    event = event, mreg = mreg, yreg = yreg,
                                    m_star = m_star, a_star = a_star, a = a)

    return(bootstrap_results)

  }


  bootstrap_4way_step <- function(data, indices, outcome, treatment, mediator, covariates, vecc = NULL,
                                  interaction = TRUE, event = NULL, m_star, a_star, a,
                                  mreg = c("linear", "logistic"),
                                  yreg = c("linear", "logistic", "loglinear", "poisson",
                                           "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) {

    data_boot <- data[indices, ]

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

    thetas <- coef$thetas

    betas <- coef1$betas

    theta_interaction <- ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")], 0)

    if (yreg == "linear") {

      if (mreg == "linear") {

        cde <- unname((thetas[treatment]+theta_interaction*m_star)*(a-a_star))

        intref <- unname(theta_interaction*
                           (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))-m_star)*(a-a_star))

        intmed <- unname(theta_interaction*betas[treatment]*(a-a_star)^2)

        pie <- unname((thetas[mediator]*betas[treatment]+
                         theta_interaction*betas[treatment]*a_star)*(a-a_star))

      }

      else if (mreg=="logistic") {

        cde <- unname((thetas[treatment]+theta_interaction*m_star)*(a-a_star))

        intref <- unname(theta_interaction*(a-a_star)*
                  (exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))/
                  (1+exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))) -
                  m_star))

        intmed <- unname(theta_interaction*(a-a_star)*
                  (exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))/
                  (1+exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))) -
                  exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))/
                  (1+exp(exp(betas['(Intercept)']+betas[treatment]*a_star+
                  sum(betas[covariates]*t(vecc)))))))

        pie <- unname((thetas[mediator]+theta_interaction*a_star)*
               (exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))/
               (1+exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))) -
               exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))/
               (1+exp(exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))))))
      }

      te <- cde + intref+ intmed + pie

      cde_prop <- cde/te

      intref_prop <- intref/te

      intmed_prop <- intmed/te

      pie_prop <- pie/te

      overall_pm <- (pie+intmed)/te

      overall_int <- (intref+intmed)/te

      overall_pe <- (intref+intmed+pie)/te

      out <- c(cde = cde, intref = intref, intmed = intmed, pie = pie,
               te = te, cde_prop = cde_prop, intref_prop = intref_prop,
               intmed_prop = intmed_prop, pie_prop = pie_prop,
               overall_pm = overall_pm, overall_int = overall_int,
               overall_pe = overall_pe)

      return(out)
    }

    else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                          "negbin", "coxph", "aft_exp", "aft_weibull")) {

      if (mreg=="linear") {
        cde_comp <- unname(exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                    theta_interaction*a*m_star- (thetas[mediator]+theta_interaction*a_star)*
                    (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                    0.5*(thetas[mediator]+theta_interaction*a_star)^2*coef$variance)-
                    exp(thetas[mediator]*m_star+theta_interaction*a_star*m_star-
                    (thetas[mediator]+theta_interaction*a_star)*(betas['(Intercept)']+
                    betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                    0.5*(thetas[mediator]+theta_interaction*a_star)^2*coef$variance))

        intref_comp <- unname(exp((thetas[treatment]+theta_interaction*
                       (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                       thetas[mediator]*coef$variance))*(a-a_star)+
                       0.5*theta_interaction^2*coef$variance*(a^2-a_star^2))-1-
                       exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                       theta_interaction*a*m_star-(thetas[mediator]+theta_interaction*a_star)*
                       (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                       0.5*(thetas[mediator]+theta_interaction*a_star)^2*coef$variance)+
                       exp(thetas[mediator]*m_star+theta_interaction*a_star*m_star-(thetas[mediator]+
                       theta_interaction*a_star)*(betas['(Intercept)']+betas[treatment]*a_star+
                       sum(betas[covariates]*t(vecc)))-0.5*(thetas[mediator]+theta_interaction*a_star)^2*
                       coef$variance))

        intmed_comp <- unname(exp((thetas[treatment]+thetas[mediator]*betas[treatment]+
                       theta_interaction*(betas['(Intercept)']+betas[treatment]*a_star+
                       betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]*
                       coef$variance))*(a-a_star)+0.5*theta_interaction^2*coef$variance*(a^2-a_star^2))-
                       exp((thetas[mediator]*betas[treatment]+theta_interaction*betas[treatment]*
                       a_star)*(a-a_star))-exp((thetas[treatment]+theta_interaction*
                       (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                       thetas[mediator]*coef$variance))*(a-a_star)+0.5*theta_interaction^2*
                       coef$variance*(a^2-a_star^2))+1)

        pie_comp <- unname(exp((thetas[mediator]*betas[treatment]+theta_interaction*
                    betas[treatment]*a_star)*(a-a_star))-1)

        tcomp <- cde_comp+intref_comp+intmed_comp+pie_comp

        total_rr <- unname(exp((thetas[treatment]+theta_interaction*
                    (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                    thetas[mediator]*coef$variance))*(a-a_star)+
                    0.5*theta_interaction^2*coef$variance*(a^2-a_star^2))*
                    exp((thetas[mediator]*betas[treatment]+theta_interaction*
                    betas[treatment]*a)*(a-a_star)))
      }

      else if (mreg=="logistic") {

        cde_comp <- unname((exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                    theta_interaction*a*m_star)*(1+exp(betas['(Intercept)']+
                    betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                    betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    theta_interaction*a_star))-exp(thetas[mediator]*m_star+
                    theta_interaction*a_star*m_star)*(1+exp(betas['(Intercept)']+
                    betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                    betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    theta_interaction*a_star))))

        intref_comp <- unname(exp(thetas[treatment]*(a-a_star))*(1+exp(betas['(Intercept)']+
                       betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a))/(1+exp(betas['(Intercept)']+
                       betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a_star))-1-exp(thetas[treatment]*(a-a_star)+
                       thetas[mediator]*m_star+theta_interaction*a*m_star)*(1+
                       exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))*
                       exp((thetas[treatment]+theta_interaction*m_star)*(a-a_star))/
                       (1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                       sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a_star))+exp(thetas[mediator]*m_star+
                       theta_interaction*a_star*m_star)*(1+exp(betas['(Intercept)']+
                       betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                       betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a_star)))

        intmed_comp <- unname(exp(thetas[treatment]*(a-a_star))*(1+exp(betas['(Intercept)']+
                       betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a))*(1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                       sum(betas[covariates]*t(vecc))))/((1+exp(betas['(Intercept)']+
                       betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a_star))*(1+exp(betas['(Intercept)']+
                       betas[treatment]*a+sum(betas[covariates]*t(vecc)))))-(1+exp(betas['(Intercept)']+
                       betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                       theta_interaction*a_star))*(1+exp(betas['(Intercept)']+
                       betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/((1+
                       exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                       thetas[mediator]+theta_interaction*a_star))*(1+exp(betas['(Intercept)']+
                       betas[treatment]*a+sum(betas[covariates]*t(vecc)))))-exp(thetas[treatment]*(a-
                       a_star))*(1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                       sum(betas[covariates]*t(vecc))+thetas[mediator]+theta_interaction*a))/(1+
                      exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                      thetas[mediator]+theta_interaction*a_star))+1)

        pie_comp <- unname((1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                    sum(betas[covariates]*t(vecc))))*(1+exp(betas['(Intercept)']+betas[treatment]*a+
                    sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    theta_interaction*a_star))/((1+exp(betas['(Intercept)']+
                    betas[treatment]*a+sum(betas[covariates]*t(vecc))))*(1+exp(betas['(Intercept)']+
                    betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    theta_interaction*a_star)))-1)

        tcomp <- cde_comp+intref_comp+intmed_comp+pie_comp

        total_rr <- unname(exp(thetas[treatment]*a)*(1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                    sum(betas[covariates]*t(vecc))))*(1+exp(betas['(Intercept)']+betas[treatment]*a+
                    sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    theta_interaction*a))/(exp(thetas[treatment]*a_star)*(1+
                    exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc))))*(1+
                    exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                    thetas[mediator]+theta_interaction*a_star))))

      }

      total_err <- total_rr-1

      cde_err <- cde_comp*(total_rr-1)/tcomp

      intmed_err <- intmed_comp*(total_rr-1)/tcomp

      intref_err <- intref_comp*(total_rr-1)/tcomp

      pie_err <- pie_comp*(total_rr-1)/tcomp

      cde_prop <- cde_comp/tcomp

      intmed_prop <- intmed_comp/tcomp

      intref_prop <- intref_comp/tcomp

      pie_prop <- pie_comp/tcomp

      overall_pm <- (pie_comp+intmed_comp)/tcomp

      overall_int <- (intref_comp+intmed_comp)/tcomp

      overall_pe <- (intref_comp+intmed_comp+pie_comp)/tcomp

      out <- c(cde_comp = cde_comp, intref_comp = intref_comp,
               intmed_comp = intmed_comp, pie_comp = pie_comp,
               total_err = total_err, cde_err = cde_err, intref_err = intref_err,
               intmed_err = intmed_err, pie_err = pie_err,
               cde_prop = cde_prop, intref_prop = intref_prop,
               intmed_prop = intmed_prop, pie_prop = pie_prop,
               overall_pm = overall_pm, overall_int = overall_int,
               overall_pe = overall_pe)

      return(out)
    }
  }


  bootstrap_4way <- function(data, nboot, outcome, treatment, mediator, covariates, vecc,
                             interaction, event, m_star, a_star, a, mreg, yreg) {

    bootstrap_results <- boot::boot(data = data, statistic = bootstrap_4way_step, R = nboot,
                                    outcome = outcome, treatment = treatment, mediator = mediator,
                                    covariates = covariates, vecc = vecc, interaction = interaction,
                                    event = event, mreg = mreg, yreg = yreg,
                                    m_star = m_star, a_star = a_star, a = a)

    return(bootstrap_results)

  }


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

  n=nrow(data)

  if (type == "delta") {

    if (nway == 3)
      effect_estimate <- format_df_delta(delta(coef = coef, vecc = vecc, m_star = m_star,
                                               a_star = a_star, a = a),
                                       conf = conf, n = n, nway = 3, yreg  = yreg)
    else if (nway == 4)
      effect_estimate <- format_df_delta(delta_4way(coef = coef, vecc = vecc, m_star = m_star,
                                                    a_star = a_star, a = a),
                                         conf = conf, n = n, nway = 4, yreg  = yreg)

    class(effect_estimate)=c("delta_out","data.frame")

  }else if (type == "bootstrap") {

    if (nway == 3)
      effect_estimate <- format_df_boot(bootstrap(data = data, outcome = outcome, treatment = treatment,
                                                  mediator = mediator, covariates = covariates,
                                                  vecc = vecc, interaction = interaction, event = event,
                                                  mreg = mreg, yreg = yreg, n = n,
                                                  nboot = nboot, conf = conf,
                                                  m_star = m_star, a_star = a_star, a = a),
                                        conf = conf, n = n, nway = 3, yreg  = yreg)
    else if (nway == 4)
      effect_estimate <- format_df_boot(bootstrap_4way(data = data, nboot = nboot, outcome = outcome,
                                                       treatment = treatment, mediator = mediator,
                                                       covariates = covariates, vecc = vecc,
                                                       interaction = interaction, event = event,
                                                       mreg = mreg, yreg = yreg,
                                                       m_star = m_star, a_star = a_star, a = a),
                                        conf = conf, n = n, nway = 4, yreg  = yreg)

    class(effect_estimate)=c("boot_out","data.frame")

    }

  return(effect_estimate)

}
