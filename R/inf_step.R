inf_step <- function(data, nboot, outcome, exposure, exposure.type, mediator,
                     covariates.pre, covariates.post, covariates.post.type, vecc,
                     EMint, MMint, EMMint, EMint.terms, MMint.terms, EMMint.terms,
                     event, mreg, yreg, m_star, a_star, a, est.method, inf.method, model) {

  if (inf.method == "delta") {

    if (length(mediator) > 1) {
      stop("Delta method inference doesn't support multiple mediator cases")
    }

    if (model != "rb") {
      stop("Delta method inference doesn't support selected model")
    } else if (length(mediator) > 1) {
      stop("Delta method inference doesn't support multiple mediator cases")
    }

    if (!(exposure.type %in% c("continuous", "binary"))) {
      stop("For the selected model, Delta method inference only supports continuous or binary exposure")
    }


    formulas <- create_formulas(model = model, outcome = outcome, exposure = exposure,
                                mediator = mediator, covariates.pre = covariates.pre,
                                covariates.post = covariates.post,
                                EMint = EMint, MMint = MMint, EMMint = EMMint,
                                EMint.terms = EMint.terms, MMint.terms = MMint.terms,
                                EMMint.terms = EMMint.terms,
                                event = event, mreg = mreg, yreg = yreg, data = data)

    regressions <- run_regressions(model = model, formulas = formulas,
                                   exposure = exposure, exposure.type = exposure.type,
                                   mediator = mediator, covariates.post = covariates.post,
                                   covariates.post.type = covariates.post.type,
                                   mreg = mreg, yreg = yreg, data = data)

    coef <- get_coef(formulas = formulas, regressions = regressions,
                     mreg = mreg, yreg = yreg, model = model)

    mlevel <- ifelse(is.factor(data[, mediator]), length(levels(data[, mediator])), 2)

    thetas <- coef$thetas

    betas <- coef$betas

    variance <- coef$variance

    vcov_block <- coef$vcov_block

    theta0 <- "x1"

    theta1 <- "x2"

    theta2 <- paste0("x", 3:(mlevel + 1))

    if (EMint == TRUE) {
      theta3 <- paste0("x", length(thetas) - ((mlevel - 2):0))
    } else {theta3 <- rep(0, mlevel - 1)}

    beta0 <- paste0("x", length(thetas) + 1 + (0:(mlevel - 2))*length(betas)/(mlevel - 1))

    beta1 <- paste0("x", length(thetas) + 2 + (0:(mlevel - 2))*length(betas)/(mlevel - 1))

    XC <- sapply(0:(mlevel-2), function(x) paste0("x", length(thetas) + 2 +
                                                    x*length(betas)/(mlevel-1) +
                                                    1:length(vecc),
                                                  "*", as.character(vecc),
                                                  collapse = "+"))


    if (yreg == "linear") {

      if (mreg == "linear") {

        cde_formula <- paste0("(", theta1, "+", theta3, "*", "m_star)*(a-a_star)")

        pnde_formula <- paste0("(", theta1, "+", theta3, "*(", beta0, "+", beta1,
                               "*a_star+", XC, "))*(a-a_star)")

        tnde_formula <- paste0("(", theta1, "+", theta3, "*(", beta0, "+", beta1,
                               "*a+", XC, "))*(a-a_star)")

        pnie_formula <- paste0(beta1, "*(", theta2, "+", theta3, "*a_star)*(a-a_star)")

        tnie_formula <- paste0(beta1, "*(", theta2, "+", theta3, "*a)*(a-a_star)")

      } else if (mreg %in% c("logistic", "multinomial")) {

        cde_formula <- paste0("(", theta1, "+", ifelse(m_star == 0, 0, theta3[m_star]), ")*(a-a_star)")

        pnde_formula <- paste0("(", theta1, "+(",
                               paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")*",
                                      theta3, collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")",
                                      collapse = "+"), "))*(a-a_star)")

        tnde_formula <- paste0("(", theta1, "+(",
                               paste0("exp(", beta0, "+", beta1, "*a+", XC, ")*",
                                      theta3, collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+", beta1, "*a+", XC, ")",
                                      collapse = "+"), "))*(a-a_star)")

        pnie_formula <- paste0("(", paste0("exp(", beta0, "+", beta1, "*a+", XC, ")*(", theta2, "+",
                                           theta3, "*a_star)", collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+", beta1, "*a+", XC, ")", collapse = "+"), ")-(",
                               paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")*(", theta2, "+",
                                      theta3, "*a_star)", collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")", collapse = "+"), ")")

        tnie_formula <-  paste0("(", paste0("exp(", beta0, "+", beta1, "*a+", XC, ")*(", theta2, "+",
                                            theta3, "*a)", collapse = "+"), ")/(1+",
                                paste0("exp(", beta0, "+", beta1, "*a+", XC, ")", collapse = "+"), ")-(",
                                paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")*(", theta2, "+",
                                       theta3, "*a)", collapse = "+"), ")/(1+",
                                paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")", collapse = "+"), ")")

      }

      te_formula <- paste0("(", tnie_formula, ")+(", pnde_formula, ")")

      pm_formula <- paste0("(", tnie_formula, ")/((", pnde_formula, ")+(", te_formula, "))")

      intref_formula <- paste0("(", pnde_formula, ")-(", cde_formula, ")")

      intmed_formula <- paste0("(", tnie_formula, ")-(", pnie_formula, ")")

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

        pnie_rr_formula <- paste0("exp((", theta2, "*", beta1, "+", theta3, "*",
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

      } else if (mreg %in% c("logistic", "multinomial")) {

        cde_rr_formula <- paste0("exp((", theta1, "+", ifelse(m_star == 0, 0, theta3[m_star]), ")*(a-a_star))")

        paste0("exp(", beta0, "+", beta1, "*a+", XC, ")", collapse = "+")

        pnde_rr_formula <- paste0("exp(", theta1, "*(a-a_star))*(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a+",beta0,
                                         "+", beta1, "*a_star+", XC, ")", collapse = "+"), ")/(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a_star+",beta0,
                                         "+", beta1, "*a_star+", XC, ")", collapse = "+"), ")")

        tnde_rr_formula <- paste0("exp(", theta1, "*(a-a_star))*(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a+",beta0,
                                         "+", beta1, "*a+", XC, ")", collapse = "+"), ")/(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a_star+",beta0,
                                         "+", beta1, "*a+", XC, ")", collapse = "+"), ")")

        pnie_rr_formula <- paste0("(1+", paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")",
                                                collapse = "+"),
                                  ")*(1+", paste0("exp(", theta2, "+", theta3, "*a_star+",beta0,
                                                  "+", beta1, "*a+", XC, ")", collapse = "+"), ")/((1+",
                                  paste0("exp(", beta0, "+", beta1, "*a+", XC, ")",
                                         collapse = "+"), ")*(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a_star+",beta0,
                                         "+", beta1, "*a_star+", XC, ")", collapse = "+"), "))")

        tnie_rr_formula <- paste0("(1+", paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")",
                                                collapse = "+"),
                                  ")*(1+", paste0("exp(", theta2, "+", theta3, "*a+",beta0,
                                                  "+", beta1, "*a+", XC, ")", collapse = "+"), ")/((1+",
                                  paste0("exp(", beta0, "+", beta1, "*a+", XC, ")",
                                         collapse = "+"), ")*(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a+",beta0,
                                         "+", beta1, "*a_star+", XC, ")", collapse = "+"), "))")

        cde_err_formula <- paste0("(exp(", ifelse(m_star == 0, 0, theta2[m_star]), ")*(exp(", theta1,
                                  "*(a-a_star)+a*", ifelse(m_star == 0, 0, theta3[m_star]),
                                  ")-exp(a_star*", ifelse(m_star == 0, 0, theta3[m_star]), "))*(1+",
                                  paste0("exp(", beta0, "+", beta1, "*a_star+", XC, ")",
                                         collapse = "+"), "))/(1+",
                                  paste0("exp(", theta2, "+", theta3, "*a_star+",beta0,
                                         "+", beta1, "*a_star+", XC, ")", collapse = "+"), ")")

      }

      te_rr_formula <- paste0("(", tnie_rr_formula, ")*(", pnde_rr_formula, ")")

      pm_formula <- paste0("((", pnde_rr_formula, ")*(", tnie_rr_formula, "-1))/(", te_rr_formula, "-1)")

      intref_err_formula <- paste0("(", pnde_rr_formula, ")-1-(", cde_err_formula," )")

      intmed_err_formula <- paste0("(", tnie_rr_formula, ")*(", pnde_rr_formula, ")-(",
                                   pnde_rr_formula, ")-(", pnie_rr_formula, ")+1")

      pie_err_formula <- paste0("(", pnie_rr_formula, ")-1")

      te_err_formula <- paste0("(", te_rr_formula, ")-1")

      cde_err_prop_formula <- paste0("(", cde_err_formula, ")/(", te_err_formula, ")")

      intmed_err_prop_formula <- paste0("(", intmed_err_formula, ")/(", te_err_formula, ")")

      intref_err_prop_formula <-paste0("(", intref_err_formula, ")/(", te_err_formula, ")")

      pie_err_prop_formula <- paste0("(", pie_err_formula, ")/(", te_err_formula, ")")

      overall_pm_formula <- paste0("((", pie_err_formula, ")+(", intmed_err_formula, "))/(",
                                   te_err_formula, ")")

      overall_int_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, "))/(",
                                    te_err_formula, ")")

      overall_pe_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, ")+(",
                                   pie_err_formula,
                                   "))/(", te_err_formula, ")")

      delta_formula <- list(cde_rr_formula = cde_rr_formula, pnde_rr_formula = pnde_rr_formula,
                            tnde_rr_formula = tnde_rr_formula, pnie_rr_formula = pnie_rr_formula,
                            tnie_rr_formula = tnie_rr_formula, te_rr_formula = te_rr_formula,
                            pm_formula = pm_formula, cde_err_formula = cde_err_formula,
                            intref_err_formula = intref_err_formula, intmed_err_formula = intmed_err_formula,
                            pie_err_formula = pie_err_formula, te_err_formula = te_err_formula,
                            cde_err_prop_formula = cde_err_prop_formula,
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

  } else if (inf.method == "bootstrap") {

    boots <- boot::boot(data = data, statistic = est_step, R = nboot,
                        outcome = outcome, exposure = exposure, exposure.type = exposure.type,
                        mediator = mediator, model = model,
                        covariates.pre = covariates.pre, covariates.post = covariates.post,
                        covariates.post.type = covariates.post.type, vecc = vecc,
                        EMint = EMint, MMint = MMint, EMMint = EMMint, EMint.terms = EMint.terms,
                        MMint.terms = MMint.terms, EMMint.terms = EMMint.terms,
                        event = event, mreg = mreg, yreg = yreg,
                        m_star = m_star, a_star = a_star, a = a, est.method = est.method)

    effect_se <- apply(boots$t, 2, sd)

  }

  return(effect_se)

}
