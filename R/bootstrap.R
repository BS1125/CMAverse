bootstrap_step <- function(data, indices, outcome, treatment, mediator, covariates, vecc,
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

  cde <- cde_function(thetas = coef$thetas, treatment = treatment,
                      mediator = mediator, interaction = interaction,
                      m_star = m_star, a_star = a_star, a = a, yreg = yreg)

  nde <- nde_function(coef$betas, coef$thetas, coef$variance, vcov_block, treatment, mediator,
                      covariates, vecc, interaction, a_star, a, mreg, yreg)

  nie <- nie_function(coef$betas, coef$thetas, treatment, mediator, covariates, vecc = vecc,
                      interaction, a_star, a, mreg, yreg)

  te <- te_function(coef$betas, coef$thetas, coef$variance, coef$vcov_block, treatment, mediator,
                    covariates, vecc, interaction, a_star, a, mreg, yreg)

  pm <- pm_function(coef$betas, coef$thetas, coef$variance, coef$vcov_block, treatment, mediator,
                    covariates, vecc, interaction, a_star, a, mreg, yreg)

  out <- c(cde = cde,
           pnde = nde$pnde, tnde = nde$tnde,
           pnie = nie$pnie, tnie = nie$tnie,
           te = te,
           pm = pm)

  return(out)

}


bootstrap <- function(data, indices, outcome, treatment, mediator, covariates, vecc,
                      interaction, event, mreg, yreg, m_star, a_star, a) {

  bootstrap_results <- boot::boot(data = data, statistic = bootstrap_step, R = nboot,
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
                            (1+exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))) - m_star))

      intmed <- unname(theta_interaction*(a-a_star)*
                         (exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))/
                            (1+exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))) -
                            exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))/
                            (1+exp(exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))))))

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
             cde_prop = cde_prop, intref_prop = intref_prop,
             intmed_prop = intmed_prop, pie_prop = pie_prop,
             overall_pm = overall_pm, overall_int = overall_int,
             overall_pe = overall_pe)

    return(out)
  }

  else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                        "negbin", "coxph", "aft_exp", "aft_weibull")) {

    if (mreg=="linear") {
      cde_comp <- unname(exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                               theta_interaction*a*m_star-
                               (thetas[mediator]+theta_interaction*a_star)*
                               (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                               0.5*(thetas[mediator]+theta_interaction*a_star)^2*
                               coef$variance)-exp(thetas[mediator]*m_star+theta_interaction*
                                                    a_star*m_star-(thetas[mediator]+theta_interaction*a_star)*
                                                    (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                                                    0.5*(thetas[mediator]+theta_interaction*a_star)^2*
                                                    coef$variance))

      intref_comp <- unname(exp((thetas[treatment]+theta_interaction*
                                   (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                                      thetas[mediator]*coef$variance))*(a-a_star)+
                                  0.5*theta_interaction^2*coef$variance*(a^2-a_star^2))-1-
                              exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                                    theta_interaction*a*m_star-
                                    (thetas[mediator]+theta_interaction*a_star)*
                                    (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                                    0.5*(thetas[mediator]+theta_interaction*a_star)^2*
                                    coef$variance)+exp(thetas[mediator]*m_star+theta_interaction*
                                                         a_star*m_star-(thetas[mediator]+theta_interaction*a_star)*
                                                         (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                                                         0.5*(thetas[mediator]+theta_interaction*a_star)^2*
                                                         coef$variance))

      intmed_comp <- unname(exp((thetas[treatment]+thetas[mediator]*betas[treatment]+
                                   theta_interaction*(betas['(Intercept)']+
                                                        betas[treatment]*a_star+betas[treatment]*a+sum(betas[covariates]*t(vecc))+
                                                        thetas[mediator]*coef$variance))*(a-a_star)+
                                  0.5*theta_interaction^2*coef$variance*(a^2-a_star^2))-
                              exp((thetas[mediator]*betas[treatment]+theta_interaction*
                                     betas[treatment]*a_star)*(a-a_star))-exp((thetas[treatment]+
                                                                                 theta_interaction*(betas['(Intercept)']+
                                                                                                      betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]*coef$variance))*
                                                                                (a-a_star)+0.5*theta_interaction^2*
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

      intmed_comp <- unname(exp(thetas[treatment]*(a-a_star))*(1+exp(betas['(Intercept)']+betas[treatment]*a+
                                                                       sum(betas[covariates]*t(vecc))+thetas[mediator]+
                                                                       theta_interaction*a))*(1+exp(betas['(Intercept)']+
                                                                                                      betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/((1+exp(betas['(Intercept)']+
                                                                                                                                                                         betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                                                                                                                                                                         theta_interaction*a_star))*(1+exp(betas['(Intercept)']+
                                                                                                                                                                                                             betas[treatment]*a+sum(betas[covariates]*t(vecc)))))-(1+exp(betas['(Intercept)']+
                                                                                                                                                                                                                                                                           betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                                                                                                                                                                                                                                                                           theta_interaction*a_star))*(1+exp(betas['(Intercept)']+
                                                                                                                                                                                                                                                                                                               betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/((1+exp(betas['(Intercept)']+
                                                                                                                                                                                                                                                                                                                                                                                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                                                                                                                                                                                                                                                                                                                                                                                  theta_interaction*a_star))*(1+exp(betas['(Intercept)']+
                                                                                                                                                                                                                                                                                                                                                                                                                      betas[treatment]*a+sum(betas[covariates]*t(vecc)))))-exp(thetas[treatment]*(a-a_star))*(1+
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      thetas[mediator]+theta_interaction*a))/(1+
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      thetas[mediator]+theta_interaction*a_star))+1)

      pie_comp <- unname((1+exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))*(1+
                                                                                                                 exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]+
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

    unname(cde_comp)

    out <- c(cde_comp = cde_comp, intref_comp = intref_comp,
             intmed_comp = intmed_comp, pie_comp = pie_comp,
             total_rr = total_rr, total_err = total_err,
             cde_err = cde_err, intref_err = intref_err,
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
                                  outcome = outcome, treatment = treatment, mediator = mediator, covariates = covariates,
                                  vecc = c(1,1,1), interaction = interaction, event = event,
                                  m_star = m_star, a_star = a_star, a = a, mreg = mreg, yreg = yreg)

  return(bootstrap_results)

}
