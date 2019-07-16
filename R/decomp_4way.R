decomp_4way_step <- function(data, outcome, treatment, mediator, covariates, vecc = NULL,
                        interaction = TRUE, event = NULL, m = 0, a_star = 1, a = 0,
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

  thetas <- coef$thetas

  betas <- coef1$betas

  if (yreg == "linear") {

    if (mreg == "linear") {

      cde <- unname((thetas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*m)*(a-a_star))

      intref <- unname(thetas[paste(treatment, mediator, sep = ":")]*
        (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))-m)*(a-a_star))

      intmed <- unname(thetas[paste(treatment, mediator, sep = ":")]*betas[treatment]*(a-a_star)^2)

      pie <- unname((thetas[mediator]*betas[treatment]+
              thetas[paste(treatment, mediator, sep = ":")]*betas[treatment]*a_star)*(a-a_star))

    }

    else if (mreg=="logistic") {

      cde <- unname((thetas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*m)*(a-a_star))

      intref <- unname(thetas[paste(treatment, mediator, sep = ":")]*(a-a_star)*
        (exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))/
           (1+exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))) - m))

      intmed <- unname(thetas[paste(treatment, mediator, sep = ":")]*(a-a_star)*
                   (exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))/
                    (1+exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc)))) -
                      exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))/
                      (1+exp(exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))))))

      pie <- unname((thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)*
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

    class(out) <- "decomp_4way"

    return(out)
  }

  else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                       "negbin", "coxph", "aft_exp", "aft_weibull")) {

    if (mreg=="linear") {
      cde_comp <- unname(exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m+
                   thetas[paste(treatment, mediator, sep = ":")]*a*m-
                  (thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)*
                  (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                   0.5*(thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)^2*
                   coef$variance)-exp(thetas[mediator]*m+thetas[paste(treatment, mediator, sep = ":")]*
                   a_star*m-(thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)*
                   (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                   0.5*(thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)^2*
                   coef$variance))

      intref_comp <- unname(exp((thetas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*
                     (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                      thetas[mediator]*coef$variance))*(a-a_star)+
                      0.5*thetas[paste(treatment, mediator, sep = ":")]^2*coef$variance*(a^2-a_star^2))-1-
                      exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m+
                      thetas[paste(treatment, mediator, sep = ":")]*a*m-
                      (thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)*
                      (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                      0.5*(thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)^2*
                      coef$variance)+exp(thetas[mediator]*m+thetas[paste(treatment, mediator, sep = ":")]*
                      a_star*m-(thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)*
                      (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                      0.5*(thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star)^2*
                      coef$variance))

      intmed_comp <- unname(exp((thetas[treatment]+thetas[mediator]*betas[treatment]+
                     thetas[paste(treatment, mediator, sep = ":")]*(betas['(Intercept)']+
                     betas[treatment]*a_star+betas[treatment]*a+sum(betas[covariates]*t(vecc))+
                     thetas[mediator]*coef$variance))*(a-a_star)+
                     0.5*thetas[paste(treatment, mediator, sep = ":")]^2*coef$variance*(a^2-a_star^2))-
                     exp((thetas[mediator]*betas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*
                     betas[treatment]*a_star)*(a-a_star))-exp((thetas[treatment]+
                     thetas[paste(treatment, mediator, sep = ":")]*(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]*coef$variance))*
                     (a-a_star)+0.5*thetas[paste(treatment, mediator, sep = ":")]^2*
                     coef$variance*(a^2-a_star^2))+1)

      pie_comp <- unname(exp((thetas[mediator]*betas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*
                          betas[treatment]*a_star)*(a-a_star))-1)

      tcomp <- cde_comp+intref_comp+intmed_comp+pie_comp

      total_rr <- unname(exp((thetas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*
               (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                thetas[mediator]*coef$variance))*(a-a_star)+
                0.5*thetas[paste(treatment, mediator, sep = ":")]^2*coef$variance*(a^2-a_star^2))*
                exp((thetas[mediator]*betas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*
                betas[treatment]*a)*(a-a_star)))
     }

    else if (mreg=="logistic") {

      cde_comp <- unname((exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m+
                  thetas[paste(treatment, mediator, sep = ":")]*a*m)*(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                  thetas[paste(treatment, mediator, sep = ":")]*a_star))-exp(thetas[mediator]*m+
                  thetas[paste(treatment, mediator, sep = ":")]*a_star*m)*(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                  thetas[paste(treatment, mediator, sep = ":")]*a_star))))

      intref_comp <- unname(exp(thetas[treatment]*(a-a_star))*(1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a))/(1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star))-1-exp(thetas[treatment]*(a-a_star)+
                     thetas[mediator]*m+thetas[paste(treatment, mediator, sep = ":")]*a*m)*(1+
                     exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))*
                     exp((thetas[treatment]+thetas[paste(treatment, mediator, sep = ":")]*m)*(a-a_star))/
                     (1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                     sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star))+exp(thetas[mediator]*m+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star*m)*(1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star)))

      intmed_comp <- unname(exp(thetas[treatment]*(a-a_star))*(1+exp(betas['(Intercept)']+betas[treatment]*a+
                     sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a))*(1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/((1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star))*(1+exp(betas['(Intercept)']+
                     betas[treatment]*a+sum(betas[covariates]*t(vecc)))))-(1+exp(betas['(Intercept)']+
                     betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star))*(1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/((1+exp(betas['(Intercept)']+
                     betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                     thetas[paste(treatment, mediator, sep = ":")]*a_star))*(1+exp(betas['(Intercept)']+
                     betas[treatment]*a+sum(betas[covariates]*t(vecc)))))-exp(thetas[treatment]*(a-a_star))*(1+
                     exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                     thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a))/(1+
                     exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                     thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star))+1)

      pie_comp <- unname((1+exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))*(1+
                   exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                   thetas[paste(treatment, mediator, sep = ":")]*a_star))/((1+exp(betas['(Intercept)']+
                   betas[treatment]*a+sum(betas[covariates]*t(vecc))))*(1+exp(betas['(Intercept)']+
                   betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                   thetas[paste(treatment, mediator, sep = ":")]*a_star)))-1)

      tcomp <- cde_comp+intref_comp+intmed_comp+pie_comp

      total_rr <- unname(exp(thetas[treatment]*a)*(1+exp(betas['(Intercept)']+betas[treatment]*a_star+
                  sum(betas[covariates]*t(vecc))))*(1+exp(betas['(Intercept)']+betas[treatment]*a+
                  sum(betas[covariates]*t(vecc))+thetas[mediator]+
                  thetas[paste(treatment, mediator, sep = ":")]*a))/(exp(thetas[treatment]*a_star)*(1+
                  exp(betas['(Intercept)']+betas[treatment]*a+sum(betas[covariates]*t(vecc))))*(1+
                  exp(betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+
                  thetas[mediator]+thetas[paste(treatment, mediator, sep = ":")]*a_star))))
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

      class(out) <- "decomp_4way"

      return(out)
  }
}

decomp_4way <- function(data, outcome, treatment, mediator, covariates, vecc = NULL,
                        interaction = TRUE, type = c("delta", "bootstrap"), nboot = 100,
                        conf = 0.95, mreg = c("linear", "logistic"),
                        yreg = c("linear", "logistic", "loglinear", "poisson",
                                 "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull"),
                        event = NULL, m = 0, a_star = 1, a = 0) {

  if (type == "delta") {
    warning("Delta method not finished")
  } else if (type == "bootstrap") {
    bootstrap_decom_4way(data = data, outcome = outcome, treatment = treatment, mediator = mediator,
                         covariates = covariates, vecc = vecc, interaction = interaction, event = event,
                         mreg = mreg, yreg = yreg, nboot = nboot, conf = conf,
                         m = m, a_star = a_star, a = a)
  }
}

