regrun <- function() {
  # Y: outcome; M: mediator; A: exposure; C: basec; L: postc
  # the variable used to calculate weights is required to be categorical
  
  # Exposure Regression For Weighting
  # for wb and msm, the exposure regression is required for calculating weights if basec is not empty, w_{a,i}=P(A=A_i)/P(A=A_i|C_i)
  # for iorw, the exposure regression is required for calculating weights, w_{a,i}=P(A=0|M_i,C_i)/P(A=A_i|M_i,C_i)
  if ((model %in% c("wb", "msm") && length(basec) > 0) | model == "iorw") {
    if (is.null(ereg)) stop("ereg is required when model is 'wb' or 'msm' with length(basec) > 0 and when model is 'iorw'")
    if (is.character(ereg)) {
      # fit glm with family = poisson() rather than family = binomial("log") for "loglinear"
      if (ereg == "loglinear" && length(unique(data[, exposure])) != 2) stop("When ereg is 'loglinear', exposure should be binary")
      if (!ereg %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop("Select character ereg from 'logistic', 'loglinear', 'multinomial', 'ordinal'")
      exposure_formula <- switch((model == "iorw") + 1, "1" = paste0(exposure, "~", paste0(basec, collapse = "+")),
                                 "2" = paste0(exposure, "~", paste0(c(mediator, basec), collapse = "+")))
      switch(ereg,
             logistic = ereg <- eval(bquote(glm(.(as.formula(exposure_formula)), family = binomial(), data = .(data)))),
             loglinear = ereg <- eval(bquote(glm(.(as.formula(exposure_formula)), family = poisson(), data = .(data)))),
             multinomial = ereg <- eval(bquote(nnet::multinom(.(as.formula(exposure_formula)), data = .(data), trace = FALSE))),
             ordinal = ereg <- eval(bquote(MASS::polr(.(as.formula(exposure_formula)), data = .(data)))))
    } else if (!(identical(class(ereg), "lm") | identical(class(ereg), c("glm", "lm")) |
                  identical(class(ereg), c("negbin", "glm", "lm")) | identical(class(ereg), c("multinom", "nnet")) |
                  identical(class(ereg), c("gam", "glm", "lm")) | identical(class(ereg), "polr"))) {
      stop("Fit ereg by lm, glm, glm.nb, gam, multinom, polr")
    }
      
  } else {
    if (!is.null(ereg)) warning("ereg is ignored when model is 'wb' or 'msm' with length(basec) = 0 or model is 'rb', 'ne' or 'gformula'")
    ereg <- NULL
  }
  
  # Mediator Regression For Weighting
  # for msm, a mediator regression for weighting is required for each mediator
  if (model == "msm") {
    if (!is.list(wmreg)) stop("wmreg should be a list")
    if (length(wmreg) != length(mediator)) stop("length(wmreg) != length(mediator)")
    for (p in 1:length(wmreg)) {
      if (is.null(wmreg[[p]])) stop(paste0("Unspecified wmreg[[", p, "]]"))
      if (is.character(wmreg[[p]])) {
        if (wmreg[[p]] == "loglinear" && length(unique(data[, mediator[p]])) != 2) stop(paste0("When wmreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!wmreg[[p]] %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop(paste0("Select character wmreg[[", p, "]] from 'logistic', 'loglinear', 'multinomial', 'ordinal'"))
        # w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
        wmreg_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec, postc), collapse = "+"))
        # regression for the denominator of w_{m_p,i}
        switch(wmreg[[p]],
               logistic = wmreg[[p]] <- eval(bquote(glm(.(as.formula(wmreg_formula)), family = binomial(), data = .(data)))),
               loglinear = wmreg[[p]] <- eval(bquote(glm(.(as.formula(wmreg_formula)), family = poisson(), data = .(data)))),
               multinomial = wmreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(wmreg_formula)), data = .(data), trace = FALSE))),
               ordinal = wmreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(wmreg_formula)), data = .(data)))))
      } else if (!(identical(class(wmreg[[p]]), "lm") | identical(class(wmreg[[p]]), c("glm", "lm")) |
                   identical(class(wmreg[[p]]), c("negbin", "glm", "lm")) | identical(class(wmreg[[p]]), c("multinom", "nnet")) |
                   identical(class(wmreg[[p]]), c("gam", "glm", "lm")) | identical(class(wmreg[[p]]), "polr"))) {
        stop(paste0("Fit wmreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      }
    }
  } else {
    if (!is.null(wmreg)) warning("wmreg is ignored when model is not 'msm'")
    wmreg <- NULL
  }
  
  # Mediator Regression
  # for rb, msm and gformula, a mediator regression is required for each mediator
  if (model %in% c("rb", "msm", "gformula")) {
    if (!is.list(mreg)) stop("mreg should be a list")
    if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
    for (p in 1:length(mreg)) {
      if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
      if (is.character(mreg[[p]])) {
        if (mreg[[p]] == "loglinear" && length(unique(data[, mediator[p]])) != 2) stop(paste0("When mreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!mreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                              "negbin", "multinomial", "ordinal")) stop(
                                paste0("Select character mreg[[", p, "]] from 'linear', 'logistic',
                                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
        # for rb, regress each mediator on A, C and previous mediators
        # for msm, regress each mediator on A and previous mediators
        # for gformula, regress each mediator on A, C, L and previous mediators
        switch(model,
               rb = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec), collapse = "+")),
               msm = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)]), collapse = "+")),
               gformula = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec, postc), collapse = "+")))
        switch(mreg[[p]],
               linear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = gaussian(), data = .(data)))),
               logistic = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = binomial(), data = .(data)))),
               loglinear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
               poisson = mreg[[p]]  <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
               quasipoisson = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = quasipoisson(), data = .(data)))),
               negbin = mreg[[p]] <- eval(bquote(MASS::glm.nb(.(as.formula(mediator_formula)), data = .(data)))),
               multinomial = mreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(mediator_formula)), data = .(data), trace = FALSE))),
               ordinal = mreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(mediator_formula)), data = .(data)))))
      } else if (!(identical(class(mreg[[p]]), "lm") | identical(class(mreg[[p]]), c("glm", "lm")) |
                   identical(class(mreg[[p]]), c("negbin", "glm", "lm")) | identical(class(mreg[[p]]), c("multinom", "nnet")) |
                   identical(class(mreg[[p]]), c("gam", "glm", "lm")) | identical(class(mreg[[p]]), "polr"))) {
        stop(paste0("Fit mreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      }
    }
  } else {
    if (!is.null(mreg)) warning("mreg is ignored when model is 'wb', 'iorw' or 'ne'")
    mreg <- NULL
  }
  
  # postc Regression
  # for gformula, a regression is required for each L
  if (model == "gformula" && length(postc) > 0) {
    if (!is.list(postcreg)) stop("postcreg should be a list")
    if (length(postcreg) != length(postc)) stop("length(postcreg) != length(postc)")
    for (p in 1:length(postcreg)) {
      if (is.null(postcreg[[p]])) stop(paste0("Unspecified postcreg[[", p, "]]"))
      if (is.character(postcreg[[p]])) {
        if (postcreg[[p]] == "loglinear" && length(unique(data[, postc[p]])) != 2) stop(paste0("When postcreg[[", p, "]] is 'loglinear', postc[[", p, "]] should be binary"))
        if (!postcreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                                  "negbin", "multinomial", "ordinal"))  stop(
                                    paste0("Select character postcreg[[", p, "]] from 'linear', 'logistic',
                                           'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
        # regress each L on A, C and previous Ls
        postc_formula <- paste0(postc[p], "~", paste(c(exposure, basec, postc[0:(p-1)]), collapse = "+"))
        switch(postcreg[[p]],
               linear = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = gaussian(), data = .(data)))),
               logistic = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = binomial(), data = .(data)))),
               loglinear = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = poisson(), data = .(data)))),
               poisson = postcreg[[p]]  <- eval(bquote(glm(.(as.formula(postc_formula)), family = poisson(), data = .(data)))),
               quasipoisson = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = quasipoisson(), data = .(data)))),
               negbin = postcreg[[p]] <- eval(bquote(MASS::glm.nb(.(as.formula(postc_formula)), data = .(data)))),
               multinomial = postcreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(postc_formula)), data = .(data), trace = FALSE))),
               ordinal = postcreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(postc_formula)), data = .(data)))))
      } else if (!(identical(class(postcreg[[p]]), "lm") | identical(class(postcreg[[p]]), c("glm", "lm")) |
                   identical(class(postcreg[[p]]), c("negbin", "glm", "lm")) | identical(class(postcreg[[p]]), c("multinom", "nnet")) |
                   identical(class(postcreg[[p]]), c("gam", "glm", "lm")) | identical(class(postcreg[[p]]), "polr"))) {
        stop(paste0("Fit postcreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
      }
    }
  } else {
    if (!is.null(postcreg)) warning("postcreg is ignored when model is not 'gformula' and when length(postc) = 0")
    postcreg <- NULL
  }
  
  # Outcome Regression
  if (is.null(yreg)) stop("yreg is required")
  if (is.character(yreg)) {
    if (yreg == "loglinear" && length(unique(data[, outcome])) != 2) stop("When yreg is 'loglinear', outcome should be binary")
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                     "negbin", "multinomial", "ordinal", "coxph", "aft_exp",
                     "aft_weibull")) stop(
                       paste0("Select character yreg from 'linear', 'logistic',
                              'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal',
                              'coxph', 'aft_exp', 'aft_weibull'"))
    if (model != "iorw") {
      out$variables$EMint <- EMint
      int.terms <- switch(EMint + 1, "1" = NULL, "2" = paste(exposure, mediator, sep = "*"))
    }
    # for rb, wb and ne, regress Y on A, M and C
    # for iorw, regress Y on A and C
    # for msm, regress Y on A and M
    # for gformula, regress Y on A, M, C and L
    switch(model,
           rb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           wb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           ne = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           iorw = outcome_formula <- paste0(outcome, "~", paste(c(exposure, basec), collapse = "+")),
           msm = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms), collapse = "+")),
           gformula = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec, postc), collapse = "+")))
    if (yreg %in% c("coxph","aft_exp","aft_weibull")) {
      if (!is.null(event)) {
        outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                                 strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
      } else outcome_formula <- paste(paste0("Surv(", outcome, ")"),
                                      strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
    }
    switch(yreg,
           linear = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                            family = gaussian(), data = .(data)))),
           logistic = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                              family = binomial(), data = .(data)))),
           loglinear = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                               family = poisson(), data = .(data)))),
           poisson = yreg  <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                              family = poisson(), data = .(data)))),
           quasipoisson = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                                  family = quasipoisson(), data = .(data)))),
           negbin = yreg <- eval(bquote(MASS::glm.nb(formula = .(as.formula(outcome_formula)),
                                                     data = .(data)))),
           multinomial = yreg <- eval(bquote(nnet::multinom(formula = .(as.formula(outcome_formula)),
                                                            data = .(data), trace = FALSE))),
           ordinal = yreg <- eval(bquote(MASS::polr(formula = .(as.formula(outcome_formula)),
                                                    data = .(data)))),
           coxph = yreg <- eval(bquote(survival::coxph(formula = .(as.formula(outcome_formula)),
                                                       data = .(data)))),
           aft_exp = yreg <- eval(bquote(survival::survreg(formula = .(as.formula(outcome_formula)),
                                                           dist = "exponential", data = .(data)))),
           aft_weibull = yreg <- eval(bquote(survival::survreg(formula = .(as.formula(outcome_formula)),
                                                               dist = "weibull", data = .(data)))))
  } else if (!(identical(class(yreg), "lm") | identical(class(yreg), c("glm", "lm")) |
               identical(class(yreg), c("negbin", "glm", "lm")) | identical(class(yreg), c("multinom", "nnet")) |
               identical(class(yreg), c("gam", "glm", "lm")) | identical(class(yreg), "polr") |
               identical(class(yreg), "coxph") | identical(class(yreg), "survreg"))) {
    stop("Fit yreg by lm, glm, glm.nb, gam, multinom, polr, coxph, survreg")
  }
  
  # for delta method inference, use survey regressions for yreg and mreg when weights are applied
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    yreg <- eval(bquote(survey::svyglm(formula = .(formula(yreg)), family = .(family(yreg)),
                               design = survey::svydesign(~ 1, data = .(data)))))
  }
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    if (inherits(mreg[[1]], "glm")) mreg[[1]] <- eval(bquote(survey::svyglm(formula = .(formula(mreg[[1]])), 
                                                                    family = .(family(mreg[[1]])),
                                                                    design = survey::svydesign(~ 1, data = .(data)))))
    if (inherits(mreg[[1]], "multinom")) mreg[[1]] <- eval(bquote(svymultinom(formula = .(formula(mreg[[1]])), 
                                                                              data = .(data))))
  }
  out <- list(yreg = yreg, ereg = ereg, mreg = mreg, wmreg = wmreg, postcreg = postcreg)
  return(out)
}