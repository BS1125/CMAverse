run_regressions <- function(formulas, data, model, exposure, mediator, postc,
                            yreg, mreg, ereg, postcreg, wmreg) {

  ###########################################Exposure Regression####################################

  # run the exposure regression

  if (model %in% c("wb", "iorw", "msm")) {

    exposure_formula <- formulas$exposure_formula

    if (is.character(ereg)) {

      if (ereg == "logistic") {

        exposure_regression <- glm(exposure_formula, data = data, family = binomial())

      } else if (ereg == "loglinear") {

        exposure_regression <- glm(exposure_formula, data = data, family = binomial("log"))

      } else if (ereg == "multinomial") {

        exposure_regression <- nnet::multinom(exposure_formula, data = data, trace = FALSE)

      } else if (ereg == "ordinal"){

        exposure_regression <- MASS::polr(exposure_formula, data = data)

      } else {

        stop("Selected model only supports logistic/loglinear/multinomial/ordinal non user-defined exposure regression objects")

      }

    } else {

      if (!((identical(class(ereg), c("glm", "lm")) &&
             family(ereg)$family %in% c("binomial", "quasibinomial"))|
            (identical(class(ereg), c("gam", "glm", "lm")) &&
             (family(ereg)$family %in% c("multinom",
                                         "binomial", "quasibinomial") |
              startsWith(family(ereg)$family,"Ordered Categorical")))|
            identical(class(ereg), c("multinom", "nnet"))|
            identical(class(ereg), "polr"))){

        stop("The selected model only applies when the exposure is categorical")

      }

      exposure_regression <- update(ereg, data = data,
                                    formula. = formula(ereg))

    }

  } else exposure_regression <- NULL

  ###########################################Mediator Regression######################################

  if (model %in% c("rb", "msm", "g-formula")) {

    mediator_formula <- formulas$mediator_formula

    mediator.nomodel.index <- which(sapply(1:length(mreg), function(x) is.character(mreg[[x]])))

    mediator.model.index <- which(sapply(1:length(mreg), function(x) !is.character(mreg[[x]])))

    mediator_regression <- list()

    for (i in mediator.model.index) {

      if (!(identical(class(mreg[[i]]), "lm") |
            (identical(class(mreg[[i]]), c("glm", "lm")))|
            (identical(class(mreg[[i]]), c("gam", "glm", "lm")))|
            identical(class(mreg[[i]]), c("multinom", "nnet"))|
            identical(class(mreg[[i]]), "polr")|
            identical(class(mreg[[i]]), "nls")|
            identical(class(mreg[[i]]), c("negbin", "glm", "lm")))) {

        stop("For mediator regression object, only lm(), glm(), glm.nb(), gam(), multinom(), polr() and nls() are supported")

      }

      mediator_regression[[i]] <- update(mreg[[i]], data = data,
                                         formula. = formula(mreg[[i]]))

    }

    for (i in mediator.nomodel.index) {

      if (mreg[[i]] == "linear") {

        mediator_regression[[i]] <- lm(mediator_formula[[i]], data = data)

      } else if (mreg[[i]] == "logistic") {

        mediator_regression[[i]] <- glm(mediator_formula[[i]], family = binomial(),
                                        data = data)

      } else if (mreg[[i]] == "loglinear") {

        mediator_regression[[i]] <- glm(mediator_formula[[i]], family = binomial("log"), data = data)

      } else if (mreg[[i]] == "poisson") {

        mediator_regression[[i]]  <- glm(mediator_formula[[i]], family = poisson(), data = data)

      } else if (mreg[[i]] == "quasipoisson") {

        mediator_regression[[i]] <- glm(mediator_formula[[i]], family = quasipoisson(), data = data)

      } else if (mreg[[i]] == "negbin") {

        mediator_regression[[i]] <- MASS::glm.nb(mediator_formula[[i]], data = data)

      } else if (mreg[[i]] == "multinomial") {

        mediator_regression[[i]] <- nnet::multinom(mediator_formula[[i]], data = data,
                                                   trace = FALSE)

      } else if (mreg[[i]] == "ordinal"){

        mediator_regression[[i]] <- MASS::polr(mediator_formula[[i]], data = data)

      } else {
        stop("Only supports linear/logistic/loglinear/poisson/quasipoisson/negbin/multinomial/ordinal non user-defined mediator regression objects")
      }
    }


    if (model == "msm") {

      for (i in 1:length(mediator_regression)) {

        if (!((identical(class(mediator_regression[[i]]), c("glm", "lm")) &&
               family(mediator_regression[[i]])$family %in% c("binomial", "quasibinomial"))|
              (identical(class(mediator_regression[[i]]), c("gam", "glm", "lm")) &&
               (family(mediator_regression[[i]])$family %in% c("multinom",
                                                               "binomial", "quasibinomial") |
                startsWith(family(mediator_regression[[i]])$family,"Ordered Categorical")))|
              identical(class(mediator_regression[[i]]), c("multinom", "nnet"))|
              identical(class(mediator_regression[[i]]), "polr"))){

          stop("The selected model only applies for categorical mediators")

        }
      }
    }
  } else mediator_regression <- NULL

  ####################################Mediator Regression For Weights##############################

  if (model == "msm") {

    wm_nom_formula <- formulas$wm_nom_formula

    wm_denom_formula <- formulas$wm_denom_formula

    wm.nomodel.index <- which(sapply(1:length(wmreg), function(x) is.character(wmreg[[x]])))

    wm.model.index <- which(sapply(1:length(wmreg), function(x) !is.character(wmreg[[x]])))

    wmdenom_regression <- list()

    for (i in wm.model.index) {

      if (!((identical(class(wmreg[[i]]), c("glm", "lm")) &&
             family(wmreg[[i]])$family %in% c("binomial", "quasibinomial"))|
            (identical(class(wmreg[[i]]), c("gam", "glm", "lm")) &&
             (family(wmreg[[i]])$family %in% c("multinom",
                                               "binomial", "quasibinomial") |
              startsWith(family(wmreg[[i]])$family,"Ordered Categorical")))|
            identical(class(wmreg[[i]]), c("multinom", "nnet"))|
            identical(class(wmreg[[i]]), "polr"))){

        stop("The selected model only applies for categorical mediators")

      }

      wmdenom_regression[[i]] <- update(wmreg[[i]], data = data,
                                        formula. = formula(wmreg[[i]]))

    }

    for (i in wm.nomodel.index) {

      if (wmreg[[i]] == "logistic") {

        wmdenom_regression[[i]] <- glm(mediator_formula[[i]], family = binomial(),
                                       data = data)

      } else if (wmreg[[i]] == "loglinear") {

        wmdenom_regression[[i]] <- glm(mediator_formula[[i]], family = binomial("log"), data = data)

      }  else if (wmreg[[i]] == "multinomial") {

        wmdenom_regression[[i]] <- nnet::multinom(mediator_formula[[i]], data = data,
                                                  trace = FALSE)

      } else if (wmreg[[i]] == "ordinal"){

        wmdenom_regression[[i]] <- MASS::polr(mediator_formula[[i]], data = data)

      } else {
        stop("Selected model only supports logistic/loglinear/multinomial/ordinal non user-defined mediator regression objects for weights")

      }
    }

    wmnom_regression <- list()

    for (i in length(mediator)) {

      wmnom_regression[[i]] <- update(wmdenom_regression[[i]],
                                      formula. = as.formula(paste0(mediator[i], "~", exposure)))


    }
  }

  #####################################Post-exposure Confounder Regression###########################

  if (model == "g-formula" && !is.null(postc)) {

    postc_formula <- formulas$postc_formula

    postc.nomodel.index <- which(sapply(1:length(postcreg), function(x) is.character(postcreg[[x]])))

    postc.model.index <- which(sapply(1:length(postcreg), function(x) !is.character(postcreg[[x]])))

    postc_regression <- list()

    for (i in postc.model.index) {

      if (!(identical(class(postcreg[[i]]), "lm") |
            identical(class(postcreg[[i]]), c("glm", "lm"))|
            identical(class(postcreg[[i]]), c("gam", "glm", "lm"))|
            identical(class(postcreg[[i]]), c("multinom", "nnet"))|
            identical(class(postcreg[[i]]), "polr")|
            identical(class(postcreg[[i]]), "nls")|
            identical(class(postcreg[[i]]), c("negbin", "glm", "lm")))){

        stop("For post-exposure confounder regression object, only lm(), glm(), glm.nb(), gam(), multinom(), polr() and nls() are supported")

      }

      postc_regression[[i]] <- update(postcreg[[i]], data = data,
                                      formula. = formula(postcreg[[i]]))

    }

    for (i in postc.nomodel.index) {

      if (postcreg[[i]] == "linear") {

        postc_regression[[i]] <- lm(postc_formula[[i]], data = data)

      } else if (postcreg[[i]] == "logistic") {

        postc_regression[[i]] <- glm(postc_formula[[i]], family = binomial(),
                                     data = data)

      } else if (postcreg[[i]] == "loglinear") {

        postc_regression[[i]] <- glm(postc_formula[[i]], family = binomial("log"), data = data)

      } else if (postcreg[[i]] == "poisson") {

        postc_regression[[i]]  <- glm(postc_formula[[i]], family = poisson(), data = data)

      } else if (postcreg[[i]] == "quasipoisson") {

        postc_regression[[i]] <- glm(postc_formula[[i]], family = quasipoisson(), data = data)

      } else if (postcreg[[i]] == "negbin") {

        postc_regression[[i]] <- MASS::glm.nb(postc_formula[[i]], data = data)

      } else if (postcreg[[i]] == "multinomial") {

        postc_regression[[i]] <- nnet::multinom(postc_formula[[i]], data = data,
                                                trace = FALSE)

      } else if (postcreg[[i]] == "ordinal"){

        postc_regression[[i]] <- MASS::polr(postc_formula[[i]], data = data)

      } else {
        stop("Only supports linear/logistic/loglinear/poisson/quasipoisson/negbin/multinomial/ordinal non user-defined post-exposure confounder regression objects")
      }
    }

  } else postc_regression <- NULL

  ###########################################Outcome Regression######################################

  outcome_formula <- formulas$outcome_formula

  if(is.character(yreg)) {

    if (yreg == "linear") {

      outcome_regression  <- lm(outcome_formula, data = data)

    } else if (yreg == "logistic") {

      outcome_regression  <- glm(outcome_formula, family = binomial(), data = data)

    } else if (yreg == "loglinear") {

      outcome_regression  <- glm(outcome_formula, family = binomial("log"), data = data)

    } else if (yreg == "poisson") {

      outcome_regression  <- glm(outcome_formula, family = poisson(), data = data)

    } else if (yreg == "quasipoisson") {

      outcome_regression  <- glm(outcome_formula, family = quasipoisson(), data = data)

    } else if (yreg == "negbin") {

      outcome_regression  <- MASS::glm.nb(outcome_formula, data = data)

    } else if (yreg == "multinomial") {

      outcome_regression  <- nnet::multinom(outcome_formula, data = data, trace = FALSE)

    } else if (yreg == "ordinal") {

      outcome_regression  <- MASS::polr(outcome_formula, data = data)

    } else if (yreg == "coxph") {

      outcome_regression <- survival::coxph(as.formula(outcome_formula), data = data)

    } else if (yreg == "aft_exp") {

      outcome_regression <- survival::survreg(as.formula(outcome_formula), dist = "exponential", data = data)

    } else if (yreg == "aft_weibull") {

      outcome_regression <- survival::survreg(as.formula(outcome_formula), dist = "weibull", data = data)

    } else {
      stop("Only supports linear/logistic/loglinear/poisson/quasipoisson/negbin/multinomial/ordinal/coxph/aft_exp/aft_weibull non user-defined outcome regression objects")
    }

  } else {

    if (!(identical(class(yreg), "lm") |
          identical(class(yreg), c("glm", "lm"))|
          identical(class(yreg), c("gam", "glm", "lm"))|
          identical(class(yreg), c("multinom", "nnet"))|
          identical(class(yreg), "polr")|
          identical(class(yreg), "nls")|
          identical(class(yreg), c("negbin", "glm", "lm"))|
          identical(class(yreg), "survreg")|
          identical(class(yreg), "coxph"))){

      stop("For outcome regression object, only lm(), glm(), glm.nb(), gam(), multinom(), polr(), nls(), coxph() and survreg() are supported")

    }

    outcome_regression <- update(yreg, data = data,
                                 formula. = formula(yreg))

  }

  ###########################################Return Regression List#############################

  if (model == "rb") {

    regressions <- list(mediator_regression = mediator_regression,
                        outcome_regression =  outcome_regression)

    regressions$outcome_regression <- update(regressions$outcome_regression,
                                             formula. = formula(regressions$outcome_regression))

    for (i in 1:length(mediator_regression)) {
      regressions$mediator_regression[[i]] <- update(regressions$mediator_regression[[i]],
                                                     formula. = formula(regressions$mediator_regression[[i]]))
    }

  } else if (model == "msm") {

    regressions <- list(exposure_regression = exposure_regression,
                        wmnom_regression = wmnom_regression,
                        wmdenom_regression = wmdenom_regression,
                        mediator_regression = mediator_regression,
                        outcome_regression =  outcome_regression)

    regressions$outcome_regression <- update(regressions$outcome_regression,
                                             formula. = formula(regressions$outcome_regression))

    regressions$exposure_regression <- update(regressions$exposure_regression,
                                              formula. = formula(regressions$exposure_regression))

    for (i in 1:length(mediator)) {

      regressions$mediator_regression[[i]] <- update(regressions$mediator_regression[[i]],
                                                     formula. = formula(regressions$mediator_regression[[i]]))

      regressions$wmnom_regression[[i]] <- update(regressions$wmnom_regression[[i]],
                                                  formula. = formula(regressions$wmnom_regression[[i]]))

      regressions$wmdenom_regression[[i]] <- update(regressions$wmdenom_regression[[i]],
                                                    formula. = formula(regressions$wmdenom_regression[[i]]))

    }

  } else if (model == "wb") {

    regressions <- list(exposure_regression = exposure_regression,
                        outcome_regression =  outcome_regression)

    regressions$outcome_regression <- update(regressions$outcome_regression,
                                             formula. = formula(regressions$outcome_regression))

    regressions$exposure_regression <- update(regressions$exposure_regression,
                                              formula. = formula(regressions$exposure_regression))

  } else if (model == "g-formula") {

    regressions <- list(postc_regression = postc_regression,
                        mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression)

    regressions$outcome_regression <- update(regressions$outcome_regression,
                                             formula. = formula(regressions$outcome_regression))

    for (i in 1:length(mediator_regression)) {
      regressions$mediator_regression[[i]] <- update(regressions$mediator_regression[[i]],
                                                     formula. = formula(regressions$mediator_regression[[i]]))
    }

    if (!is.null(postc_regression)) {

      for (i in 1:length(postc_regression)) {
        regressions$postc_regression[[i]] <- update(regressions$postc_regression[[i]],
                                                    formula. = formula(regressions$postc_regression[[i]]))
      }
    }

  } else if (model == "iorw") {

    regressions <- list(exposure_regression = exposure_regression,
                        outcome_regression = outcome_regression)

    regressions$exposure_regression <- update(regressions$exposure_regression,
                                              formula. = formula(regressions$exposure_regression))

    regressions$outcome_regression <- update(regressions$outcome_regression,
                                             formula. = formula(regressions$outcome_regression))

  } else if (model == "ne") {

    regressions <- list(outcome_regression =  outcome_regression)

    regressions$outcome_regression <- update(regressions$outcome_regression,
                                             formula. = formula(regressions$outcome_regression))

  }

  return(regressions)

}
