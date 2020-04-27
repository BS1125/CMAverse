run_regressions <- function(formulas, data, model, exposure, mediator, postc,
                            yreg, mreg, ereg, postcreg, wmreg) {

  require(dplyr)
  require(survival)

  ###########################################Exposure Regression####################################

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

        stop("Selected model only supports logistic/multinomial/ordinal non user-defined exposure regression model")

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

        stop("Selected model doesn't support this user-defined exposure regression object")

      }

      exposure_regression <- update(ereg, data = data,
                                    formula. = formula(ereg))

    }

    if (model == "msm") {

      pa <- left_join(select(data, exposure),
                      count(data, !!as.name(exposure)),
                      by = exposure)[, "n"]/nrow(data)

      if (inherits(exposure_regression, "glm") &&
          family(exposure_regression)$family %in% c("binomial", "quasibinomial")) {

        wadenom.prob <- predict(exposure_regression, newdata = data,
                                type = "response")

        class <- as.numeric(as.factor(data[, exposure])) - 1

        wadenom <-  wadenom.prob ^class *
          (1 - wadenom.prob)^(1 - class)

      } else if ((inherits(exposure_regression, "gam") &&
                  (family(exposure_regression)$family == "multinom" |
                   startsWith(family(exposure_regression)$family, "Ordered Categorical")))|
                 inherits(exposure_regression, "multinom")|inherits(exposure_regression, "polr")) {

        wadenom.prob <- predict(exposure_regression, newdata = data,
                                type = ifelse(inherits(exposure_regression, "multinom")|
                                                inherits(exposure_regression, "polr"),
                                              "probs","response"))

        class <- as.numeric(as.factor(data[, exposure]))

        wadenom <- sapply(1:nrow(data),FUN = function(i) wadenom.prob[i, class[i]])

      }

      wa <- pa/wadenom

    } else if (model == "iorw") {

      if (inherits(exposure_regression, "glm") &&
          family(exposure_regression)$family %in% c("binomial", "quasibinomial")) {

        wadenom.prob <- predict(exposure_regression, newdata = data,
                                type = "response")

        class <- as.numeric(as.factor(data[, exposure])) - 1

        wa <-  (1 - wadenom.prob) /
          (wadenom.prob ^ class *
             (1 - wadenom.prob) ^ (1 - class))


      } else if ((inherits(exposure_regression, "gam") &&
                  (family(exposure_regression)$family == "multinom" |
                   startsWith(family(exposure_regression)$family, "Ordered Categorical")))|
                 inherits(exposure_regression, "multinom")|inherits(exposure_regression, "polr")) {

        wadenom.prob <- predict(exposure_regression, newdata = data,
                                type = ifelse(inherits(exposure_regression, "multinom")|
                                                inherits(exposure_regression, "polr"),
                                              "probs","response"))

        class <- as.numeric(as.factor(data[, exposure]))

        wa <- prob[, 1] / sapply(1:nrow(data),FUN = function(i) wadenom.prob[i, class[i]])

      }

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
            (identical(class(mreg[[i]]), c("glm", "lm")) &&
             family(mreg[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                             "binomial", "quasibinomial","poisson", "quasipoisson"))|
            (identical(class(mreg[[i]]), c("gam", "glm", "lm")) &&
             (family(mreg[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi","multinom",
                                              "binomial", "quasibinomial","poisson", "quasipoisson") |
              startsWith(family(mreg[[i]])$family,"Negative Binomial")|
              startsWith(family(mreg[[i]])$family,"Ordered Categorical")))|
            identical(class(mreg[[i]]), c("multinom", "nnet"))|
            identical(class(mreg[[i]]), "polr")|identical(class(mreg[[i]]), "nls")|
            identical(class(mreg[[i]]), c("negbin", "glm", "lm")))){

        stop("Unsupported user-defined mediator regression object")

      }

      mediator_regression[[i]] <- update(mreg[[i]], data = data,
                                         formula. = formula(mreg[[i]]))

    }

    for (i in mediator.nomodel.index) {

      if (mreg[[i]] == "linear") {

        mediator_regression[[i]] <- glm(mediator_formula[[i]], family = gaussian(), data = data)

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
        stop("Unsupported non user-defined mediator regression model")
      }
    }

    if (model == "msm") {

      for (i in length(mediator_regression)) {

        if (!((identical(class(mediator_regression[[i]]), c("glm", "lm")) &&
               family(mediator_regression[[i]])$family %in% c("binomial", "quasibinomial"))|
              (identical(class(mediator_regression[[i]]), c("gam", "glm", "lm")) &&
               (family(mediator_regression[[i]])$family %in% c("multinom",
                                                               "binomial", "quasibinomial") |
                startsWith(family(mediator_regression[[i]])$family,"Ordered Categorical")))|
              identical(class(mediator_regression[[i]]), c("multinom", "nnet"))|
              identical(class(mediator_regression[[i]]), "polr"))){

          stop("Selected model doesn't support this mediator regression object")

        }

        mediator_regression[[i]]$call[["weights"]] <- wa

        mediator_regression[[i]] <- update(mediator_regression[[i]],
                                           formula. = formula(mediator_regression[[i]]))

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

        stop("Selected model doesn't support this user-defined mediator regression object for weights")

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
        stop("Selected model doesn't support this non user-defined mediator regression model for weights")
      }
    }


    wmnom_regression <- list()

    for (i in length(mediator)) {

      wmnom_regression[[i]] <- update(wmdenom_regression[[i]],
                                      formula. = as.formula(paste0(mediator[i], "~", exposure)))


    }

    wmnom <- rep(1, nrow(data))

    wmdenom <- rep(1, nrow(data))

    for (i in 1:length(wmdenom_regression)) {

      if (inherits(wmdenom_regression[[i]], "glm") &&
          family(wmdenom_regression[[i]])$family %in% c("binomial", "quasibinomial")) {

        wmdenom.prob <- predict(wmdenom_regression[[i]], newdata = data,
                                type = "response")

        wmnom.prob <- predict(wmnom_regression[[i]], newdata = data,
                              type = "response")

        class <- as.numeric(as.factor(data[, mediator[i]])) - 1

        wmdenom <-  wmdenom * wmdenom.prob ^ class *
          (1 - wmdenom.prob)^(1 - class)

        wmnom <-  wmnom * wmnom.prob ^ class *
          (1 - wmnom.prob)^(1 - class)

      } else if ((inherits(wmdenom_regression[[i]], "gam") &&
                  (family(wmdenom_regression[[i]])$family == "multinom" |
                   startsWith(family(wmdenom_regression[[i]])$family, "Ordered Categorical")))|
                 inherits(wmdenom_regression[[i]], "multinom")|inherits(wmdenom_regression[[i]], "polr")) {

        wmdenom.prob <- predict(wmdenom_regression[[i]], newdata = data,
                                type = ifelse(inherits(wmdenom_regression[[i]], "multinom")|
                                                inherits(wmdenom_regression[[i]], "polr"),
                                              "probs","response"))

        wmnom.prob <- predict(wmnom_regression[[i]], newdata = data,
                              type = ifelse(inherits(wmnom_regression[[i]], "multinom")|
                                              inherits(wmnom_regression[[i]], "polr"),
                                            "probs","response"))

        class <- as.numeric(as.factor(data[, mediator[i]]))

        wmdenom <- wmdenom * sapply(1:nrow(data), FUN = function(i) wmdenom.prob[i, class[i]])

        wmnom <- wmnom * sapply(1:nrow(data), FUN = function(i) wmnom.prob[i, class[i]])

      }

    }

    wm <- wmnom / wmdenom

    wy <- wa * wm

  }

  #####################################Post-exposure Confounder Regression###########################

  if (model == "g-formula" && !is.null(postc)) {

    postc_formula <- formulas$postc_formula

    postc.nomodel.index <- which(sapply(1:length(postcreg), function(x) is.character(postcreg[[x]])))

    postc.model.index <- which(sapply(1:length(postcreg), function(x) !is.character(postcreg[[x]])))

    postc_regression <- list()

    for (i in postc.model.index) {

      if (!(identical(class(postcreg[[i]]), "lm") |
            (identical(class(postcreg[[i]]), c("glm", "lm")) &&
             family(postcreg[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                 "binomial", "quasibinomial","poisson", "quasipoisson"))|
            (identical(class(postcreg[[i]]), c("gam", "glm", "lm")) &&
             (family(postcreg[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi","multinom",
                                                  "binomial", "quasibinomial","poisson", "quasipoisson") |
              startsWith(family(postcreg[[i]])$family,"Negative Binomial")|
              startsWith(family(postcreg[[i]])$family,"Ordered Categorical")))|
            identical(class(postcreg[[i]]), c("multinom", "nnet"))|
            identical(class(postcreg[[i]]), "polr")|identical(class(postcreg[[i]]), "nls")|
            identical(class(postcreg[[i]]), c("negbin", "glm", "lm")))){

        stop("Unsupported user-defined post-exposure counfounder regression object")

      }

      postc_regression[[i]] <- update(postcreg[[i]], data = data,
                                      formula. = formula(postcreg[[i]]))

    }

    for (i in postc.nomodel.index) {

      if (postcreg[[i]] == "linear") {

        postc_regression[[i]] <- glm(postc_formula[[i]], family = gaussian(), data = data)

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
        stop("Unsupported non user-defined post-exposure confounder regression model")
      }
    }

  } else postc_regression <- NULL


  ###########################################Outcome Regression######################################

  outcome_formula <- formulas$outcome_formula

  if(is.character(yreg)) {

    if (yreg == "linear") {

      outcome_regression  <- glm(outcome_formula, family = gaussian(), data = data)

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

      outcome_regression <- survreg(as.formula(outcome_formula), dist = "exponential", data = data)

    } else if (yreg == "aft_weibull") {

      outcome_regression <- survreg(as.formula(outcome_formula), dist = "weibull", data = data)

    } else {
      stop("Unsupported non user-defined outcome regression")
    }

  } else {

    if (!(identical(class(yreg), "lm") |
          (identical(class(yreg), c("glm", "lm")) &&
           family(yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                      "binomial", "quasibinomial","poisson", "quasipoisson"))|
          (identical(class(yreg), c("gam", "glm", "lm")) &&
           (family(yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi","multinom",
                                       "binomial", "quasibinomial","poisson", "quasipoisson") |
            startsWith(family(yreg)$family,"Negative Binomial")|
            startsWith(family(yreg)$family,"Ordered Categorical")))|
          identical(class(yreg), c("multinom", "nnet"))|
          identical(class(yreg), "polr")|identical(class(yreg), "nls")|
          identical(class(yreg), c("negbin", "glm", "lm"))|
          identical(class(yreg), "survreg")|
          identical(class(yreg), "coxph"))){

      stop("Unsupported user-defined outcome regression object")

    }

    outcome_regression <- update(yreg, data = data,
                                 formula. = formula(yreg))

  }

  if (model == "msm") {

    outcome_regression$call[["weights"]] <- wy

    outcome_regression <- update(outcome_regression,
                                 formula. = formula(outcome_regression))

  } else if (model == "iorw") {

    tot_regression <- outcome_regression

    dir_regression <- outcome_regression

    dir_regression$call[["weights"]] <- wa

    dir_regression <- update(dir_regression,
                             formula. = formula(dir_regression))

  }


  ###########################################Return Regression List#############################

  if (model %in% c("rb", "msm")) {

    regressions <- list(mediator_regression = mediator_regression,
                        outcome_regression =  outcome_regression)

  } else if (model == "wb") {

    regressions <- list(exposure_regression = exposure_regression,
                        outcome_regression =  outcome_regression)

  } else if (model == "g-formula") {

    regressions <- list(postc_regression = postc_regression,
                        mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression)

  } else if (model == "iorw") {

    regressions <- list(tot_outcome_regression = tot_outcome_regression,
                        dir_outcome_regression = dir_outcome_regression)

  } else if (model == "ne") {

    regressions <- list(outcome_regression =  outcome_regression)

  }

  return(regressions)

}
