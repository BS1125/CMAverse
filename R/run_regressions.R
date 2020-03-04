run_regressions <- function(model, formulas, exposure, exposure.type, mediator,
                            covariates.post, covariates.post.type,
                            yreg, mreg, data) {

  if (model %in% c("wb", "iorw", "ne")) mreg = NULL

  outcome_formula <- formulas$outcome_formula

  mediator_formula <- formulas$mediator_formula

  if (!is.null(mreg)) {

    mediator.nomodel.index <- which(sapply(1:length(mreg), function(x) is.character(mreg[[x]])))

    mediator.model.index <- which(sapply(1:length(mreg), function(x) !is.character(mreg[[x]])))

    mediator_regression <- list()

    for (i in mediator.model.index) {
      if (!(inherits(mreg[[i]],"glm")|inherits(mreg[[i]],"lm")|inherits(mreg[[i]],"multinom")|
            inherits(mreg[[i]],"survreg")|inherits(mreg[[i]],"coxph"))) {
        stop("Only glm/lm/multinom/survreg/coxph user-defined mediator models are supported")
      }

      mediator_regression[[i]] <- mreg[[i]]

      mediator_regression[[i]]$call$data <- data

      mediator_regression[[i]] <- update(mediator_regression[[i]])

    }

    for (i in 1:length(mediator.nomodel.index)) {

      if (!(mreg[[i]] %in% c("linear", "logistic", "multinomial"))) {
        stop("Unsupported non user-defined mediator regression model")
      }

      if (length(mediator_formula) > 1) {

        if (mreg[[mediator.nomodel.index[i]]] == "linear") {
          mediator_regression[[mediator.nomodel.index[i]]] <-
            lm(mediator_formula[[i]], data = data)
        } else if (mreg[[i]] == "logistic") {
          mediator_regression[[mediator.nomodel.index[i]]] <-
            glm(mediator_formula[[i]],
                family = binomial("logit"), data = data)
        } else if (mreg[[i]] == "multinomial") {
          mediator_regression[[mediator.nomodel.index[i]]] <-
            nnet::multinom(mediator_formula[[i]], data = data, trace = FALSE)
        }

      } else if (length(mediator_formula) == 1) {

        if (mreg[[mediator.nomodel.index[i]]] == "linear") {
          mediator_regression[[mediator.nomodel.index[i]]] <-
            lm(mediator_formula[[1]], data = data)
        } else if (mreg[[i]] == "logistic") {
          mediator_regression[[mediator.nomodel.index[i]]] <-
            glm(mediator_formula[[1]],
                family = binomial("logit"), data = data)
        } else if (mreg[[i]] == "multinomial") {
          mediator_regression[[mediator.nomodel.index[i]]] <-
            nnet::multinom(mediator_formula[[1]], data = data, trace = FALSE)
        }


      }
    }

  } else {mediator_regression <- NULL}

  if(is.character(yreg)) {

    if (!(yreg %in% c("linear", "logistic", "loglinear", "poisson",
                      "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull"))) {
      stop("Unsupported non user-defined outcome regression model")
    }

    if (yreg == "linear") {
      outcome_regression  <- glm(outcome_formula, family = gaussian(), data = data)
    }

    if (yreg == "logistic") {
      outcome_regression  <- glm(outcome_formula, family = binomial(), data = data)
    }

    if (yreg == "loglinear") {
      outcome_regression  <- glm(outcome_formula, family = binomial("log"), data = data)
    }

    if (yreg == "poisson") {
      outcome_regression  <- glm(outcome_formula, family = poisson(), data = data)
    }

    if (yreg == "quasipoisson") {
      outcome_regression  <- glm(outcome_formula, family = quasipoisson(), data = data)
    }

    if (yreg == "negbin") {
      outcome_regression  <- glm.nb(outcome_formula, data = data)
    }

    if (yreg == "coxph") {
      outcome_regression <- survival::coxph(as.formula(outcome_formula), data = data)
    }

    if (yreg == "aft_exp") {
      outcome_regression <- survreg(as.formula(outcome_formula), dist = "exponential", data = data)
    }

    if (yreg == "aft_weibull") {
      outcome_regression <- survreg(as.formula(outcome_formula), dist = "weibull", data = data)
    }
  } else {

    if (!(inherits(yreg,"glm")|inherits(yreg,"lm")|
          inherits(yreg,"survreg")|inherits(yreg,"coxph"))) {
      stop("Only glm/lm/survreg/coxph user-defined outcome models are supported")
    }

    outcome_regression <- update(yreg, data = data)

  }

  #######################################Weighting-based Approach##############################

  if (model == "wb") {

    exposure_formula <- formulas$exposure_formula

    if (exposure.type == "continuous") {

      exposure_regression <- lm(exposure_formula, data = data)

    } else if (exposure.type == "binary") {

      exposure_regression <- glm(exposure_formula, data = data, family = binomial())

    }

    regressions <- list(exposure_regression = exposure_regression,
                        outcome_regression =  outcome_regression)

    #######################################Marginal Structural Model#################################

  } else if (model == "msm") {

    wa_denom_formula <- formulas$wa_denom_formula

    wz_nom_formula <- formulas$wz_nom_formula

    wz_denom_formula <- formulas$wz_denom_formula

    cde_outcome_formula <- formulas$cde_outcome_formula

    pa <- left_join(select(data, exposure),
                    count(data, !!as.name(exposure)),
                    by = exposure)[, "n"]/nrow(data)

    if (exposure.type == "continuous") {

      model.wa.denom <- lm(wa_denom_formula, data = data)

      wa.denom <- dnorm(data[, exposure],
                        mean = model.wa.denom$fitted.values,
                        sd = summary(model.wa.denom)$sigma)

    } else if (exposure.type == "binary") {

      model.wa.denom <- glm(wa_denom_formula, family = binomial(link = "logit"), data = data)

      wa.denom <- predict(model.wa.denom, newdata = data,
                          type = "response")

    }

    wa <- pa / wa.denom

    model.wz.nom <- mediator_regression

    model.wz.denom <- mediator_regression

    wz.nom <- rep(1, nrow(data))

    wz.denom <- rep(1, nrow(data))

    for (i in 1:length(mediator_regression)) {

      mediator_regression[[i]]$call[["weights"]] <- wa
      mediator_regression[[i]] <- update(mediator_regression[[i]])

      model.wz.nom[[i]]$call$formula <- wz_nom_formula[[i]]
      model.wz.nom[[i]] <- update(model.wz.nom[[i]])

      model.wz.denom[[i]]$call$formula <- wz_denom_formula[[i]]
      model.wz.denom[[i]] <- update(model.wz.denom[[i]])

      if (inherits(mediator_regression[[i]], "glm")|inherits(mediator_regression[[i]], "survreg")|
          inherits(mediator_regression[[i]], "coxph")) {

        type <- ifelse(inherits(mediator_regression[[i]], "coxph"), "risk", "response")

        wz.nom <- wz.nom * predict(model.wz.nom[[i]], newdata = data, type = type) ^ data[, mediator[i]] *
          (1 - predict(model.wz.nom[[i]], newdata = data, type = type))^(1 - data[, mediator[i]])

        wz.denom <- wz.denom * predict(model.wz.denom[[i]], newdata = data, type = type) ^ data[, mediator[i]] *
          (1 - predict(model.wz.denom[[i]], newdata = data, type = type))^(1 - data[, mediator[i]])

      } else if (inherits(mediator_regression[[i]], "lm")){

        wz.nom <- wz.nom * dnorm(data[, mediator[i]],
                                 mean = model.wz.nom[[i]]$fitted.values,
                                 sd = summary(model.wz.nom[[i]])$sigma)

        wz.denom <- wz.denom * dnorm(data[, mediator[i]],
                                     mean = model.wz.denom[[i]]$fitted.values,
                                     sd = summary(model.wz.denom[[i]])$sigma)

      } else if (inherits(mediator_regression[[i]], "multinom")) {

        class <- as.numeric(data$M_cat)

        prob.wz.nom <- predict(model.wz.nom[[i]],
                               newdata = data, type = "prob")

        wz.nom <- wz.nom * sapply(1:nrow(data),FUN = function(i) prob.wz.nom[i, class[i]])

        prob.wz.denom <- predict(model.wz.denom[[i]],
                                 newdata = data, type = "prob")

        wz.denom <- wz.denom * sapply(1:nrow(data),FUN = function(i) prob.wz.denom[i, class[i]])

      }

    }

    wz <- wz.nom / wz.denom

    outcome_regression$call[["weights"]] <- wa * wz

    outcome_regression <- update(outcome_regression)

    cde_outcome_regression <- update(outcome_regression, formula. = cde_outcome_formula)

    regressions <- list(mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression,
                        cde_outcome_regression = cde_outcome_regression)

  } else if (model == "iorw") {

    ####################################Inverse Odds Ratio Weighting Approach########################

    exposure_formula <- formulas$exposure_formula

    if (exposure.type == "continuous") {

      exposure_regression <- lm(exposure_formula, data = data)

      w <- dnorm(rep(0, nrow(data)),
                 mean = exposure_regression$fitted.values,
                 sd = summary(exposure_regression)$sigma) / dnorm(data[, exposure],
                                                                  mean = exposure_regression$fitted.values,
                                                                  sd = summary(exposure_regression)$sigma)

    } else if (exposure.type == "binary") {

      exposure_regression <- glm(exposure_formula, data = data, family = binomial())

      w <- (1 - predict(exposure_regression, newdata = data, type = "response")) /
        (predict(exposure_regression, newdata = data, type = "response") ^ (as.numeric(data[, exposure])-1) *
           (1 - predict(exposure_regression, newdata = data, type = "response")) ^ (1 - (as.numeric(data[, exposure])-1)))


    }

    tot_outcome_regression <- outcome_regression

    outcome_regression$call[["weights"]] <- w

    dir_outcome_regression <- update(outcome_regression)

    regressions <- list(tot_outcome_regression = tot_outcome_regression,
                        dir_outcome_regression = dir_outcome_regression)

  } else if (model == "g-formula") {

    ########################################G-formula Approach#######################################

    if (!is.null(covariates.post)){

      for (i in 1:length(covariates.post)) {
        if (!(covariates.post.type[i] %in% c("continuous", "binary", "categorical"))) {
          stop("Unsupported types for post-exposure covariates")
        }
      }

      postcovar_formula <- formulas$postcovar_formula

      postcovar_regression <- lapply(1:length(postcovar_formula), FUN = function(x) {
        if (covariates.post.type[x] == "continuous") {
          lm(postcovar_formula[[x]], data = data)
        } else if (covariates.post.type[x] == "binary") {
          glm(postcovar_formula[[x]],
              family = binomial("logit"), data = data)
        } else if (covariates.post.type[x] == "categorical") {
          nnet::multinom(postcovar_formula[[x]], data = data, trace = FALSE)
        }})

    } else postcovar_regression <- NULL

    regressions <- list(postcovar_regression = postcovar_regression,
                        mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression)

  } else {

    regressions <- list(mediator_regression = mediator_regression,
                        outcome_regression =  outcome_regression)

  }

  return(regressions)

}
