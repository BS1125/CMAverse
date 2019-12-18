run_regressions <- function(model, formulas = NULL, exposure.type, covariates.post.type,
                            yreg, mreg, data = NULL) {

  outcome_formula <- formulas$outcome_formula

  mediator_formula <- formulas$mediator_formula

  if (yreg == "linear") {
      outcome_regression  <- lm(outcome_formula, data = data)
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

  if (!is.null(mreg) && length(mediator_formula) == 1) {

    if (mreg == "linear") {
      mediator_regression <- lm(mediator_formula[[1]], data = data)
    } else if (mreg == "binary") {
      mediator_regression <- glm(mediator_formula[[1]], family = binomial("logit"),data = data)
    }
   } else if (!is.null(mreg) && length(mediator_formula) > 1){
    mediator_regression <- lapply(1:length(mediator_formula), FUN = function(x) {
      if (mreg[x] == "linear") {
        lm(mediator_formula[[x]], data = data)
      } else if (mreg[x] == "logistic") {
        glm(mediator_formula[[x]],
            family = binomial("logit"), data = data)}})
  } else mediator_regression <- NULL

#######################################Marginal Structural Model#################################

  if (model == "msm") {

    wa_denom_formula <-  formulas$wa_denom_formula

    wz_nom_formula <-  formulas$wz_nom_formula

    wz_denom_formula <-  formulas$wz_denom_formula

    cde_outcome_formula <-  formulas$cde_outcome_formula

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

    if (mreg == "linear") {

      model.wz.nom <- lm(wz_nom_formula, data = data)

      model.wz.denom <- lm(wz_denom_formula, data = data)

      wz.nom <- dnorm(data[, mediator],
                      mean = model.wz.nom$fitted.values,
                      sd = summary(model.wz.nom)$sigma)

      wz.denom <- dnorm(data[, mediator],
                        mean = model.wz.denom$fitted.values,
                        sd = summary(model.wz.denom)$sigma)

    } else if (mreg == "binary") {

      model.wz.nom <- glm(wz_nom_formula, family = binomial(link = "logit"), data = data)

      model.wz.denom <- glm(wz_denom_formula, family = binomial(link = "logit"), data = data)

      wz.nom <- predict(model.wz.nom, newdata = data, type = "response") * data[, mediator] +
        (1 - predict(model.wz.nom, newdata = data, type = "response"))^(1 - data[, mediator])

      wz.denom <- predict(model.wz.denom, newdata = data, type = "response") * data[, mediator] +
        (1 - predict(model.wz.denom, newdata = data, type = "response"))^(1 - data[, mediator])

    }

    wa <- pa / wa.denom

    wz <- wz.nom / wz.denom

    outcome_regression$call[["weights"]] <- wa * wz

    mediator_regression$call[["weights"]] <- wa

    outcome_regression <- update(outcome_regression)

    mediator_regression <- update(mediator_regression)

    outcome_regression$call$formula <- cde_outcome_formula

    cde_outcome_regression <- update(outcome_regression)

    regressions <- list(mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression,
                        cde_outcome_regression = cde_outcome_regression)

  }

####################################Inverse Odds Ratio Weighting Approach########################

  if (model == "iorw") {

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
        predict(exposure_regression, newdata = data, type = "response") * data[, exposure] +
        (1 - predict(exposure_regression, newdata = data, type = "response")) * (1 - data[, exposure])


    }

    tot_outcome_regression <- outcome_regression

    outcome_regression$call[["weights"]] <- w

    dir_outcome_regression <- update(outcome_regression)

    regressions <- list(tot_outcome_regression = tot_outcome_regression,
                        dir_outcome_regression = dir_outcome_regression)
  }

########################################G-formula Approach#######################################

  if (model == "g-formula") {

    post_covar_formula <- formulas$post_covar_formula

    post_covar_regression <- lapply(1:length(post_covar_formula), FUN = function(x) {
      if (covariates.post.type[x] == "continuous") {
        lm(post_covar_formula[[x]], data = data)
      } else if (covariates.post.type[x] == "binary") {
        glm(post_covar_formula[[x]],
            family = binomial("logit"), data = data)}})

    regressions <- list(post_covar_regression = post_covar_regression,
                        mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression)

  }

  return(regressions)

}
