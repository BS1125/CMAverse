causal_mediation <- function(data = NULL, outcome = NULL, exposure = NULL, exposure.type = "binary",
                             mediator = NULL, Mintertwined = FALSE, Mjoint = TRUE,
                             covariates.pre = NULL, covariates.post = NULL, vecc = NULL,
                             EMint = FALSE, MMint = FALSE, EMMint = FALSE,
                             EMint.terms = NULL, MMint.terms = NULL, EMMint.terms = NULL,
                             event = NULL, mreg = "linear", yreg = "linear",
                             m_star = NULL, a_star = 0, a = 1,
                             model = "standard", est.method = "delta",
                             nboot = 1000, conf = 0.95, nsims = 1000, nrep = 5) {

  require(dplyr)

  if (is.null(exposure) | is.null(outcome) | is.null(mediator)) {
    stop("Unspecified exposure, mediator, or outcome")
  }

  if (!(model %in% c("ne", "msm")) && length(mreg) < length(mediator)) {
    stop("Unspecified mediator model")
  }

  if (!all(mreg %in% c("linear", "logistic"))) {
    stop("Unsupported mediator regression model")
  }

  if (!(yreg %in% c("linear", "logistic", "loglinear", "poisson",
      "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull"))) {
    stop("Unsupported outcome regression model")
  }

  if (!(model %in% c("rb", "ne", "msm"))) {
    stop("Unsupported causal mediation model")
  }

  if (is.null(covariates.post) == FALSE && !(model %in% c("ne", "msm"))) {
    stop("Unsupported causal mediation model for post-exposure covariates.")
  }

  if (is.null(m_star) == TRUE) {
    m_star = rep(0, length(mediator))
  } else if (is.null(m_star) == FALSE && length(m_star) != length(mediator)) {
      stop("m_star doesn't have equal length with mediator")
  }

  if (EMint == FALSE) EMint.terms = NULL
  if (MMint == FALSE) MMint.terms = NULL

  if (EMMint == FALSE) {
    EMMint.terms = NULL
  } else {
    EMint = TRUE
    MMint = TRUE
  }

  if (length(mediator) == 1 && MMint == TRUE) {
    stop("Mediator-mediator interaction doesn't apply for a single mediator")
  }

  if (model %in% c("wb", "iorw", "ne")) mreg = NULL

  formulas <- create_formulas(exposure = exposure, mediator = mediator, outcome = outcome,
                              covariates.pre = covariates.pre, covariates.post = covariates.post,
                              EMint = EMint, MMint = MMint, EMMint = EMMint,
                              EMint.terms = EMint.terms, MMint.terms = MMint.terms,
                              EMMint.terms = EMMint.terms, event = event, mreg = mreg,
                              yreg = yreg, model = model)

  regressions <- run_regressions(model = model, formulas = formulas,
                                 mreg = mreg, yreg = yreg, data = data)

##############################################################################################################
##############################################Regression-based Model##########################################
##############################################################################################################

  if (model == "rb") {

    if (!(exposure.type %in% c("continuous", "binary"))) {
      stop("The regression-based model only supports binary and continuous exposure")
    }

    if ((is.null(covariates.pre) & !is.null(vecc))) {

      warning("Incompatible arguments")

    } else if (!is.null(covariates.pre) & is.null(vecc) & is.numeric(data[,covariates.pre])) {

      vecc <- colMeans(as.data.frame(data[, covariates.pre]))

    } else if (is.numeric(data[,covariates.pre]) == FALSE) {
      stop("The standard model only supports numeric covariates.pre ")
    }

    if (est.method %in% c("delta", "bootstrap")) {

      if (length(mediator) > 1) {
        stop("When length(mediator) > 1, the regression-based model only supports simulation-based estimation approach")
      }


      coef <- get_coef(regressions = regressions, mreg = mreg, yreg = yreg)

      n <- nrow(data)

      betas <- coef$betas

      thetas <- coef$thetas

      variance <- coef$variance

      vcov_thetas <- coef$vcov_thetas

      vcov_betas <- coef$vcov_betas

      vcov_block <- coef$vcov_block

      effect_estimates <- bootstrap_step(data = data, indices = c(1:nrow(data)), outcome = outcome,
                                         exposure = exposure, mediator = mediator, covariates.pre = covariates.pre,
                                         event = event, model = model, vecc = vecc,
                                         EMint = EMint, EMint.terms = EMint.terms,
                                         m_star = m_star, a_star = a_star, a = a,
                                         mreg = mreg, yreg = yreg)

####################################################Delta method################################################

      if (est.method == "delta") {

      effect_se <- delta_method(mediator = mediator, thetas = thetas, betas = betas,
                                variance = variance, vcov_block = vcov_block, EMint = EMint,
                                mreg = mreg, yreg, vecc = vecc, m_star = m_star,
                                a = a, a_star = a_star)

      out <- list(effect = effect_estimates, se = effect_se)

      }

#################################################Bootstrapping##################################################

      if (est.method == "bootstrap") {

      out <- boot::boot(data = data, statistic = bootstrap_step, R = nboot,
                      outcome = outcome, exposure = exposure,
                      mediator = mediator, model = model,
                      covariates.pre = covariates.pre, vecc = vecc, EMint = EMint,
                      EMint.terms = EMint.terms,
                      event = event, mreg = mreg, yreg = yreg,
                      m_star = m_star, a_star = a_star, a = a)

      effect_se <- apply(out$t, 2, sd)

      names(effect_se) <- paste0(names(effect_estimates), "_se")

      out <- list(effect = effect_estimates, se = effect_se)

    }

}
#############################################Simulation-based approach#########################################

      if (est.method == "simulation") {

        if (length(mediator) > 1 && MMint == TRUE) {
           stop("The regression-based model doesn't support mediator-mediator interaction")}

        out <- simulation_approach(data = data, formulas = formulas, exposure = exposure,
                                   mediator = mediator, outcome = outcome,
                                   covariates.pre = covariates.pre, covariates.post = covariates.post,
                                   mreg = mreg, yreg = yreg, model = model, nsims = nsims)

      }
    }

##########################################################################################################
#############################################Weighting-based Approach#####################################
##########################################################################################################

  if (model == "wb") {

    if (Mintertwined == TRUE && Mjoint == FALSE && length(mediator) > 1) {
      stop("Weighting-based model can't be used to calculate individual mediated effects for intertwined mediators")
    }

    if (exposure.type != "binary") {
      stop("The weighting-based approach only supports binary exposure")
    }

    effect_estimates <- bootstrap_step(data = data, indices = 1:nrow(data),
                                       outcome = outcome, exposure = exposure,
                                       model = "wb",
                                       mediator = mediator, covariates = covariates,
                                       EMint = EMint, MMint = MMint, EMMint = EMMint,
                                       EMint.terms = EMint.terms,
                                       MMint.terms = MMint.terms, EMMint.terms = EMMint.terms,
                                       event = event, mreg = NULL, yreg = yreg,
                                       m_star = m_star, a_star = a_star, a = a)

    out <- boot::boot(data = data, statistic = bootstrap_step, R = nboot,
                      outcome = outcome, exposure = exposure,
                      mediator = mediator, EMint = EMint, MMint = MMint, EMMint = EMMint,
                      EMint.terms = EMint.terms,
                      MMint.terms = MMint.terms, EMMint.terms = EMMint.terms,
                      covariates = covariates, model = "wb",
                      event = event, mreg = NULL, yreg = yreg,
                      m_star = m_star, a_star = a_star, a = a)

    se <- apply(out$t, 2, sd)

    out <- effect_estimates

    for (i in 1:length(effect_estimates))
      out[paste0(names(effect_estimates)[i], "_se")] <- se[i]

  }


################################################################################################################
###############################################Natural Effect Model#############################################
################################################################################################################

  if (model == "ne") {

    if (Mintertwined == TRUE) {
      stop("Natural effect model doesn't support intertwined mediators")
    }

    if (!(yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson"))) {
      stop("Unsupported outcome regression model for the natural effect model")
    }

    expData <- medflex::neImpute(regressions$outcome_regression, data = data,
                                 nMed = length(mediator), nrep = nrep)

    if (EMint == TRUE) {

      medflex_formula <- paste(c(paste(outcome, paste0(exposure, "0"), sep = "~"), paste0(exposure, "1"),
                         paste(paste0(exposure, "0"), paste0(exposure, "1"), sep = "*"), covariates),
                         collapse = "+")

    } else if (EMint == FALSE) {

      medflex_formula <- paste(c(paste(outcome, paste0(exposure, "0"), sep = "~"), paste0(exposure, "1"),
                                 covariates),
                               collapse = "+")

    }

    assign("medflex_formula", medflex_formula, envir = globalenv())

    family <- regressions$outcome_regression$family

    out <- summary(medflex::neEffdecomp(medflex::neModel(as.formula(medflex_formula),
                          family = family,
                          expData = expData, se = "robust")))$coefficients[, c("Estimate", "Std. Error")]

    if (yreg != "linear") {

      out[, "Estimate"] <- exp(out[, "Estimate"])

      for (i in 1:length(out[, "Std. Error"]))
        out[, "Std. Error"][i] <- msm::deltamethod(as.formula("~exp(x1)"), out[, "Estimate"][i],
                                                   out[, "Std. Error"][i])

    }

    rm(medflex_formula, pos = ".GlobalEnv")

  }

################################################################################################################
###############################################Marginal Structural Model########################################
################################################################################################################

if (model == "msm") {

  if (length(mediator) > 1) {
    stop("Marginal Structural Model only applys when length(mediator) = 1")
  }

  if (!(exposure.type %in% c("continuous", "binary"))) {
    stop("Marginal Structural Model only applys for binary or continuous exposure")
  }

  formulas <- create_formulas(exposure = exposure, mediator = mediator, outcome = outcome,
                              covariates.pre = covariates.pre, covariates.post = covariates.post,
                              EMint = EMint, MMint = MMint, EMMint = EMMint,
                              EMint.terms = EMint.terms, MMint.terms = MMint.terms,
                              EMMint.terms = EMMint.terms, event = event, mreg = mreg,
                              yreg = yreg, model = model)

  out <- simulation_approach(data = data, formulas = formulas, exposure = exposure,
                             mediator = mediator, outcome = outcome,
                             covariates.pre = covariates.pre, covariates.post = covariates.post,
                             mreg = mreg, yreg = yreg, model = model, nsims = nsims)

}


##################################################################################################
#######################################Inverse OR Weighting Approach##############################
##################################################################################################

if (model == "iorw") {

  out <- boot::boot(data = data, statistic = bootstrap_step, R = nboot,
                    outcome = outcome, exposure = exposure,
                    mediator = mediator, model = model,
                    covariates.pre = covariates.pre, vecc = vecc, EMint = EMint,
                    EMint.terms = EMint.terms,
                    event = event, mreg = mreg, yreg = yreg,
                    m_star = m_star, a_star = a_star, a = a)

  effect_se <- apply(out$t, 2, sd)

  names(effect_se) <- paste0(names(effect_estimates), "_se")

  out <- list(effect = effect_estimates, se = effect_se)

}

################################################################################################################
#################################################g-formula Approach##############################
################################################################################################################

  if (model == "g-formula") {

    out <- simulation_approach(data = data, formulas = formulas, exposure = exposure,
                               mediator = mediator, outcome = outcome,
                               covariates.pre = covariates.pre, covariates.post = covariates.post,
                               mreg = mreg, yreg = yreg, model = model, nsims = nsims)

  }

  return(out)

}



