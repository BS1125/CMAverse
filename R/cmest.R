cmest <- function(data = NULL, model = "rb",
                  outcome = NULL, event = NULL,
                  exposure = NULL, mediator = NULL, EMint = FALSE,
                  prec = NULL, postc = NULL,
                  yreg = "linear", mreg = "linear", ereg = NULL, postcreg = NULL, wmreg = NULL,
                  astar = 0, a = 1, mval = NULL, yref = NULL, precval = NULL,
                  estimation = "paramfunc", inference = "delta",
                  nboot = 200, nrep = 5) {

  require(dplyr)

  ##########################################Argument Restrictions####################################

  if (is.null(data)) {
    stop("Unspecified dataset")
  }

  if (is.null(model)) {
    stop("Unspecified model")
  } else if (!(model %in% c("rb", "wb", "ne", "msm", "iorw", "g-formula"))) {
    stop("Unsupported causal mediation model")
  }

  if (is.null(exposure) | is.null(outcome) | is.null(mediator)) {
    stop("Unspecified exposure, mediator, or outcome")
  }

  if (is.null(event) && ((yreg %in% c("coxph","aft_exp","aft_weibull"))|(inherits(yreg, "coxph"))|
                         (inherits(yreg, "survreg")))) {
    stop("Unspecified event")
  }

  if (!is.null(postc)|!is.null(postcreg)) {

    if (model != "g-formula") {
      stop("Specified model doesn't support post-exposure confounders")
    }

    if (is.null(postc)|is.null(postcreg)) {
      stop("Unspecified postc or postcreg")
    } else if (length(postc) != length(postcreg)) {
      stop("length(postc) != length(postcreg)")
    }
  }

  if (model %in% c("rb", "msm", "g-formula")) {
    if (is.null(mreg)) {
      stop("Unspecified mreg")
    } else if (length(mreg) != length(mediator)) {
      stop("Unspecified mreg")
    } else if (model == "msm" && length(wmreg) != length(mediator)) {
      stop("Unspecified wmreg")
    }
  }

  if (is.null(yreg)) {
    stop("Unspecified yreg")
  }

  if (is.null(ereg)) {
    if (model %in% c("wb", "iorw", "msm")) {
      stop("Unspecified ereg")
    }
  }

  if (is.null(wmreg)) {
    if (model == "msm") {
      stop("Unspecified wmreg")
    }
  }

  if (is.null(a)|is.null(astar)) {
    stop("Unspecified a or astar")
  }

  if (is.null(mval)) {
    if (model %in% c("rb", "wb", "msm", "g-formula")) {
      mval <- as.list(rep(0, length(mediator)))
    }
  } else if ((length(mval) != length(mediator))) {
    stop("length(mval) != length(mediator)")
  }


  if (model == "rb" && estimation == "paramfunc") {

    if (is.null(prec)) {
      vecc <- NULL
    } else {

      vecc <- c()

      if (is.null(precval)) {

        precval <- as.list(rep(NA, length(prec)))

        for (i in 1:length(prec)) {

          if (is.numeric(data[, prec[i]])) {

            vecc <- c(vecc, mean(data[, prec[i]],na.rm=TRUE))

          } else {

            vecc <- c(vecc, unname(colMeans(as.matrix(model.matrix(as.formula(paste0("~", prec[i])),
                                                                   data = data)[,-1]),na.rm = TRUE)))
          }
        }
      } else {

        if (length(precval) != length(prec)) {
          stop("length(precval) != length(prec)")
        }

        for (i in 1:length(prec)) {

          if (is.numeric(data[, prec[i]])) {

            vecc <- c(vecc, precval[[i]])

          } else {

            vecc <- c(vecc, as.numeric(which(levels(as.factor(data[, prec[i]])) == precval[[i]]))[-1])

          }
        }
      }
    }
  } else vecc <- NULL


  ######################################Effect estimation and inference###########################


  if (model != "ne") {

    n <- nrow(data)

    formulas <- create_formulas(model = model,
                                outcome = outcome, event = event,
                                exposure = exposure, mediator = mediator, EMint = EMint,
                                prec = prec, postc = postc,
                                yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

    regressions <- run_regressions(formulas = formulas, data = data, model = model,
                                   exposure = exposure, mediator = mediator, postc = postc,
                                   yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                                   wmreg = wmreg)

    effect_estimate <- est_step(data = data, indices = 1:n, model = model,
                                outcome = outcome, event = event, exposure = exposure,
                                mediator = mediator, EMint = EMint,
                                prec = prec, postc = postc,
                                yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                                astar = astar, a = a, mval = mval, yref = yref, vecc = vecc,
                                estimation = estimation)

    effect_se <- inf_step(data = data, nboot = nboot, model = model,
                          outcome = outcome, event = event, exposure = exposure,
                          mediator = mediator, EMint = EMint,
                          prec = prec, postc = postc,
                          yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                          astar = astar, a = a, mval = mval, yref = yref, vecc = vecc,
                          estimation = estimation, inference = inference)

    out <- list(effect_estimate = effect_estimate,
                effect_se = effect_se)

  } else if (model == "ne") {

    formulas <- create_formulas(model = model,
                                outcome = outcome, event = event,
                                exposure = exposure, mediator = mediator, EMint = EMint,
                                prec = prec, postc = postc,
                                yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

    regressions <- run_regressions(formulas = formulas, data = data, model = model,
                                   exposure = exposure, mediator = mediator, postc = postc,
                                   yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

    outcome_regression <- regressions$outcome_regression

    expData <- medflex::neImpute(outcome_regression, data = data,
                                 nMed = length(mediator), nrep = nrep)

    if (EMint == TRUE) {

      medflex_formula <- paste(c(paste(outcome, paste0(exposure, "0"), sep = "~"), paste0(exposure, "1"),
                                 paste(paste0(exposure, "0"), paste0(exposure, "1"), sep = "*"), prec),
                               collapse = "+")

    } else if (EMint == FALSE) {

      medflex_formula <- paste(c(paste(outcome, paste0(exposure, "0"), sep = "~"), paste0(exposure, "1"),
                                 prec),
                               collapse = "+")

    }

    assign("medflex_formula", medflex_formula, envir = globalenv())

    assign("outcome_formula", formula(outcome_regression), envir = globalenv())

    family <- outcome_regression$family

    se <- ifelse(inference == "delta", "robust", "bootstrap")

    results <- summary(medflex::neEffdecomp(medflex::neModel(formula = as.formula(medflex_formula),
                                                             family = family,
                                                             expData = expData, se = se),
                                            xRef = c(astar, a)))$coefficients[, c("Estimate", "Std. Error")]

    if (!((inherits(regressions$outcome_regression, "glm")|
           inherits(regressions$outcome_regression, "lm"))&&
          family(regressions$outcome_regression)$family == "gaussian")) {

      for (i in 1:length(results[, "Std. Error"]))
        results[, "Std. Error"][i] <- msm::deltamethod(as.formula("~exp(x1)"), results[, "Estimate"][i],
                                                       (results[, "Std. Error"][i])^2)

      results[, "Estimate"] <- exp(results[, "Estimate"])

    }

    rm(medflex_formula, pos = ".GlobalEnv")
    rm(outcome_formula, pos = ".GlobalEnv")

    effect_estimate <- results[, "Estimate"]
    effect_se <- results[, "Std. Error"]

    out <- list(effect_estimate = effect_estimate,
                effect_se = effect_se)

  }

  out$data <- data

  out$model <- model

  out$variables <- list(outcome = outcome, event = event,
                     exposure = exposure, mediator = mediator,
                     prec = prec, postc = postc)

  out$reg <- list(yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

  out$ref <- list(astar = astar, a = a, mval = mval, yref = yref, vecc = vecc)

  out$method <- c(estimation = estimation, inference = inference)

  if(inference == "bootstrap") out$nboot <- nboot

  out$regressions <- regressions

  class(out) <- "cmest"

  return(out)

}


