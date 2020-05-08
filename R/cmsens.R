cmsens <- function(cmest_out = NULL, sens = "uc",
                   MEvariable = NULL, MEvariable.type = NULL, measurement.error = NULL,
                   lambda = c(0.5, 1, 1.5, 2), B = 100) {

  data <- cmest_out$data

  model <- cmest_out$model

  outcome <- cmest_out$variables$outcome

  event <- cmest_out$variables$event

  exposure <- cmest_out$variables$exposure

  mediator <- cmest_out$variables$mediator

  prec <- cmest_out$variables$prec

  postc <- cmest_out$variables$postc

  yreg <- cmest_out$reg$yreg

  mreg <- cmest_out$reg$mreg

  ereg <- cmest_out$reg$ereg

  postcreg <- cmest_out$reg$postcreg

  wmreg <- cmest_out$reg$wmreg

  astar <- cmest_out$ref$astar

  a <- cmest_out$ref$a

  mval <- cmest_out$ref$mval

  yref <- cmest_out$ref$yref

  vecc <- cmest_out$ref$vecc

  estimation <- cmest_out$method["estimation"]

  inference <- cmest_out$method["inference"]

  if(inference == "bootstrap") nboot <- cmest_out$nboot

  regressions <- cmest_out$regressions

  outcome_regression <- regressions$outcome_regression

  effect_estimate <- cmest_out$effect_estimate

  effect_se <- cmest_out$effect_se

  if (!sens %in% c("uc", "me")) {
    stop("Only supports sensitivity analysis for unmeasured confounding and measurement error")
  }


  #################################################################################################
  ############################Sensitiviti Analysis for Unmeasure Confounding#######################
  #################################################################################################

  if (sens == "uc") {

    if (model %in% c("rb", "wb", "msm", "g-formula")) {

      index_biased <- 1:6

    } else if (model == "iorw") {

      index_biased <- 1:3

    } else if (model == "ne") {

      index_biased <- 1:length(effect_estimate)

    }

    est <- c()

    evalue_pe <- c()

    evalue_ci_lo <- c()

    evalue_ci_hi <- c()

    if (identical(class(outcome_regression), "lm") |
        (identical(class(outcome_regression), c("glm", "lm")) &&
         (family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi")))|
        (identical(class(outcome_regression), c("gam", "glm", "lm")) &&
         (family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                   "gaulss", "gevlss")|
          startsWith(family(outcome_regression)$family, "Tweedie")|
          startsWith(family(outcome_regression)$family, "Beta regression")|
          startsWith(family(outcome_regression)$family, "Scaled t")))|
        identical(class(outcome_regression), "nls")) {

      for (i in index_biased) {

        d <- unname(effect_estimate[i] / sd(data[, outcome]))

        sd <- unname(effect_se[i] / sd(data[, outcome]))

        est <- c(est, exp(0.91 * d))

        evalues <- EValue::evalues.RR(est = exp(0.91 * d), lo = exp(0.91 * d - 1.78 * sd),
                                      hi = exp(0.91 * d + 1.78 * sd))

        evalue_pe <- c(evalue_pe, unname(evalues["E-values","point"]))

        evalue_ci_lo <- c(evalue_ci_lo, unname(evalues["E-values","lower"]))

        evalue_ci_hi <- c(evalue_ci_hi, unname(evalues["E-values","upper"]))

      }
    } else {

      for (i in index_biased) {

        est <- c(est, effect_estimate[i])

        evalues <- EValue::evalues.RR(est = effect_estimate[i],
                                      lo = effect_estimate[i] - 1.96 * effect_se[i],
                                      hi = effect_estimate[i] + 1.96 * effect_se[i])

        evalue_pe <- c(evalue_pe, unname(evalues["E-values","point"]))

        evalue_ci_lo <- c(evalue_ci_lo, unname(evalues["E-values","lower"]))

        evalue_ci_hi <- c(evalue_ci_hi, unname(evalues["E-values","upper"]))

      }

    }

    out <- cbind(Point = evalue_pe, CIlower = evalue_ci_lo, CIupper = evalue_ci_hi)

    rownames(out) <- names(effect_estimate)[index_biased]

    class(out) <- "cmsens.uc"

    #################################################################################################
    ############################Sensitiviti Analysis for Unmeasure Confounding#######################
    #################################################################################################

  } else if (sens == "me") {

    out <- list(rep(NA, length(measurement.error)))

    for (i in 1:length(measurement.error)) {

      n <- nrow(data)

      effect_estimate <- est_step(data = data, indices = 1:n, model = model,
                                  outcome = outcome, event = event, exposure = exposure,
                                  mediator = mediator, EMint = EMint,
                                  prec = prec, postc = postc,
                                  yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                                  astar = astar, a = a, mval = mval, yref = yref, vecc = vecc,
                                  estimation = estimation,
                                  ME = TRUE, MEvariable = MEvariable,
                                  MEvariable.type = MEvariable.type,
                                  measurement.error = measurement.error[[i]], lambda = lambda, B = B)

      effect_se <- inf_step(data = data, nboot = nboot, model = model,
                            outcome = outcome, event = event, exposure = exposure,
                            mediator = mediator, EMint = EMint,
                            prec = prec, postc = postc,
                            yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                            astar = astar, a = a, mval = mval, yref = yref, vecc = vecc,
                            estimation = estimation, inference = inference,
                            ME = TRUE, MEvariable = MEvariable,
                            MEvariable.type = MEvariable.type,
                            measurement.error = measurement.error[[i]], lambda = lambda, B = B)

      out[[i]] <- list(data = data,
                       effect_estimate = effect_estimate,
                            effect_se = effect_se)

      class(out[[i]]) <- "cmest"

      out[[i]] <- summary(out[[i]])

    }

    if (MEvariable.type == "continuous") {

      relia_ratio <- 1 - measurement.error/sd(data[, MEvariable])

      names(out) <- paste("Sigma = ", measurement.error, ", Reliability Ratio = ",
                               round(relia_ratio, 4), sep = "")

    } else if (MEvariable.type == "categorical") {

      name_vec <- c()

      for (i in 1:length(measurement.error)) {

        name_vec <- c(name_vec, paste("Misclassification Matrix = matrix(c(",
                                      paste(as.vector(measurement.error[[i]]), sep = "", collapse = ","), "), nrow = ",
                                      length(unique(data[, MEvariable])), ")", sep = ""))
      }

      names(out) <- name_vec

    }

    out <- list(cmsens = out)

    out$cmest <- cmest_out

    out$ME <- list(MEvariable = MEvariable, MEvariable.type = MEvariable.type,
                   measurement.error = measurement.error, lambda = lambda, B = B)

    class(out) <- "cmsens.me"

  }

  return(out)

}
