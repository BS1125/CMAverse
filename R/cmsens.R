cmsens <- function(cmest_out = NULL, sens = "uc",
                   MEvariable = NULL, MEvariable.type = NULL, measurement.error = NULL,
                   data = NULL, model = "rb",
                   outcome = NULL, event = NULL,
                   exposure = NULL, mediator = NULL, EMint = FALSE,
                   prec = NULL, postc = NULL,
                   yreg = "linear", mreg = "linear", ereg = NULL, postcreg = NULL, wmreg = NULL,
                   astar = 0, a = 1, mval = NULL, yref = NULL, precval = NULL,
                   estimation = "paramfunc", inference = "delta",
                   nboot = 200, nrep = 5) {


  require(EValue)
  require(simex)

  if (!sens %in% c("uc", "me")) {
    stop("Only supports sensitivity analysis for unmeasured confounding and measurement error")
  }

  if (sens == "uc") {


    if (model %in% c("rb", "wb", "msm", "g-formula")) {

        index_biased <- 1:6

      } else if (model == "iorw") {

        index_biased <- 1:3

      } else if (model == "ne") {

        index_biased <- 1:length(cmest_out$effect_estimate)

      }

    est <- c()

    evalue_pe <- c()

    evalue_ci <- c()

    if (yreg == "linear"|
        (identical(class(yreg), "lm") |
         (identical(class(yreg), c("glm", "lm")) &&
          family(yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))|
         (identical(class(yreg), c("gam", "glm", "lm")) &&
          (family(yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi")))|
         identical(class(yreg), "nls"))) {

      for (i in index_biased) {

        d <- unname(cmest_out$effect_estimate[i] / sd(data[, outcome]))

        sd <- unname(cmest_out$effect_se[i] / sd(data[, outcome]))

        est <- c(est, exp(0.91 * d))

        evalues <- evalues.RR(est = exp(0.91 * d), lo = exp(0.91 * d - 1.78 * sd),
                                hi = exp(0.91 * d + 1.78 * sd))

        evalue_pe <- c(evalue_pe, unname(evalues["E-values","point"]))

        evalue_ci <- c(evalue_ci, unname(evalues["E-values",c("lower", "upper")][which(!is.na(evalues["E-values",c("lower", "upper")]))]))

        }

    } else {

      for (i in index_biased) {

        est <- c(est, cmest_out$effect_estimate[i])

        evalues <- evalues.RR(est = cmest_out$effect_estimate[i],
                              lo = cmest_out$effect_estimate[i] - 1.96 * cmest_out$effect_se[i],
                              hi = cmest_out$effect_estimate[i] + 1.96 * cmest_out$effect_se[i])

        evalue_pe <- c(evalue_pe, unname(evalues["E-values","point"]))

        evalue_ci <- c(evalue_ci, unname(evalues["E-values",c("lower", "upper")][which(!is.na(evalues["E-values",c("lower", "upper")]))]))

      }

    }

    sens_out <- cbind(Point = evalue_pe, CI = evalue_ci)

    rownames(sens_out) <- names(cmest_out$effect_estimate)[index_biased]

  } else if (sens == "me") {

    formulas <- create_formulas(data = data, model = model,
                                  outcome = outcome, event = event,
                                  exposure = exposure, mediator = mediator, EMint = EMint,
                                  prec = prec, postc = postc,
                                  yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)


    if (MEvariable.type == "continuous") {

      regressions_naive <- run_regressions(formulas = formulas, data = data, model = model,
                                     exposure = exposure, mediator = mediator, postc = postc,
                                     yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

       if (model == "rb") {

         if (MEvariable %in% c(exposure, prec)) {

          sens_out <- list(rep(NA, length(measurement.error)))

          for (j in 1:length(measurement.error)) {

            outcome_regression_naive <- regressions_naive$outcome_regression

            outcome_regression_simex <- simex(model = outcome_regression_naive,
                                              SIMEXvariable = MEvariable,
                                              measurement.error = measurement.error[j],
                                              asymptotic = FALSE)

            mediator_regression_simex <- list()

            for (i in 1:length(mediator)) {

              mediator_regression_naive <- regressions_naive$mediator_regression[[i]]

              mediator_regression_simex[[i]] <- simex(model = mediator_regression_naive,
                                              SIMEXvariable = MEvariable,
                                              measurement.error = measurement.error[j],
                                              asymptotic = FALSE)

            }

            reg.simex <- list(outcome_regression = outcome_regression_simex,
                              mediator_regression = mediator_regression_simex)

            sens_out[[j]] <- summary(cmest(data = data, model = model,
                               outcome = outcome, event = event, exposure = exposure,
                               mediator = mediator, EMint = EMint,
                               prec = prec, postc = postc,
                               yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                               wmreg = wmreg, reg.simex = reg.simex,
                               astar = astar, a = a, mval = mval, yref = yref, precval = precval,
                               estimation = estimation, inference = inference,
                               nboot = nboot, nrep = nrep))

          }

         } else if (MEvariable %in% c(mediator)) {

           sens_out <- list(rep(NA, length(measurement.error)))

           for (j in 1:length(measurement.error)) {

             outcome_regression_naive <- regressions_naive$outcome_regression

             outcome_regression_simex <- simex(model = outcome_regression_naive,
                                               SIMEXvariable = MEvariable,
                                               measurement.error = measurement.error[j],
                                               asymptotic = FALSE)

             mediator_regression_simex <- list()

             for (i in 1:length(mediator)) {

               mediator_regression_naive <- regressions_naive$mediator_regression[[i]]

               if (MEvariable %in% names(attr(mediator_regression_naive$terms,"dataClasses"))) {

                 mediator_regression_simex[[i]] <- simex(model = mediator_regression_naive,
                                                       SIMEXvariable = MEvariable,
                                                       measurement.error = measurement.error[j],
                                                       asymptotic = FALSE)

               } else mediator_regression_simex[[i]] <- NULL

             }

             reg.simex <- list(outcome_regression = outcome_regression_simex,
                               mediator_regression = mediator_regression_simex)

             sens_out[[j]] <- summary(cmest(data = data, model = model,
                                            outcome = outcome, event = event, exposure = exposure,
                                            mediator = mediator, EMint = EMint,
                                            prec = prec, postc = postc,
                                            yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                                            wmreg = wmreg, reg.simex = reg.simex,
                                            astar = astar, a = a, mval = mval, yref = yref, precval = precval,
                                            estimation = estimation, inference = inference,
                                            nboot = nboot, nrep = nrep))

           }

         } else if (MEvariable == outcome) {

           sens_out <- list(rep(NA, length(measurement.error)))

           for (j in 1:length(measurement.error)) {

             outcome_regression_naive <- regressions_naive$outcome_regression

             outcome_regression_simex <- simex(model = outcome_regression_naive,
                                               SIMEXvariable = MEvariable,
                                               measurement.error = measurement.error[j],
                                               asymptotic = FALSE)

             reg.simex <- list(outcome_regression = outcome_regression_simex)

             sens_out[[j]] <- summary(cmest(data = data, model = model,
                                            outcome = outcome, event = event, exposure = exposure,
                                            mediator = mediator, EMint = EMint,
                                            prec = prec, postc = postc,
                                            yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                                            wmreg = wmreg, reg.simex = reg.simex,
                                            astar = astar, a = a, mval = mval, yref = yref, precval = precval,
                                            estimation = estimation, inference = inference,
                                            nboot = nboot, nrep = nrep))

           }



         }

         names(sens_out) <- paste("measurement.error = ", measurement.error)

         class(sens_out) <- "cmsens.me"
      }

    }

  }

  return(sens_out)

}
