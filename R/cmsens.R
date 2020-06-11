cmsens <- function(cmest_out = NULL, sens = "uc", MEmethod = "simex",
                   MEvariable = NULL, MEvariable.type = NULL, measurement.error = NULL,
                   lambda = c(0.5, 1, 1.5, 2), B = 100) {

  data <- cmest_out$data

  model <- cmest_out$model

  outcome <- cmest_out$variables$outcome

  event <- cmest_out$variables$event

  exposure <- cmest_out$variables$exposure

  mediator <- cmest_out$variables$mediator

  EMint <- cmest_out$reg$EMint

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
  ############################Sensitivity Analysis for Unmeasured Confounding#######################
  #################################################################################################

  if (sens == "uc") {

    if (model %in% c("rb", "wb", "msm", "g-formula")) {

      out_index <- 1:6

    } else if (model == "iorw") {

      out_index <- 1:3

    } else if (model == "ne") {

      out_index <- 1:length(effect_estimate)

    }

    effect_estimate <- effect_estimate[out_index]

    effect_se <- effect_se[out_index]

    out <- c()

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

      outcome_se <- sd(data[, outcome], na.rm = TRUE)

      d <- unname(effect_estimate / outcome_se)

      sd <- unname(effect_se / outcome_se)

      for (i in out_index) {

        evalues <- EValue::evalues.RR(est = exp(0.91 * d[i]), lo = exp(0.91 * d[i] - 1.78 * sd[i]),
                                      hi = exp(0.91 * d[i] + 1.78 * sd[i]))

        evalues <- c(evalues[1, ], evalues[2, ])

        names(evalues) <- c("estRR", "lowerRR", "upperRR", "Evalue.estRR", "Evalue.lowerRR", "Evalue.upperRR")

        out <- rbind(out, evalues)

      }

    } else {

      summ <- summary(cmest_out)

      ci_lo <- summ[, "95% CIL"]

      ci_up <-summ[, "95% CIU"]

      for (i in out_index) {

        evalues <- EValue::evalues.RR(est = effect_estimate[i],
                                      lo = ci_lo[i],
                                      hi = ci_up[i])

        evalues <- c(evalues[1, ], evalues[2, ])

        names(evalues) <- c("estRR", "lowerRR", "upperRR", "Evalue.estRR", "Evalue.lowerRR", "Evalue.upperRR")

        out <- rbind(out, evalues)

      }

    }

    class(out) <- "cmsens.uc"

    #################################################################################################
    ############################Sensitiviti Analysis for Measurement Error###########################
    #################################################################################################

  } else if (sens == "me") {

    if (MEmethod == "simex") {

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

        relia_ratio <- 1 - measurement.error/sd(data[, MEvariable], na.rm = TRUE)

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

    } else if (MEmethod == "rc") {

      if (MEvariable.type == "continuous" && MEvariable != outcome &&
          ((model == "rb"&& !MEvariable %in% c(mediator))|(model == "wb"&& MEvariable != exposure)|
           (model == "msm" && !MEvariable %in% c(exposure, mediator))|
           (model == "g-formula" && !MEvariable %in% c(mediator, postc))|
           (model == "iorw" && MEvariable != exposure))) {

        n <- nrow(data)

        for (i in 1:length(measurement.error)) {



          if (model == "rb") {

            if (MEvariable %in% c(exposure, mediator, prec)) {

              MEvariable.index <- which(c(exposure, mediator, prec) == MEvariable)

              data_rc <- model.matrix(as.formula(paste0("~", paste(c(MEvariable,
                                                                     c(exposure, mediator, prec)[-MEvariable.index]), collapse = "+"))),
                                      data = data[, c(exposure, mediator, prec)])[, -1]

              mean <- colMeans(data_rc)

              W <- data_rc - rbind(mean)[rep(1, n),]

              covmat <- cov(data_rc[, MEvariable], data_rc)

              covmat[1] <- covmat[1] - measurement.error[i]^2

              EX <- rep(mean[1], n) + W %*% t(covmat %*% solve(cov(data_rc)))
            }



          }


        }





      } else {
        stop("Regression calibration only supports a continuous independent variable measured with error")
      }



    } else {
      stop("MEmethod !%in% c('simex', 'rc')")
    }




  }

  return(out)

}
