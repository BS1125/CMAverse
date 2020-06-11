multimp <- function(args_mice = list(m = 5), data = NULL, model = "rb",
                    outcome = NULL, event = NULL,
                    exposure = NULL, mediator = NULL, EMint = FALSE,
                    prec = NULL, postc = NULL,
                    yreg = "linear", mreg = "linear", ereg = NULL, postcreg = NULL, wmreg = NULL,
                    astar = 0, a = 1, mval = NULL, yref = NULL, precval = NULL,
                    estimation = "paramfunc", inference = "delta",
                    nboot = 200, nrep = 5) {

  m <- args_mice$m

  if (inference == "delta") {

    data_imp <- complete(mice(data, m = m), action = "all")

    est_imp <- do.call(rbind, lapply(1:m, function(x)
      cmest(data = data_imp[[x]], model = model,
            outcome = outcome, event = event,
            exposure = exposure, mediator = mediator, EMint = EMint,
            prec = prec, postc = postc,
            yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
            astar = astar, a = a, mval = mval, yref = yref, precval = precval,
            estimation = estimation, inference = inference)$effect_estimate))

    effect_estimate <- colMeans(est_imp)

    var_imp <- do.call(rbind, lapply(1:m, function(x)
      cmest(data = data_imp[[x]], model = model,
            outcome = outcome, event = event,
            exposure = exposure, mediator = mediator, EMint = EMint,
            prec = prec, postc = postc,
            yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
            astar = astar, a = a, mval = mval, yref = yref, precval = precval,
            estimation = estimation, inference = inference)$effect_se^2))

    var_within <- colMeans(var_imp)

    var_between <- colMeans((est_imp - t(matrix(rep(effect_estimate, m), ncol = m)))^2)/(m - 1)

    effect_se <- sqrt(var_within + var_between * (1 + 1/m))

  } else if (inference == "bootstrap") {

    imp_step <- function(data, indices, m, model, outcome, event, exposure, mediator, EMint,
                         prec, postc, yreg, mreg, ereg, postcreg, wmreg, astar, a, mval, yref,
                         precval, estimation, inference) {

      data_boot <- data[indices, ]

      data_imp <- complete(mice(data_boot, m = m), action = "all")

      out <- colMeans(do.call(rbind, lapply(1:m, function(x)
        cmest(data = data_imp[[x]], model = model,
              outcome = outcome, event = event,
              exposure = exposure, mediator = mediator, EMint = EMint,
              prec = prec, postc = postc,
              yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
              astar = astar, a = a, mval = mval, yref = yref, precval = precval,
              estimation = estimation, inference = inference)$effect_estimate)))

      return(out)

    }

    effect_estimate <- imp_step(data = data, indices = 1:nrow(data), m = m, model = model,
                                outcome = outcome, event = event,
                                exposure = exposure, mediator = mediator, EMint = EMint,
                                prec = prec, postc = postc,
                                yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                                astar = astar, a = a, mval = mval, yref = yref, precval = precval,
                                estimation = estimation, inference = inference)

    boots <- boot::boot(data = data, statistic =imp_step, R = nboot, m = m,
                        model = model,
                        outcome = outcome, event = event,
                        exposure = exposure, mediator = mediator, EMint = EMint,
                        prec = prec, postc = postc,
                        yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                        astar = astar, a = a, mval = mval, yref = yref, precval = precval,
                        estimation = estimation, inference = inference)

    effect_se <- apply(boots$t, 2, sd)


  }

  out <- list(effect_estimate = effect_estimate, effect_se = effect_se)

  class(out) <- c("cmest", "cmest_multimp")

  return(out)

}
