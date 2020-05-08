est_step <- function(data, indices, model, regressions,
                     outcome, event, exposure, mediator, EMint, prec, postc,
                     astar, a, mval, yref, vecc,
                     estimation,
                     ME = FALSE, MEvariable = NULL, MEvariable.type = NULL,
                     measurement.error = NULL, lambda = c(0.5, 1, 1.5, 2), B = 100) {

  data_boot <- data[indices, ]

  outcome_regression <- update(regressions$outcome_regression, data = data_boot)

  if ((inherits(outcome_regression, "gam") &&
       (family(outcome_regression)$family == "multinom" |
        startsWith(family(outcome_regression)$family, "Ordered Categorical")))|
      (inherits(outcome_regression, "multinom") && length(outcome_regression$lev) > 2)|
      inherits(outcome_regression, "polr")) {

    if (is.null(yref)) {
      warnings("Unspecified yref, 1 is used")
      yref = 1
    } else {
      yref = which(levels(as.factor(data_boot[, outcome])) == yref)
    }

  }

  ###################################################################################################
  #################################Closed-form Parameter Function Estimation#########################
  ###################################################################################################

  if (estimation == "paramfunc") {

    if (!(model %in% c("rb", "iorw"))) {
      stop("Closed-form parameter function estimation doesn't support the selected model")
    } else if (length(mediator) > 1&&model == "rb") {
      stop("For the selected model, closed-form parameter function estimation only supports a single mediator")
    }

    if (!(yreg %in% c("linear", "logistic", "loglinear", "poisson",
                      "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull"))|
        (!(mreg[[1]] %in% c("linear", "logistic", "multinominal")) && model == "rb")) {

      stop("Closed-form parameter function estimation only
           supports yreg in {linear, logistic, loglinear, poisson, quasipoisson,
           negbin, coxph, aft_exp, aft_weibull} and mreg in {linear, logistic, multinominal}")
    }

    if (is.factor(data_boot[, exposure]) | is.character(data_boot[, exposure])) {
      a <- as.numeric(levels(as.factor(data_boot[, exposure])) == a)[-1]
      astar <- as.numeric(levels(as.factor(data_boot[, exposure])) == astar)[-1]
    }

    if (model == "rb"){

      if (is.factor(data_boot[, mediator])|is.character(data_boot[, mediator])) {
        mstar <- as.numeric(levels(as.factor(data_boot[, mediator])) == mval[[1]])[-1]
      } else {mstar <- mval[[1]]}

      elevel <- ifelse(is.character(data_boot[, exposure])|is.factor(data_boot[, exposure]),
                       length(levels(data_boot[, exposure])), 2)

      mlevel <- ifelse(is.character(data_boot[, mediator])|is.factor(data_boot[, mediator]),
                       length(levels(data_boot[, mediator])), 2)

      mediator_regression <- update(regressions$mediator_regression[[1]], data = data_boot)

      if (ME && MEvariable %in% all.vars(formula(outcome_regression))) {

        outcome_SIMEXreg <- simex_reg(reg = outcome_regression, data = data_boot,
                                      MEvariable = MEvariable, MEvariable.type = MEvariable.type,
                                      measurement.error = measurement.error, lambda = lambda,
                                      B = B)

        thetas <- outcome_SIMEXreg$SIMEXcoef

      } else thetas <- coef(outcome_regression)

      if (ME && MEvariable %in% all.vars(formula(mediator_regression))) {

        mediator_SIMEXreg <- simex_reg(reg = mediator_regression, data = data_boot,
                                       MEvariable = MEvariable, MEvariable.type = MEvariable.type,
                                       measurement.error = measurement.error, lambda = lambda,
                                       B = B)

        betas <- mediator_SIMEXreg$SIMEXcoef

        if (mreg == "linear") { variance <- mediator_SIMEXreg$SIMEXsigma^2
        } else variance <- NULL

      } else {

        betas  <- as.vector(t(coef(mediator_regression)))

        if (mreg == "linear") { variance <- sigma(mediator_regression)^2
        } else variance <- NULL

      }

      theta0 <- thetas[1]

      theta1 <- thetas[2:elevel]

      theta2 <- thetas[(elevel+1):(elevel + mlevel - 1)]

      if (EMint == TRUE) {
        theta3 <- t(matrix(thetas[length(thetas)-(((elevel-1)*(mlevel-1)-1):0)],
                           ncol = mlevel - 1))
      } else {theta3 <- t(matrix(rep(0, (elevel-1)*(mlevel-1)),
                                 ncol = mlevel - 1))}

      beta0 <- betas[1+(0:(mlevel-2))*length(betas)/(mlevel-1)]

      beta1 <- t(matrix(betas[rowSums(expand.grid(2:elevel,(0:(mlevel-2))*
                                                    length(betas)/(mlevel-1)))],
                        ncol = mlevel - 1))

      covariatesTerm <- sapply(0:(mlevel-2), function(x) sum(betas[elevel + 1:length(vecc) + x*length(betas)/(mlevel-1)]*vecc))

      if (yreg == "linear") {

        if (mreg == "linear") {

          cde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                          (sum(theta3 * a) - sum(theta3 * astar)) * mstar)

          pnde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                           (sum(theta3 * a) - sum(theta3 * astar)) *
                           (beta0 + sum(beta1 * astar) + covariatesTerm))

          tnde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                           (sum(theta3 * a) - sum(theta3 * astar)) *
                           (beta0 + sum(beta1 * a) + covariatesTerm))

          pnie <- unname((theta2 + sum(theta3  * astar)) *
                           (sum(beta1 * a) - sum(beta1 * astar)))

          tnie <- unname((theta2 + sum(theta3  * a)) *
                           (sum(beta1 * a) - sum(beta1 * astar)))

        } else if (mreg %in% c("logistic", "multinomial")) {

          cde <- unname((sum(theta1 * a) - sum(theta1 * astar) +
                           ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a ) -
                           ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar)))

          pnde <- unname((sum(theta1 * a) - sum(theta1 * astar) +
                            (sum((theta3 %*% a - theta3 %*% astar) *
                                   exp(beta0 + beta1 %*% astar + covariatesTerm)) /
                               (1 + sum(exp(beta0 + beta1 %*% astar +
                                              covariatesTerm))))))

          tnde <- unname((sum(theta1 * a) - sum(theta1 * astar) +
                            (sum((theta3 %*% a - theta3 %*% astar) *
                                   exp(beta0 + beta1 %*% a + covariatesTerm)) /
                               (1 + sum(exp(beta0 + beta1 %*% a +
                                              covariatesTerm))))))

          pnie <- unname(sum((theta2+theta3 %*% astar)*exp(beta0 + beta1 %*% a + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) -
                           sum((theta2+theta3 %*% astar)*exp(beta0 + beta1 %*% astar + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))))

          tnie <- unname(sum((theta2+theta3 %*% a)*exp(beta0 + beta1 %*% a + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) -
                           sum((theta2+theta3 %*% a)*exp(beta0 + beta1 %*% astar + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))))

        }

        te <- pnde + tnie

        pm <- tnie / (pnde + te)

        intref <- pnde - cde

        intmed <- tnie - pnie

        pie <- pnie

        cde_prop <- cde/te

        intref_prop <- intref/te

        intmed_prop <- intmed/te

        pie_prop <- pie/te

        overall_pm <- (pie+intmed)/te

        overall_int <- (intref+intmed)/te

        overall_pe <- (intref+intmed+pie)/te

        out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie, te = te, pm = pm,
                 intref = intref, intmed = intmed, pie = pie,
                 cde_prop = cde_prop, intref_prop = intref_prop,
                 intmed_prop = intmed_prop, pie_prop = pie_prop,
                 overall_pm = overall_pm, overall_int = overall_int, overall_pe = overall_pe)

      } else if (!yreg == "linear") {

        if (mreg == "linear") {

          cde_rr <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                 (sum(theta3 * a) - sum(theta3 * astar)) * mstar))

          pnde_rr <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                  (sum(theta3 * a) - sum(theta3 * astar)) *
                                  (beta0 + sum(beta1 * astar) +
                                     covariatesTerm + theta2  * variance) +
                                  0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar))))

          tnde_rr <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                  (sum(theta3 * a) - sum(theta3 * astar)) *
                                  (beta0 + sum(beta1 * a) +
                                     covariatesTerm + theta2  * variance) +
                                  0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar))))

          pnie_rr <- unname(exp(theta2 * (sum(beta1 * a) - sum(beta1 * astar)) +
                                  sum(theta3  * astar) * (sum(beta1 * a) - sum(beta1 * astar))))

          tnie_rr <- unname(exp(theta2 * (sum(beta1 * a) - sum(beta1 * astar)) +
                                  sum(theta3  * a) * (sum(beta1 * a) - sum(beta1 * astar))))

          cde_err <- unname((exp(sum(theta1 * a) - sum(theta1 * astar) +
                                   sum(theta3  * a) * mstar) -
                               exp(sum(theta3 * astar) * mstar)) *
                              exp(theta2 * mstar - (theta2 + sum(theta3 * astar)) *
                                    (beta0 + sum(beta1 * astar) + covariatesTerm) -
                                    0.5 * (theta2 + sum(theta3 * astar)) ^ 2 * variance))

        } else if (mreg %in% c("logistic", "multinomial")) {

          cde_rr <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                 ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a) -
                                 ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar)))


          pnde_rr <- unname((exp(sum(theta1 * a) - sum(theta1 * astar)) *
                               (1 + sum(exp(theta2 + theta3 %*% a +
                                              beta0 + beta1 %*% astar + covariatesTerm)))) /
                              (1 + sum(exp(theta2 + theta3 %*% astar +
                                             beta0 + beta1 %*% astar + covariatesTerm))))

          tnde_rr <- unname((exp(sum(theta1 * a) - sum(theta1 * astar)) *
                               (1 + sum(exp(theta2 + theta3 %*% a +
                                              beta0 + beta1 %*% a + covariatesTerm)))) /
                              (1 + sum(exp(theta2 + theta3 %*% astar +
                                             beta0 + beta1 %*% a + covariatesTerm))))

          pnie_rr <- unname(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                               (1 + sum(exp(theta2 + theta3 %*% astar +
                                              beta0 + beta1 %*% a + covariatesTerm)))) /
                              ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                                 (1 + sum(exp(theta2 + theta3 %*% astar +
                                                beta0 + beta1 %*% astar + covariatesTerm)))))

          tnie_rr <- unname(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                               (1 + sum(exp(theta2 + theta3 %*% a +
                                              beta0 + beta1 %*% a + covariatesTerm)))) /
                              ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                                 (1 + sum(exp(theta2 + theta3 %*% a +
                                                beta0 + beta1 %*% astar + covariatesTerm)))))

          cde_err <- unname(exp(sum(theta2*mstar)) *
                              (exp(sum(theta1 * a) - sum(theta1 * astar) +
                                     ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a)) -
                                 exp(ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar))) *
                              (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) /
                              (1+ sum(exp(theta2 + theta3 %*% astar + beta0 +
                                            beta1 %*% astar + covariatesTerm))))

        }

        intref_err <- pnde_rr - 1 - cde_err

        intmed_err <- tnie_rr * pnde_rr - pnde_rr - pnie_rr + 1

        pie_err <- pnie_rr - 1

        te_rr <- tnie_rr * pnde_rr

        pm <- (pnde_rr * (tnie_rr - 1)) / (te_rr - 1)

        te_err <- te_rr - 1

        cde_err_prop <- cde_err/te_err

        intmed_err_prop <- intmed_err/te_err

        intref_err_prop <- intref_err/te_err

        pie_err_prop <- pie_err/te_err

        overall_pm <- (pie_err+intmed_err)/te_err

        overall_int <- (intref_err+intmed_err)/te_err

        overall_pe <- (intref_err+intmed_err+pie_err)/te_err

        out <- c(cde_rr = cde_rr, pnde_rr = pnde_rr, tnde_rr = tnde_rr, pnie_rr = pnie_rr,
                 tnie_rr = tnie_rr, te_rr = te_rr, pm = pm,
                 cde_err = cde_err, intref_err = intref_err,
                 intmed_err = intmed_err, pie_err = pie_err,
                 cde_err_prop = cde_err_prop, intref_err_prop = intref_err_prop,
                 intmed_err_prop = intmed_err_prop, pie_err_prop = pie_err_prop,
                 overall_pm = overall_pm, overall_int = overall_int,
                 overall_pe = overall_pe)


      }

    } else if (model == "iorw") {

      tot_regression <- update(outcome_regression, data = data_boot)

      dir_regression <- update(outcome_regression, data = data_boot)

      exposure_regression <- update(regressions$exposure_regression, data = data_boot)

      exposure_regression_pred <- exposure_regression

      if (ME && MEvariable %in% all.vars(formula(exposure_regression))) {

        exposure_regression_pred <- simex_reg(reg = exposure_regression_pred,
                                              data = data_boot, MEvariable = MEvariable,
                                              MEvariable.type = MEvariable.type,
                                              measurement.error = measurement.error,
                                              lambda = lambda, B = B)

      }

      if ((inherits(exposure_regression, "glm") &&
           family(exposure_regression)$family %in% c("binomial", "quasibinomial"))|
          (identical(class(exposure_regression), c("multinom", "nnet"))&&
           length(exposure_regression$lev)==2)) {

        category <- as.numeric(as.factor(data_boot[, exposure])) - 1

        wadenom.prob <- predict(exposure_regression_pred, newdata = data_boot,
                                type = ifelse(inherits(exposure_regression, "multinom"),
                                              "probs", "response"))

        wadenom <-  wadenom.prob ^ category *
          (1 - wadenom.prob)^(1 - category)

        wanom <- 1 - wadenom.prob

      } else {

        data_pred <- data_boot[, exposure, drop = FALSE]

        data_pred[, exposure] <- as.factor(data_pred[, exposure])

        category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = data_pred)

        wadenom.prob <- predict(exposure_regression_pred, newdata = data_boot,
                                type = ifelse(inherits(exposure_regression, "multinom")|
                                                inherits(exposure_regression, "polr"),
                                              "probs", "response"))

        wadenom <- rowSums(category * wadenom.prob)

        wanom <- wadenom.prob[, 1]

      }

      wa <- wanom/wadenom

      dir_regression$call[["weights"]] <- wa

      dir_regression <- update(dir_regression)


      if (ME && MEvariable %in% all.vars(formula(outcome_regression))) {

        dir_coef <- simex_reg(reg = dir_regression,
                              data = data_boot, MEvariable = MEvariable,
                              MEvariable.type = MEvariable.type,
                              measurement.error = measurement.error,
                              lambda = lambda, B = B)$SIMEXcoef

        tot_coef <- simex_reg(reg = tot_regression,
                              data = data_boot, MEvariable = MEvariable,
                              MEvariable.type = MEvariable.type,
                              measurement.error = measurement.error,
                              lambda = lambda, B = B)$SIMEXcoef

      } else {

        dir_coef <- coef(dir_regression)

        tot_coef <- coef(tot_regression)

      }


      if (is.factor(data_boot[, exposure]) | is.character(data_boot[, exposure])) {

        dir_coef <- unname(dir_coef[2+(0:(length(unique(data[, exposure]))-2))])

        tot_coef <- unname(tot_coef[2+(0:(length(unique(data[, exposure]))-2))])

      } else {

        dir_coef <- unname(coef(dir_regression)[2])

        tot_coef <- unname(coef(tot_regression)[2])

      }

      if (yreg == "linear") {

        dir <- unname(dir_coef %*% a - dir_coef %*% astar)

        dir <- unname(tot_coef %*% a - tot_coef %*% astar)

        ind <- tot - dir

        out <- c(tot = tot, dir = dir, ind = ind)

      } else {

        RRdir <- unname(exp(dir_coef %*% a - dir_coef %*% astar))

        RRtot <- unname((tot_coef %*% a - tot_coef %*% astar))

        RRind <- RRtot / RRdir

        out <- c(RRtot = RRtot, RRdir = RRdir, RRind = RRind)

      }
    }

    ###################################################################################################
    ###############################Direct Counterfactual Imputation Estimation#########################
    ###################################################################################################

  } else if (estimation == "imputation") {

    if (model %in% c("rb", "msm", "g-formula")) {

      if (model == "msm") {

        # calculate P(A=ai)/P(A=ai|C=ci)

        exposure_regression <- update(regressions$exposure_regression, data = data_boot)

        exposure_regression_pred <- exposure_regression

        if (ME && MEvariable %in% all.vars(formula(exposure_regression))) {

          exposure_regression_pred <- simex_reg(reg = exposure_regression,
                                                data = data_boot, MEvariable = MEvariable,
                                                MEvariable.type = MEvariable.type,
                                                measurement.error = measurement.error,
                                                lambda = lambda, B = B)

        }

        wanom <- left_join(select(data_boot, exposure),
                           count(data_boot, !!as.name(exposure)),
                           by = exposure)[, "n"]/nrow(data_boot)

        if ((inherits(exposure_regression, "glm") &&
             family(exposure_regression)$family %in% c("binomial", "quasibinomial"))|
            (inherits(exposure_regression, "multinom")&&
             length(exposure_regression$lev)==2)) {

          category <- as.numeric(as.factor(data_boot[, exposure])) - 1

          wadenom.prob <- predict(exposure_regression, newdata = data_boot,
                                  type = ifelse(inherits(exposure_regression, "multinom"),
                                                "probs", "response"))

          wadenom <-  wadenom.prob ^ category *
            (1 - wadenom.prob)^(1 - category)

        } else {

          data_pred <- data_boot[, exposure, drop = FALSE]

          data_pred[, exposure] <- as.factor(data_pred[, exposure])

          category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = data_pred)

          wadenom.prob <- predict(exposure_regression, newdata = data_boot,
                                  type = ifelse(inherits(exposure_regression, "multinom")|
                                                  inherits(exposure_regression, "polr"),
                                                "probs", "response"))

          wadenom <- rowSums(category * wadenom.prob)

        }

        wa <- wanom/wadenom

        # calculate P(Mp=Mp,i|A=ai)/P(Mp=Mp,i|A=ai,C=ci,L=Li)

        wmnom_regression <- wmnom_regression_pred <- wmdenom_regression <-
          wmdenom_regression_pred <- list()

         for (i in 1:length(mediator)) {

            wmnom_regression[[i]] <- update(regressions$wmnom_regression[[i]], data = data_boot)

            if (ME && MEvariable %in% all.vars(formula(wmnom_regression[[i]]))) {

              wmnom_regression_pred[[i]] <- simex_reg(reg = wmnom_regression[[i]],
                                                      data = data_boot, MEvariable = MEvariable,
                                                      MEvariable.type = MEvariable.type,
                                                      measurement.error = measurement.error,
                                                      lambda = lambda, B = B)

            } else wmnom_regression_pred[[i]] <- wmnom_regression[[i]]

            wmdenom_regression[[i]] <- update(regressions$wmdenom_regression[[i]], data = data_boot)

            if (ME && MEvariable %in% all.vars(formula(wmdenom_regression[[i]]))) {

              wmdenom_regression_pred[[i]] <- simex_reg(reg = wmdenom_regression[[i]],
                                                        data = data_boot, MEvariable = MEvariable,
                                                        MEvariable.type = MEvariable.type,
                                                        measurement.error = measurement.error,
                                                        lambda = lambda, B = B)

            } else wmdenom_regression_pred[[i]] <- wmdenom_regression[[i]]
          }


        wmnom <- rep(1, nrow(data_boot))

        wmdenom <- rep(1, nrow(data_boot))

        for (i in 1:length(wmdenom_regression)) {

          if ((inherits(wmdenom_regression[[i]], "glm") &&
               family(wmdenom_regression[[i]])$family %in% c("binomial", "quasibinomial"))|
              (identical(class(wmdenom_regression[[i]]), c("multinom", "nnet"))&&
               length(wmdenom_regression[[i]]$lev) == 2)) {

            category <- as.numeric(as.factor(data_boot[, mediator[i]])) - 1

            wmdenom.prob <- predict(wmdenom_regression_pred[[i]], newdata = data_boot,
                                    type = ifelse(inherits(wmdenom_regression[[i]], "multinom"),
                                                  "probs", "response"))

            wmnom.prob <- predict(wmnom_regression_pred[[i]], newdata = data_boot,
                                  type = ifelse(inherits(wmdenom_regression[[i]], "multinom"),
                                                "probs", "response"))

            wmdenom <-  wmdenom * wmdenom.prob ^ category *
              (1 - wmdenom.prob)^(1 - category)

            wmnom <-  wmnom * wmnom.prob ^ category *
              (1 - wmnom.prob)^(1 - category)

          } else {

            data_pred <- data_boot[, mediator[i], drop = FALSE]

            data_pred[, mediator[i]] <- as.factor(data_pred[, mediator[i]])

            category <- model.matrix(as.formula(paste("~0+", mediator[i], sep = "")), data = data_pred)

            wmdenom.prob <- predict(wmdenom_regression_pred[[i]], newdata = data_boot,
                                    type = ifelse(inherits(wmdenom_regression[[i]], "multinom")|
                                                    inherits(wmdenom_regression[[i]], "polr"),
                                                  "probs","response"))

            wmnom.prob <- predict(wmnom_regression_pred[[i]], newdata = data_boot,
                                  type = ifelse(inherits(wmnom_regression[[i]], "multinom")|
                                                  inherits(wmnom_regression[[i]], "polr"),
                                                "probs","response"))

            wmdenom <- wmdenom * rowSums(category * wmdenom.prob)

            wmnom <- wmnom * rowSums(category * wmnom.prob)

          }

        }

        wm <- wmnom / wmdenom

        wy <- wa * wm


        mediator_regression <- regressions$mediator_regression

        for (i in 1:length(mediator_regression)) {

          mediator_regression[[i]] <- update(mediator_regression[[i]],
                                             data = data_boot, weights = wa)

        }

        outcome_regression <- update(outcome_regression, weights = wy)

      }


      if (model == "g-formula" && !is.null(postc)) {

        postcdesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                                    data_boot[,prec])

        postcdesign_astar <- data.frame(c(rep(astar,nrow(data_boot))),
                                        data_boot[,prec])

        if (is.factor(data_boot[, exposure])) {

          postcdesign_a[, exposure] <- factor(postcdesign_a[, exposure], levels = levels(data_boot[, exposure]))

          postcdesign_astar[, exposure] <- factor(postcdesign_astar[, exposure], levels = levels(data_boot[, exposure]))

        }

        postc_regression <- postc_regression_pred <- list()

        for (i in 1:length(postc_regression)) {

          postc_regression[[i]] <- update(regressions$postc_regression[[i]], data = data_boot)

            if (ME && MEvariable %in% all.vars(formula(postc_regression[[i]]))) {

              postc_regression_pred[[i]] <- simex_reg(reg = postc_regression[[i]],
                                                      data = data_boot, MEvariable = MEvariable,
                                                      MEvariable.type = MEvariable.type,
                                                      measurement.error = measurement.error,
                                                      lambda = lambda, B = B)

            } else postc_regression_pred[[i]] <- postc_regression[[i]]

            }


        postc_a <- data.frame(matrix(nrow = nrow(data_boot), ncol = length(postc)))

        postc_astar <- data.frame(matrix(nrow = nrow(data_boot), ncol = length(postc)))

        colnames(postc_a) <- colnames(postc_astar) <-postc

        for (i in 1:length(postc)) {

          postcdesign_a <- cbind(postcdesign_a, postc_a[, i-1, drop = FALSE])

          postcdesign_astar <- cbind(postcdesign_astar, postc_astar[, i-1, drop = FALSE])

          if ((inherits(postc_regression[[i]], "glm") &&
               (family(postc_regression[[i]])$family %in% c("binomial", "quasibinomial")))|
              (identical(class(postc_regression[[i]]), c("multinom", "nnet"))&&
               length(postc_regression[[i]]$lev) == 2)) {

            prob_a <- predict(postc_regression_pred[[i]], newdata = postcdesign_a,
                              type = ifelse(inherits(postc_regression[[i]], "multinom"),
                                            "probs", "response"))

            pred_a <- rbinom(nrow(postcdesign_a), size = 1, prob = prob_a)

            prob_astar <- predict(postc_regression_pred[[i]], newdata = postcdesign_astar,
                                  type = ifelse(inherits(postc_regression[[i]], "multinom"),
                                                "probs", "response"))

            pred_astar <- rbinom(nrow(postcdesign_astar), size = 1, prob = prob_astar)

            if (is.numeric(data_boot[, postc[i]])) {

              mid_a <- as.numeric(levels(as.factor(data_boot[, postc[i]]))[pred_a + 1])

              mid_astar <- as.numeric(levels(as.factor(data_boot[, postc[i]]))[pred_astar + 1])

            } else if (is.factor(data_boot[, postc[i]])) {

              mid_a <- factor(levels(data_boot[, postc[i]])[pred_a + 1],
                              levels = levels(data_boot[, postc[i]]))

              mid_astar <- factor(levels(data_boot[, postc[i]])[pred_astar + 1],
                                  levels = levels(data_boot[, postc[i]]))

            } else if (is.character(data_boot[, postc[i]])) {

              mid_a <- levels(as.factor(data_boot[, postc[i]]))[pred_a + 1]

              mid_astar <- levels(as.factor(data_boot[, postc[i]]))[pred_astar + 1]

            }

          } else if ((identical(class(postc_regression[[i]]), c("gam", "glm", "lm")) &&
                      (family(postc_regression[[i]])$family %in% c("multinom") |
                       startsWith(family(postc_regression[[i]])$family,"Ordered Categorical")))|
                     (identical(class(postc_regression[[i]]), c("multinom", "nnet")) &&
                      length(postc_regression[[i]]$lev) > 2)|
                     identical(class(postc_regression[[i]]), "polr")) {

            prob_a <- predict(postc_regression_pred[[i]], newdata = postcdesign_a,
                              type = ifelse(inherits(postc_regression[[i]], "gam"), "response", "probs"))

            pred_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

            prob_astar <- predict(postc_regression_pred[[i]], newdata = postcdesign_astar,
                                  type = ifelse(inherits(postc_regression[[i]], "gam"), "response", "probs"))

            pred_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

            if (is.numeric(data_boot[, postc[i]])) {

              mid_a <- as.numeric(levels(as.factor(data_boot[, postc[i]]))[pred_a])

              mid_astar <- as.numeric(levels(as.factor(data_boot[, postc[i]]))[pred_astar])

            } else if (is.factor(data_boot[, postc[i]])) {

              mid_a <- factor(levels(data_boot[, postc[i]])[pred_a],
                              levels = levels(data_boot[, postc[i]]))

              mid_astar <- factor(levels(data_boot[, postc[i]])[pred_astar],
                                  levels = levels(data_boot[, postc[i]]))

            } else if (is.character(data_boot[, postc[i]])) {

              mid_a <- levels(as.factor(data_boot[, postc[i]]))[pred_a]

              mid_astar <- levels(as.factor(data_boot[, postc[i]]))[pred_astar]

            }

          } else {

            mid_a <- predict(postc_regression_pred[[i]], newdata = postcdesign_a, type = "response")

            mid_astar <- predict(postc_regression_pred[[i]], newdata = postcdesign_astar, type = "response")

          }

          postc_a[, i] <- mid_a

          postc_astar[, i] <- mid_astar

        }


      } else {

        postc_a <- postc_astar <- data.frame()[1:nrow(data_boot), ]

      }

      mdesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                              data_boot[,prec],
                              postc_a)

      mdesign_astar <- data.frame(c(rep(astar,nrow(data_boot))),
                                  data_boot[,prec],
                                  postc_astar)

      colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, prec, postc)


      if (is.factor(data_boot[, exposure])) {

        mdesign_a[, exposure] <- factor(mdesign_a[, exposure], levels = levels(data_boot[, exposure]))

        mdesign_astar[, exposure] <- factor(mdesign_astar[, exposure], levels = levels(data_boot[, exposure]))

      }

      mediator_regression <- mediator_regression_pred <- list()

      for (i in 1:length(mediator)) {

        mediator_regression[[i]] <- update(regressions$mediator_regression[[i]], data = data_boot)

          if (ME && MEvariable %in% all.vars(formula(mediator_regression[[i]]))) {

            mediator_regression_pred[[i]] <- simex_reg(reg = mediator_regression[[i]],
                                                       data = data_boot, MEvariable = MEvariable,
                                                       MEvariable.type = MEvariable.type,
                                                       measurement.error = measurement.error,
                                                       lambda = lambda, B = B)

          } else mediator_regression_pred[[i]] <- mediator_regression[[i]]

          }


      m_a <- data.frame(matrix(nrow = nrow(data_boot), ncol = length(mediator)))

      m_astar <- data.frame(matrix(nrow = nrow(data_boot), ncol = length(mediator)))

      colnames(m_a) <- colnames(m_astar) <- mediator

      for (i in 1:length(mediator)) {

        mdesign_a <- cbind(mdesign_a, m_a[, i-1, drop = FALSE])

        mdesign_astar <- cbind(mdesign_astar, m_astar[, i-1, drop = FALSE])

        if ((inherits(mediator_regression[[i]], "glm") &&
             (family(mediator_regression[[i]])$family %in% c("binomial", "quasibinomial")))|
            (identical(class(mediator_regression[[i]]), c("multinom", "nnet"))&&
             length(mediator_regression[[i]]$lev) == 2)) {

          prob_a <- predict(mediator_regression_pred[[i]], newdata = mdesign_a,
                            type = ifelse(inherits(mediator_regression[[i]], "multinom"),
                                          "probs", "response"))

          pred_a <- rbinom(nrow(mdesign_a), size = 1, prob = prob_a)

          prob_astar <- predict(mediator_regression_pred[[i]], newdata = mdesign_astar,
                                type = ifelse(inherits(mediator_regression[[i]], "multinom"),
                                              "probs", "response"))

          pred_astar <- rbinom(nrow(mdesign_astar), size = 1, prob = prob_astar)

          if (is.numeric(data_boot[, mediator[i]])) {

            mid_a <- as.numeric(levels(as.factor(data_boot[, mediator[i]]))[pred_a + 1])

            mid_astar <- as.numeric(levels(as.factor(data_boot[, mediator[i]]))[pred_astar + 1])

          } else if (is.factor(data_boot[, mediator[i]])) {

            mid_a <- factor(levels(data_boot[, mediator[i]])[pred_a + 1],
                            levels = levels(data_boot[, mediator[i]]))

            mid_astar <- factor(levels(data_boot[, mediator[i]])[pred_astar + 1],
                                levels = levels(data_boot[, mediator[i]]))

          } else if (is.character(data_boot[, mediator[i]])) {

            mid_a <- levels(as.factor(data_boot[, mediator[i]]))[pred_a + 1]

            mid_astar <- levels(as.factor(data_boot[, mediator[i]]))[pred_astar + 1]

          }

        } else if ((identical(class(mediator_regression[[i]]), c("gam", "glm", "lm")) &&
                    (family(mediator_regression[[i]])$family %in% c("multinom") |
                     startsWith(family(mediator_regression[[i]])$family,"Ordered Categorical")))|
                   (identical(class(mediator_regression[[i]]), c("multinom", "nnet")) &&
                    length(mediator_regression[[i]]$lev) > 2)|
                   identical(class(mediator_regression[[i]]), "polr")) {

          prob_a <- predict(mediator_regression_pred[[i]], newdata = mdesign_a,
                            type = ifelse(inherits(mediator_regression[[i]], "gam"), "response", "probs"))

          pred_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

          prob_astar <- predict(mediator_regression_pred[[i]], newdata = mdesign_astar,
                                type = ifelse(inherits(mediator_regression[[i]], "gam"), "response", "probs"))

          pred_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

          if (is.numeric(data_boot[, mediator[i]])) {

            mid_a <- as.numeric(levels(as.factor(data_boot[, mediator[i]]))[pred_a])

            mid_astar <- as.numeric(levels(as.factor(data_boot[, mediator[i]]))[pred_astar])

          } else if (is.factor(data_boot[, mediator[i]])) {

            mid_a <- factor(levels(data_boot[, mediator[i]])[pred_a],
                            levels = levels(data_boot[, mediator[i]]))

            mid_astar <- factor(levels(data_boot[, mediator[i]])[pred_astar],
                                levels = levels(data_boot[, mediator[i]]))

          } else if (is.character(data_boot[, mediator[i]])) {

            mid_a <- levels(as.factor(data_boot[, mediator[i]]))[pred_a]

            mid_astar <- levels(as.factor(data_boot[, mediator[i]]))[pred_astar]

          }

        } else {

          mid_a <- predict(mediator_regression_pred[[i]], newdata = mdesign_a, type = "response")

          mid_astar <- predict(mediator_regression_pred[[i]], newdata = mdesign_astar, type = "response")

        }

        m_a[, i] <- mid_a

        m_astar[, i] <- mid_astar

      }

      if (model == "g-formula" && !is.null(prec)) {

        m_a <- m_a[sample(1:nrow(data_boot), replace = FALSE), ]

        m_astar <- m_astar[sample(1:nrow(data_boot), replace = FALSE), ]

      }


      outcome_regression_pred <- outcome_regression

      if (ME && MEvariable %in% all.vars(formula(outcome_regression))) {

        outcome_regression_pred <- simex_reg(reg = outcome_regression,
                                             data = data_boot, MEvariable = MEvariable,
                                             MEvariable.type = MEvariable.type,
                                             measurement.error = measurement.error,
                                             lambda = lambda, B = B)


      }

      ydesign0m <- data.frame(rep(astar, nrow(data_boot)), as.data.frame(mval)[rep(1,nrow(data_boot)),],
                              data_boot[,prec], postc_astar)

      ydesign1m <- data.frame(rep(a, nrow(data_boot)), as.data.frame(mval)[rep(1,nrow(data_boot)),],
                              data_boot[,prec], postc_a)

      ydesign00 <- data.frame(rep(astar, nrow(data_boot)), m_astar,
                              data_boot[,prec], postc_astar)

      ydesign01 <- data.frame(rep(astar, nrow(data_boot)), m_a,
                              data_boot[,prec], postc_astar)

      ydesign10 <- data.frame(rep(a, nrow(data_boot)), m_astar,
                              data_boot[,prec], postc_a)

      ydesign11 <- data.frame(rep(a, nrow(data_boot)), m_a,
                              data_boot[,prec], postc_a)

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, prec,
                                                        postc)

      if (is.factor(data_boot[, exposure])) {

        ydesign0m[, exposure] <- factor(ydesign0m[, exposure], levels = levels(data_boot[, exposure]))

        ydesign1m[, exposure] <- factor(ydesign1m[, exposure], levels = levels(data_boot[, exposure]))

        ydesign00[, exposure] <- factor(ydesign00[, exposure], levels = levels(data_boot[, exposure]))

        ydesign10[, exposure] <- factor(ydesign10[, exposure], levels = levels(data_boot[, exposure]))

        ydesign01[, exposure] <- factor(ydesign01[, exposure], levels = levels(data_boot[, exposure]))

        ydesign11[, exposure] <- factor(ydesign11[, exposure], levels = levels(data_boot[, exposure]))

      }

      for (i in 1:length(mediator)) {

        if (is.factor(data_boot[, mediator[i]])) {
          ydesign0m[,1+i] <- factor(ydesign0m[,1+i], levels = levels(data_boot[, mediator[i]]))
          ydesign1m[,1+i] <- factor(ydesign1m[,1+i], levels = levels(data_boot[, mediator[i]]))
        }
      }

      if ((inherits(outcome_regression, "gam") &&
           (family(outcome_regression)$family == "multinom" |
            startsWith(family(outcome_regression)$family, "Ordered Categorical")))|
          (inherits(outcome_regression, "multinom") && length(outcome_regression$lev) > 2)|
          inherits(outcome_regression, "polr")) {

        type = ifelse(inherits(outcome_regression, "gam"), "response", "probs")

        EY0m <- mean(predict(outcome_regression_pred, newdata =  ydesign0m, type = type)[, yref], na.rm = TRUE)

        EY1m <- mean(predict(outcome_regression_pred, newdata =  ydesign1m, type = type)[, yref], na.rm = TRUE)

        EY00 <- mean(predict(outcome_regression_pred, newdata =  ydesign00, type = type)[, yref], na.rm = TRUE)

        EY01 <- mean(predict(outcome_regression_pred, newdata =  ydesign01, type = type)[, yref], na.rm = TRUE)

        EY10 <- mean(predict(outcome_regression_pred, newdata =  ydesign10, type = type)[, yref], na.rm = TRUE)

        EY11 <- mean(predict(outcome_regression_pred, newdata =  ydesign11, type = type)[, yref], na.rm = TRUE)

      } else {

        type <- ifelse(inherits(outcome_regression, "coxph"), "risk",
                       ifelse(inherits(outcome_regression, "multinom"), "probs", "response"))

        EY0m <- mean(predict(outcome_regression_pred, newdata =  ydesign0m, type = type), na.rm = TRUE)

        EY1m <- mean(predict(outcome_regression_pred, newdata =  ydesign1m, type = type), na.rm = TRUE)

        EY00 <- mean(predict(outcome_regression_pred, newdata =  ydesign00, type = type), na.rm = TRUE)

        EY01 <- mean(predict(outcome_regression_pred, newdata =  ydesign01, type = type), na.rm = TRUE)

        EY10 <- mean(predict(outcome_regression_pred, newdata =  ydesign10, type = type), na.rm = TRUE)

        EY11 <- mean(predict(outcome_regression_pred, newdata =  ydesign11, type = type), na.rm = TRUE)

      }

    } else if (model == "wb") {

      exposure_regression <- update(regressions$exposure_regression, data = data_boot)

      exposure_regression_pred <- exposure_regression

      if (ME && MEvariable %in% all.vars(formula(exposure_regression))) {

        exposure_regression_pred <- simex_reg(reg = exposure_regression_pred,
                                              data = data_boot, MEvariable = MEvariable,
                                              MEvariable.type = MEvariable.type,
                                              measurement.error = measurement.error,
                                              lambda = lambda, B = B)

      }

      outcome_regression_pred <- outcome_regression

      if (ME && MEvariable %in% all.vars(formula(outcome_regression))) {

        outcome_regression_pred <- simex_reg(reg = outcome_regression_pred,
                                             data = data_boot, MEvariable = MEvariable,
                                             MEvariable.type = MEvariable.type,
                                             measurement.error = measurement.error,
                                             lambda = lambda, B = B)

      }


      index_astar <- which(data_boot[, exposure] == astar)

      index_a <- which(data_boot[, exposure] == a)

      astar_lev <- which(levels(as.factor(data_boot[, exposure])) == astar)

      a_lev <- which(levels(as.factor(data_boot[, exposure])) == a)

      wnom <- left_join(select(data_boot, exposure),
                        count(data_boot, !!as.name(exposure)),
                        by = exposure)[, "n"]/nrow(data_boot)

      if ((inherits(exposure_regression, "glm") &&
           family(exposure_regression)$family %in% c("binomial", "quasibinomial"))|
          (inherits(exposure_regression, "multinom")&&
           length(exposure_regression$lev) == 2)) {

        category <- as.numeric(as.factor(data_boot[, exposure])) - 1

        wdenom.prob <- predict(exposure_regression, newdata = data_boot,
                               type = ifelse(inherits(exposure_regression, "multinom"),
                                             "probs", "response"))

        wdenom0 <- wdenom.prob^(astar_lev - 1)*(1 - wdenom.prob)^(2 - astar_lev)

        wdenom1 <- wdenom.prob^(a_lev - 1)*(1 - wdenom.prob)^(2 - a_lev)

      } else {

        data_pred <- data_boot[, exposure, drop = FALSE]

        data_pred[, exposure] <- as.factor(data_pred[, exposure])

        category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = data_pred)

        wdenom.prob <- predict(exposure_regression, newdata = data_boot,
                               type = ifelse(inherits(exposure_regression, "multinom")|
                                               inherits(exposure_regression, "polr"),
                                             "probs", "response"))

        wdenom0 <- wdenom.prob[, astar_lev]

        wdenom1 <- wdenom.prob[, a_lev]

      }

      ydesign0m <- data.frame(rep(astar, length(index_astar)), as.data.frame(mval)[rep(1,length(index_astar)),],
                              data_boot[index_astar,prec])

      ydesign1m <- data.frame(rep(a, length(index_a)), as.data.frame(mval)[rep(1,length(index_a)),],
                              data_boot[index_a,prec])

      ydesign01 <- data.frame(rep(astar, length(index_a)), data_boot[index_a, mediator],
                              data_boot[index_a, prec])

      ydesign10 <- data.frame(rep(a, length(index_astar)), data_boot[index_astar, mediator],
                              data_boot[index_astar, prec])

      colnames(ydesign0m) <- colnames(ydesign1m) <-
        colnames(ydesign10) <- colnames(ydesign01) <- c(exposure, mediator, prec)

      if (is.factor(data_boot[, exposure])) {

        ydesign0m[, exposure] <- factor(ydesign0m[, exposure], levels = levels(data_boot[, exposure]))

        ydesign1m[, exposure] <- factor(ydesign1m[, exposure], levels = levels(data_boot[, exposure]))

        ydesign10[, exposure] <- factor(ydesign10[, exposure], levels = levels(data_boot[, exposure]))

        ydesign01[, exposure] <- factor(ydesign01[, exposure], levels = levels(data_boot[, exposure]))

      }

      for (i in 1:length(mediator)) {

        if (is.factor(data_boot[, mediator[i]])) {
          ydesign0m[,1+i] <- factor(ydesign0m[,1+i], levels = levels(data_boot[, mediator[i]]))
          ydesign1m[,1+i] <- factor(ydesign1m[,1+i], levels = levels(data_boot[, mediator[i]]))
        }
      }

      EY0m <- weighted.mean(predict(outcome_regression, newdata = ydesign0m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = (wnom/wdenom0)[index_astar],
                            na.rm = TRUE)

      EY1m <- weighted.mean(predict(outcome_regression, newdata = ydesign1m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = (wnom/wdenom1)[index_a], na.rm = TRUE)

      EY00 <- weighted.mean(data_boot[index_astar, outcome], w = (wnom/wdenom0)[index_astar], na.rm = TRUE)

      EY11 <- weighted.mean(data_boot[index_a, outcome], w = (wnom/wdenom1)[index_a], na.rm = TRUE)

      EY01 <- weighted.mean(predict(outcome_regression, newdata = ydesign01,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = (wnom/wdenom1)[index_a], na.rm = TRUE)

      EY10 <- weighted.mean(predict(outcome_regression, newdata = ydesign10,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = (wnom/wdenom0)[index_astar], na.rm = TRUE)

    }


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

      cde <- EY1m - EY0m

      pnde <- EY10 - EY00

      tnde <- EY11 - EY01

      pnie <- EY01 - EY00

      tnie <- EY11 - EY10

      te <- tnie + pnde

      pm <- tnie / (pnde + te)

      intref <- pnde - cde

      intmed <- tnie - pnie

      pie <- pnie

      cde_prop <- cde/te

      intref_prop <- intref/te

      intmed_prop <- intmed/te

      pie_prop <- pie/te

      overall_pm <- (pie + intmed)/te

      overall_int <- (intref + intmed)/te

      overall_pe <- (intref + intmed + pie)/te

      out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie,
               te = te, pm = pm,
               intref = intref, intmed = intmed, pie = pie,
               cde_prop = cde_prop, intref_prop = intref_prop,
               intmed_prop = intmed_prop, pie_prop = pie_prop,
               overall_pm = overall_pm, overall_int = overall_int, overall_pe = overall_pe)

      if (!is.null(postc) && model %in% c("msm", "g-formula")) {

        names(out) <- c("cde","rpnde","rtnde","rpnie","rtnie","rte","pm",
                        "intref","intmed","pie",
                        "cde_prop","intref_prop","intmed_prop","pie_prop",
                        "overall_pm", "overall_int", "overall_pe")

      }

    } else {

      cde_rr <- EY1m/EY0m

      pnde_rr <- EY10/EY00

      tnde_rr <- EY11/EY01

      pnie_rr <- EY01/EY00

      tnie_rr <- EY11/EY10

      te_rr <- tnie_rr * pnde_rr

      pm <- (pnde_rr * (tnie_rr - 1)) / (te_rr - 1)

      cde_err <- (EY1m-EY0m)/EY00

      intref_err <- pnde_rr - 1 - cde_err

      intmed_err <- tnie_rr * pnde_rr - pnde_rr - pnie_rr + 1

      pie_err <- pnie_rr - 1

      te_err <- te_rr - 1

      cde_err_prop <- cde_err/te_err

      intmed_err_prop <- intmed_err/te_err

      intref_err_prop <- intref_err/te_err

      pie_err_prop <- pie_err/te_err

      overall_pm <- (pie_err+intmed_err)/te_err

      overall_int <- (intref_err+intmed_err)/te_err

      overall_pe <- (intref_err+intmed_err+pie_err)/te_err

      out <- c(cde_rr = cde_rr, pnde_rr = pnde_rr, tnde_rr = tnde_rr, pnie_rr = pnie_rr,
               tnie_rr = tnie_rr, te_rr = te_rr, pm = pm,
               cde_err = cde_err, intref_err = intref_err,
               intmed_err = intmed_err, pie_err = pie_err,
               cde_err_prop = cde_err_prop, intref_err_prop = intref_err_prop,
               intmed_err_prop = intmed_err_prop, pie_err_prop = pie_err_prop,
               overall_pm = overall_pm, overall_int = overall_int,
               overall_pe = overall_pe)

      if (!is.null(postc) && model %in% c("msm", "g-formula")) {

        names(out) <- c("cde_rr","rpnde_rr","rtnde_rr","rpnie_rr","rtnie_rr","rte_rr","pm",
                        "cde_err","intref_err","intmed_err","pie_err",
                        "cde_err_prop","intref_err_prop","intmed_err_prop","pie_err_prop",
                        "overall_pm", "overall_int", "overall_pe")

      }

    }
  }

  return(out)

}
