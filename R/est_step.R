est_step <- function(data, indices, model,
                     outcome, event, exposure, mediator, EMint, prec, postc,
                     yreg, mreg, ereg, postcreg, wmreg, reg.simex,
                     astar, a, mval, yref, vecc,
                     estimation) {

  data_boot <- data[indices, ]

  formulas <- create_formulas(data = data_boot, model = model,
                              outcome = outcome, event = event,
                              exposure = exposure, mediator = mediator, EMint = EMint,
                              prec = prec, postc = postc,
                              yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

  regressions <- run_regressions(formulas = formulas, data = data_boot, model = model,
                                 exposure = exposure, mediator = mediator, postc = postc,
                                 yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

  outcome_regression <- regressions$outcome_regression

  if ((identical(class(outcome_regression), c("gam", "glm", "lm")) &&
       (family(outcome_regression)$family %in% c("multinom") |
        startsWith(family(outcome_regression)$family,"Ordered Categorical")))|
      identical(class(outcome_regression), c("multinom", "nnet"))|
      identical(class(outcome_regression), "polr")) {

    if (is.null(yref)) {
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



    if (model == "rb"){

      if (is.character(data_boot[, exposure])|is.factor(data_boot[, exposure])) {
        a <- as.numeric(levels(as.factor(data_boot[, exposure])) == a)[-1]
        astar <- as.numeric(levels(as.factor(data_boot[, exposure])) == astar)[-1]
      }

      if (is.character(data_boot[, mediator])|is.factor(data_boot[, mediator])) {
        mstar <- as.numeric(levels(as.factor(data_boot[, mediator])) == mval[[1]])[-1]
      } else {mstar <- mval[[1]]}

      coef <- get_coef(formulas = formulas, regressions = regressions, model = model,
                     yreg = yreg, mreg = mreg, data = data)

      elevel <- ifelse(is.character(data_boot[, exposure])|is.factor(data_boot[, exposure]),
                       length(levels(data_boot[, exposure])), 2)

      mlevel <- ifelse(is.character(data_boot[, mediator])|is.factor(data_boot[, mediator]),
                       length(levels(data_boot[, mediator])), 2)

      thetas <- coef$thetas

      betas <- coef$betas

      variance <- coef$variance

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

      regressions_coef <- regressions

      if (!is.null(reg.simex$tot_outcome_regression)) regressions_coef$tot_outcome_regression <- reg.simex$tot_outcome_regression
      if (!is.null(reg.simex$dir_outcome_regression)) regressions_coef$dir_outcome_regression <- reg.simex$dir_outcome_regression

      coef <- get_coef(formulas = formulas, regressions = regressions_coef, model = model,
                       yreg = yreg, mreg = mreg, data = data)

      if (yreg == "linear") {

        dir <- unname(ifelse(sum(a) == 0, 0, coef[1:(elevel - 1)] %*% a) -
                        ifelse(sum(astar) == 0, 0, coef[1:(elevel - 1)] %*% astar))

        tot <- unname(ifelse(sum(a) == 0, 0, coef[elevel:(2*(elevel - 1))] %*% a) -
                        ifelse(sum(astar) == 0, 0, coef[elevel:(2*(elevel - 1))] %*% astar))

        ind <- tot - dir

        out <- c(tot = tot, dir = dir, ind = ind)

      } else {

        RRdir <- unname(exp(ifelse(sum(a) == 0, 0, coef[1:(elevel - 1)] %*% a) -
                              ifelse(sum(astar) == 0, 0, coef[1:(elevel - 1)] %*% astar)))

        RRtot <- unname(exp(ifelse(sum(a) == 0, 0, coef[elevel:(2*(elevel - 1))] %*% a) -
                              ifelse(sum(astar) == 0, 0, coef[elevel:(2*(elevel - 1))] %*% astar)))



        RRind <- RRtot / RRdir

        out <- c(RRtot = RRtot, RRdir = RRdir, RRind = RRind)

      }

    }

    ###################################################################################################
    ###############################Direct Counterfactual Imputation Estimation#########################
    ###################################################################################################

  } else if (estimation == "imputation") {

    if (model %in% c("rb", "msm", "g-formula")) {

      if (model == "g-formula") {

        if (!is.null(postc)){

          postc_regression <- regressions$postc_regression

          postc_regression_pred <- postc_regression

          for (i in 1:length(postc_regression)) {

           if (!is.null(reg.simex$postc_regression[[i]])) postc_regression_pred[[i]] <- reg.simex$postc_regression[[i]]

          }

          postcdesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                                      data_boot[,prec])

          postcdesign_astar <- data.frame(c(rep(astar,nrow(data_boot))),
                                          data_boot[,prec])

          colnames(postcdesign_a) <- colnames(postcdesign_astar) <-
            c(exposure, prec)

          if (is.factor(data_boot[, exposure])) {

            postcdesign_a[, exposure] <- factor(postcdesign_a[, exposure], levels = levels(data_boot[, exposure]))

            postcdesign_astar[, exposure] <- factor(postcdesign_astar[, exposure], levels = levels(data_boot[, exposure]))

          }

          if (identical(class(postc_regression[[1]]), "lm") |
              (identical(class(postc_regression[[1]]), c("glm", "lm")) &&
               family(postc_regression[[1]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                           "poisson", "quasipoisson"))|
              (identical(class(postc_regression[[1]]), c("gam", "glm", "lm")) &&
               (family(postc_regression[[1]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                            "poisson", "quasipoisson") |
                startsWith(family(postc_regression[[1]])$family,"Negative Binomial")))|
              identical(class(postc_regression[[1]]), "nls")|
              identical(class(postc_regression[[1]]), c("negbin", "glm", "lm"))) {

            postc_a <- predict(postc_regression_pred[[1]], newdata = postcdesign_a, type = "response")

            postc_astar <- predict(postc_regression_pred[[1]], newdata = postcdesign_astar, type = "response")

          } else if ((identical(class(postc_regression[[1]]), c("glm", "lm")) &&
                      family(postc_regression[[1]])$family %in% c("binomial", "quasibinomial"))|
                     (identical(class(postc_regression[[1]]), c("gam", "glm", "lm")) &&
                      (family(postc_regression[[1]])$family %in% c("binomial", "quasibinomial")))) {

            prob_a <- predict(postc_regression_pred[[1]], newdata = postcdesign_a, type = "response")

            pred_a <- rbinom(nrow(postcdesign_a), size=1, prob=prob_a)

            prob_astar <- predict(postc_regression_pred[[1]], newdata = postcdesign_astar, type = "response")

            pred_astar <- rbinom(nrow(postcdesign_astar), size=1, prob=prob_astar)

            if (is.numeric(data_boot[, postc[1]])) {

              postc_a <- as.numeric(levels(as.factor(data_boot[, postc[1]]))[pred_a + 1])

              postc_astar <- as.numeric(levels(as.factor(data_boot[, postc[1]]))[pred_star + 1])

            } else if (is.factor(data_boot[, postc[1]])) {

              postc_a <- factor(levels(as.factor(data_boot[, postc[1]]))[pred_a + 1],
                                levels = levels(as.factor(data_boot[, postc[1]])))

              postc_astar <- factor(levels(as.factor(data_boot[, postc[1]]))[pred_star + 1],
                                    levels = levels(as.factor(data_boot[, postc[1]])))

            }

          } else if ((identical(class(postc_regression[[1]]), c("gam", "glm", "lm")) &&
                      (family(postc_regression[[1]])$family %in% c("multinom") |
                       startsWith(family(postc_regression[[1]])$family,"Ordered Categorical")))|
                     identical(class(postc_regression[[1]]), c("multinom", "nnet"))|
                     identical(class(postc_regression[[1]]), "polr")) {

            prob_a <- predict(postc_regression_pred[[1]], newdata = postcdesign_a,
                              type = ifelse(inherits(postc_regression[[1]], "gam"), "response", "probs"))

            pred_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

            prob_astar <- predict(postc_regression_pred[[1]], newdata = postcdesign_astar,
                                  type = ifelse(inherits(postc_regression[[1]], "gam"), "response", "probs"))

            pred_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

            if (is.numeric(data_boot[, postc[1]])) {

              postc_a <- as.numeric(levels(as.factor(data_boot[, postc[1]]))[pred_a])

              postc_astar <- as.numeric(levels(as.factor(data_boot[, postc[1]]))[pred_astar])

            } else if (is.factor(data_boot[, postc[1]])) {

              postc_a <- factor(levels(as.factor(data_boot[, postc[1]]))[pred_a],
                                levels = levels(as.factor(data_boot[, postc[1]])))

              postc_astar <- factor(levels(as.factor(data_boot[, postc[1]]))[pred_astar],
                                    levels = levels(as.factor(data_boot[, postc[1]])))

            }
          }

          postc_a <- data_frame(postc_a)

          postc_astar <- data_frame(postc_astar)

          if (length(postc) > 1) {

            for (i in 2:length(postc)) {

              postcdesign_a <- data_frame(postcdesign_a, postc_a[, i-1])

              colnames(postcdesign_a) <- c(colnames(postcdesign_a), postc[i-1])

              postcdesign_astar <- data_frame(postcdesign_astar, postc_astar[, i-1])

              colnames(postcdesign_astar) <- c(colnames(postcdesign_astar), postc[i-1])

              if (identical(class(postc_regression[[i]]), "lm") |
                  (identical(class(postc_regression[[i]]), c("glm", "lm")) &&
                   family(postc_regression[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                               "poisson", "quasipoisson"))|
                  (identical(class(postc_regression[[i]]), c("gam", "glm", "lm")) &&
                   (family(postc_regression[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                                "poisson", "quasipoisson") |
                    startsWith(family(postc_regression[[i]])$family,"Negative Binomial")))|
                  identical(class(postc_regression[[i]]), "nls")|
                  identical(class(postc_regression[[i]]), c("negbin", "glm", "lm"))) {

                mid_a <- predict(postc_regression_pred[[i]], newdata = postcdesign_a, type = "response")

                mid_astar <- predict(postc_regression_pred[[i]], newdata = postcdesign_astar, type = "response")

              } else if ((identical(class(postc_regression[[i]]), c("glm", "lm")) &&
                          family(postc_regression[[i]])$family %in% c("binomial", "quasibinomial"))|
                         (identical(class(postc_regression[[i]]), c("gam", "glm", "lm")) &&
                          (family(postc_regression[[i]])$family %in% c("binomial", "quasibinomial")))) {

                prob_a <- predict(postc_regression_pred[[i]], newdata = postcdesign_a, type = "response")

                pred_a <- rbinom(nrow(postcdesign_a), size=1, prob=prob_a)

                prob_astar <- predict(postc_regression_pred[[i]], newdata = postcdesign_astar, type = "response")

                pred_astar <- rbinom(nrow(postcdesign_astar), size=1, prob=prob_astar)

                if (is.numeric(data_boot[, postc[i]])) {

                  mid_a <- as.numeric(levels(as.factor(data_boot[, postc[i]]))[pred_a + 1])

                  mid_astar <- as.numeric(levels(as.factor(data_boot[, postc[i]]))[pred_star + 1])

                } else if (is.factor(data_boot[, postc[i]])) {

                  mid_a <- factor(levels(as.factor(data_boot[, postc[i]]))[pred_a + 1],
                                  levels = levels(as.factor(data_boot[, postc[i]])))

                  mid_astar <- factor(levels(as.factor(data_boot[, postc[i]]))[pred_star + 1],
                                      levels = levels(as.factor(data_boot[, postc[i]])))

                }

              } else if ((identical(class(postc_regression[[i]]), c("gam", "glm", "lm")) &&
                          (family(postc_regression[[i]])$family %in% c("multinom") |
                           startsWith(family(postc_regression[[i]])$family,"Ordered Categorical")))|
                         identical(class(postc_regression[[i]]), c("multinom", "nnet"))|
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

                  mid_a <- factor(levels(as.factor(data_boot[, postc[i]]))[pred_a],
                                  levels = levels(as.factor(data_boot[, postc[i]])))

                  mid_astar <- factor(levels(as.factor(data_boot[, postc[i]]))[pred_astar],
                                      levels = levels(as.factor(data_boot[, postc[i]])))

                }
              }

              postc_a <- data_frame(postc_a, mid_a)

              postc_astar <- data_frame(postc_a, mid_astar)

            }
          }

        } else {
          postc_astar <- postc_a <-data.frame()[1:nrow(data_boot), ]
        }


        mdesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                                data_boot[,prec],
                                postc_a)

        mdesign_astar <- data.frame(c(rep(astar,nrow(data_boot))),
                                    data_boot[,prec],
                                    postc_astar)

        colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, prec, postc)

      } else {

        mdesign_a <- data.frame(rep(a,nrow(data_boot)),data_boot[,prec])

        mdesign_astar <- data.frame(rep(astar,nrow(data_boot)),data_boot[,prec])

        colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, prec)

      }


      mediator_regression <- regressions$mediator_regression

      mediator_regression_pred <- mediator_regression

      for (i in 1:length(mediator_regression)) {

        if (!is.null(reg.simex$mediator_regression[[i]])) mediator_regression_pred[[i]] <- reg.simex$mediator_regression[[i]]

      }

      if (is.factor(data_boot[, exposure])) {

        mdesign_a[, exposure] <- factor(mdesign_a[, exposure], levels = levels(data_boot[, exposure]))

        mdesign_astar[, exposure] <- factor(mdesign_astar[, exposure], levels = levels(data_boot[, exposure]))

      }

      if (identical(class(mediator_regression[[1]]), "lm") |
          (identical(class(mediator_regression[[1]]), c("glm", "lm")) &&
           family(mediator_regression[[1]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                          "poisson", "quasipoisson"))|
          (identical(class(mediator_regression[[1]]), c("gam", "glm", "lm")) &&
           (family(mediator_regression[[1]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                           "poisson", "quasipoisson") |
            startsWith(family(mediator_regression[[1]])$family,"Negative Binomial")))|
          identical(class(mediator_regression[[1]]), "nls")|
          identical(class(mediator_regression[[1]]), c("negbin", "glm", "lm"))) {

        m_a <- predict(mediator_regression_pred[[1]], newdata = mdesign_a, type = "response")

        m_astar <- predict(mediator_regression_pred[[1]], newdata = mdesign_astar, type = "response")

      } else if ((identical(class(mediator_regression[[1]]), c("glm", "lm")) &&
                  family(mediator_regression[[1]])$family %in% c("binomial", "quasibinomial"))|
                 (identical(class(mediator_regression[[1]]), c("gam", "glm", "lm")) &&
                  (family(mediator_regression[[1]])$family %in% c("binomial", "quasibinomial")))) {

        prob_a <- predict(mediator_regression_pred[[1]], newdata = mdesign_a, type = "response")

        pred_a <- rbinom(nrow(mdesign_a), size=1, prob=prob_a)

        prob_astar <- predict(mediator_regression_pred[[1]], newdata = mdesign_astar, type = "response")

        pred_astar <- rbinom(nrow(mdesign_astar), size=1, prob=prob_astar)

        if (is.numeric(data_boot[, mediator[1]])) {

          m_a <- as.numeric(levels(as.factor(data_boot[, mediator[1]]))[pred_a + 1])

          m_astar <- as.numeric(levels(as.factor(data_boot[, mediator[1]]))[pred_star + 1])

        } else if (is.factor(data_boot[, mediator[1]])) {

          m_a <- factor(levels(as.factor(data_boot[, mediator[1]]))[pred_a + 1],
                        levels = levels(as.factor(data_boot[, mediator[1]])))

          m_astar <- factor(levels(as.factor(data_boot[, mediator[1]]))[pred_star + 1],
                            levels = levels(as.factor(data_boot[, mediator[1]])))

        }

      } else if ((identical(class(mediator_regression[[1]]), c("gam", "glm", "lm")) &&
                  (family(mediator_regression[[1]])$family %in% c("multinom") |
                   startsWith(family(mediator_regression[[1]])$family,"Ordered Categorical")))|
                 identical(class(mediator_regression[[1]]), c("multinom", "nnet"))|
                 identical(class(mediator_regression[[1]]), "polr")) {

        prob_a <- predict(mediator_regression_pred[[1]], newdata = mdesign_a,
                          type = ifelse(inherits(mediator_regression[[1]], "gam"), "response", "probs"))

        pred_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

        prob_astar <- predict(mediator_regression_pred[[1]], newdata = mdesign_astar,
                              type = ifelse(inherits(mediator_regression[[1]], "gam"), "response", "probs"))

        pred_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

        if (is.numeric(data_boot[, mediator[1]])) {

          m_a <- as.numeric(levels(as.factor(data_boot[, mediator[1]]))[pred_a])

          m_astar <- as.numeric(levels(as.factor(data_boot[, mediator[1]]))[pred_astar])

        } else if (is.factor(data_boot[, mediator[1]])) {

          m_a <- factor(levels(as.factor(data_boot[, mediator[1]]))[pred_a],
                        levels = levels(as.factor(data_boot[, mediator[1]])))

          m_astar <- factor(levels(as.factor(data_boot[, mediator[1]]))[pred_astar],
                            levels = levels(as.factor(data_boot[, mediator[1]])))

        }
      }

      m_a <- data_frame(m_a)

      m_astar <- data_frame(m_astar)

      if (length(mediator) > 1) {

        for (i in 2:length(mediator)) {

          mdesign_a <- data_frame(mdesign_a, m_a[, i-1])

          colnames(mdesign_a) <- c(colnames(mdesign_a), mediator[i-1])

          mdesign_astar <- data_frame(mdesign_astar, m_astar[, i-1])

          colnames(mdesign_astar) <- c(colnames(mdesign_astar), mediator[i-1])

          if (identical(class(mediator_regression[[i]]), "lm") |
              (identical(class(mediator_regression[[i]]), c("glm", "lm")) &&
               family(mediator_regression[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                              "poisson", "quasipoisson"))|
              (identical(class(mediator_regression[[i]]), c("gam", "glm", "lm")) &&
               (family(mediator_regression[[i]])$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                               "poisson", "quasipoisson") |
                startsWith(family(mediator_regression[[i]])$family,"Negative Binomial")))|
              identical(class(mediator_regression[[i]]), "nls")|
              identical(class(mediator_regression[[i]]), c("negbin", "glm", "lm"))) {

            mid_a <- predict(mediator_regression_pred[[i]], newdata = mdesign_a, type = "response")

            mid_astar <- predict(mediator_regression_pred[[i]], newdata = mdesign_astar, type = "response")

          } else if ((identical(class(mediator_regression[[i]]), c("glm", "lm")) &&
                      family(mediator_regression[[i]])$family %in% c("binomial", "quasibinomial"))|
                     (identical(class(mediator_regression[[i]]), c("gam", "glm", "lm")) &&
                      (family(mediator_regression[[i]])$family %in% c("binomial", "quasibinomial")))) {

            prob_a <- predict(mediator_regression_pred[[i]], newdata = mdesign_a, type = "response")

            pred_a <- rbinom(nrow(mdesign_a), size=1, prob=prob_a)

            prob_astar <- predict(mediator_regression_pred[[i]], newdata = mdesign_astar, type = "response")

            pred_astar <- rbinom(nrow(mdesign_astar), size=1, prob=prob_astar)

            if (is.numeric(data_boot[, mediator[i]])) {

              mid_a <- as.numeric(levels(as.factor(data_boot[, mediator[i]]))[pred_a + 1])

              mid_astar <- as.numeric(levels(as.factor(data_boot[, mediator[i]]))[pred_star + 1])

            } else if (is.factor(data_boot[, mediator[i]])) {

              mid_a <- factor(levels(as.factor(data_boot[, mediator[i]]))[pred_a + 1],
                              levels = levels(as.factor(data_boot[, mediator[i]])))

              mid_astar <- factor(levels(as.factor(data_boot[, mediator[i]]))[pred_star + 1],
                                  levels = levels(as.factor(data_boot[, mediator[i]])))

            }

          } else if ((identical(class(mediator_regression[[i]]), c("gam", "glm", "lm")) &&
                      (family(mediator_regression[[i]])$family %in% c("multinom") |
                       startsWith(family(mediator_regression[[i]])$family,"Ordered Categorical")))|
                     identical(class(mediator_regression[[i]]), c("multinom", "nnet"))|
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

              mid_a <- factor(levels(as.factor(data_boot[, mediator[i]]))[pred_a],
                              levels = levels(as.factor(data_boot[, mediator[i]])))

              mid_astar <- factor(levels(as.factor(data_boot[, mediator[i]]))[pred_astar],
                                  levels = levels(as.factor(data_boot[, mediator[i]])))

            }
          }

          m_a <- data_frame(m_a, mid_a)

          m_astar <- data_frame(m_a, mid_astar)

        }
      }


      if (model != "g-formula") {

        ydesign0m <- data.frame(rep(astar, nrow(data_boot)), as.data.frame(mval)[rep(1,nrow(data_boot)),],
                                data_boot[,prec])

        ydesign1m <- data.frame(rep(a, nrow(data_boot)), as.data.frame(mval)[rep(1,nrow(data_boot)),],
                                data_boot[,prec])

        ydesign00 <- data.frame(rep(astar, nrow(data_boot)), m_astar, data_boot[,prec])

        ydesign01 <- data.frame(rep(astar, nrow(data_boot)), m_a, data_boot[,prec])

        ydesign10 <- data.frame(rep(a, nrow(data_boot)), m_astar, data_boot[,prec])

        ydesign11 <- data.frame(rep(a, nrow(data_boot)), m_a, data_boot[,prec])

        colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <-
          colnames(ydesign01) <- colnames(ydesign10) <- colnames(ydesign11) <-
          c(exposure, mediator, prec)

      } else {

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

      }


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

      if (!is.null(reg.simex$outcome_regression)) {
        outcome_regression_pred <- reg.simex$outcome_regression
      } else {outcome_regression_pred <- outcome_regression}

      if (identical(class(outcome_regression), "lm") |
          (identical(class(outcome_regression), c("glm", "lm")) &&
           family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                    "binomial", "quasibinomial","poisson", "quasipoisson"))|
          (identical(class(outcome_regression), c("gam", "glm", "lm")) &&
           (family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                     "binomial", "quasibinomial","poisson", "quasipoisson") |
            startsWith(family(outcome_regression)$family,"Negative Binomial")))|
          identical(class(outcome_regression), "nls")|
          identical(class(outcome_regression), c("negbin", "glm", "lm"))|
          identical(class(outcome_regression), "survreg")|
          identical(class(outcome_regression), "coxph")) {

        type <- ifelse(inherits(outcome_regression, "coxph"), "risk", "response")

        EY0m <- mean(predict(outcome_regression_pred, newdata =  ydesign0m, type = type), na.rm = TRUE)

        EY1m <- mean(predict(outcome_regression_pred, newdata =  ydesign1m, type = type), na.rm = TRUE)

        EY00 <- mean(predict(outcome_regression_pred, newdata =  ydesign00, type = type), na.rm = TRUE)

        EY01 <- mean(predict(outcome_regression_pred, newdata =  ydesign01, type = type), na.rm = TRUE)

        EY10 <- mean(predict(outcome_regression_pred, newdata =  ydesign10, type = type), na.rm = TRUE)

        EY11 <- mean(predict(outcome_regression_pred, newdata =  ydesign11, type = type), na.rm = TRUE)


      } else if ((identical(class(outcome_regression), c("gam", "glm", "lm")) &&
                  (family(outcome_regression)$family %in% c("multinom") |
                   startsWith(family(outcome_regression)$family,"Ordered Categorical")))|
                 identical(class(outcome_regression), c("multinom", "nnet"))|
                 identical(class(outcome_regression), "polr")) {

        type = ifelse(inherits(outcome_regression, "gam"), "response", "probs")

        EY0m <- mean(predict(outcome_regression_pred, newdata =  ydesign0m, type = type)[, yref], na.rm = TRUE)

        EY1m <- mean(predict(outcome_regression_pred, newdata =  ydesign1m, type = type)[, yref], na.rm = TRUE)

        EY00 <- mean(predict(outcome_regression_pred, newdata =  ydesign00, type = type)[, yref], na.rm = TRUE)

        EY01 <- mean(predict(outcome_regression_pred, newdata =  ydesign01, type = type)[, yref], na.rm = TRUE)

        EY10 <- mean(predict(outcome_regression_pred, newdata =  ydesign10, type = type)[, yref], na.rm = TRUE)

        EY11 <- mean(predict(outcome_regression_pred, newdata =  ydesign11, type = type)[, yref], na.rm = TRUE)

      }

    } else if (model == "wb") {

      exposure_regression <- regressions$exposure_regression

      if (!is.null(reg.simex$exposure_regression)) {
       exposure_regression_pred <- reg.simex$exposure_regression
      } else {exposure_regression_pred <- exposure_regression}

      outcome_regression <- regressions$outcome_regression

      if (!(identical(class(outcome_regression), "lm") |
           (identical(class(outcome_regression), c("glm", "lm")) &&
            family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                     "binomial", "quasibinomial","poisson", "quasipoisson"))|
           (identical(class(outcome_regression), c("gam", "glm", "lm")) &&
            (family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi",
                                                      "binomial", "quasibinomial","poisson", "quasipoisson") |
             startsWith(family(outcome_regression)$family,"Negative Binomial")))|
           identical(class(outcome_regression), "nls")|
           identical(class(outcome_regression), c("negbin", "glm", "lm")))) {

        stop("The selected model doesn't support this outcome regression")

      } else if ((identical(class(outcome_regression), c("glm", "lm")) &&
                  family(outcome_regression)$family %in% c("binomial", "quasibinomial"))|
                 (identical(class(outcome_regression), c("gam", "glm", "lm")) &&
                  (family(outcome_regression)$family %in% c("binomial", "quasibinomial")))){

        data_boot[, outcome] <- as.factor(data_boot[, outcome])

        levels(data_boot[, outcome]) <- 0:1

        data_boot[, outcome] <- as.numeric(data_boot[, outcome]) - 1

        outcome_regression <- update(outcome_regression, formula. = formula(outcome_regression),
                                     data = data_boot)

      }

      index_astar <- which(data_boot[, exposure] == astar)

      index_a <- which(data_boot[, exposure] == a)

      wnom <- left_join(select(data_boot, exposure),
                         count(data_boot, !!as.name(exposure)),
                         by = exposure)[, "n"]/nrow(data_boot)


      if (inherits(exposure_regression, "glm") &&
          family(exposure_regression)$family %in% c("binomial", "quasibinomial")) {

        wdenom.prob <- predict(exposure_regression_pred, newdata = data_boot,
                                type = "response")

        class <- as.numeric(as.factor(data_boot[, exposure])) - 1

        wdenom <-  wdenom.prob ^ class * (1 - wdenom.prob) ^ (1 - class)


      } else if ((inherits(exposure_regression, "gam") &&
                  (family(exposure_regression)$family == "multinom" |
                   startsWith(family(exposure_regression)$family, "Ordered Categorical")))|
                 inherits(exposure_regression, "multinom")|inherits(exposure_regression, "polr")) {

        wdenom.prob <- predict(exposure_regression_pred, newdata = data_boot,
                                type = ifelse(inherits(exposure_regression, "multinom")|
                                                inherits(exposure_regression, "polr"),
                                              "probs","response"))

        class <- as.numeric(as.factor(data_boot[, exposure]))

        wdenom <- sapply(1:nrow(data_boot),FUN = function(i) wdenom.prob[i, class[i]])

      }

      w <- wnom/wdenom

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
                            w = w[index_astar],
                            na.rm = TRUE)

      EY1m <- weighted.mean(predict(outcome_regression, newdata = ydesign1m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a], na.rm = TRUE)

      EY00 <- weighted.mean(data_boot[index_astar, outcome], w = w[index_astar], na.rm = TRUE)

      EY11 <- weighted.mean(data_boot[index_a, outcome], w = w[index_a], na.rm = TRUE)

      EY01 <- weighted.mean(predict(outcome_regression, newdata = ydesign01,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a], na.rm = TRUE)

      EY10 <- weighted.mean(predict(outcome_regression, newdata = ydesign10,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_astar], na.rm = TRUE)

    }


    if (identical(class(outcome_regression), "lm") |
        (identical(class(outcome_regression), c("glm", "lm")) &&
         family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))|
        (identical(class(outcome_regression), c("gam", "glm", "lm")) &&
         (family(outcome_regression)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi")))|
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
