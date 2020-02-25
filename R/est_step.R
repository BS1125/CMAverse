est_step <- function(data, indices, outcome, exposure, exposure.type, mediator,
                     covariates.pre, covariates.post, covariates.post.type, vecc,
                     EMint, MMint, EMMint, EMint.terms,
                     MMint.terms, EMMint.terms,
                     event, mreg, yreg, m_star, a_star, a, est.method, model) {

  data_boot <- data[indices, ]

  formulas <- create_formulas(model = model, outcome = outcome, exposure = exposure,
                              mediator = mediator, covariates.pre = covariates.pre,
                              covariates.post = covariates.post,
                              EMint = EMint, MMint = MMint, EMMint = EMMint,
                              EMint.terms = EMint.terms, MMint.terms = MMint.terms,
                              EMMint.terms = EMMint.terms,
                              event = event, mreg = mreg, yreg = yreg)

  regressions <- run_regressions(model = model, formulas = formulas,
                                 exposure = exposure,  exposure.type = exposure.type,
                                 mediator = mediator, covariates.post = covariates.post,
                                 covariates.post.type = covariates.post.type,
                                 mreg = mreg, yreg = yreg, data = data_boot)

  if (est.method == "paramfunc") {

    if (!(model %in% c("rb", "iorw"))) {
      stop("Closed-form parameter function estimation doesn't support the selected model")
    } else if (length(mediator) > 1) {
      stop("Closed-form parameter function estimation doesn't support multiple mediator cases")
    }

    coef <- get_coef(formulas = formulas, regressions = regressions,
                     mreg = mreg, yreg = yreg, model = model)

    if (model == "rb"){

      if (!(exposure.type %in% c("continuous", "binary"))) {
        stop("For the selected model, closed-form parameter function estimation only supports continuous or binary exposure")
      }

      mlevel <- ifelse(is.factor(data_boot[, mediator]), length(levels(data_boot[, mediator])), 2)

      ############################################DE and IE#############################################

      thetas <- coef$thetas

      betas <- coef$betas

      variance <- coef$variance

      theta0 <- thetas[1]

      theta1 <- thetas[2]

      theta2 <- thetas[3:(mlevel + 1)]

      if (EMint == TRUE) {
        theta3 <- thetas[length(thetas)-((mlevel-2):0)]
      } else {theta3 <- rep(0, mlevel-1)}

      beta0 <- betas[1+(0:(mlevel-2))*length(betas)/(mlevel-1)]

      beta1 <- betas[2+(0:(mlevel-2))*length(betas)/(mlevel-1)]

      covariatesTerm <- sapply(0:(mlevel-2), function(x) sum(betas[2 + 1:length(vecc) + x*length(betas)/(mlevel-1)]*vecc))

      if (yreg == "linear") {

        if (mreg == "linear") {

          cde <- unname((theta1 + theta3 * m_star) * (a-a_star))

          pnde <- unname((theta1 + theta3 * (beta0 + beta1 * a_star +
                                               covariatesTerm)) * (a - a_star))

          tnde <- unname((theta1 + theta3 * (beta0 + beta1 * a +
                                               covariatesTerm)) * (a - a_star))

          pnie <- unname((theta2 + theta3  * a_star) * beta1 * (a - a_star))

          tnie <- unname((theta2 + theta3  * a) * beta1 * (a - a_star))

        } else if (mreg %in% c("logistic", "multinomial")) {

          cde <- unname((theta1 + ifelse(m_star == 0, 0, theta3[m_star])) * (a-a_star))

          pnde <- unname((theta1 +
                            (sum(theta3 * exp(beta0 + beta1 * a_star +
                                                covariatesTerm)) /
                               (1 + sum(exp(beta0 + beta1 * a_star +
                                              covariatesTerm))))) * (a - a_star))
          tnde <- unname((theta1  +
                            (sum(theta3 * exp(beta0 + beta1 * a +
                                                covariatesTerm)) /
                               (1 + sum(exp(beta0 + beta1 * a +
                                              covariatesTerm))))) * (a - a_star))

          pnie <- unname(sum((theta2+theta3*a_star)*exp(beta0 + beta1 * a + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 * a + covariatesTerm))) -
                           sum((theta2+theta3*a_star)*exp(beta0 + beta1 * a_star + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 * a_star + covariatesTerm))))

          tnie <- unname(sum((theta2+theta3*a)*exp(beta0 + beta1 * a + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 * a + covariatesTerm))) -
                           sum((theta2+theta3*a)*exp(beta0 + beta1 * a_star + covariatesTerm)) /
                           (1 + sum(exp(beta0 + beta1 * a_star + covariatesTerm))))

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

      } else if (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                             "negbin", "coxph", "aft_exp", "aft_weibull")) {

        if (mreg == "linear") {

          cde_rr <- unname(exp((theta1 + theta3 * m_star) * (a - a_star)))

          pnde_rr <- unname(exp((theta1 + theta3 * (beta0 + beta1 * a_star +
                                                      covariatesTerm + theta2  * variance)) * (a - a_star) +
                                  0.5 * theta3 ^ 2 * variance * (a ^ 2 - a_star ^ 2)))

          tnde_rr <- unname(exp((theta1 + theta3 * (beta0 + beta1 * a +
                                                      covariatesTerm + theta2  * variance)) * (a - a_star) +
                                  0.5 * theta3 ^ 2 * variance * (a ^ 2 - a_star ^ 2)))

          pnie_rr <- unname(exp((theta2 + theta3  * a_star) * beta1 * (a - a_star)))

          tnie_rr <- unname(exp((theta2 + theta3  * a) * beta1 * (a - a_star)))

          cde_err <- unname((exp(theta1 * (a - a_star) + theta3 * a * m_star) -
                               exp(theta3 * a_star * m_star)) *
                              exp(theta2 * m_star - (theta2 + theta3 * a_star) *
                                    (beta0 + beta1 * a_star + covariatesTerm) -
                                    0.5 * (theta2 + theta3 * a_star) ^ 2 * variance))

        } else if (mreg %in% c("logistic", "multinomial")) {

          cde_rr <- unname(exp((theta1 + ifelse(m_star == 0, 0, theta3[m_star])) * (a-a_star)))

          pnde_rr <- unname((exp(theta1 * (a - a_star)) *
                               (1 + sum(exp(theta2 + theta3 * a +
                                              beta0 + beta1 * a_star + covariatesTerm)))) /
                              (1 + sum(exp(theta2 + theta3 * a_star + beta0 + beta1 * a_star + covariatesTerm))))

          tnde_rr <- unname((exp(theta1 * (a - a_star)) *
                               (1 + sum(exp(theta2 + theta3 * a +
                                              beta0 + beta1 * a + covariatesTerm)))) /
                              (1 + sum(exp(theta2 + theta3 * a_star + beta0 + beta1 * a + covariatesTerm))))

          pnie_rr <- unname(((1 + sum(exp(beta0 + beta1 * a_star + covariatesTerm))) *
                               (1 + sum(exp(theta2 + theta3 * a_star + beta0 +
                                              beta1 * a + covariatesTerm)))) /
                              ((1 + sum(exp(beta0 + beta1 * a + covariatesTerm))) *
                                 (1+ sum(exp(theta2 + theta3 * a_star + beta0 +
                                               beta1 * a_star + covariatesTerm)))))

          tnie_rr <- unname(((1 + sum(exp(beta0 + beta1 * a_star + covariatesTerm))) *
                               (1 + sum(exp(theta2 + theta3 * a + beta0 +
                                              beta1 * a + covariatesTerm)))) /
                              ((1 + sum(exp(beta0 + beta1 * a + covariatesTerm))) *
                                 (1+ sum(exp(theta2 + theta3 * a + beta0 +
                                               beta1 * a_star + covariatesTerm)))))

          cde_err <- unname(exp(ifelse(m_star == 0, 0, theta2[m_star])) *
                              (exp(theta1 * (a - a_star) + a * ifelse(m_star == 0, 0, theta3[m_star])) -
                                 exp(a_star * ifelse(m_star == 0, 0, theta3[m_star]))) *
                              (1 + sum(exp(beta0 + beta1 * a_star + covariatesTerm))) /
                              (1+ sum(exp(theta2 + theta3 * a_star + beta0 +
                                            beta1 * a_star + covariatesTerm))))

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
                 intmed_err = intmed_err, pie_err = pie_err, te_err = te_err,
                 cde_err_prop = cde_err_prop, intref_err_prop = intref_err_prop,
                 intmed_err_prop = intmed_err_prop, pie_err_prop = pie_err_prop,
                 overall_pm = overall_pm, overall_int = overall_int,
                 overall_pe = overall_pe)


      }

    } else if (model == "iorw") {

      if (yreg == "linear") {

        dir <- unname(coef["dir_coef"])

        tot <- unname(coef["tot_coef"])

        ind <- tot - dir

        out <- c(tot = tot, dir = dir, ind = ind)

      } else {

        ORdir <- unname(exp(coef["dir_coef"]))

        ORtot <- unname(exp(coef["tot_coef"]))

        ORind <- ORtot / ORdir

        out <- c(ORtot = ORtot, ORdir = ORdir, ORind = ORind)

      }

    }

  } else if (est.method == "imputation") {

    EY0m_sim <- EY1m_sim <- EY00_sim <- EY01_sim <- EY10_sim <- EY11_sim <-c()

    if (model %in% c("rb", "msm")) {

      mediator_regression <- regressions$mediator_regression

      outcome_regression <- regressions$outcome_regression

      if (model == "msm") cde_outcome_regression <- regressions$cde_outcome_regression

      mdesign_a <- data.frame(rep(a,nrow(data_boot)),data_boot[,covariates.pre])

      mdesign_a_star <- data.frame(rep(a_star,nrow(data_boot)),data_boot[,covariates.pre])

      colnames(mdesign_a) <- colnames(mdesign_a_star) <- c(exposure, covariates.pre)

      if (is.factor(data[, exposure])) {
        mdesign_a_star[,1] <- as.factor(mdesign_a_star[,1])
        mdesign_a[,1] <- as.factor(mdesign_a[,1])
      }

      m_a <- do.call(data.frame, lapply(1:length(mediator_regression),
                                        FUN = function(x) {
                                          if (inherits(mediator_regression[[x]], "lm")|
                                              inherits(mediator_regression[[x]], "glm")|
                                              inherits(mediator_regression[[x]], "survreg")|
                                              inherits(mediator_regression[[x]], "coxph")) {
                                            predict(mediator_regression[[x]],
                                                    newdata = mdesign_a,
                                                    type = ifelse(inherits(mediator_regression[[x]],
                                                                           "coxph"), "risk", "response"))

                                          } else if (inherits(mediator_regression[[x]], "multinom")) {
                                            predict(mediator_regression[[x]], newdata = mdesign_a, type = "class")
                                          }}))

      m_a_star <- do.call(data.frame, lapply(1:length(mediator_regression),
                                             FUN = function(x) {
                                               if (inherits(mediator_regression[[x]], "lm")|
                                                   inherits(mediator_regression[[x]], "glm")|
                                                   inherits(mediator_regression[[x]], "survreg")|
                                                   inherits(mediator_regression[[x]], "coxph")) {
                                                 predict(mediator_regression[[x]],
                                                         newdata = mdesign_a_star,
                                                         type = ifelse(inherits(mediator_regression[[x]],
                                                                                "coxph"), "risk", "response"))

                                               } else if (inherits(mediator_regression[[x]], "multinom")) {
                                                 predict(mediator_regression[[x]], newdata = mdesign_a_star, type = "class")
                                               }}))

      ydesign0m <- data.frame(rep(a_star, nrow(data_boot)), t(m_star)[rep(1,nrow(data_boot)),],
                              data_boot[,covariates.pre])

      ydesign1m <- data.frame(rep(a, nrow(data_boot)), t(m_star)[rep(1,nrow(data_boot)),],
                              data_boot[,covariates.pre])

      ydesign00 <- data.frame(rep(a_star, nrow(data_boot)), m_a_star, data_boot[,covariates.pre])

      ydesign01 <- data.frame(rep(a_star, nrow(data_boot)), m_a, data_boot[,covariates.pre])

      ydesign10 <- data.frame(rep(a, nrow(data_boot)), m_a_star, data_boot[,covariates.pre])

      ydesign11 <- data.frame(rep(a, nrow(data_boot)), m_a, data_boot[,covariates.pre])

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, covariates.pre)

      if (is.factor(data[, exposure])) {
        ydesign0m[,1] <- as.factor(ydesign0m[,1])
        ydesign1m[,1] <- as.factor(ydesign1m[,1])
        ydesign01[,1] <- as.factor(ydesign01[,1])
        ydesign10[,1] <- as.factor(ydesign10[,1])
        ydesign00[,1] <- as.factor(ydesign00[,1])
        ydesign11[,1] <- as.factor(ydesign11[,1])
      }

      for (i in 1:length(mediator)) {

        if (is.factor(data[, mediator[i]])) {
          ydesign0m[,1+i] <- as.factor(ydesign0m[,1+i])
          ydesign1m[,1+i] <- as.factor(ydesign1m[,1+i])
        }
      }

      if (model == "rb") {

        EY0m <- mean(predict(outcome_regression, newdata =  ydesign0m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

        EY1m <- mean(predict(outcome_regression, newdata =  ydesign1m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))


      } else if (model == "msm") {

        EY0m <- mean(predict(cde_outcome_regression, newdata =  ydesign0m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

        EY1m <- mean(predict(cde_outcome_regression, newdata =  ydesign1m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      }

      EY00 <- mean(predict(outcome_regression, newdata =  ydesign00,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY01 <- mean(predict(outcome_regression, newdata =  ydesign01,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY10 <- mean(predict(outcome_regression, newdata =  ydesign10,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY11 <- mean(predict(outcome_regression, newdata =  ydesign11,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

    } else if (model == "g-formula") {

      outcome_regression <- regressions$outcome_regression

      mediator_regression <- regressions$mediator_regression

      postcovar_regression <- regressions$postcovar_regression

      postcovardesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                                      data_boot[,covariates.pre])

      postcovardesign_a_star <- data.frame(c(rep(a_star,nrow(data_boot))),
                                           data_boot[,covariates.pre])

      mdesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                              data_boot[,covariates.pre],
                              data_boot[,covariates.post])

      mdesign_a_star <- data.frame(c(rep(a_star,nrow(data_boot))),
                                   data_boot[,covariates.pre],
                                   data_boot[,covariates.post])

      colnames(postcovardesign_a) <- colnames(postcovardesign_a_star) <-
        c(exposure, covariates.pre)

      colnames(mdesign_a) <- colnames(mdesign_a_star) <-
        c(exposure, covariates.pre, covariates.post)

      if (is.factor(data[, exposure])) {
        mdesign_a_star[,1] <- as.factor(mdesign_a_star[,1])
        mdesign_a[,1] <- as.factor(mdesign_a[,1])
        postcovardesign_a_star[,1] <- as.factor(postcovardesign_a_star[,1])
        postcovardesign_a[,1] <- as.factor(postcovardesign_a[,1])
      }

      m_a_star <- do.call(data.frame, lapply(1:length(mediator_regression),
                                             FUN = function(x) {
                                               if (inherits(mediator_regression[[x]], "lm")|
                                                   inherits(mediator_regression[[x]], "glm")|
                                                   inherits(mediator_regression[[x]], "survreg")|
                                                   inherits(mediator_regression[[x]], "coxph")) {
                                                 predict(mediator_regression[[x]],
                                                         newdata = mdesign_a_star,
                                                         type = ifelse(inherits(mediator_regression[[x]],
                                                                                "coxph"), "risk", "response"))

                                               } else if (inherits(mediator_regression[[x]], "multinom")) {
                                                 predict(mediator_regression[[x]], newdata = mdesign_a_star, type = "class")
                                               }}))

      m_a <- do.call(data.frame, lapply(1:length(mediator_regression),
                                        FUN = function(x) {
                                          if (inherits(mediator_regression[[x]], "lm")|
                                              inherits(mediator_regression[[x]], "glm")|
                                              inherits(mediator_regression[[x]], "survreg")|
                                              inherits(mediator_regression[[x]], "coxph")) {
                                            predict(mediator_regression[[x]],
                                                    newdata = mdesign_a,
                                                    type = ifelse(inherits(mediator_regression[[x]],
                                                                           "coxph"), "risk", "response"))

                                          } else if (inherits(mediator_regression[[x]], "multinom")) {
                                            predict(mediator_regression[[x]], newdata = mdesign_a, type = "class")
                                          }}))

      if (!is.null(covariates.post)){
        postcovar_a_star <- do.call(data.frame, lapply(1:length(postcovar_regression),
                                                       FUN = function(x) {
                                                         if (inherits(postcovar_regression[[x]], "lm")|
                                                             inherits(postcovar_regression[[x]], "glm")|
                                                             inherits(postcovar_regression[[x]], "survreg")|
                                                             inherits(postcovar_regression[[x]], "coxph")) {
                                                           predict(postcovar_regression[[x]],
                                                                   newdata = postcovardesign_a_star,
                                                                   type = ifelse(inherits(postcovar_regression[[x]],
                                                                                          "coxph"), "risk", "response"))

                                                         } else if (inherits(postcovar_regression[[x]], "multinom")) {
                                                           predict(postcovar_regression[[x]], newdata = postcovardesign_a_star, type = "class")
                                                         }}))

        postcovar_a <- do.call(data.frame, lapply(1:length(postcovar_regression),
                                                  FUN = function(x) {
                                                    if (inherits(postcovar_regression[[x]], "lm")|
                                                        inherits(postcovar_regression[[x]], "glm")|
                                                        inherits(postcovar_regression[[x]], "survreg")|
                                                        inherits(postcovar_regression[[x]], "coxph")) {
                                                      predict(postcovar_regression[[x]],
                                                              newdata = postcovardesign_a,
                                                              type = ifelse(inherits(postcovar_regression[[x]],
                                                                                     "coxph"), "risk", "response"))

                                                    } else if (inherits(postcovar_regression[[x]], "multinom")) {
                                                      predict(postcovar_regression[[x]], newdata = postcovardesign_a, type = "class")
                                                    }}))
      } else {
        postcovar_a_star <- postcovar_a <-data.frame()[1:nrow(data_boot), ]
      }

      ydesign0m <- data.frame(rep(a_star, nrow(data_boot)), t(m_star)[rep(1,nrow(data_boot)),],
                              data_boot[,covariates.pre], postcovar_a_star)

      ydesign1m <- data.frame(rep(a, nrow(data_boot)), t(m_star)[rep(1,nrow(data_boot)),],
                              data_boot[,covariates.pre], postcovar_a)

      ydesign00 <- data.frame(rep(a_star, nrow(data_boot)), m_a_star,
                              data_boot[,covariates.pre], postcovar_a_star)

      ydesign01 <- data.frame(rep(a_star, nrow(data_boot)), m_a,
                              data_boot[,covariates.pre], postcovar_a_star)

      ydesign10 <- data.frame(rep(a, nrow(data_boot)), m_a_star,
                              data_boot[,covariates.pre], postcovar_a)

      ydesign11 <- data.frame(rep(a, nrow(data_boot)), m_a,
                              data_boot[,covariates.pre], postcovar_a)

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, covariates.pre,
                                                        covariates.post)

      if (is.factor(data[, exposure])) {
        ydesign0m[,1] <- as.factor(ydesign0m[,1])
        ydesign1m[,1] <- as.factor(ydesign1m[,1])
        ydesign01[,1] <- as.factor(ydesign01[,1])
        ydesign10[,1] <- as.factor(ydesign10[,1])
        ydesign00[,1] <- as.factor(ydesign00[,1])
        ydesign11[,1] <- as.factor(ydesign11[,1])
      }

      for (i in 1:length(mediator)) {

        if (is.factor(data[, mediator[i]])) {
          ydesign0m[,1+i] <- as.factor(ydesign0m[,1+i])
          ydesign1m[,1+i] <- as.factor(ydesign1m[,1+i])
        }
      }

      EY0m <- mean(predict(outcome_regression, newdata =  ydesign0m,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY1m <- mean(predict(outcome_regression, newdata =  ydesign1m,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY00 <- mean(predict(outcome_regression, newdata =  ydesign00,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY01 <- mean(predict(outcome_regression, newdata =  ydesign01,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY10 <- mean(predict(outcome_regression, newdata =  ydesign10,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

      EY11 <- mean(predict(outcome_regression, newdata =  ydesign11,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")))

    } else if (model == "wb") {

      index_a_star <- which(data_boot[, exposure] == a_star)

      index_a <- which(data_boot[, exposure] == a)

      w.nom <- left_join(select(data_boot, exposure),
                         count(data_boot, !!as.name(exposure)),
                         by = exposure)[, "n"]/nrow(data_boot)

      exposure_regression <- regressions$exposure_regression

      outcome_regression <- regressions$outcome_regression

      if (exposure.type == "continuous") {

        w.denom <- dnorm(data_boot[, exposure],
                         mean = exposure_regression$fitted.values,
                         sd = summary(exposure_regression)$sigma)

      } else if (exposure.type == "binary") {

        w.denom <- (predict(exposure_regression, type = "response", data = data_boot)) ^
          (as.numeric(data_boot[, exposure])-1) * (1 - predict(exposure_regression, type = "response", data = data_boot)) ^
          (1 - (as.numeric(data_boot[, exposure])-1))

      }

      w <- w.nom/w.denom

      ydesign0m <- data.frame(rep(a_star, length(index_a_star)), t(m_star)[rep(1,length(index_a_star)),],
                              data_boot[index_a_star,covariates.pre])

      ydesign1m <- data.frame(rep(a, length(index_a)), t(m_star)[rep(1,length(index_a)),],
                              data_boot[index_a,covariates.pre])

      ydesign01 <- data.frame(rep(a_star, length(index_a)), data_boot[index_a, mediator],
                              data_boot[index_a, covariates.pre])

      ydesign10 <- data.frame(rep(a, length(index_a_star)), data_boot[index_a_star, mediator],
                              data_boot[index_a_star, covariates.pre])

      colnames(ydesign0m) <- colnames(ydesign1m) <-
        colnames(ydesign10) <- colnames(ydesign01) <- c(exposure, mediator, covariates.pre)


      for (i in 1:length(mediator)) {

        if (is.factor(data[, mediator[i]])) {
          ydesign0m[,1+i] <- as.factor(ydesign0m[,1+i])
          ydesign1m[,1+i] <- as.factor(ydesign1m[,1+i])
        }
      }

      if (is.factor(data[, exposure])) {
        ydesign0m[,1] <- as.factor(ydesign0m[,1])
        ydesign1m[,1] <- as.factor(ydesign1m[,1])
        ydesign01[,1] <- as.factor(ydesign01[,1])
        ydesign10[,1] <- as.factor(ydesign10[,1])
      }

      EY0m <- weighted.mean(predict(outcome_regression, newdata = ydesign0m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a_star])

      EY1m <- weighted.mean(predict(outcome_regression, newdata = ydesign1m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a])

      EY00 <- weighted.mean(data_boot[index_a_star, outcome], w = w[index_a_star])

      EY11 <- weighted.mean(data_boot[index_a, outcome], w = w[index_a])

      EY01 <- weighted.mean(predict(outcome_regression, newdata = ydesign01,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a])

      EY10 <- weighted.mean(predict(outcome_regression, newdata = ydesign10,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a_star])

    }


    if (yreg == "linear"|(!is.character(yreg)&&family(yreg)$Family == "gaussian")) {

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

    } else if (yreg != "linear"|(!is.character(yreg)&&family(yreg)$Family != "gaussian")) {

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
               intmed_err = intmed_err, pie_err = pie_err, te_err = te_err,
               cde_err_prop = cde_err_prop, intref_err_prop = intref_err_prop,
               intmed_err_prop = intmed_err_prop, pie_err_prop = pie_err_prop,
               overall_pm = overall_pm, overall_int = overall_int,
               overall_pe = overall_pe)

    }
  }

  return(out)

}
