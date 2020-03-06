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
                              event = event, mreg = mreg, yreg = yreg, data = data_boot)

  regressions <- run_regressions(model = model, formulas = formulas,
                                 exposure = exposure,  exposure.type = exposure.type,
                                 mediator = mediator, covariates.post = covariates.post,
                                 covariates.post.type = covariates.post.type,
                                 mreg = mreg, yreg = yreg, data = data_boot)

  if (est.method == "paramfunc") {

    if (!(model %in% c("rb", "iorw"))) {
      stop("Closed-form parameter function estimation doesn't support the selected model")
    } else if (length(mediator) > 1&&model != "iorw") {
      stop("For the selected model, closed-form parameter function estimation doesn't support multiple mediator cases")
    }

    mediator_regression <- regressions$mediator_regression
    outcome_regression <- regressions$outcome_regression

    if (model == "rb" && !(((inherits(mediator_regression[[1]], "glm") |
                             inherits(mediator_regression[[1]], "lm"))&&
           family(mediator_regression[[1]])$family %in% c("gaussian", "binomial")) |
        inherits(mediator_regression[[1]], "multinom"))) {
      stop("Closed-form parameter function estimation doesn't support the selected mediator regression model.")
    }

    coef <- get_coef(formulas = formulas, regressions = regressions,
                     mreg = mreg, yreg = yreg, model = model)

    if (model == "rb"){

      if (!(exposure.type %in% c("continuous", "binary"))) {
        stop("For the selected model, closed-form parameter function estimation only supports continuous or binary exposure")
      }

      mlevel <- ifelse(is.factor(data_boot[, mediator]), length(levels(data_boot[, mediator])), 2)

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

      if ((inherits(outcome_regression, "glm")|inherits(outcome_regression, "lm"))&&
          family(outcome_regression)$family == "gaussian") {

        if ((inherits(mediator_regression[[1]], "glm")|inherits(mediator_regression[[1]], "lm"))&&
            family(mediator_regression[[1]])$family == "gaussian") {

          cde <- unname((theta1 + theta3 * m_star) * (a-a_star))

          pnde <- unname((theta1 + theta3 * (beta0 + beta1 * a_star +
                                               covariatesTerm)) * (a - a_star))

          tnde <- unname((theta1 + theta3 * (beta0 + beta1 * a +
                                               covariatesTerm)) * (a - a_star))

          pnie <- unname((theta2 + theta3  * a_star) * beta1 * (a - a_star))

          tnie <- unname((theta2 + theta3  * a) * beta1 * (a - a_star))

        } else if ((inherits(mediator_regression[[1]], "glm")&&
                   family(mediator_regression[[1]])$family == "binomial") |
                   inherits(mediator_regression[[1]], "multinom")) {

          cde <- unname((theta1 + sum(theta3*m_star)) * (a-a_star))

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

      } else if (!((inherits(outcome_regression, "glm")|inherits(outcome_regression, "lm"))&&
                 family(outcome_regression)$family == "gaussian")) {

        if ((inherits(mediator_regression[[1]], "glm")|inherits(mediator_regression[[1]], "lm"))&&
              family(mediator_regression[[1]])$family == "gaussian") {

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

        } else if ((inherits(mediator_regression[[1]], "glm")&&
                    family(mediator_regression[[1]])$family == "binomial") |
                   inherits(mediator_regression[[1]], "multinom")) {

          cde_rr <- unname(exp((theta1 + sum(theta3*m_star)) * (a-a_star)))

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

          cde_err <- unname(exp(sum(theta2*m_star)) *
                              (exp(theta1 * (a - a_star) + a * sum(theta3*m_star)) -
                                 exp(a_star * sum(theta3*m_star))) *
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

      if ((inherits(outcome_regression, "glm")|inherits(outcome_regression, "lm"))&&
          family(regressions$outcome_regression)$family == "gaussian") {

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

    if (model %in% c("rb", "msm")) {

      mediator_regression <- regressions$mediator_regression

      outcome_regression <- regressions$outcome_regression

      if (model == "msm") cde_outcome_regression <- regressions$cde_outcome_regression

      mdesign_a <- data.frame(rep(a,nrow(data_boot)),data_boot[,covariates.pre])

      mdesign_a_star <- data.frame(rep(a_star,nrow(data_boot)),data_boot[,covariates.pre])

      colnames(mdesign_a) <- colnames(mdesign_a_star) <- c(exposure, covariates.pre)

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
                                            predict(mediator_regression[[x]], newdata = mdesign_a, type = "prob")[, -1]
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
                                                 predict(mediator_regression[[x]], newdata = mdesign_a_star, type = "prob")[, -1]
                                               }}))

      ydesign0m <- data.frame(rep(a_star, nrow(data_boot)), t(m_star)[rep(1,nrow(data_boot)),],
                              data_boot[,covariates.pre])

      ydesign1m <- data.frame(rep(a, nrow(data_boot)), t(m_star)[rep(1,nrow(data_boot)),],
                              data_boot[,covariates.pre])

      ydesign00 <- data.frame(rep(a_star, nrow(data_boot)), m_a_star, data_boot[,covariates.pre])

      ydesign01 <- data.frame(rep(a_star, nrow(data_boot)), m_a, data_boot[,covariates.pre])

      ydesign10 <- data.frame(rep(a, nrow(data_boot)), m_a_star, data_boot[,covariates.pre])

      ydesign11 <- data.frame(rep(a, nrow(data_boot)), m_a, data_boot[,covariates.pre])

      mname <- c()

      for (i in 1:length(mediator)) {

        if (inherits(mediator_regression[[i]], "multinom")) {
          mname <- c(mname, paste0(mediator[i], 1:(length(levels(data[, mediator[i]])) - 1)))
        } else mname <- c(mname, mediator[i])

      }

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <-
        colnames(ydesign01) <- colnames(ydesign10) <- colnames(ydesign11) <-
        c(exposure, mname, covariates.pre)

      if (model == "rb") {

        EY0m <- mean(predict(outcome_regression, newdata =  ydesign0m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                     na.rm = TRUE)

        EY1m <- mean(predict(outcome_regression, newdata =  ydesign1m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                     na.rm = TRUE)


      } else if (model == "msm") {

        EY0m <- mean(predict(cde_outcome_regression, newdata =  ydesign0m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                     na.rm = TRUE)

        EY1m <- mean(predict(cde_outcome_regression, newdata =  ydesign1m,
                             type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                     na.rm = TRUE)

      }

      EY00 <- mean(predict(outcome_regression, newdata =  ydesign00,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY01 <- mean(predict(outcome_regression, newdata =  ydesign01,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY10 <- mean(predict(outcome_regression, newdata =  ydesign10,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY11 <- mean(predict(outcome_regression, newdata =  ydesign11,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

    } else if (model == "g-formula") {

      outcome_regression <- regressions$outcome_regression

      mediator_regression <- regressions$mediator_regression

      if (!is.null(covariates.post)){

        postcovar_regression <- regressions$postcovar_regression

        postcovardesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                                        data_boot[,covariates.pre])

        postcovardesign_a_star <- data.frame(c(rep(a_star,nrow(data_boot))),
                                             data_boot[,covariates.pre])

        colnames(postcovardesign_a) <- colnames(postcovardesign_a_star) <-
          c(exposure, covariates.pre)

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
                                                           predict(postcovar_regression[[x]], newdata = postcovardesign_a_star, type = "prob")[, -1]
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
                                                      predict(postcovar_regression[[x]], newdata = postcovardesign_a, type = "prob")[ ,-1]
                                                    }}))
      } else {
        postcovar_a_star <- postcovar_a <-data.frame()[1:nrow(data_boot), ]
      }


      mdesign_a <- data.frame(c(rep(a,nrow(data_boot))),
                              data_boot[,covariates.pre],
                              postcovar_a)

      mdesign_a_star <- data.frame(c(rep(a_star,nrow(data_boot))),
                                   data_boot[,covariates.pre],
                                   postcovar_a_star)

      pcname <- c()

      if (!is.null(covariates.post)){

        for (i in 1:length(covariates.post)) {

          if (inherits(postcovar_regression[[i]], "multinom")) {
            pcname <- c(pcname, paste0(covariates.post[i], 1:(length(levels(data[, covariates.post[i]])) - 1)))
          } else pcname <- c(pcname, covariates.post[i])

        }
      }

      colnames(mdesign_a) <- colnames(mdesign_a_star) <-
        c(exposure, covariates.pre, pcname)

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
                                                 predict(mediator_regression[[x]], newdata = mdesign_a_star, type = "prob")[, -1]
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
                                            predict(mediator_regression[[x]], newdata = mdesign_a, type = "prob")[, -1]
                                          }}))



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

      mname <- c()

      for (i in 1:length(mediator)) {

        if (inherits(mediator_regression[[i]], "multinom")) {
          mname <- c(mname, paste0(mediator[i], 1:(length(levels(data[, mediator[i]])) - 1)))
        } else mname <- c(mname, mediator[i])

      }

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mname, covariates.pre,
                                                        pcname)

      EY0m <- mean(predict(outcome_regression, newdata =  ydesign0m,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY1m <- mean(predict(outcome_regression, newdata =  ydesign1m,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY00 <- mean(predict(outcome_regression, newdata =  ydesign00,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY01 <- mean(predict(outcome_regression, newdata =  ydesign01,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY10 <- mean(predict(outcome_regression, newdata =  ydesign10,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

      EY11 <- mean(predict(outcome_regression, newdata =  ydesign11,
                           type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                   na.rm = TRUE)

    } else if (model == "wb") {

      index_a_star <- which(data_boot[, exposure] == a_star)

      index_a <- which(data_boot[, exposure] == a)

      w.nom <- left_join(select(data_boot, exposure),
                         count(data_boot, !!as.name(exposure)),
                         by = exposure)[, "n"]/nrow(data_boot)

      exposure_regression <- update(regressions$exposure_regression,
                                    formula. = regressions$exposure_regression$formula,
                                    data = data_boot)

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

      EY0m <- weighted.mean(predict(outcome_regression, newdata = ydesign0m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a_star],
                            na.rm = TRUE)

      EY1m <- weighted.mean(predict(outcome_regression, newdata = ydesign1m,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a], na.rm = TRUE)

      EY00 <- weighted.mean(data_boot[index_a_star, outcome], w = w[index_a_star], na.rm = TRUE)

      EY11 <- weighted.mean(data_boot[index_a, outcome], w = w[index_a], na.rm = TRUE)

      EY01 <- weighted.mean(predict(outcome_regression, newdata = ydesign01,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a], na.rm = TRUE)

      EY10 <- weighted.mean(predict(outcome_regression, newdata = ydesign10,
                                    type = ifelse(inherits(outcome_regression, "coxph"), "risk", "response")),
                            w = w[index_a_star], na.rm = TRUE)

    }


    if (family(regressions$outcome_regression)$family == "gaussian") {

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

    } else if (family(regressions$outcome_regression)$family != "gaussian") {

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
