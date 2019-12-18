bootstrap_step <- function(data, indices, outcome, exposure, mediator, covariates, vecc, model,
                           EMint = FALSE, MMint = FALSE, EMMint = FALSE, EMint.terms = NULL,
                           MMint.terms = NULL, EMMint.terms = NULL, weights,
                           event, mreg, yreg, m_star, a_star, a) {

  data_boot <- data[indices, ]

  formulas <- create_formulas(model = model, outcome = outcome, exposure = exposure,
                              mediator = mediator, covariates = covariates,
                              EMint = EMint, MMint = MMint, EMMint = EMMint,
                              EMint.terms = EMint.terms, MMint.terms = MMint.terms,
                              EMMint.terms = EMMint.terms,
                              event = event, mreg = mreg, yreg = yreg)

  regressions <- run_regressions(model = model, formulas = formulas,
                                 mreg = mreg, yreg = yreg, data = data_boot)

###########################################Regression-based Approach#############################

  if (model == "rb") {

    coef <- get_coef(regressions = regressions, mreg = mreg, yreg = yreg)

    thetas <- coef$thetas

    betas <- unlist(coef$betas)

    variance <- unlist(coef$variance)

    covariatesTerm <- ifelse(!is.null(vecc), sum(betas[3:(2+length(vecc))]*t(vecc)), 0)

    if (EMint) { EMintTerm <- thetas[stringr::str_replace_all(EMint.terms,
                                                              pattern = "[*]", replacement = ":")]
    } else if (EMint == FALSE) EMintTerm <- 0

    if (yreg == "linear") {

      cde <- tnde <- pnde <- unname(thetas[exposure] * (a - a_star))

      tnie <- pnie <-  ifelse(mreg == "linear", betas[exposure] * thetas[mediator] *
                                (a - a_star), (exp(betas[1] + betas[exposure] * a +
                                                     covariatesTerm) /
                                                 (1 + exp(betas[1] + betas[exposure] *
                                                            a + covariatesTerm)) -
           exp(betas[1] + betas[exposure] * a_star + covariatesTerm) /
           (1 + exp(betas[1] + betas[exposure] * a_star + covariatesTerm))) *
          thetas[mediator])

      if (EMint == TRUE) {

        cde <- cde + EMintTerm * m_star * (a - a_star)

        pnde <- pnde + ifelse(mreg == "linear", EMintTerm *
                                (betas[1] + betas[exposure] * a_star +
                                   covariatesTerm) * (a - a_star), EMintTerm *
                                (exp(betas[1] + betas[exposure] * a_star +
                                  covariatesTerm) / (1 + exp(betas[1] +
                                  betas[exposure] * a_star + covariatesTerm))) * (a - a_star))

        tnde <- tnde + ifelse(mreg == "linear", EMintTerm *
                                (betas[1] + betas[exposure] * a +
                                   covariatesTerm) * (a - a_star), EMintTerm *
                                (exp(betas[1] + betas[exposure] * a +
                                       covariatesTerm) / (1 + exp(betas[1] +
                                                                    betas[exposure] * a + covariatesTerm))) * (a - a_star))

        tnie <- tnie + ifelse(mreg == "linear", EMintTerm * betas[exposure] *
                                a * (a - a_star), EMintTerm * (exp(betas[1] + betas[exposure] * a +
                                  covariatesTerm) / (1 + exp(betas[1] + betas[exposure] * a +
                                                               covariatesTerm)) -
                              exp(betas[1] + betas[exposure] * a_star +
                                    covariatesTerm) / (1 + exp(betas[1] + betas[exposure] * a_star +
                                                                                   covariatesTerm))) * a)

        pnie <- pnie + ifelse(mreg == "linear", EMintTerm * betas[exposure] *
                                a_star * (a - a_star), EMintTerm * (exp(betas[1] + betas[exposure] * a +
                                                                     covariatesTerm) / (1 + exp(betas[1] + betas[exposure] * a +
                                                                                                  covariatesTerm)) -
                                                                 exp(betas[1] + betas[exposure] * a_star +
                                                                       covariatesTerm) / (1 + exp(betas[1] + betas[exposure] * a_star +
                                                                                                    covariatesTerm))) * a_star)
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

      out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie,
               intref = intref, intmed = intmed, pie = pie,
               te = te, pm = pm, cde_prop = cde_prop, intref_prop = intref_prop,
               intmed_prop = intmed_prop, pie_prop = pie_prop,
               overall_pm = overall_pm, overall_int = overall_int, overall_pe = overall_pe)

      } else if (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                    "negbin", "coxph", "aft_exp", "aft_weibull")) {

        if (mreg == "linear") {

          cde_rr <- unname(exp((thetas[exposure] + EMintTerm * m_star) * (a - a_star)))

          pnde_rr <- unname(exp((thetas[exposure] + EMintTerm * (betas[1] + betas[exposure] * a_star +
                                                                   covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                                  0.5 * EMintTerm ^ 2 * variance * (a ^ 2 - a_star ^ 2)))

          tnde_rr <- unname(exp((thetas[exposure] + EMintTerm * (betas[1] + betas[exposure] * a +
                                                                   covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                                  0.5 * EMintTerm ^ 2 * variance * (a ^ 2 - a_star ^ 2)))

          pnie_rr <- unname(exp((thetas[mediator] * betas[exposure] +
                                   EMintTerm * betas[exposure] * a_star) * (a - a_star)))

          tnie_rr <- unname(exp((thetas[mediator] * betas[exposure] +
                                   EMintTerm * betas[exposure] * a) * (a - a_star)))

          cde_err <- unname(exp(thetas[exposure]*(a-a_star)+thetas[mediator]*m_star+
                                  EMintTerm*a*m_star- (thetas[mediator]+EMintTerm*a_star)*
                                  (betas[1]+betas[exposure]*a_star+
                                     covariatesTerm)-
                                  0.5*(thetas[mediator]+EMintTerm*a_star)^2*variance)-
                              exp(thetas[mediator]*m_star+EMintTerm*a_star*m_star-
                                    (thetas[mediator]+EMintTerm*a_star)*(betas[1]+
                                                                           betas[exposure]*a_star+covariatesTerm)-
                                    0.5*(thetas[mediator]+EMintTerm*a_star)^2*variance))

        }

        if (mreg == "logistic") {

          cde_rr <- unname(exp((thetas[exposure] + EMintTerm*m_star) * (a - a_star)))

          pnde_rr <- unname((exp(thetas[exposure] * (a - a_star)) * (1 + exp(thetas[mediator] +
                                                                               EMintTerm * a + betas[1] + betas[exposure] * a_star +
                                                                               covariatesTerm))) /(1 + exp(thetas[mediator] + EMintTerm * a_star +
                                                                                                             betas[1] + betas[exposure] * a_star + covariatesTerm)))

          tnde_rr <- unname((exp(thetas[exposure] * (a - a_star)) * (1 + exp(thetas[mediator] +
                                                                               EMintTerm * a + betas[1] + betas[exposure] * a + covariatesTerm))) /
                              (1 + exp(thetas[mediator] + EMintTerm * a_star +  betas[1] +
                                         betas[exposure] * a + covariatesTerm)))

          pnie_rr <- unname(((1 + exp(betas[1] + betas[exposure] * a_star + covariatesTerm)) *
                               (1 + exp(thetas[mediator] + EMintTerm * a_star + betas[1] +
                                          betas[exposure] * a + covariatesTerm))) /
                              ((1 + exp(betas[1] + betas[exposure] * a + covariatesTerm)) *
                                 (1+ exp(thetas[mediator] + EMintTerm * a_star + betas[1] +
                                           betas[exposure] * a_star + covariatesTerm))))

          tnie_rr <- unname(((1 + exp(betas[1] + betas[exposure] * a_star + covariatesTerm)) *
                               (1 + exp(thetas[mediator] + EMintTerm * a + betas[1] + betas[exposure] * a +
                                          covariatesTerm))) / ((1 + exp(betas[1] + betas[exposure] * a + covariatesTerm)) *
                                                                 (1 + exp(thetas[mediator] + EMintTerm * a + betas[1] + betas[exposure] * a_star +
                                                                            covariatesTerm))))

          cde_err <- unname((exp(thetas[exposure]*(a-a_star)+thetas[mediator]*m_star+
                                   EMintTerm*a*m_star)*(1+exp(betas[1]+
                                                                betas[exposure]*a_star+covariatesTerm))/(1+exp(betas[1]+
                                                                                                                      betas[exposure]*a_star+covariatesTerm+thetas[mediator]+
                                                                                                                      EMintTerm*a_star))-exp(thetas[mediator]*m_star+
                                                                                                                                               EMintTerm*a_star*m_star)*(1+exp(betas[1]+
                                                                                                                                                                                 betas[exposure]*a_star+covariatesTerm))/(1+exp(betas[1]+
                                                                                                                                                                                                                                       betas[exposure]*a_star+covariatesTerm+thetas[mediator]+
                                                                                                                                                                                                                                       EMintTerm*a_star))))

        }

        intref_err <- pnde_rr - 1 - cde_err

        intmed_err <- tnie_rr * pnde_rr - pnde_rr - pnie_rr + 1

        pie_err <- pnie_rr - 1

        te_rr <- tnie_rr * pnde_rr

        pm <- (pnde_rr * (tnie_rr - 1)) / (pnde_rr * tnie_rr - 1)

        total_err <- te_rr - 1

        cde_err_prop <- cde_err/total_err

        intmed_err_prop <- intmed_err/total_err

        intref_err_prop <- intref_err/total_err

        pie_err_prop <- pie_err/total_err

        overall_pm <- (pie_err+intmed_err)/total_err

        overall_int <- (intref_err+intmed_err)/total_err

        overall_pe <- (intref_err+intmed_err+pie_err)/total_err

        out <- c(cde_rr = cde_rr, pnde_rr = pnde_rr, tnde_rr = tnde_rr, pnie_rr = pnie_rr,
                 tnie_rr = tnie_rr, te_rr = te_rr, pm = pm,
                 cde_err = cde_err, intref_err = intref_err,
                 intmed_err = intmed_err, pie_err = pie_err, total_err = total_err,
                 cde_err_prop = cde_err_prop, intref_err_prop = intref_err_prop,
                 intmed_err_prop = intmed_err_prop, pie_err_prop = pie_err_prop,
                 overall_pm = overall_pm, overall_int = overall_int,
                 overall_pe = overall_pe)


      }
    }

###########################################Weighting-based Approach##############################

  if (model == "wb") {

    prob = predict.glm(exposure_regression, type = "response",
                       newdata = setNames(as.data.frame(data_boot[, covariates]), nm = covariates))

    w <- 1/dbinom(x = data_boot[, exposure], size = 1, prob = prob)

    data_cde <- data_boot[1,]
    data_cde[, exposure] <- a
    data_cde[, mediator] <- m_star
    EY1m <- as.numeric(predict(regressions$outcome_regression, newdata = data_cde,
                               type = ifelse(yreg != "coxph", "response", "risk")))

    data_cde[, exposure] <- a_star
    EY0m <- as.numeric(predict(regressions$outcome_regression, newdata = data_cde,
                               type = ifelse(yreg != "coxph", "response", "risk")))

    Yfitted <- predict(regressions$outcome_regression, newdata = data_boot,
                       type = ifelse(yreg != "coxph", "response", "risk"))

    data_cf <- data_boot
    data_cf[, exposure] <- 1 - data_boot[, exposure]
    Ycf <- predict(regressions$outcome_regression, newdata = data_cf,
                   type = ifelse(yreg != "coxph", "response", "risk"))

    EY00 <- weighted.mean(x = Yfitted[which(data_boot[, exposure] == a_star)],
                          w = w[which(data_boot[, exposure] == a_star)])
    EY11 <- weighted.mean(x = Yfitted[which(data_boot[, exposure] == a)],
                          w = w[which(data_boot[, exposure] == a)])
    EY10 <- weighted.mean(x = Ycf[which(data_boot[, exposure] == a_star)],
                          w = w[which(data_boot[, exposure] == a_star)])
    EY01 <- weighted.mean(x = Ycf[which(data_boot[, exposure] == a)],
                          w = w[which(data_boot[, exposure] == a)])

    if (yreg == "linear") {

      cde <- EY1m - EY0m

      pnde <- EY10 - EY00

      tnie <- EY11 - EY10

      tnde <- EY11 - EY01

      pnie <- EY01 - EY00

      te <- tnie + pnde

      pm <- tnie / (pnde + te)

      out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie, te = te, pm = pm)

    } else if (yreg != "linear") {

      ORcde <- EY1m / EY0m

      ORpnde <- EY10 / EY00

      ORtnie <- EY11 / EY10

      ORtnde <- EY11 / EY01

      ORpnie <- EY01 / EY00

      ORte <- ORtnie * ORpnde

      pm <- (ORpnde * (ORtnie - 1)) / (ORpnde * ORtnie - 1)

      out <- c(ORcde = ORcde, ORpnde = ORpnde, ORtnde = ORtnde, ORpnie = ORpnie, ORtnie = ORtnie,
               ORte = ORte, pm = pm)

      }
    }

####################################Inverse Odds Ratio Weighting Approach########################

  if (model == "iorw") {

      tot_outcome_regression <-  regressions$tot_outcome_regression
      dir_outcome_regression <-  regressions$dir_outcome_regression

      if (yreg == "linear") {

        dir <- unname(coef(dir_outcome_regression)[exposure])

        tot <- unname(coef(tot_outcome_regression)[exposure])

        ind <- tot - dir

        out <- c(tot = tot, dir = dir, ind = ind)

      } else {

        ORdir <- unname(exp(coef(dir_outcome_regression)[exposure]))

        ORtot <- unname(exp(coef(tot_outcome_regression)[exposure]))

        ORind <- ORtot / ORdir

        out <- c(ORtot = ORtot, ORdir = ORdir, ORind = ORind)

      }

    }

  return(out)

}
