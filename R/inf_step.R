inf_step <- function(nboot, data, model,
                     outcome, event, exposure, mediator, EMint, prec, postc,
                     yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                     astar, a, mval, yref, vecc,
                     estimation, inference,
                     ME = FALSE, MEvariable = NULL, MEvariable.type = NULL,
                     measurement.error = NULL, lambda = c(0.5, 1, 1.5, 2), B = 100) {



  if (inference == "delta") {

    formulas <- create_formulas(model = model,
                                outcome = outcome, event = event,
                                exposure = exposure, mediator = mediator, EMint = EMint,
                                prec = prec, postc = postc,
                                yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg)

    regressions <- run_regressions(formulas = formulas, data = data, model = model,
                                   exposure = exposure, mediator = mediator, postc = postc,
                                   yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                                   wmreg = wmreg)

    if (model != "rb") {
      stop("Delta method inference doesn't support selected model")
    } else if (model == "rb" && length(mediator) > 1) {
      stop("For the selected model, delta method inference only supports a single mediator")
    }

    if (estimation == "imputation") {
      stop("Delta method inference doesn't support direct imputation estimation")
    }

    if (is.character(data[, exposure])|is.factor(data[, exposure])) {
      a <- as.numeric(levels(as.factor(data[, exposure])) == a)[-1]
      astar <- as.numeric(levels(as.factor(data[, exposure])) == astar)[-1]
    }

    if (is.character(data[, mediator])|is.factor(data[, mediator])) {
      mstar <- as.numeric(levels(as.factor(data[, mediator])) == mval[[1]])[-1]
    } else {mstar <- mval[[1]]}

    elevel <- ifelse(is.character(data[, exposure])|is.factor(data[, exposure]),
                     length(levels(data[, exposure])), 2)

    mlevel <- ifelse(is.character(data[, mediator])|is.factor(data[, mediator]),
                     length(levels(data[, mediator])), 2)

    outcome_regression <- regressions$outcome_regression

    mediator_regression <- regressions$mediator_regression[[1]]

    if (ME && MEvariable %in% all.vars(formula(outcome_regression))) {

      outcome_SIMEXreg <- simex_reg(reg = outcome_regression, data = data,
                                    MEvariable = MEvariable, MEvariable.type = MEvariable.type,
                                    measurement.error = measurement.error, lambda = lambda,
                                    B = B)

      thetas <- outcome_SIMEXreg$SIMEXcoef

      vcov_thetas <- outcome_SIMEXreg$SIMEXvar

    } else {

      thetas <- coef(outcome_regression)

      vcov_thetas <- vcov(outcome_regression)

    }

    if (ME && MEvariable %in% all.vars(formula(mediator_regression))) {

      mediator_SIMEXreg <- simex_reg(reg = mediator_regression, data = data,
                                     MEvariable = MEvariable, MEvariable.type = MEvariable.type,
                                     measurement.error = measurement.error, lambda = lambda,
                                     B = B)

      betas <- mediator_SIMEXreg$SIMEXcoef

      if (mreg == "linear") { variance <- mediator_SIMEXreg$SIMEXsigma^2
      } else variance <- NULL

      vcov_betas <- mediator_SIMEXreg$SIMEXvar

    } else {

      betas  <- as.vector(t(coef(mediator_regression)))

      if (mreg == "linear") { variance <- sigma(mediator_regression)^2
      } else variance <- NULL

      vcov_betas <- vcov(mediator_regression)

    }

    vcov_block <- Matrix::bdiag(vcov_thetas, vcov_betas)

    theta0 <- "x1"

    theta1 <- paste0("x", 2:elevel)

    theta2 <- paste0("x", (elevel + 1):(elevel + mlevel - 1))

    if (EMint == TRUE) {
      theta3 <- t(matrix(paste0("x", length(thetas)-(((elevel-1)*(mlevel-1)-1):0)),
                         ncol = mlevel - 1))

    } else {theta3 <- t(matrix(rep(0, (elevel-1)*(mlevel-1)),
                               ncol = mlevel - 1))}

    beta0 <- paste0("x", length(thetas) + 1 + (0:(mlevel - 2))*length(betas)/(mlevel - 1))

    beta1 <- t(matrix(paste0("x", length(thetas) + rowSums(expand.grid(2:elevel,(0:(mlevel-2))*
                                                                         length(betas)/(mlevel-1)))),
                      ncol = mlevel - 1))

    XC <- sapply(0:(mlevel-2), function(x) paste0("x", length(thetas) + elevel +
                                                    x*length(betas)/(mlevel-1) +
                                                    1:length(vecc),
                                                  "*", paste0("(", vecc, ")"),
                                                  collapse = "+"))

    if (yreg == "linear") {

      if (mreg == "linear") {

        cde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                              paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                              paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                              paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*mstar")

        pnde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                               paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                               beta0, "+(", paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                               ")+", XC, ")")

        tnde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                               paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                               beta0, "+(", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"),
                               ")+", XC, ")")

        pnie_formula <- paste0("(", theta2, "+(", paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                               "))*((", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))")

        tnie_formula <- paste0("(", theta2, "+(", paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"),
                               "))*((", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))")

      } else if (mreg %in% c("logistic", "multinomial")) {

        cde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                              paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+(",
                              ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", a, ")"),
                                                               sep = "*", collapse = "+")), ")-(",
                              ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", astar, ")"),
                                                               sep = "*", collapse = "+")),
                              ")")

        pnde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                               paste0("((",sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                          sep = "*", collapse = "+")),
                                      ")-(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                            sep = "*", collapse = "+")),
                                      "))*exp(", beta0, "+(", sapply(1:(mlevel - 1),
                                                                     FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                             sep = "*", collapse = "+")),
                                      ")", "+", XC, ")", collapse = "+"), "))/(1+",
                               paste0("exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      ")+", XC, ")", collapse = "+"), ")")

        tnde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                               paste0("((",sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                          sep = "*", collapse = "+")),
                                      ")-(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                            sep = "*", collapse = "+")),
                                      "))*exp(", beta0, "+(", sapply(1:(mlevel - 1),
                                                                     FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                             sep = "*", collapse = "+")),
                                      ")", "+", XC, ")", collapse = "+"), "))/(1+",
                               paste0("exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      ")+", XC, ")", collapse = "+"), ")")

        pnie_formula <- paste0("(", paste0("(", theta2, "+(",
                                           sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                          sep = "*", collapse = "+")),
                                           "))*exp(", beta0, "+(",
                                           sapply(1:(mlevel - 1),
                                                  FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                          sep = "*", collapse = "+")),
                                           ")", "+", XC, ")", collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      ")+", XC, ")", collapse = "+"), ")-(",
                               paste0("(", theta2, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ],
                                                                                     paste0("(", astar, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      "))*exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"), sep = "*", collapse = "+")),
                                      ")", "+", XC, ")", collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      ")+", XC, ")", collapse = "+"), ")")

        tnie_formula <- paste0("(", paste0("(", theta2, "+(",
                                           sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                          sep = "*", collapse = "+")),
                                           "))*exp(", beta0, "+(",
                                           sapply(1:(mlevel - 1),
                                                  FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                          sep = "*", collapse = "+")),
                                           ")", "+", XC, ")", collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      ")+", XC, ")", collapse = "+"), ")-(",
                               paste0("(", theta2, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"), sep = "*", collapse = "+")),
                                      "))*exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"), sep = "*", collapse = "+")),
                                      ")", "+", XC, ")", collapse = "+"), ")/(1+",
                               paste0("exp(", beta0, "+(",
                                      sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                     sep = "*", collapse = "+")),
                                      ")+", XC, ")", collapse = "+"), ")")
      }

      te_formula <- paste0("(", tnie_formula, ")+(", pnde_formula, ")")

      pm_formula <- paste0("(", tnie_formula, ")/(", te_formula, ")")

      intref_formula <- paste0("(", pnde_formula, ")-(", cde_formula, ")")

      intmed_formula <- paste0("(", tnie_formula, ")-(", pnie_formula, ")")

      pie_formula <- pnie_formula

      cde_prop_formula <- paste0("(", cde_formula, ")/(", te_formula, ")")

      intref_prop_formula <- paste0("(", intref_formula, ")/(", te_formula, ")")

      intmed_prop_formula <- paste0("(", intmed_formula, ")/(", te_formula, ")")

      pie_prop_formula <- paste0("(", pie_formula, ")/(", te_formula, ")")

      overall_pm_formula <- paste0("((", pnie_formula, ")+(", intmed_formula, "))/(", te_formula, ")")

      overall_int_formula <- paste0("((", intref_formula, ")+(", intmed_formula, "))/(", te_formula, ")")

      overall_pe_formula <- paste0("((", intref_formula, ")+(", intmed_formula, ")+(", pnie_formula,
                                   "))/(", te_formula, ")")

      delta_formula <- list(cde_formula = cde_formula, pnde_formula = pnde_formula,
                            tnde_formula = tnde_formula, pnie_formula = pnie_formula,
                            tnie_formula = tnie_formula, te_formula = te_formula, pm_formula = pm_formula,
                            intref_formula = intref_formula, intmed_formula = intmed_formula, pie_formula = pie_formula,
                            cde_prop_formula = cde_prop_formula, intref_prop_formula = intref_prop_formula,
                            intmed_prop_formula = intmed_prop_formula, pie_prop_formula = pie_prop_formula,
                            overall_pm_formula = overall_pm_formula, overall_int_formula = overall_int_formula,
                            overall_pe_formula = overall_pe_formula)

    } else if  (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                            "negbin", "coxph", "aft_exp", "aft_weibull")) {

      if (mreg == "linear") {

        cde_rr_formula <- paste0("exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                 paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                                 paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                 paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*mstar)")


        pnde_rr_formula <- paste0("exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                                  paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                                  beta0, "+", paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                  "+", XC, "+", theta2, "*variance)+0.5*variance*((",
                                  paste(paste0(theta3, "^2"), paste0("(", a, ")"), sep = "*", collapse = "+"),
                                  ")-(", paste(paste0(theta3, "^2"), paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                  ")))")

        tnde_rr_formula <- paste0("exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                                  paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                                  beta0, "+", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"),
                                  "+", XC, "+", theta2, "*variance)+0.5*variance*((",
                                  paste(paste0(theta3, "^2"), paste0("(", a, ")"), sep = "*", collapse = "+"),
                                  ")-(", paste(paste0(theta3, "^2"), paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                  ")))")

        pnie_rr_formula <- paste0("exp(", theta2, "*((", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))+(",
                                  paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")*((",
                                  paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")))")

        tnie_rr_formula <- paste0("exp(", theta2, "*((", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))+(",
                                  paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")*((",
                                  paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")))")

        cde_err_formula <- paste0("(exp(((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))+(",
                                  paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")*mstar)-exp((",
                                  paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                  ")*mstar))*exp(", theta2, "*mstar-(", theta2, "+(",
                                  paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                                  beta0, "+(", paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                  ")+", XC, ")-0.5*(", theta2, "+(",
                                  paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))^2*variance)")

      } else if (mreg %in% c("logistic", "multinomial")) {

        cde_rr_formula <- paste0("exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                 paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+(",
                                 ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", a, ")"),
                                                                  sep = "*", collapse = "+"))
                                 , ")-(", ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", astar, ")"),
                                                                           sep = "*", collapse = "+")),
                                 "))")

        pnde_rr_formula <- paste0("(exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                            sep = "*", collapse = "+")),

                                         ")+", XC, ")", collapse = "+"), "))/(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                            sep = "*", collapse = "+")),

                                         ")+", XC, ")", collapse = "+"), ")")

        tnde_rr_formula <- paste0("(exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                                            sep = "*", collapse = "+")),

                                         ")+", XC, ")", collapse = "+"), "))/(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                                            sep = "*", collapse = "+")),

                                         ")+", XC, ")", collapse = "+"), ")")

        pnie_rr_formula <- paste0("((1+", paste0("exp(", beta0, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                                     sep = "*", collapse = "+")),
                                                 ")+", XC, ")", collapse = "+"), ")*(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                                            sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), "))/((1+",
                                  paste0("exp(", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), ")*(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                            sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), "))")

        tnie_rr_formula <- paste0("((1+", paste0("exp(", beta0, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                                     sep = "*", collapse = "+")),
                                                 ")+", XC, ")", collapse = "+"), ")*(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                                            sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), "))/((1+",
                                  paste0("exp(", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", a, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), ")*(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", a, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                            sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), "))")

        cde_err_formula <- paste0("(exp((", ifelse(sum(mstar) == 0, 0, theta2[which(mstar == 1)]), "))*(1+",
                                  paste0("exp(", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), ")*(exp((",
                                  paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                  paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+(",
                                  ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", a, ")"),
                                                                   sep = "*", collapse = "+")), "))-exp((",
                                  ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", astar, ")"),
                                                                   sep = "*", collapse = "+")), "))))/(1+",
                                  paste0("exp(", theta2, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(theta3[x, ], paste0("(", astar, ")"),
                                                                                                              sep = "*", collapse = "+")),
                                         ")+", beta0, "+(" , sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
                                                                                                            sep = "*", collapse = "+")),
                                         ")+", XC, ")", collapse = "+"), ")")

      }

      te_rr_formula <- paste0("(", tnie_rr_formula, ")*(", pnde_rr_formula, ")")

      pm_formula <- paste0("((", pnde_rr_formula, ")*(", tnie_rr_formula, "-1))/(", te_rr_formula, "-1)")

      intref_err_formula <- paste0("(", pnde_rr_formula, ")-1-(", cde_err_formula," )")

      intmed_err_formula <- paste0("(", tnie_rr_formula, ")*(", pnde_rr_formula, ")-(",
                                   pnde_rr_formula, ")-(", pnie_rr_formula, ")+1")

      pie_err_formula <- paste0("(", pnie_rr_formula, ")-1")

      te_err_formula <- paste0("(", te_rr_formula, ")-1")

      cde_err_prop_formula <- paste0("(", cde_err_formula, ")/(", te_err_formula, ")")

      intmed_err_prop_formula <- paste0("(", intmed_err_formula, ")/(", te_err_formula, ")")

      intref_err_prop_formula <-paste0("(", intref_err_formula, ")/(", te_err_formula, ")")

      pie_err_prop_formula <- paste0("(", pie_err_formula, ")/(", te_err_formula, ")")

      overall_pm_formula <- paste0("((", pie_err_formula, ")+(", intmed_err_formula, "))/(",
                                   te_err_formula, ")")

      overall_int_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, "))/(",
                                    te_err_formula, ")")

      overall_pe_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, ")+(",
                                   pie_err_formula,
                                   "))/(", te_err_formula, ")")

      delta_formula <- list(cde_rr_formula = cde_rr_formula, pnde_rr_formula = pnde_rr_formula,
                            tnde_rr_formula = tnde_rr_formula, pnie_rr_formula = pnie_rr_formula,
                            tnie_rr_formula = tnie_rr_formula, te_rr_formula = te_rr_formula,
                            pm_formula = pm_formula, cde_err_formula = cde_err_formula,
                            intref_err_formula = intref_err_formula, intmed_err_formula = intmed_err_formula,
                            pie_err_formula = pie_err_formula,
                            cde_err_prop_formula = cde_err_prop_formula,
                            intref_err_prop_formula = intref_err_prop_formula,
                            intmed_err_prop_formula = intmed_err_prop_formula, pie_err_prop_formula = pie_err_prop_formula,
                            overall_pm_formula = overall_pm_formula, overall_int_formula = overall_int_formula,
                            overall_pe_formula = overall_pe_formula)

    }

    effect_se <- c()

    for (formula in names(delta_formula)) {

      delta_formula[[formula]] <- as.formula(paste0("~", stringr::str_replace_all(
        delta_formula[[formula]],
        pattern = c("\\bmstar\\b" = paste0("(", mstar, ")"),
                    "\\bvariance\\b" = as.character(variance)))))

      effect_se <- c(effect_se, msm::deltamethod(delta_formula[[formula]], c(thetas, betas), vcov_block))

    }

  } else if (inference == "bootstrap") {

    boots <- boot::boot(data = data, statistic = est_step, R = nboot,
                        model = model,
                        outcome = outcome, event = event, exposure = exposure,
                        mediator = mediator, EMint = EMint,
                        prec = prec, postc = postc,
                        yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg, wmreg = wmreg,
                        astar = astar, a = a, mval = mval, yref = yref, vecc = vecc,
                        estimation = estimation,
                        ME = ME, MEvariable = MEvariable, MEvariable.type = MEvariable.type,
                        measurement.error = measurement.error, lambda = lambda, B = B)

    effect_se <- apply(boots$t, 2, sd)

  } else if (inference == "none") {

    effect_se <- NULL

  } else {

    stop("Unsupported inference method")

    }

  return(effect_se)

}
