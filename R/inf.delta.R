inf.delta <- function(data = NULL, yreg = NULL, mreg = NULL) {

  # for categorical exposure, create indicator vectors for a and astar
  if (is.factor(data[, exposure]) | is.character(data[, exposure])) {
    a_lev<- levels(droplevels(as.factor(data[, exposure])))
    a <- as.numeric(a_lev == a)[-1]
    astar <- as.numeric(a_lev == astar)[-1]
    elevel <- length(a_lev)
    rm(a_lev)
  } else if (is.numeric(data[, exposure]) | is.logical(data[, exposure])) {
    elevel <- 2
  }

  # create covariate values to calculate conditional causal effects
  vecc <- c()
  if (length(prec) != 0) for (i in 1:length(prec)) vecc <- c(vecc, precval[[i]])

  # for categorical mediator, create an indicator vector for mstar
  if (is.factor(data[, mediator]) | is.character(data[, mediator])) {
    m_lev <- levels(droplevels(as.factor(data[, mediator])))
    mstar <- as.numeric(m_lev == mval[[1]])[-1]
    mlevel <- length(m_lev)
    rm(m_lev)
  } else if (is.numeric(data[, mediator]) | is.logical(data[, mediator])) {
    mstar <- mval[[1]]
    mlevel <- 2
  }

  # coefficients for yreg
  thetas <- coef(yreg)
  # coefficients for mreg
  betas  <- as.vector(t(coef(mreg)))

  # variance-covariance matrix of thetas
  vcov_thetas <- vcov(yreg)
  # variance-covariance matrix of betas
  vcov_betas <- vcov(mreg)
  # stack the two diagonally
  vcov_block <- bdiag(vcov_thetas, vcov_betas)

  theta0 <- "x1"
  theta1 <- paste0("x", 2:elevel)
  theta2 <- paste0("x", (elevel + 1):(elevel + mlevel - 1))
  if (EMint == TRUE) {
    theta3 <- t(matrix(paste0("x", length(thetas)-(((elevel-1)*(mlevel-1)-1):0)),
                       ncol = mlevel - 1))
  } else {theta3 <- t(matrix(rep(0, (elevel-1)*(mlevel-1)),
                             ncol = mlevel - 1))}
  beta0 <- paste0("x", length(thetas) + 1 + (0:(mlevel - 2))*length(betas)/(mlevel - 1))
  beta1 <- t(matrix(paste0("x", length(thetas) + rowSums(expand.grid(2:elevel,(0:(mlevel-2))*length(betas)/(mlevel-1)))),
                    ncol = mlevel - 1))
  XC <- sapply(0:(mlevel-2), function(x) paste0("x", length(thetas) + elevel + x*length(betas)/(mlevel-1) +
                                                  1:length(vecc), "*", paste0("(", vecc, ")"), collapse = "+"))

  if ((is_lm_yreg | is_glm_yreg ) && family_yreg$family == "gaussian" && family_yreg$link == "identity") {

    if ((is_lm_mreg | is_glm_mreg) && family_mreg[[1]]$family == "gaussian" && family_mreg[[1]]$link == "identity") {

      # linear Y with linear M
      cde_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                            paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                            paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                            paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*", mstar)

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

    } else {

      # linear Y with categorical M
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

    if (full) {
      pm_formula <- paste0("(", tnie_formula, ")/(", te_formula, ")")
      intref_formula <- paste0("(", pnde_formula, ")-(", cde_formula, ")")
      intmed_formula <- paste0("(", tnie_formula, ")-(", pnie_formula, ")")
      cde_prop_formula <- paste0("(", cde_formula, ")/(", te_formula, ")")
      intref_prop_formula <- paste0("(", intref_formula, ")/(", te_formula, ")")
      intmed_prop_formula <- paste0("(", intmed_formula, ")/(", te_formula, ")")
      pnie_prop_formula <- paste0("(", pnie_formula, ")/(", te_formula, ")")
      int_formula <- paste0("((", intref_formula, ")+(", intmed_formula, "))/(", te_formula, ")")
      pe_formula <- paste0("((", intref_formula, ")+(", intmed_formula, ")+(", pnie_formula,
                                   "))/(", te_formula, ")")
      delta_formula <- list(cde_formula = cde_formula, pnde_formula = pnde_formula,
                            tnde_formula = tnde_formula, pnie_formula = pnie_formula,
                            tnie_formula = tnie_formula, te_formula = te_formula, pm_formula = pm_formula,
                            intref_formula = intref_formula, intmed_formula = intmed_formula, 
                            cde_prop_formula = cde_prop_formula, intref_prop_formula = intref_prop_formula,
                            intmed_prop_formula = intmed_prop_formula, pnie_prop_formula = pnie_prop_formula,
                            int_formula = int_formula,
                            pe_formula = pe_formula)
    } else delta_formula <- list(cde_formula = cde_formula, pnde_formula = pnde_formula,
                                 tnde_formula = tnde_formula, pnie_formula = pnie_formula,
                                 tnie_formula = tnie_formula, te_formula = te_formula)

  } else {

    if ((is_lm_mreg | is_glm_mreg) && family_mreg[[1]]$family == "gaussian" && family_mreg[[1]]$link == "identity") {

      # nonlinear Y with linear M
      variance <- sigma(mreg)^2

      cde_logrr_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                               paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*", mstar)

      pnde_logrr_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                                paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                                beta0, "+", paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                "+", XC, "+", theta2, "*", variance, ")+0.5*", variance, "*((",
                                paste(paste0(theta3, "^2"), paste0("(", a, ")"), sep = "*", collapse = "+"),
                                ")-(", paste(paste0(theta3, "^2"), paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                "))")

      tnde_logrr_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+((",
                                paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                                beta0, "+", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"),
                                "+", XC, "+", theta2, "*", variance, ")+0.5*", variance, "*((",
                                paste(paste0(theta3, "^2"), paste0("(", a, ")"), sep = "*", collapse = "+"),
                                ")-(", paste(paste0(theta3, "^2"), paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                "))")

      pnie_logrr_formula <- paste0(theta2, "*((", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))+(",
                                paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")*((",
                                paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))")

      tnie_logrr_formula <- paste0(theta2, "*((", paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))+(",
                                paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")*((",
                                paste(beta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))")

      if (full) cde_err_formula <- paste0("exp(((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                                          paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))+(",
                                          paste(theta3, paste0("(", a, ")"), sep = "*", collapse = "+"), ")*", mstar, "-exp((",
                                          paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                          ")*", mstar, "))*exp(", theta2, "*", mstar, "-(", theta2, "+(",
                                          paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))*(",
                                          beta0, "+(", paste(beta1, paste0("(", astar, ")"), sep = "*", collapse = "+"),
                                          ")+", XC, ")-0.5*(", theta2, "+(",
                                          paste(theta3, paste0("(", astar, ")"), sep = "*", collapse = "+"), "))^2*", variance, ")")

    } else {

      # nonlinear Y with categorical M
      cde_logrr_formula <- paste0("(", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
                               paste(theta1, paste0("(", astar, ")"), sep = "*", collapse = "+"), ")+(",
                               ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", a, ")"),
                                                                sep = "*", collapse = "+"))
                               , ")-(", ifelse(sum(mstar) == 0, 0, paste(theta3[which(mstar == 1),], paste0("(", astar, ")"),
                                                                         sep = "*", collapse = "+")),")")

      pnde_logrr_formula <- paste0("log((exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
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
                                       ")+", XC, ")", collapse = "+"), "))")

      tnde_logrr_formula <- paste0("log((exp((", paste(theta1, paste0("(", a, ")"), sep = "*", collapse = "+"), ")-(",
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
                                       ")+", XC, ")", collapse = "+"), "))")

      pnie_logrr_formula <- paste0("log(((1+", paste0("exp(", beta0, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
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
                                       ")+", XC, ")", collapse = "+"), ")))")

      tnie_logrr_formula <- paste0("log(((1+", paste0("exp(", beta0, "+(", sapply(1:(mlevel - 1), FUN = function(x) paste(beta1[x, ], paste0("(", astar, ")"),
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
                                       ")+", XC, ")", collapse = "+"), ")))")

      if (full) cde_err_formula <- paste0("(exp((", ifelse(sum(mstar) == 0, 0, theta2[which(mstar == 1)]), "))*(1+",
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

    te_logrr_formula <- paste0("(", tnie_logrr_formula, ") + (", pnde_logrr_formula, ")")

    if (full) {
      pm_formula <- paste0("exp(", pnde_logrr_formula, ")*(exp(", tnie_logrr_formula, ")-1)/(exp(", te_logrr_formula, ")-1)")
      intref_err_formula <- paste0("exp(", pnde_logrr_formula, ")-1-(", cde_err_formula," )")
      intmed_err_formula <- paste0("exp(", tnie_logrr_formula, ")*exp(", pnde_logrr_formula, ")-exp(",
                                   pnde_logrr_formula, ")-exp(", pnie_logrr_formula, ")+1")
      pnie_err_formula <- paste0("exp(", pnie_logrr_formula, ")-1")
      te_err_formula <- paste0("exp(", te_logrr_formula, ")-1")
      cde_err_prop_formula <- paste0("(", cde_err_formula, ")/(", te_err_formula, ")")
      intmed_err_prop_formula <- paste0("(", intmed_err_formula, ")/(", te_err_formula, ")")
      intref_err_prop_formula <-paste0("(", intref_err_formula, ")/(", te_err_formula, ")")
      pnie_err_prop_formula <- paste0("(", pnie_err_formula, ")/(", te_err_formula, ")")
      int_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, "))/(",
                                    te_err_formula, ")")
      pe_formula <- paste0("((", intref_err_formula, ")+(", intmed_err_formula, ")+(",
                                   pnie_err_formula,
                                   "))/(", te_err_formula, ")")
      delta_formula <- list(cde_logrr_formula = cde_logrr_formula, pnde_logrr_formula = pnde_logrr_formula,
                            tnde_logrr_formula = tnde_logrr_formula, pnie_logrr_formula = pnie_logrr_formula,
                            tnie_logrr_formula = tnie_logrr_formula, te_logrr_formula = te_logrr_formula,
                            pm_formula = pm_formula, cde_err_formula = cde_err_formula,
                            intref_err_formula = intref_err_formula, intmed_err_formula = intmed_err_formula,
                            pnie_err_formula = pnie_err_formula,
                            cde_err_prop_formula = cde_err_prop_formula,
                            intref_err_prop_formula = intref_err_prop_formula,
                            intmed_err_prop_formula = intmed_err_prop_formula, pnie_err_prop_formula = pnie_err_prop_formula,
                            int_formula = int_formula,
                            pe_formula = pe_formula)
    } else delta_formula <- list(cde_logrr_formula = cde_logrr_formula, pnde_logrr_formula = pnde_logrr_formula,
                                 tnde_logrr_formula = tnde_logrr_formula, pnie_logrr_formula = pnie_logrr_formula,
                                 tnie_logrr_formula = tnie_logrr_formula, te_logrr_formula = te_logrr_formula)

  }

  effect_se <- c()
  for (formula in names(delta_formula)) {
    delta_formula[[formula]] <- as.formula(paste0("~", delta_formula[[formula]]))
    effect_se <- c(effect_se, deltamethod(delta_formula[[formula]], c(thetas, betas), vcov_block))
  }

  out <- effect_se
  return(out)

  }
