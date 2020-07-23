est.rb <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {

  data_orig <- data
  if (is.null(indices)) indices <- 1:n
  # resample data
  data <- data[indices, ]

  # for case control study
  # method 1: weight subjects with y=1 by yprevalence/p(y=1) and weight subjects with y=0 by (1-yprevalence)/p(y=0)
  # method 2: fit yreg with all data and fit other regs on data from controls
  # use method 1 when yprevalence is provided
  # when yprevalence is not provided but the outcome is rare, use method 2
  if (casecontrol && !is.null(yprevalence)) {

    # method 1 for a case control design
    prob1 <- mean(data[, outcome] == y_case)
    w4casecon <- ifelse(data[, outcome] == y_case, yprevalence / prob1, (1 - yprevalence) / (1 - prob1))

    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon
    if (is.null(weights_yreg)) weights_yreg <- w4casecon
    # update yreg
    if (inference == "delta") {
      call_yreg$design <- eval(bquote(survey::svydesign(~1, weights = ~.(weights_yreg), data = .(data))))
    } else {
      call_yreg$weights <- weights_yreg
      call_yreg$data <- data
    }
    yreg <- eval.parent(call_yreg)

    for (p in 1:length(mediator)) {
      # weights for mreg[[p]]
      if (!is.null(weights_mreg[[p]])) weights_mreg <- weights_mreg[[p]][indices] * w4casecon
      if (is.null(weights_mreg[[p]])) weights_mreg <- w4casecon
      # update mreg[[p]]
      if (inference == "delta" && !is_svymultinom_mreg[p]) {
        call_mreg[[p]]$design <- eval(bquote(survey::svydesign(~1, weights = ~.(weights_mreg), data = .(data))))
      } else {
        call_mreg[[p]]$weights <- weights_mreg
        call_mreg[[p]]$data <- data
      }
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }

    rm(prob1, w4casecon)

  } else if (casecontrol && yrare) {

    # method 2 for a case control design
    # data from controls
    control_indices <- which(data[, outcome] == y_control)

    # update yreg
    if (inference == "delta"&& !is.null(weights_yreg)) {
      call_yreg$design <- eval(bquote(survey::svydesign(~1, weights = ~.(weights_yreg[indices]), data = .(data))))
    } else {
      call_yreg$weights <- weights_yreg[indices]
      call_yreg$data <- data
    }
    yreg <- eval.parent(call_yreg)

    # update mreg
    for (p in 1:length(mediator)) {
      # update mreg[[p]]
      if (inference == "delta" && !is.null(weights_mreg[[p]]) && !is_svymultinom_mreg[p]) {
        call_mreg[[p]]$design <-
          eval(bquote(survey::svydesign(~1, weights = ~.(weights_mreg[[p]][indices][control_indices]), data = .(data[control_indices, ]))))
      } else {
        call_mreg[[p]]$weights <- weights_mreg[[p]][indices][control_indices]
        call_mreg[[p]]$data <- data[control_indices, ]
      }
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }

    rm(control_indices)

  } else {

    # not a case control design
    # update yreg
    if (inference == "delta" && !is.null(weights_yreg)) {
      call_yreg$design <- eval(bquote(survey::svydesign(~1, weights = ~.(weights_yreg[indices]), data = .(data))))
    } else {
      call_yreg$weights <- weights_yreg[indices]
      call_yreg$data <- data
    }
    yreg <- eval.parent(call_yreg)

    # update mreg
    for (p in 1:length(mediator)) {
      # update mreg[[p]]
      if (inference == "delta" && !is.null(weights_mreg[[p]]) && !is_svymultinom_mreg[p]) {
        call_mreg[[p]]$design <-
          eval(bquote(survey::svydesign(~1, weights = ~.(weights_mreg[[p]][indices]), data = .(data))))
      } else {
        call_mreg[[p]]$weights <- weights_mreg[[p]][indices]
        call_mreg[[p]]$data <- data
      }
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }

  }

  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yreg <- yreg
    out$reg.output$mreg <- mreg
  }

  ###################################################################################################
  #################################Closed-form Parameter Function Estimation#########################
  ###################################################################################################
  if (estimation == "paramfunc") {

    if (length(mediator) != 1) stop("estimation can be paramfunc for model = 'rb' when length(mediator) = 1")

    mreg <- mreg[[1]]

    # for categorical exposure, create indicator vectors for a and astar
    if (is.factor(data[, exposure]) | is.character(data[, exposure])) {
      a_lev<- levels(droplevels(as.factor(data[, exposure])))
      a <- as.numeric(a_lev == a)[-1]
      astar <- as.numeric(a_lev == astar)[-1]
      elevel <- length(a_lev)
      rm(a_lev)
    } else if (is.numeric(data[, exposure]) | is.logical(data[, exposure])) {
      elevel <- 2
    } else stop("The exposure variable should be numeric, logical, factor or character")

    # create covariate values to calculate conditional causal effects
    vecc <- c()
    if (length(prec) != 0) {
      for (i in 1:length(prec)) {
        if (is.factor(data[, prec[i]]) | is.character(data[, prec[i]])) {
          # extract conditional values of levels existing in the new data set
          c_lev_orig <- levels(droplevels(as.factor(data_orig[, prec[i]])))
          c_lev_new <- levels(droplevels(as.factor(data[, prec[i]])))
          c_lev_index <- which(c_lev_orig %in% c_lev_new)
          vecc <- c(vecc, c(NA, precval[[i]])[c_lev_index][-1])
          rm(c_lev_orig, c_lev_new, c_lev_index)
        } else if (is.numeric(data[, prec[i]]) | is.logical(data[, prec[i]])) {
          vecc <- c(vecc, precval[[i]])
        } else stop(paste0("The prec[", i, "] variable should be numeric, logical, factor or character"))
      }
    }

    # for categorical mediator, create an indicator vector for mstar
    if (is.factor(data[, mediator]) | is.character(data[, mediator])) {
      m_lev <- levels(droplevels(as.factor(data[, mediator])))
      mstar <- as.numeric(m_lev == mval[[1]])[-1]
      mlevel <- length(m_lev)
      rm(m_lev)
    } else if (is.numeric(data[, mediator]) | is.logical(data[, mediator])) {
      mstar <- mval[[1]]
      mlevel <- 2
    } else stop("The mediator variable should be numeric, logical, factor or character")

    # coefficients for yreg
    thetas <- coef(yreg)
    # coefficients for mreg
    betas  <- as.vector(t(coef(mreg)))

    # intercept coefficient for yreg
    theta0 <- thetas[1]
    # exposure coefficient for yreg
    theta1 <- thetas[2:elevel]
    # mediator coefficient for yreg
    theta2 <- thetas[(elevel+1):(elevel + mlevel - 1)]
    # exposure-mediator interaction coefficient for yreg
    switch(as.character(EMint),
           "TRUE" = theta3 <- t(matrix(thetas[length(thetas) - (((elevel-1)*(mlevel-1)-1):0)],
                                       ncol = mlevel - 1)),
           "FALSE" = theta3 <- t(matrix(rep(0, (elevel-1)*(mlevel-1)), ncol = mlevel - 1)))

    # intercept coefficient for mreg
    beta0 <- betas[1+(0:(mlevel-2))*length(betas)/(mlevel-1)]
    # exposure coefficient for mreg
    beta1 <- t(matrix(betas[rowSums(expand.grid(2:elevel, (0:(mlevel-2)) *
                                                  length(betas)/(mlevel-1)))], ncol = mlevel - 1))
    # t(vecc)%*%beta'2 for mreg
    covariatesTerm <- sapply(0:(mlevel-2), function(x)
      ifelse(length(prec) == 0, 0, sum(betas[elevel + 1:length(vecc) +
                                               x * length(betas)/(mlevel-1)] * vecc)))

    # closed-form parameter function estimation
    if ((is_lm_yreg | is_glm_yreg) && family_yreg$family == "gaussian" && family_yreg$link == "identity") {
      if ((is_lm_mreg[1] | (is_glm_mreg[1])) && family_mreg[[1]]$family == "gaussian" && family_mreg[[1]]$link == "identity") {

        # linear Y with linear M
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

      } else {

        # linear Y with categorical M
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

      if (full) {
        pm <- tnie / te
        intref <- pnde - cde
        intmed <- tnie - pnie
        cde_prop <- cde/te
        intref_prop <- intref/te
        intmed_prop <- intmed/te
        pnie_prop <- pnie/te
        overall_pm <- (pnie+intmed)/te
        overall_int <- (intref+intmed)/te
        overall_pe <- (intref+intmed+pnie)/te
        est <- c(cde, pnde, tnde, pnie, tnie, te, pm, intref, intmed, cde_prop, intref_prop, 
                 intmed_prop, pnie_prop, overall_pm, overall_int, overall_pe)
      } else est <- c(cde, pnde, tnde, pnie, tnie, te)

    } else {

      if ((is_lm_mreg[1] | (is_glm_mreg[1])) && family_mreg[[1]]$family == "gaussian" && family_mreg[[1]]$link == "identity") {

        # nonlinear Y with linear M
        variance <- sigma(mreg)^2
        logRRcde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                             (sum(theta3 * a) - sum(theta3 * astar)) * mstar)
        logRRpnde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                              (sum(theta3 * a) - sum(theta3 * astar)) *
                              (beta0 + sum(beta1 * astar) +
                                 covariatesTerm + theta2  * variance) +
                              0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar)))
        logRRtnde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                              (sum(theta3 * a) - sum(theta3 * astar)) *
                              (beta0 + sum(beta1 * a) +
                                 covariatesTerm + theta2  * variance) +
                              0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar)))
        logRRpnie <- unname(theta2 * (sum(beta1 * a) - sum(beta1 * astar)) +
                              sum(theta3  * astar) * (sum(beta1 * a) - sum(beta1 * astar)))
        logRRtnie <- unname(theta2 * (sum(beta1 * a) - sum(beta1 * astar)) +
                              sum(theta3  * a) * (sum(beta1 * a) - sum(beta1 * astar)))
        if (full) ERRcde <- unname((exp(sum(theta1 * a) - sum(theta1 * astar) +
                                          sum(theta3  * a) * mstar) -
                                      exp(sum(theta3 * astar) * mstar)) *
                                     exp(theta2 * mstar - (theta2 + sum(theta3 * astar)) *
                                           (beta0 + sum(beta1 * astar) + covariatesTerm) -
                                           0.5 * (theta2 + sum(theta3 * astar)) ^ 2 * variance))

      } else {

        # nonlinear Y with categorical M
        logRRcde <- unname(sum(theta1 * a) - sum(theta1 * astar) +
                             ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a) -
                             ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar))
        logRRpnde <- unname(log((exp(sum(theta1 * a) - sum(theta1 * astar)) *
                                   (1 + sum(exp(theta2 + theta3 %*% a +
                                                  beta0 + beta1 %*% astar + covariatesTerm)))) /
                                  (1 + sum(exp(theta2 + theta3 %*% astar +
                                                 beta0 + beta1 %*% astar + covariatesTerm)))))
        logRRtnde <- unname(log((exp(sum(theta1 * a) - sum(theta1 * astar)) *
                                   (1 + sum(exp(theta2 + theta3 %*% a +
                                                  beta0 + beta1 %*% a + covariatesTerm)))) /
                                  (1 + sum(exp(theta2 + theta3 %*% astar +
                                                 beta0 + beta1 %*% a + covariatesTerm)))))
        logRRpnie <- unname(log(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                                   (1 + sum(exp(theta2 + theta3 %*% astar +
                                                  beta0 + beta1 %*% a + covariatesTerm)))) /
                                  ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                                     (1 + sum(exp(theta2 + theta3 %*% astar +
                                                    beta0 + beta1 %*% astar + covariatesTerm))))))
        logRRtnie <- unname(log(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                                   (1 + sum(exp(theta2 + theta3 %*% a +
                                                  beta0 + beta1 %*% a + covariatesTerm)))) /
                                  ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                                     (1 + sum(exp(theta2 + theta3 %*% a +
                                                    beta0 + beta1 %*% astar + covariatesTerm))))))
        if (full) ERRcde <- unname(exp(sum(theta2*mstar)) *
                                     (exp(sum(theta1 * a) - sum(theta1 * astar) +
                                            ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a)) -
                                        exp(ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar))) *
                                     (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) /
                                     (1+ sum(exp(theta2 + theta3 %*% astar + beta0 +
                                                   beta1 %*% astar + covariatesTerm))))

      }

      logRRte <- logRRtnie + logRRpnde

      if (full) {
        pm <- (exp(logRRpnde) * (exp(logRRtnie) - 1)) / (exp(logRRte) - 1)
        ERRintref <- exp(logRRpnde) - 1 - ERRcde
        ERRintmed <- exp(logRRtnie) * exp(logRRpnde) - exp(logRRpnde) - exp(logRRpnie) + 1
        ERRpnie <- exp(logRRpnie) - 1
        ERRte <- exp(logRRte) - 1
        ERRcde_prop <- ERRcde/ERRte
        ERRintmed_prop <- ERRintmed/ERRte
        ERRintref_prop <- ERRintref/ERRte
        ERRpnie_prop <- ERRpnie/ERRte
        overall_pm <- (ERRpnie+ERRintmed)/ERRte
        overall_int <- (ERRintref+ERRintmed)/ERRte
        overall_pe <- (ERRintref+ERRintmed+ERRpnie)/ERRte
        est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, pm,
                 ERRcde, ERRintref, ERRintmed, ERRpnie,
                 ERRcde_prop, ERRintref_prop, ERRintmed_prop, ERRpnie_prop,
                 overall_pm, overall_int, overall_pe)
      } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte)

    }

    ###################################################################################################
    ###############################Direct Counterfactual Imputation Estimation#########################
    ###################################################################################################
  } else if (estimation == "imputation") {

    # the index of the reference level for a categorical outcome
    if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                         startsWith(family_yreg$family, "Ordered Categorical"))) |
        is_multinom_yreg | is_polr_yreg) {
      y_lev <- levels(droplevels(as.factor(data[, outcome])))
      yref_index <- switch((yref %in% y_lev) + 1, "1" = NULL, "2" = which(y_lev == yref))
      rm(y_lev)
    }

    # simulate A
    if (is.factor(data[, exposure])) {
      a_sim <- factor(c(rep(a, n)), levels = a_lev)
      astar_sim <- factor(c(rep(astar, n)), levels = a_lev)
    } else {
      a_sim <- c(rep(a, n))
      astar_sim <- c(rep(astar, n))
    }

    # simulate C
    prec_sim <- data[, prec]

    # design matrices for simulating mediator[1]
    mdesign_a <- data.frame(a_sim, prec_sim)
    mdesign_astar <- data.frame(astar_sim, prec_sim)
    colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, prec)

    m_a <- m_astar <- data.frame(matrix(nrow = n, ncol = length(mediator)))
    colnames(m_a) <- colnames(m_astar) <- mediator

    # simulating mediator[p]
    for (p in 1:length(mediator)) {

      # design matrices for simulating mediator[p]
      mdesign_a <- cbind(mdesign_a, m_a[, p - 1, drop = FALSE])
      mdesign_astar <- cbind(mdesign_astar, m_astar[, p - 1, drop = FALSE])

      # predict mediator[p]
      type <- ifelse(is_multinom_mreg[p] | is_polr_mreg[p], "probs", "response")
      mpred_a <- predict(mreg[[p]], newdata = mdesign_a, type = type)
      mpred_astar <- predict(mreg[[p]], newdata = mdesign_astar, type = type)

      # categorical M
      if ((is_glm_mreg[p] && ((family_mreg[[p]]$family %in% c("binomial", "multinom")) |
                              startsWith(family_mreg[[p]]$family, "Ordered Categorical")))|
          is_multinom_mreg[p] | is_polr_mreg[p]) {

        prob_a <- as.matrix(mpred_a)
        prob_astar <- as.matrix(mpred_astar)

        if (dim(prob_a)[2] == 1) {
          # simulate mediator[p] for exposure=a
          msim_a <- rbinom(n, size = 1, prob = prob_a[, 1]) + 1
          # simulate mediator[p] for exposure=astar
          msim_astar <- rbinom(n, size = 1, prob = prob_astar[, 1]) + 1
        } else {
          msim_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))
          msim_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))
        }

        m_lev <- levels(droplevels(as.factor(model.frame(mreg[[p]])[, mediator[p]])))
        # mid_a: simulated mediator[p] for exposure = a
        # mid_astar: simulated mediator[p] for exposure = astar
        if (is.factor(data[, mediator[p]])) {
          mid_a <- factor(m_lev[msim_a], levels = m_lev)
          mid_astar <- factor(m_lev[msim_astar], levels = m_lev)
        } else if (is.character(data[, mediator[p]])) {
          mid_a <- m_lev[msim_a]
          mid_astar <- m_lev[msim_astar]
        } else if (is.numeric(data[, mediator[p]])) {
          mid_a <- as.numeric(m_lev[msim_a])
          mid_astar <- as.numeric(m_lev[msim_astar])
        } else if (is.logical(data[, mediator[p]])) {
          mid_a <- as.logical(m_lev[msim_a])
          mid_astar <- as.logical(m_lev[msim_astar])
        } else stop("The mediator[", p, "] variable should be numeric, logical, factor or character")

        rm(prob_a, prob_astar, msim_a, msim_astar, m_lev)

        # linear M
      } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "gaussian") {

        error <- rnorm(n, mean = 0, sd = sigma(mreg[[p]]))
        mid_a <- mpred_a + error
        mid_astar <- mpred_astar + error

        rm(error)

        # gamma M
      } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "Gamma") {

        shape_mreg <- gamma.shape(mreg[[p]])$alpha
        mid_a <- rgamma(n, shape = shape_mreg, scale = mpred_a/shape_mreg)
        mid_astar <- rgamma(n, shape = shape_mreg, scale = mpred_astar/shape_mreg)

        rm(shape_mreg)

        # inverse gaussian M
      } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "inverse.gaussian") {

        lambda <- 1/summary(mreg[[p]])$dispersion
        mid_a <- SuppDists::rinvGauss(n, nu = mpred_a, lambda = lambda)
        mid_astar <- SuppDists::rinvGauss(n, nu = mpred_astar, lambda = lambda)

        rm(lambda)

        # poisson M
      } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "poisson") {

        mid_a <- rpois(n, lambda = mpred_a)
        mid_astar <- rpois(n, lambda = mpred_astar)

        # negative binomial M
      } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && startsWith(family_reg[[p]]$family, "Negative Binomial")) {

        theta <- summary(mreg[[p]])$theta
        mid_a <- MASS::rnegbin(n, mu = mpred_a, theta = theta)
        mid_astar <- MASS::rpois(n, mu = mpred_astar, theta = theta)

        rm(theta)

      } else stop(paste0("Unsupported mreg[[", p, "]]"))

      m_a[, p] <- mid_a
      m_astar[, p] <- mid_astar

    }

    rm(mdesign_a, mdesign_astar, type, mpred_a, mpred_astar, mid_a, mid_astar)

    # simulate mstar for cde
    mstar_sim <- do.call(cbind, lapply(1:length(mediator), function(x)
      if (is.factor(data[, mediator[x]])) {
        data.frame(factor(rep(mval[[x]], n), levels = levels(data[, mediator[x]])))
      } else data.frame(rep(mval[[x]], n))))

    # design matrices for outcome simulation
    ydesign0m <- data.frame(astar_sim, mstar_sim, prec_sim)
    ydesign1m <- data.frame(a_sim, mstar_sim, prec_sim)
    ydesign00 <- data.frame(astar_sim, m_astar, prec_sim)
    ydesign01 <- data.frame(astar_sim, m_a, prec_sim)
    ydesign10 <- data.frame(a_sim, m_astar, prec_sim)
    ydesign11 <- data.frame(a_sim, m_a, prec_sim)
    rm(a_sim, astar_sim, m_a, m_astar, mstar_sim, prec_sim)
    colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
      colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, prec)

    # predict Y
    type <- ifelse(is_coxph_yreg, "risk", ifelse(is_multinom_yreg | is_polr_yreg, "probs", "response"))
    EY0m_pred <- as.matrix(predict(yreg, newdata =  ydesign0m, type = type))
    EY1m_pred <- as.matrix(predict(yreg, newdata =  ydesign1m, type = type))
    EY00_pred <- as.matrix(predict(yreg, newdata =  ydesign00, type = type))
    EY01_pred <- as.matrix(predict(yreg, newdata =  ydesign01, type = type))
    EY10_pred <- as.matrix(predict(yreg, newdata =  ydesign10, type = type))
    EY11_pred <- as.matrix(predict(yreg, newdata =  ydesign11, type = type))
    rm(type, ydesign0m, ydesign1m, ydesign00, ydesign01, ydesign10, ydesign11)

    # weights of yreg
    weightsEY <- model.frame(yreg)$'(weights)'
    if (is.null(weightsEY)) weightsEY <- rep(1, n)

    # categorical Y
    if ((is_glm_yreg && ((family_yreg$family %in% c("binomial", "quasibinomial", "multinom")) |
                         startsWith(family_yreg$family, "Ordered Categorical")))|
        is_multinom_yreg | is_polr_yreg) {

      if (!is.null(yref_index)) {
        if (dim(EY0m_pred)[2] == 1) {
          EY0m <- weighted.mean(cbind(1 - EY0m_pred, EY0m_pred)[, yref_index], na.rm = TRUE, w = weightsEY)
          EY1m <- weighted.mean(cbind(1 - EY1m_pred, EY1m_pred)[, yref_index], na.rm = TRUE, w = weightsEY)
          EY00 <- weighted.mean(cbind(1 - EY00_pred, EY00_pred)[, yref_index], na.rm = TRUE, w = weightsEY)
          EY01 <- weighted.mean(cbind(1 - EY01_pred, EY01_pred)[, yref_index], na.rm = TRUE, w = weightsEY)
          EY10 <- weighted.mean(cbind(1 - EY10_pred, EY10_pred)[, yref_index], na.rm = TRUE, w = weightsEY)
          EY11 <- weighted.mean(cbind(1 - EY11_pred, EY11_pred)[, yref_index], na.rm = TRUE, w = weightsEY)
        } else {
          EY0m <- weighted.mean(EY0m_pred[, yref_index], na.rm = TRUE, w = weightsEY)
          EY1m <- weighted.mean(EY1m_pred[, yref_index], na.rm = TRUE, w = weightsEY)
          EY00 <- weighted.mean(EY00_pred[, yref_index], na.rm = TRUE, w = weightsEY)
          EY01 <- weighted.mean(EY01_pred[, yref_index], na.rm = TRUE, w = weightsEY)
          EY10 <- weighted.mean(EY10_pred[, yref_index], na.rm = TRUE, w = weightsEY)
          EY11 <- weighted.mean(EY11_pred[, yref_index], na.rm = TRUE, w = weightsEY)
        }
      } else EY0m <- EY1m <- EY00 <- EY01 <- EY10 <- EY11 <- 0

    } else {

      # non-categorical Y
      EY0m <- weighted.mean(EY0m_pred, na.rm = TRUE, w = weightsEY)
      EY1m <- weighted.mean(EY1m_pred, na.rm = TRUE, w = weightsEY)
      EY00 <- weighted.mean(EY00_pred, na.rm = TRUE, w = weightsEY)
      EY01 <- weighted.mean(EY01_pred, na.rm = TRUE, w = weightsEY)
      EY10 <- weighted.mean(EY10_pred, na.rm = TRUE, w = weightsEY)
      EY11 <- weighted.mean(EY11_pred, na.rm = TRUE, w = weightsEY)

    }

    rm(weightsEY, EY0m_pred, EY1m_pred, EY00_pred, EY01_pred, EY10_pred, EY11_pred)

    # output causal effects in additive scale for continuous Y
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi", "gaulss", "gevlss") |
         startsWith(family_yreg$family, "Tweedie") |
         startsWith(family_yreg$family, "Beta regression") |
         startsWith(family_yreg$family, "Scaled t"))) {

      cde <- EY1m - EY0m
      pnde <- EY10 - EY00
      tnde <- EY11 - EY01
      pnie <- EY01 - EY00
      tnie <- EY11 - EY10
      te <- tnie + pnde

      if (full) {
        pm <- tnie / te
        intref <- pnde - cde
        intmed <- tnie - pnie
        cde_prop <- cde/te
        intref_prop <- intref/te
        intmed_prop <- intmed/te
        pnie_prop <- pnie/te
        overall_pm <- (pnie + intmed)/te
        overall_int <- (intref + intmed)/te
        overall_pe <- (intref + intmed + pnie)/te
        est <- c(cde, pnde, tnde, pnie, tnie, te, pm, intref, intmed, cde_prop, intref_prop, 
                 intmed_prop, pnie_prop, overall_pm, overall_int, overall_pe)

      } else est <- c(cde, pnde, tnde, pnie, tnie, te)

    } else if (((is_lm_yreg | is_glm_yreg) &&
                (family_yreg$family %in% c("binomial", "quasibinomial", "multinom", "poisson", "quasipoisson", "ziplss") |
                 startsWith(family_yreg$family, "Negative Binomial") |
                 startsWith(family_yreg$family, "Zero inflated Poisson") |
                 startsWith(family_yreg$family, "Ordered Categorical"))) |
               is_multinom_yreg | is_polr_yreg | is_survreg_yreg | is_coxph_yreg){

      # output causal effects in ratio scale for non-continuous Y
      logRRcde <- log(EY1m) - log(EY0m)
      logRRpnde <- log(EY10) - log(EY00)
      logRRtnde <- log(EY11) - log(EY01)
      logRRpnie <- log(EY01) - log(EY00)
      logRRtnie <- log(EY11) - log(EY10)
      logRRte <- logRRtnie + logRRpnde

      if (full) {
        pm <- (exp(logRRpnde) * (exp(logRRtnie) - 1)) / (exp(logRRte) - 1)
        ERRcde <- (EY1m-EY0m)/EY00
        ERRintref <- exp(logRRpnde) - 1 - ERRcde
        ERRintmed <- exp(logRRtnie) * exp(logRRpnde) - exp(logRRpnde) - exp(logRRpnie) + 1
        ERRpnie <- exp(logRRpnie) - 1
        ERRte <- exp(logRRte) - 1
        ERRcde_prop <- ERRcde/ERRte
        ERRintmed_prop <- ERRintmed/ERRte
        ERRintref_prop <- ERRintref/ERRte
        ERRpnie_prop <- ERRpnie/ERRte
        overall_pm <- (ERRpnie + ERRintmed)/ERRte
        overall_int <- (ERRintref + ERRintmed)/ERRte
        overall_pe <- (ERRintref + ERRintmed + ERRpnie)/ERRte
        est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, pm,
                 ERRcde, ERRintref, ERRintmed, ERRpnie,
                 ERRcde_prop, ERRintref_prop, ERRintmed_prop, ERRpnie_prop,
                 overall_pm, overall_int, overall_pe)
      } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte)

    } else stop("Unsupported yreg")

  }

  # progress bar
  if (inference == "bootstrap") {
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("progbar", envir = env), curVal + 1)
  }

  if (outReg) out$est <- est
  if (!outReg) out <- est

  return(out)

}

