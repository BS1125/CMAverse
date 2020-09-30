est.rb <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {
  data_orig <- data
  if (is.null(indices)) indices <- 1:n
  # resample data
  data <- data[indices, ]
  
  # for case control study
  # method 1: weight subjects with y=1 by yprevalence/p(y=1) and weight subjects with y=0 by (1-yprevalence)/p(y=0)
  # method 2: fit yreg with all data and fit other regs on data among controls
  # use method 1 when yprevalence is provided
  # when yprevalence is not provided but the outcome is rare, use method 2
  if (casecontrol && !is.null(yprevalence)) {
    # method 1 for a case control design
    prob1 <- mean(data[, outcome] == y_case, na.rm = TRUE)
    w4casecon <- as.vector(ifelse(data[, outcome] == y_case, yprevalence / prob1, (1 - yprevalence) / (1 - prob1)))
    # update yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon
    if (is.null(weights_yreg)) weights_yreg <- w4casecon
    if (is_svyglm_yreg && !(inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) {
      call_design <- getCall(yreg$survey.design)
      call_design$weights <- weights_yreg
      call_design$data <- data
      call_yreg$design <- eval.parent(call_design)
    } else {
      call_yreg$weights <- weights_yreg
      call_yreg$data <- data
    }
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      if (!is.null(weights_mreg[[p]])) weights_mreg[[p]] <- weights_mreg[[p]][indices] * w4casecon
      if (is.null(weights_mreg[[p]])) weights_mreg[[p]] <- w4casecon
      if (is_svyglm_mreg[p] && !(inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) {
        call_design <- getCall(mreg[[p]]$survey.design)
        call_design$weights <- weights_mreg[[p]]
        call_design$data <- data
        call_mreg[[p]]$design <- eval.parent(call_design)
      } else {
        call_mreg[[p]]$weights <- weights_mreg[[p]]
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
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      call_mreg[[p]]$weights <- weights_mreg[[p]][indices][control_indices]
      call_mreg[[p]]$data <- data[control_indices, ]
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    rm(control_indices)
  } else {
    # not a case control design
    # update yreg
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      call_mreg[[p]]$weights <- weights_mreg[[p]][indices]
      call_mreg[[p]]$data <- data
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
    mreg <- mreg[[1]]
    
    # for categorical exposure, create indicator vectors for a and astar
    if (is.factor(data[, exposure]) | is.character(data[, exposure])) {
      a_lev <- levels(droplevels(as.factor(data[, exposure])))
      a <- as.numeric(a_lev == a)[-1]
      astar <- as.numeric(a_lev == astar)[-1]
      elevel <- length(a_lev)
      rm(a_lev)
    } else if (is.numeric(data[, exposure])) elevel <- 2
    
    # create covariate values to calculate conditional causal effects
    vecc <- c()
    if (length(basec) != 0) {
      for (i in 1:length(basec)) {
        if (is.factor(data[, basec[i]]) | is.character(data[, basec[i]])) {
          # extract conditional values of levels existing in the new data set
          c_lev_orig <- levels(droplevels(as.factor(data_orig[, basec[i]])))
          c_lev_new <- levels(droplevels(as.factor(data[, basec[i]])))
          c_lev_index <- which(c_lev_orig %in% c_lev_new)
          vecc <- c(vecc, c(NA, basecval[[i]])[c_lev_index][-1])
          rm(c_lev_orig, c_lev_new, c_lev_index)
        } else if (is.numeric(data[, basec[i]])) vecc <- c(vecc, basecval[[i]])
      }
    }
    
    # for categorical mediator, create an indicator vector for mstar
    if (is.factor(data[, mediator]) | is.character(data[, mediator])) {
      m_lev <- levels(droplevels(as.factor(data[, mediator])))
      mstar <- as.numeric(m_lev == mval[[1]])[-1]
      mlevel <- length(m_lev)
      rm(m_lev)
    } else if (is.numeric(data[, mediator])) {
      mstar <- mval[[1]]
      mlevel <- 2
    }
    
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
           "TRUE" = theta3 <- t(matrix(thetas[length(thetas) - (((elevel-1)*(mlevel-1)-1):0)], ncol = mlevel - 1)),
           "FALSE" = theta3 <- t(matrix(rep(0, (elevel-1)*(mlevel-1)), ncol = mlevel - 1)))
    
    # intercept coefficient for mreg
    beta0 <- betas[1+(0:(mlevel-2))*length(betas)/(mlevel-1)]
    # exposure coefficient for mreg
    beta1 <- t(matrix(betas[rowSums(expand.grid(2:elevel, (0:(mlevel-2)) * length(betas)/(mlevel-1)))], ncol = mlevel - 1))
    # t(vecc)%*%beta'2 for mreg
    covariatesTerm <- sapply(0:(mlevel-2), function(x)
      ifelse(length(basec) == 0, 0, sum(betas[elevel + 1:length(vecc) + x * length(betas)/(mlevel-1)] * vecc)))
    
    # closed-form parameter function estimation
    if ((is_lm_yreg | is_glm_yreg) && family_yreg$family == "gaussian") {
      if ((is_lm_mreg[1] | (is_glm_mreg[1])) && family_mreg[[1]]$family == "gaussian") {
        # linear Y with linear M
        cde <- sum(theta1 * a) - sum(theta1 * astar) + (sum(theta3 * a) - sum(theta3 * astar)) * mstar
        pnde <- sum(theta1 * a) - sum(theta1 * astar) + (sum(theta3 * a) - sum(theta3 * astar)) *
          (beta0 + sum(beta1 * astar) + covariatesTerm)
        tnde <- sum(theta1 * a) - sum(theta1 * astar) + (sum(theta3 * a) - sum(theta3 * astar)) *
          (beta0 + sum(beta1 * a) + covariatesTerm)
        pnie <- (theta2 + sum(theta3  * astar)) * (sum(beta1 * a) - sum(beta1 * astar))
        tnie <- (theta2 + sum(theta3  * a)) * (sum(beta1 * a) - sum(beta1 * astar))
      } else {
        # linear Y with categorical M
        cde <- (sum(theta1 * a) - sum(theta1 * astar) +
                  ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a ) -
                  ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar))
        pnde <- (sum(theta1 * a) - sum(theta1 * astar) +
                   (sum((theta3 %*% a - theta3 %*% astar) *
                          exp(beta0 + beta1 %*% astar + covariatesTerm)) /
                      (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm)))))
        tnde <- (sum(theta1 * a) - sum(theta1 * astar) +
                   (sum((theta3 %*% a - theta3 %*% astar) *
                          exp(beta0 + beta1 %*% a + covariatesTerm)) /
                      (1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm)))))
        pnie <- sum((theta2+theta3 %*% astar)*exp(beta0 + beta1 %*% a + covariatesTerm)) /
          (1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) -
          sum((theta2+theta3 %*% astar)*exp(beta0 + beta1 %*% astar + covariatesTerm)) /
          (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm)))
        tnie <- sum((theta2+theta3 %*% a)*exp(beta0 + beta1 %*% a + covariatesTerm)) /
          (1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) -
          sum((theta2+theta3 %*% a)*exp(beta0 + beta1 %*% astar + covariatesTerm)) /
          (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm)))
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
        int <- (intref+intmed)/te
        pe <- (intref+intmed+pnie)/te
        est <- unname(c(cde, pnde, tnde, pnie, tnie, te, intref, intmed, cde_prop, intref_prop, 
                        intmed_prop, pnie_prop, pm, int, pe))
      } else est <- unname(c(cde, pnde, tnde, pnie, tnie, te))
    } else {
      if ((is_lm_mreg[1] | (is_glm_mreg[1])) && family_mreg[[1]]$family == "gaussian") {
        # nonlinear Y with linear M
        variance <- sigma(mreg)^2
        logRRcde <- sum(theta1 * a) - sum(theta1 * astar) + (sum(theta3 * a) - sum(theta3 * astar)) * mstar
        logRRpnde <- sum(theta1 * a) - sum(theta1 * astar) + (sum(theta3 * a) - sum(theta3 * astar)) *
          (beta0 + sum(beta1 * astar) + covariatesTerm + theta2  * variance) +
          0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar))
        logRRtnde <- sum(theta1 * a) - sum(theta1 * astar) + (sum(theta3 * a) - sum(theta3 * astar)) *
          (beta0 + sum(beta1 * a) + covariatesTerm + theta2  * variance) +
          0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar))
        logRRpnie <- theta2 * (sum(beta1 * a) - sum(beta1 * astar)) + sum(theta3  * astar) * (sum(beta1 * a) - sum(beta1 * astar))
        logRRtnie <- theta2 * (sum(beta1 * a) - sum(beta1 * astar)) + sum(theta3  * a) * (sum(beta1 * a) - sum(beta1 * astar))
        if (full) ERRcde <- (exp(sum(theta1 * a) - sum(theta1 * astar) + sum(theta3  * a) * mstar) -
                               exp(sum(theta3 * astar) * mstar)) * exp(theta2 * mstar - (theta2 + sum(theta3 * astar)) *
                                                                         (beta0 + sum(beta1 * astar) + covariatesTerm) -
                                                                         0.5 * (theta2 + sum(theta3 * astar)) ^ 2 * variance)
      } else {
        # nonlinear Y with categorical M
        logRRcde <- sum(theta1 * a) - sum(theta1 * astar) +
          ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a) -
          ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar)
        logRRpnde <- log((exp(sum(theta1 * a) - sum(theta1 * astar)) * (1 + sum(exp(theta2 + theta3 %*% a +
                                                                                      beta0 + beta1 %*% astar + covariatesTerm)))) /
                           (1 + sum(exp(theta2 + theta3 %*% astar + beta0 + beta1 %*% astar + covariatesTerm))))
        logRRtnde <- log((exp(sum(theta1 * a) - sum(theta1 * astar)) * (1 + sum(exp(theta2 + theta3 %*% a +
                                                                                      beta0 + beta1 %*% a + covariatesTerm)))) /
                           (1 + sum(exp(theta2 + theta3 %*% astar + beta0 + beta1 %*% a + covariatesTerm))))
        logRRpnie <- log(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                            (1 + sum(exp(theta2 + theta3 %*% astar + beta0 + beta1 %*% a + covariatesTerm)))) /
                           ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                              (1 + sum(exp(theta2 + theta3 %*% astar + beta0 + beta1 %*% astar + covariatesTerm)))))
        logRRtnie <- log(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                            (1 + sum(exp(theta2 + theta3 %*% a + beta0 + beta1 %*% a + covariatesTerm)))) /
                           ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                              (1 + sum(exp(theta2 + theta3 %*% a + beta0 + beta1 %*% astar + covariatesTerm)))))
        if (full) ERRcde <- exp(sum(theta2*mstar)) * (exp(sum(theta1 * a) - sum(theta1 * astar) +
                                                            ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a)) -
                                                        exp(ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar))) *
          (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) /
          (1+ sum(exp(theta2 + theta3 %*% astar + beta0 + beta1 %*% astar + covariatesTerm)))
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
        int <- (ERRintref+ERRintmed)/ERRte
        pe <- (ERRintref+ERRintmed+ERRpnie)/ERRte
        est <- unname(c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, 
                        ERRcde, ERRintref, ERRintmed, ERRpnie,
                        ERRcde_prop, ERRintref_prop, ERRintmed_prop, ERRpnie_prop,
                        pm, int, pe))
      } else est <- unname(c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte))
    }
    
    ###################################################################################################
    ###############################Direct Counterfactual Imputation Estimation#########################
    ###################################################################################################
  } else if (estimation == "imputation") {
    # the index of the reference level for a categorical outcome
    if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                         startsWith(family_yreg$family, "Ordered Categorical"))) |
        is_multinom_yreg | is_polr_yreg) {
      yref_index <- switch((yref %in% y_lev) + 1, "1" = NULL, "2" = which(y_lev == yref))
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
    basec_sim <- data[, basec]
    
    # design matrices for simulating mediator[p]
    mdesign_a <- data.frame(a_sim, basec_sim)
    mdesign_astar <- data.frame(astar_sim, basec_sim)
    colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, basec)
    m_a <- m_astar <- data.frame(matrix(nrow = n, ncol = length(mediator)))
    colnames(m_a) <- colnames(m_astar) <- mediator
    # simulating mediator[p]
    for (p in 1:length(mediator)) {
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
        m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
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
        } 
        rm(prob_a, prob_astar, msim_a, msim_astar, m_lev)
        # linear M
      } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "gaussian") {
        error <- rnorm(n, mean = 0, sd = sigma(mreg[[p]]))
        mid_a <- mpred_a + error
        mid_astar <- mpred_astar + error
        rm(error)
        # gamma M
      } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "Gamma") {
        shape_mreg <- MASS::gamma.shape(mreg[[p]])$alpha
        mid_a <- rgamma(n, shape = shape_mreg, scale = mpred_a/shape_mreg)
        mid_astar <- rgamma(n, shape = shape_mreg, scale = mpred_astar/shape_mreg)
        rm(shape_mreg)
        # inverse gaussian M
      } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "inverse.gaussian") {
        lambda <- 1/summary(mreg[[p]])$dispersion
        mid_a <- SuppDists::rinvGauss(n, nu = mpred_a, lambda = lambda)
        mid_astar <- SuppDists::rinvGauss(n, nu = mpred_astar, lambda = lambda)
        rm(lambda)
        # poisson M
      } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "poisson") {
        mid_a <- rpois(n, lambda = mpred_a)
        mid_astar <- rpois(n, lambda = mpred_astar)
        # quasipoisson M
      } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "quasipoisson") {
        phi <- summary(mreg[[p]])$dispersion
        mid_a <- rqpois(n, lambda = mpred_a, phi = phi)
        mid_astar <- rqpois(n, lambda = mpred_astar, phi = phi)
        rm(phi)
        # negative binomial M
      } else if ( is_glm_mreg[p] && startsWith(family_mreg[[p]]$family, "Negative Binomial")) {
        theta <- summary(mreg[[p]])$theta
        mid_a <- MASS::rnegbin(n, mu = mpred_a, theta = theta)
        mid_astar <- MASS::rnegbin(n, mu = mpred_astar, theta = theta)
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
    ydesign0m <- data.frame(astar_sim, mstar_sim, basec_sim)
    ydesign1m <- data.frame(a_sim, mstar_sim, basec_sim)
    ydesign00 <- data.frame(astar_sim, m_astar, basec_sim)
    ydesign01 <- data.frame(astar_sim, m_a, basec_sim)
    ydesign10 <- data.frame(a_sim, m_astar, basec_sim)
    ydesign11 <- data.frame(a_sim, m_a, basec_sim)
    rm(a_sim, astar_sim, m_a, m_astar, mstar_sim, basec_sim)
    colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
      colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, basec)
    
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
    weightsEY <- as.vector(model.frame(yreg)$'(weights)')
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
    
    # output causal effects on the difference scale for continuous Y
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
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
        int <- (intref + intmed)/te
        pe <- (intref + intmed + pnie)/te
        est <- c(cde, pnde, tnde, pnie, tnie, te, intref, intmed, cde_prop, intref_prop, 
                 intmed_prop, pnie_prop, pm, int, pe)
      } else est <- c(cde, pnde, tnde, pnie, tnie, te)
    } else {
      # output causal effects on the ratio scale for non-continuous Y
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
        int <- (ERRintref + ERRintmed)/ERRte
        pe <- (ERRintref + ERRintmed + ERRpnie)/ERRte
        est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, 
                 ERRcde, ERRintref, ERRintmed, ERRpnie,
                 ERRcde_prop, ERRintref_prop, ERRintmed_prop, ERRpnie_prop,
                 pm, int, pe)
      } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte)
    }
  }
  
  # progress bar
  if (!multimp) {
    if (inference == "bootstrap") {
      curVal <- get("counter", envir = env)
      assign("counter", curVal + 1, envir = env)
      setTxtProgressBar(get("progbar", envir = env), curVal + 1)
    }
  }
  if (outReg) out$est <- est
  if (!outReg) out <- est
  return(out)
}

