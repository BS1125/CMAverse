est.gformula <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {
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
    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon
    if (is.null(weights_yreg)) weights_yreg <- w4casecon
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      if (!is.null(weights_mreg[[p]])) weights_mreg[[p]] <- weights_mreg[[p]][indices] * w4casecon
      if (is.null(weights_mreg[[p]])) weights_mreg[[p]] <- w4casecon
      call_mreg[[p]]$weights <- weights_mreg[[p]]
      call_mreg[[p]]$data <- data
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    if (length(postc) != 0) {
      # update postcreg
      for (p in 1:length(postc)) {
        if (!is.null(weights_postcreg[[p]])) weights_postcreg[[p]] <- weights_postcreg[[p]][indices] * w4casecon
        if (is.null(weights_postcreg[[p]])) weights_postcreg[[p]] <- w4casecon
        call_postcreg[[p]]$weights <- weights_postcreg[[p]]
        call_postcreg[[p]]$data <- data
        if (outReg && (inherits(postcreg[[p]], "rcreg") | inherits(postcreg[[p]], "simexreg"))) call_postcreg[[p]]$variance <- TRUE
        postcreg[[p]] <- eval.parent(call_postcreg[[p]])
      }
    }
    rm(prob1, w4casecon)
  } else if (casecontrol && yrare) {
    # method 2 for a case control design
    # data from controls
    control_indices <- which(data[, outcome] == y_control)
    # update yreg
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      call_mreg[[p]]$weights <- weights_mreg[[p]][indices][control_indices]
      call_mreg[[p]]$data <- data[control_indices, ]
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    # update postcreg
    if (length(postc) != 0) {
      for (p in 1:length(postc)) {
        # update postcreg[[p]]
        call_postcreg[[p]]$weights <- weights_postcreg[[p]][indices][control_indices]
        call_postcreg[[p]]$data <- data[control_indices, ]
        if (outReg && (inherits(postcreg[[p]], "rcreg") | inherits(postcreg[[p]], "simexreg"))) call_postcreg[[p]]$variance <- TRUE
        postcreg[[p]] <- eval.parent(call_postcreg[[p]])
      }
    }
    rm(control_indices)
  } else {
    # not a case control design
    # update yreg
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # update mreg
    for (p in 1:length(mediator)) {
      call_mreg[[p]]$weights <- weights_mreg[[p]][indices]
      call_mreg[[p]]$data <- data
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    # update postcreg
    if (length(postc) != 0) {
      for (p in 1:length(postc)) {
        call_postcreg[[p]]$weights <- weights_postcreg[[p]][indices]
        call_postcreg[[p]]$data <- data
        if (outReg && (inherits(postcreg[[p]], "rcreg") | inherits(postcreg[[p]], "simexreg"))) call_postcreg[[p]]$variance <- TRUE
        postcreg[[p]] <- eval.parent(call_postcreg[[p]])
      }
    }
  }
  
  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yreg <- yreg
    out$reg.output$mreg <- mreg
    if (length(postc) != 0) out$reg.output$postcreg <- postcreg
  }
  
  # the index of the reference level for a categorical outcome
  if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                       startsWith(family_yreg$family, "Ordered Categorical"))) |
      is_multinom_yreg | is_polr_yreg) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    yval_index <- switch((yval %in% y_lev) + 1, "1" = NULL, "2" = which(y_lev == yval))
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
  basec_sim <- data[, basec]
  
  # design matrices for simulating postc[p]
  postcdesign_a <- data.frame(a_sim, basec_sim)
  postcdesign_astar <- data.frame(astar_sim, basec_sim)
  colnames(postcdesign_a) <- colnames(postcdesign_astar) <- c(exposure, basec)
  postc_a <- postc_astar <- data.frame(matrix(nrow = n, ncol = length(postc)))
  colnames(postc_a) <- colnames(postc_astar) <- postc
  
  if (length(postc) != 0) {
    # simulating postc[p]
    for (p in 1:length(postc)) {
      # predict postc[p]
      type <- ifelse(is_multinom_postcreg[p] | is_polr_postcreg[p], "probs", "response")
      postcpred_a <- predict(postcreg[[p]], newdata = postcdesign_a, type = type)
      postcpred_astar <- predict(postcreg[[p]], newdata = postcdesign_astar, type = type)
      full_index <- which(rowSums(is.na(postcdesign_a))==0)
      n_full <- length(full_index)
      # categorical L
      if ((is_glm_postcreg[p] && ((family_postcreg[[p]]$family %in% c("binomial", "multinom")) |
                                  startsWith(family_postcreg[[p]]$family, "Ordered Categorical")))|
          is_multinom_postcreg[p] | is_polr_postcreg[p]) {
        l_lev <- levels(droplevels(as.factor(data[, postc[p]])))
        prob_a <- as.matrix(postcpred_a)
        prob_astar <- as.matrix(postcpred_astar)
        if (dim(prob_a)[2] == 1) {
          mid_a <- l_lev[rbinom(n_full, size = 1, prob = prob_a[full_index, 1]) + 1]
          mid_astar <- l_lev[rbinom(n_full, size = 1, prob = prob_astar[full_index, 1]) + 1]
        } else {
          mid_a <- l_lev[apply(prob_a[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
          mid_astar <- l_lev[apply(prob_astar[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
        }
        
        if (is.numeric(data[, postc[p]])) {
          mid_a <- as.numeric(l_lev[mid_a])
          mid_astar <- as.numeric(l_lev[mid_astar])
        } 
        
        rm(prob_a, prob_astar, l_lev)
        # linear L
      } else if ((is_lm_postcreg[p] | is_glm_postcreg[p]) && family_postcreg[[p]]$family == "gaussian") {
        error <- rnorm(n_full, mean = 0, sd = sigma(postcreg[[p]]))
        mid_a <- postcpred_a[full_index] + error
        mid_astar <- postcpred_astar[full_index] + error
        rm(error)
        # gamma L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "Gamma") {
        shape_postcreg <- MASS::gamma.shape(postcreg[[p]])$alpha
        mid_a <- rgamma(n_full, shape = shape_postcreg, scale = postcpred_a[full_index]/shape_postcreg)
        mid_astar <- rgamma(n_full, shape = shape_postcreg, scale = postcpred_astar[full_index]/shape_postcreg)
        rm(shape_postcreg)
        # inverse gaussian L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "inverse.gaussian") {
        lambda <- 1/summary(postcreg[[p]])$dispersion
        mid_a <- SuppDists::rinvGauss(n_full, nu = postcpred_a[full_index], lambda = lambda)
        mid_astar <- SuppDists::rinvGauss(n_full, nu = postcpred_astar[full_index], lambda = lambda)
        rm(lambda)
        # poisson L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "poisson") {
        mid_a <- rpois(n_full, lambda = postcpred_a[full_index])
        mid_astar <- rpois(n_full, lambda = postcpred_astar[full_index])
        # quasipoisson L
      } else if (is_glm_postcreg[p] && family_postcreg[[p]]$family == "quasipoisson") {
        phi <- summary(postcreg[[p]])$dispersion
        mid_a <- rqpois(n_full, lambda = postcpred_a[full_index], phi = phi)
        mid_astar <- rqpois(n_full, lambda = postcpred_astar[full_index], phi = phi)
        rm(phi)
        # negative binomial L
      } else if (is_glm_postcreg[p] && startsWith(family_postcreg[[p]]$family, "Negative Binomial")) {
        theta <- summary(postcreg[[p]])$theta
        mid_a <- MASS::rnegbin(n_full, mu = postcpred_a[full_index], theta = theta)
        mid_astar <- MASS::rnegbin(n_full, mu = postcpred_astar[full_index], theta = theta)
        rm(theta)
      } else stop(paste0("Unsupported postcreg[[", p, "]]"))
      
      postc_a[full_index, p] <- mid_a
      postc_astar[full_index, p] <- mid_astar
      
      if (is.factor(data[, postc[p]])) {
        l_lev <- levels(droplevels(as.factor(data[, postc[p]])))
        postc_a[, p] <- factor(postc_a[, p], levels = l_lev)
        postc_astar[, p] <- factor(postc_astar[, p], levels = l_lev)
      }
    }
    
    rm(postcdesign_a, postcdesign_astar, type, postcpred_a, postcpred_astar, mid_a, mid_astar, full_index, n_full)
  }
  
  # design matrices for simulating mediator[p]
  mdesign_a <- data.frame(a_sim, basec_sim, postc_a)
  mdesign_astar <- data.frame(astar_sim, basec_sim, postc_astar)
  colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, basec, postc)
  m_a <- m_astar <- data.frame(matrix(nrow = n, ncol = length(mediator)))
  colnames(m_a) <- colnames(m_astar) <- mediator
  
  # simulating mediator[p]
  for (p in 1:length(mediator)) {
    # predict mediator[p]
    type <- ifelse(is_multinom_mreg[p] | is_polr_mreg[p], "probs", "response")
    mpred_a <- predict(mreg[[p]], newdata = mdesign_a, type = type)
    mpred_astar <- predict(mreg[[p]], newdata = mdesign_astar, type = type)
    full_index <- which(rowSums(is.na(mdesign_a))==0)
    n_full <- length(full_index)
    # categorical M
    if ((is_glm_mreg[p] && ((family_mreg[[p]]$family %in% c("binomial", "multinom")) |
                            startsWith(family_mreg[[p]]$family, "Ordered Categorical")))|
        is_multinom_mreg[p] | is_polr_mreg[p]) {
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      prob_a <- as.matrix(mpred_a)
      prob_astar <- as.matrix(mpred_astar)
      if (dim(prob_a)[2] == 1) {
        # simulate mediator[p] for exposure=a
        mid_a <- m_lev[rbinom(n_full, size = 1, prob = prob_a[full_index, 1]) + 1]
        # simulate mediator[p] for exposure=astar
        mid_astar <- m_lev[rbinom(n_full, size = 1, prob = prob_astar[full_index, 1]) + 1]
      } else {
        mid_a <- m_lev[apply(prob_a[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
        mid_astar <- m_lev[apply(prob_astar[full_index,], 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))]
      }
      
      if (is.numeric(data[, mediator[p]])) {
        mid_a <- as.numeric(mid_a)
        mid_astar <- as.numeric(mid_astar)
      }
      
      rm(prob_a, prob_astar, m_lev)
      # linear M
    } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "gaussian") {
      error <- rnorm(n_full, mean = 0, sd = sigma(mreg[[p]]))
      mid_a <- mpred_a[full_index] + error
      mid_astar <- mpred_astar[full_index] + error
      rm(error)
      # gamma M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "Gamma") {
      shape_mreg <- MASS::gamma.shape(mreg[[p]])$alpha
      mid_a <- rgamma(n_full, shape = shape_mreg, scale = mpred_a[full_index]/shape_mreg)
      mid_astar <- rgamma(n_full, shape = shape_mreg, scale = mpred_astar[full_index]/shape_mreg)
      rm(shape_mreg)
      # inverse gaussian M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "inverse.gaussian") {
      lambda <- 1/summary(mreg[[p]])$dispersion
      mid_a <- SuppDists::rinvGauss(n_full, nu = mpred_a[full_index], lambda = lambda)
      mid_astar <- SuppDists::rinvGauss(n_full, nu = mpred_astar[full_index], lambda = lambda)
      rm(lambda)
      # poisson M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "poisson") {
      mid_a <- rpois(n_full, lambda = mpred_a[full_index])
      mid_astar <- rpois(n_full, lambda = mpred_astar[full_index])
      # quasipoisson M
    } else if (is_glm_mreg[p] && family_mreg[[p]]$family == "quasipoisson") {
      phi <- summary(mreg[[p]])$dispersion
      mid_a <- rqpois(n_full, lambda = mpred_a[full_index], phi = phi)
      mid_astar <- rqpois(n_full, lambda = mpred_astar[full_index], phi = phi)
      rm(phi)
      # negative binomial M
    } else if ( is_glm_mreg[p] && startsWith(family_mreg[[p]]$family, "Negative Binomial")) {
      theta <- summary(mreg[[p]])$theta
      mid_a <- MASS::rnegbin(n_full, mu = mpred_a[full_index], theta = theta)
      mid_astar <- MASS::rnegbin(n_full, mu = mpred_astar[full_index], theta = theta)
      rm(theta)
    } else stop(paste0("Unsupported mreg[[", p, "]]"))
    
    # randomly shuffle values of simulated mediator[p] if postc is not empty
    if (length(postc) != 0) {
      m_a[full_index, p] <- sample(mid_a, replace = FALSE)
      m_astar[full_index, p] <- sample(mid_astar, replace = FALSE)
    } else {
      m_a[full_index, p] <- mid_a
      m_astar[full_index, p] <- mid_astar
    }
    
    if (is.factor(data[, mediator[p]])) {
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      m_a[, p] <- factor(m_a[, p], levels = m_lev)
      m_astar[, p] <- factor(m_astar[, p], levels = m_lev)
    }
    
  }
  rm(mdesign_a, mdesign_astar, type, mpred_a, mpred_astar, mid_a, mid_astar, full_index, n_full)
  
  # simulate mstar for cde
  mstar_sim <- do.call(cbind, lapply(1:length(mediator), function(x)
    if (is.factor(data[, mediator[x]])) {
      data.frame(factor(rep(mval[[x]], n), levels = levels(data[, mediator[x]])))
    } else data.frame(rep(mval[[x]], n))))
  
  # design matrices for outcome simulation
  ydesign0m <- data.frame(astar_sim, mstar_sim, basec_sim, postc_astar)
  ydesign1m <- data.frame(a_sim, mstar_sim, basec_sim, postc_a)
  ydesign00 <- data.frame(astar_sim, m_astar, basec_sim, postc_astar)
  ydesign01 <- data.frame(astar_sim, m_a, basec_sim, postc_astar)
  ydesign10 <- data.frame(a_sim, m_astar, basec_sim, postc_a)
  ydesign11 <- data.frame(a_sim, m_a, basec_sim, postc_a)
  rm(a_sim, astar_sim, m_a, m_astar, mstar_sim, basec_sim, postc_a, postc_astar)
  colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
    colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, basec, postc)
  
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
  weightsEY <- as.vector(model.frame(yreg, na.action = na.pass)$'(weights)')
  if (is.null(weightsEY)) weightsEY <- rep(1, n)
  
  # categorical Y
  if ((is_glm_yreg && ((family_yreg$family %in% c("binomial", "quasibinomial", "multinom")) |
                       startsWith(family_yreg$family, "Ordered Categorical")))|
      is_multinom_yreg | is_polr_yreg) {
    if (!is.null(yval_index)) {
      if (dim(EY0m_pred)[2] == 1) {
        EY0m <- weighted_mean(cbind(1 - EY0m_pred, EY0m_pred)[, yval_index], w = weightsEY)
        EY1m <- weighted_mean(cbind(1 - EY1m_pred, EY1m_pred)[, yval_index], w = weightsEY)
        EY00 <- weighted_mean(cbind(1 - EY00_pred, EY00_pred)[, yval_index], w = weightsEY)
        EY01 <- weighted_mean(cbind(1 - EY01_pred, EY01_pred)[, yval_index], w = weightsEY)
        EY10 <- weighted_mean(cbind(1 - EY10_pred, EY10_pred)[, yval_index], w = weightsEY)
        EY11 <- weighted_mean(cbind(1 - EY11_pred, EY11_pred)[, yval_index], w = weightsEY)
      } else {
        EY0m <- weighted_mean(EY0m_pred[, yval_index], w = weightsEY)
        EY1m <- weighted_mean(EY1m_pred[, yval_index], w = weightsEY)
        EY00 <- weighted_mean(EY00_pred[, yval_index], w = weightsEY)
        EY01 <- weighted_mean(EY01_pred[, yval_index], w = weightsEY)
        EY10 <- weighted_mean(EY10_pred[, yval_index], w = weightsEY)
        EY11 <- weighted_mean(EY11_pred[, yval_index], w = weightsEY)
      }
    } else EY0m <- EY1m <- EY00 <- EY01 <- EY10 <- EY11 <- 0
  } else {
    # non-categorical Y
    EY0m <- weighted_mean(EY0m_pred, w = weightsEY)
    EY1m <- weighted_mean(EY1m_pred, w = weightsEY)
    EY00 <- weighted_mean(EY00_pred, w = weightsEY)
    EY01 <- weighted_mean(EY01_pred, w = weightsEY)
    EY10 <- weighted_mean(EY10_pred, w = weightsEY)
    EY11 <- weighted_mean(EY11_pred, w = weightsEY)
  }
  rm(weightsEY, EY0m_pred, EY1m_pred, EY00_pred, EY01_pred, EY10_pred, EY11_pred)
  
  # output causal effects in additive scale for continuous Y
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
      if (EMint) {
        intref <- pnde - cde
        intmed <- tnie - pnie
        cde_prop <- cde/te
        intref_prop <- intref/te
        intmed_prop <- intmed/te
        pnie_prop <- pnie/te
        int <- (intref + intmed)/te
        pe <- (intref + intmed + pnie)/te
        est <- c(cde, pnde, tnde, pnie, tnie, te, intref, intmed, 
                 cde_prop, intref_prop, intmed_prop, pnie_prop, pm, int, pe)
      } else est <- c(cde, pnde, tnde, pnie, tnie, te, pm)
    } else est <- c(cde, pnde, tnde, pnie, tnie, te)
  } else {
    # output causal effects in ratio scale for non-continuous Y
    
    ## output effects on the odds ratio scale for logistic regressions
    if (is_glm_yreg && family_yreg$family %in% c("binomial", "quasibinomial") &&
        yreg$family$link == "logit") {
      logRRcde <- log(EY1m/(1-EY1m)) - log(EY0m/(1-EY0m))
      logRRpnde <- log(EY10/(1-EY10)) - log(EY00/(1-EY00))
      logRRtnde <- log(EY11/(1-EY11)) - log(EY01/(1-EY01))
      logRRpnie <- log(EY01/(1-EY01)) - log(EY00/(1-EY00))
      logRRtnie <- log(EY11/(1-EY11)) - log(EY10/(1-EY10))
      ## otherwise on the risk ratio scale
    } else {
      logRRcde <- log(EY1m) - log(EY0m)
      logRRpnde <- log(EY10) - log(EY00)
      logRRtnde <- log(EY11) - log(EY01)
      logRRpnie <- log(EY01) - log(EY00)
      logRRtnie <- log(EY11) - log(EY10)
    }
    
    logRRte <- logRRtnie + logRRpnde
    if (full) {
      pm <- (exp(logRRpnde) * (exp(logRRtnie) - 1)) / (exp(logRRte) - 1)
      if (EMint) {
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
      } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, pm)
    } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte)
    
  } 
  
  # progress bar
  if (!multimp) {
    curVal <- get("counter", envir = env)
    assign("counter", curVal + 1, envir = env)
    setTxtProgressBar(get("progbar", envir = env), curVal + 1)
  }
  if (outReg) out$est <- est
  if (!outReg) out <- est
  return(out)
}

