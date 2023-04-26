est.msm <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {
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
    if (length(basec) != 0) {
      # weights for ereg
      if (!is.null(weights_ereg)) weights_ereg <- weights_ereg[indices] * w4casecon
      if (is.null(weights_ereg)) weights_ereg <- w4casecon
      # update ereg
      call_ereg$weights <- weights_ereg
      call_ereg$data <- data
      if (outReg && (inherits(ereg, "rcreg") | inherits(ereg, "simexreg"))) call_ereg$variance <- TRUE
      ereg <- eval.parent(call_ereg)
      # calculate w_{a,i}=P(A=A_i)/P(A=A_i|C_i)
      wanom <- left_join(select(data, exposure), count(data, !!as.name(exposure)), by = exposure)[, "n"]/n
      wadenom_prob <- as.matrix(predict(ereg, newdata = data,
                                        type = ifelse(is_multinom_ereg | is_polr_ereg, "probs", "response")))
      a_lev <- levels(droplevels(as.factor(data[, exposure])))
      wa_data <- data[, exposure, drop = FALSE]
      wa_data[, exposure] <- factor(wa_data[, exposure], levels = a_lev)
      if (dim(wadenom_prob)[2] == 1) {
        category <- as.numeric(wa_data[, 1]) - 1
        wadenom <-  wadenom_prob[, 1] ^ category * (1 - wadenom_prob[, 1]) ^ (1 - category)
      } else {
        category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)
        wadenom <- rowSums(category * wadenom_prob)
      }
      wa <- as.vector(wanom / wadenom)
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)
    
    # calculate w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/
    # P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
    wmdenom <- wmnom <- matrix(nrow = n, ncol = length(mediator))
    for (p in 1:length(mediator)) {
      # weights for wmdenomreg[[p]] and wmnomreg[[p]]
      if (!is.null(weights_wmnomreg[[p]])) weights_wmnomreg[[p]] <- weights_wmnomreg[[p]][indices] * w4casecon
      if (is.null(weights_wmnomreg[[p]])) weights_wmnomreg[[p]] <- w4casecon      
      if (!is.null(weights_wmdenomreg[[p]])) weights_wmdenomreg[[p]] <- weights_wmdenomreg[[p]][indices] * w4casecon
      if (is.null(weights_wmdenomreg[[p]])) weights_wmdenomreg[[p]] <- w4casecon
      # update wmdenomreg[[p]] and wmnomreg[[p]]
      call_wmdenomreg[[p]]$weights <- weights_wmdenomreg[[p]]
      call_wmdenomreg[[p]]$data <- data
      if (outReg && (inherits(wmdenomreg[[p]], "rcreg") | inherits(wmdenomreg[[p]], "simexreg"))) call_wmdenomreg[[p]]$variance <- TRUE
      call_wmnomreg[[p]]$weights <- weights_wmnomreg[[p]]
      call_wmnomreg[[p]]$data <- data
      if (outReg && (inherits(wmnomreg[[p]], "rcreg") | inherits(wmnomreg[[p]], "simexreg"))) call_wmnomreg[[p]]$variance <- TRUE
      wmdenomreg[[p]] <- eval.parent(call_wmdenomreg[[p]])
      wmnomreg[[p]] <- eval.parent(call_wmnomreg[[p]])
      
      # calculate w_{m_p,i}
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      wm_data <- data[, mediator[p], drop = FALSE]
      wm_data[, mediator[p]] <- factor(wm_data[, mediator[p]], levels = m_lev)
      wmdenom_prob <- as.matrix(predict(wmdenomreg[[p]], newdata = data, 
                                        type = ifelse(is_multinom_wmdenomreg[p] | is_polr_wmdenomreg[p], "probs", "response")))
      wmnom_prob <- as.matrix(predict(wmnomreg[[p]], newdata = data, 
                                      type = ifelse(is_multinom_wmnomreg[p] | is_polr_wmnomreg[p], "probs", "response")))
      if (dim(wmdenom_prob)[2] == 1) {
        category <- as.numeric(wm_data[, 1]) - 1
        wmdenom[, p] <-  wmdenom_prob[, 1] ^ category * (1 - wmdenom_prob[, 1]) ^ (1 - category)
        wmnom[, p] <-  wmnom_prob[, 1] ^ category * (1 - wmnom_prob[, 1]) ^ (1 - category)
      } else {
        category <- model.matrix(as.formula(paste("~0+", mediator[p], sep = "")), data = wm_data)
        wmdenom[, p] <- rowSums(category * wmdenom_prob)
        wmnom[, p] <- rowSums(category * wmnom_prob)
      }
    }
    wm <- wmnom / wmdenom
    rm(m_lev, wm_data, wmdenom_prob, wmnom_prob, category, wmdenom, wmnom)
    
    # weights for yreg
    w_yreg <- wa
    for (p in 1:ncol(wm)) w_yreg <- w_yreg * wm[, p]
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon * w_yreg
    if (is.null(weights_yreg)) weights_yreg <- w4casecon * w_yreg
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # weights for mreg[[p]]
    for (p in 1:length(mediator)) {
      if (!is.null(weights_mreg[[p]])) weights_mreg[[p]] <- weights_mreg[[p]][indices] * w4casecon * wa
      if (is.null(weights_mreg[[p]])) weights_mreg[[p]] <- w4casecon * wa
      # update mreg[[p]]
      call_mreg[[p]]$weights <- weights_mreg[[p]]
      call_mreg[[p]]$data <- data
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    rm(prob1, w4casecon)
    
  } else if (casecontrol && yrare) {
    # method 2 for a case control design
    # data among controls
    control_indices <- which(data[, outcome] == y_control)
    
    if (length(basec) != 0) {
      # update ereg
      call_ereg$weights <- weights_ereg[indices][control_indices]
      call_ereg$data <- data[control_indices, ]
      if (outReg && (inherits(ereg, "rcreg") | inherits(ereg, "simexreg"))) call_ereg$variance <- TRUE
      ereg <- eval.parent(call_ereg)
      # calculate w_{a,i}=P(A=A_i)/P(A=A_i|C_i)
      wanom <- left_join(select(data, exposure), count(data, !!as.name(exposure)), by = exposure)[, "n"]/n
      wadenom_prob <- as.matrix(predict(ereg, newdata = data,
                                        type = ifelse(is_multinom_ereg | is_polr_ereg, "probs", "response")))
      a_lev <- levels(droplevels(as.factor(data[, exposure])))
      wa_data <- data[, exposure, drop = FALSE]
      wa_data[, exposure] <- factor(wa_data[, exposure], levels = a_lev)
      if (dim(wadenom_prob)[2] == 1) {
        category <- as.numeric(wa_data[, 1]) - 1
        wadenom <-  wadenom_prob[, 1] ^ category * (1 - wadenom_prob[, 1]) ^ (1 - category)
      } else {
        category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)
        wadenom <- rowSums(category * wadenom_prob)
      }
      wa <- as.vector(wanom / wadenom)
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)
    
    # calculate w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/
    # P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
    wmdenom <- wmnom <- matrix(nrow = n, ncol = length(mediator))
    for (p in 1:length(mediator)) {
      # update wmdenomreg[[p]] and wmnomreg[[p]]
      call_wmdenomreg[[p]]$weights <- weights_wmdenomreg[[p]][indices][control_indices]
      call_wmdenomreg[[p]]$data <- data[control_indices, ]
      if (outReg && (inherits(wmdenomreg[[p]], "rcreg") | inherits(wmdenomreg[[p]], "simexreg"))) call_wmdenomreg[[p]]$variance <- TRUE
      call_wmnomreg[[p]]$weights <- weights_wmnomreg[[p]][indices][control_indices]
      call_wmnomreg[[p]]$data <- data[control_indices, ]
      if (outReg && (inherits(wmnomreg[[p]], "rcreg") | inherits(wmnomreg[[p]], "simexreg"))) call_wmnomreg[[p]]$variance <- TRUE
      wmdenomreg[[p]] <- eval.parent(call_wmdenomreg[[p]])
      wmnomreg[[p]] <- eval.parent(call_wmnomreg[[p]])
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      wm_data <- data[, mediator[p], drop = FALSE]
      wm_data[, mediator[p]] <- factor(wm_data[, mediator[p]], levels = m_lev)
      wmdenom_prob <- as.matrix(predict(wmdenomreg[[p]], newdata = data, 
                                        type = ifelse(is_multinom_wmdenomreg[p] | is_polr_wmdenomreg[p], "probs", "response")))
      wmnom_prob <- as.matrix(predict(wmnomreg[[p]], newdata = data, 
                                      type = ifelse(is_multinom_wmnomreg[p] | is_polr_wmnomreg[p], "probs", "response")))
      if (dim(wmdenom_prob)[2] == 1) {
        category <- as.numeric(wm_data[, 1]) - 1
        wmdenom[, p] <-  wmdenom_prob[, 1] ^ category * (1 - wmdenom_prob[, 1]) ^ (1 - category)
        wmnom[, p] <-  wmnom_prob[, 1] ^ category * (1 - wmnom_prob[, 1]) ^ (1 - category)
      } else {
        category <- model.matrix(as.formula(paste("~0+", mediator[p], sep = "")), data = wm_data)
        wmdenom[, p] <- rowSums(category * wmdenom_prob)
        wmnom[, p] <- rowSums(category * wmnom_prob)
      }
    }
    wm <- wmnom / wmdenom
    rm(m_lev, wm_data, wmdenom_prob, wmnom_prob, category, wmdenom, wmnom)
    
    # weights for yreg
    w_yreg <- wa
    for (p in 1:ncol(wm)) w_yreg <- w_yreg * wm[, p]
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w_yreg
    if (is.null(weights_yreg)) weights_yreg <- w_yreg
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # weights for mreg[[p]]
    for (p in 1:length(mediator)) {
      if (!is.null(weights_mreg[[p]])) weights_mreg[[p]] <- (weights_mreg[[p]][indices] * wa)[control_indices]
      if (is.null(weights_mreg[[p]])) weights_mreg[[p]] <- wa[control_indices]
      # update mreg[[p]]
      call_mreg[[p]]$weights <- weights_mreg[[p]]
      call_mreg[[p]]$data <- data[control_indices, ]
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
    rm(control_indices)
    
  } else {
    # not a case control design
    if (length(basec) != 0) {
      # update ereg
      call_ereg$weights <- weights_ereg[indices]
      call_ereg$data <- data
      if (outReg && (inherits(ereg, "rcreg") | inherits(ereg, "simexreg"))) call_ereg$variance <- TRUE
      ereg <- eval.parent(call_ereg)
      # calculate w_{a,i}=P(A=A_i)/P(A=A_i|C_i)
      wanom <- left_join(select(data, exposure), count(data, !!as.name(exposure)), by = exposure)[, "n"]/n
      wadenom_prob <- as.matrix(predict(ereg, newdata = data,
                                        type = ifelse(is_multinom_ereg | is_polr_ereg, "probs", "response")))
      a_lev <- levels(droplevels(as.factor(data[, exposure])))
      wa_data <- data[, exposure, drop = FALSE]
      wa_data[, exposure] <- factor(wa_data[, exposure], levels = a_lev)
      if (dim(wadenom_prob)[2] == 1) {
        category <- as.numeric(wa_data[, 1]) - 1
        wadenom <-  wadenom_prob[, 1] ^ category * (1 - wadenom_prob[, 1]) ^ (1 - category)
      } else {
        category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)
        wadenom <- rowSums(category * wadenom_prob)
      }
      wa <- as.vector(wanom / wadenom)
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)
    
    # calculate w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/
    # P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
    wmdenom <- wmnom <- matrix(nrow = n, ncol = length(mediator))
    for (p in 1:length(mediator)) {
      # update wmdenomreg[[p]] and wmnomreg[[p]]
      call_wmdenomreg[[p]]$weights <- weights_wmdenomreg[[p]][indices]
      call_wmdenomreg[[p]]$data <- data
      if (outReg && (inherits(wmdenomreg[[p]], "rcreg") | inherits(wmdenomreg[[p]], "simexreg"))) call_wmdenomreg[[p]]$variance <- TRUE
      call_wmnomreg[[p]]$weights <- weights_wmnomreg[[p]][indices]
      call_wmnomreg[[p]]$data <- data
      if (outReg && (inherits(wmnomreg[[p]], "rcreg") | inherits(wmnomreg[[p]], "simexreg"))) call_wmnomreg[[p]]$variance <- TRUE
      wmdenomreg[[p]] <- eval.parent(call_wmdenomreg[[p]])
      wmnomreg[[p]] <- eval.parent(call_wmnomreg[[p]])
      m_lev <- levels(droplevels(as.factor(data[, mediator[p]])))
      wm_data <- data[, mediator[p], drop = FALSE]
      wm_data[, mediator[p]] <- factor(wm_data[, mediator[p]], levels = m_lev)
      wmdenom_prob <- as.matrix(predict(wmdenomreg[[p]], newdata = data, 
                                        type = ifelse(is_multinom_wmdenomreg[p] | is_polr_wmdenomreg[p], "probs", "response")))
      wmnom_prob <- as.matrix(predict(wmnomreg[[p]], newdata = data, 
                                      type = ifelse(is_multinom_wmnomreg[p] | is_polr_wmnomreg[p], "probs", "response")))
      if (dim(wmdenom_prob)[2] == 1) {
        category <- as.numeric(wm_data[, 1]) - 1
        wmdenom[, p] <-  wmdenom_prob[, 1] ^ category * (1 - wmdenom_prob[, 1]) ^ (1 - category)
        wmnom[, p] <-  wmnom_prob[, 1] ^ category * (1 - wmnom_prob[, 1]) ^ (1 - category)
      } else {
        category <- model.matrix(as.formula(paste("~0+", mediator[p], sep = "")), data = wm_data)
        wmdenom[, p] <- rowSums(category * wmdenom_prob)
        wmnom[, p] <- rowSums(category * wmnom_prob)
      }
    }
    wm <- wmnom / wmdenom
    rm(m_lev, wm_data, wmdenom_prob, wmnom_prob, category, wmdenom, wmnom)
    
    # weights for outcome regression
    w_yreg <- wa
    for (p in 1:ncol(wm)) w_yreg <- w_yreg * wm[, p]
    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w_yreg
    if (is.null(weights_yreg)) weights_yreg <- w_yreg
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
    # weights for mreg[[p]]
    for (p in 1:length(mediator)) {
      if (!is.null(weights_mreg[[p]])) weights_mreg[[p]] <- weights_mreg[[p]][indices] * wa
      if (is.null(weights_mreg[[p]])) weights_mreg[[p]] <- wa
      # update mreg[[p]]
      call_mreg[[p]]$weights <- weights_mreg[[p]]
      call_mreg[[p]]$data <- data
      if (outReg && (inherits(mreg[[p]], "rcreg") | inherits(mreg[[p]], "simexreg"))) call_mreg[[p]]$variance <- TRUE
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }
  }
  
  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yreg <- yreg
    out$reg.output$mreg <- mreg
    out$reg.output$wmdenomreg <- wmdenomreg
    out$reg.output$wmnomreg <- wmnomreg
    if (length(basec) != 0) out$reg.output$ereg <- ereg
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
  
  # design matrices for simulating mediator[p]
  mdesign_a <- data.frame(a_sim)
  mdesign_astar <- data.frame(astar_sim)
  colnames(mdesign_a) <- colnames(mdesign_astar) <- exposure
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
    } else stop(paste0("Unsupported mreg[[", p, "]]"))
    
    m_a[full_index, p] <- mid_a
    m_astar[full_index, p] <- mid_astar
    
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
  ydesign0m <- data.frame(astar_sim, mstar_sim)
  ydesign1m <- data.frame(a_sim, mstar_sim)
  ydesign00 <- data.frame(astar_sim, m_astar)
  ydesign01 <- data.frame(astar_sim, m_a)
  ydesign10 <- data.frame(a_sim, m_astar)
  ydesign11 <- data.frame(a_sim, m_a)
  rm(a_sim, astar_sim, m_a, m_astar, mstar_sim)
  colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
    colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator)
  
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
      if (EMint) {
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
      } else est <- c(cde, pnde, tnde, pnie, tnie, te, pm)
    } else est <- c(cde, pnde, tnde, pnie, tnie, te)
  } else {
    # output causal effects on the ratio scale for non-continuous Y
    
    ## output effects on the odds ratio scale for logistic regressions
    if (is_glm_yreg && family_yreg$family %in% c("binomial", "quasibinomial") &&
        family_yreg$link == "logit") {
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

