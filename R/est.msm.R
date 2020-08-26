est.msm <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {

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

    if (length(basec) != 0) {
      # weights for ereg
      if (!is.null(weights_ereg)) weights_ereg <- weights_ereg[indices] * w4casecon
      if (is.null(weights_ereg)) weights_ereg <- w4casecon
      # update ereg
      call_ereg$weights <- weights_ereg
      call_ereg$data <- data
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
      wa <- wanom / wadenom
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)

    # calculate w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/
    # P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
    wmdenom <- wmnom <- matrix(nrow = n, ncol = length(mediator))
    wmdenomreg <- wmnomreg <- list()
    for (p in 1:length(mediator)) {
      call_wmdenomreg <- call_wmreg[[p]]
      call_wmnomreg <- call_mreg[[p]]
      # weights for wmdenomreg[[p]] and wmnomreg[[p]]
      if (!is.null(weights_wmreg[[p]])) weights_wmdenomreg <- weights_wmreg[[p]][indices] * w4casecon
      if (is.null(weights_wmreg[[p]])) weights_wmdenomreg <- w4casecon
      if (!is.null(weights_mreg[[p]])) weights_wmnomreg <- weights_mreg[[p]][indices] * w4casecon
      if (is.null(weights_mreg[[p]])) weights_wmnomreg <- w4casecon
      # update wmdenomreg[[p]] and wmnomreg[[p]]
      call_wmdenomreg$weights <- weights_wmdenomreg
      call_wmdenomreg$data <- data
      call_wmnomreg$weights <- weights_wmnomreg
      call_wmnomreg$data <- data
      wmdenomreg[[p]] <- eval.parent(call_wmdenomreg)
      wmnomreg[[p]] <- eval.parent(call_wmnomreg)
      m_lev <- levels(droplevels(as.factor(data_orig[, mediator[p]])))
      wm_data <- data[, mediator[p], drop = FALSE]
      wm_data[, mediator[p]] <- factor(wm_data[, mediator[p]], levels = m_lev)
      type <- ifelse(is_multinom_wmreg[p] | is_polr_wmreg[p], "probs", "response")
      wmdenom_prob <- as.matrix(predict(wmdenomreg[[p]], newdata = data, type = type))
      wmnom_prob <- as.matrix(predict(wmnomreg[[p]], newdata = data, type = type))
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
    rm(call_wmdenomreg, call_wmnomreg, m_lev, wm_data, type, wmdenom_prob, wmnom_prob,
       category, wmdenom, wmnom)
    # weights for each mediator regression
    w_mreg <- replicate(length(mediator), wa)
    # weights for outcome regression
    w_yreg <- wa
    for (p in 1:ncol(wm)) w_yreg <- w_yreg * wm[, p]

    for (p in 1:length(mediator)) {
      # weights for mreg[[p]]
      if (!is.null(weights_mreg[[p]])) weights_mreg <- weights_mreg[[p]][indices] * w4casecon * w_mreg[, p]
      if (is.null(weights_mreg[[p]])) weights_mreg <- w4casecon * w_mreg[, p]
      # update mreg[[p]]
      call_mreg[[p]]$weights <- weights_mreg
      call_mreg[[p]]$data <- data
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }

    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon * w_yreg
    if (is.null(weights_yreg)) weights_yreg <- w4casecon * w_yreg
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)

    rm(prob1, w4casecon)

  } else if (casecontrol && yrare) {

    # method 2 for a case control design
    # data from controls
    control_indices <- which(data[, outcome] == y_control)

    if (length(basec) != 0) {
      # update ereg
      call_ereg$weights <- weights_ereg[indices][control_indices]
      call_ereg$data <- data[control_indices, ]
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
      wa <- wanom / wadenom
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)

    # calculate w_{m_p,i}=P(M_p=M_{p,i}|A=A_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})/
    # P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
    wmdenom <- wmnom <- matrix(nrow = n, ncol = length(mediator))
    wmdenomreg <- wmnomreg <- list()
    for (p in 1:length(mediator)) {
      call_wmdenomreg <- call_wmreg[[p]]
      call_wmnomreg <- call_mreg[[p]]
      # update wmdenomreg[[p]] and wmnomreg[[p]]
      call_wmdenomreg$weights <- weights_wmreg[[p]][indices][control_indices]
      call_wmdenomreg$data <- data[control_indices, ]
      call_wmnomreg$weights <- weights_mreg[[p]][indices][control_indices]
      call_wmnomreg$data <- data[control_indices, ]
      wmdenomreg[[p]] <- eval.parent(call_wmdenomreg)
      wmnomreg[[p]] <- eval.parent(call_wmnomreg)
      m_lev <- levels(droplevels(as.factor(data_orig[, mediator[p]])))
      wm_data <- data[, mediator[p], drop = FALSE]
      wm_data[, mediator[p]] <- factor(wm_data[, mediator[p]], levels = m_lev)
      type <- ifelse(is_multinom_wmreg[p] | is_polr_wmreg[p], "probs", "response")
      wmdenom_prob <- as.matrix(predict(wmdenomreg[[p]], newdata = data, type = type))
      wmnom_prob <- as.matrix(predict(wmnomreg[[p]], newdata = data, type = type))
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
    rm(call_wmdenomreg, call_wmnomreg, m_lev, wm_data, type, wmdenom_prob, wmnom_prob,
       category, wmdenom, wmnom)
    # weights for each mediator regression
    w_mreg <- replicate(length(mediator), wa)
    # weights for outcome regression
    w_yreg <- wa
    for (p in 1:ncol(wm)) w_yreg <- w_yreg * wm[, p]

    for (p in 1:length(mediator)) {
      # weights for mreg[[p]]
      if (!is.null(weights_mreg[[p]])) w.mid <- (weights_mreg[[p]][indices] * w_mreg[, p])[control_indices]
      if (is.null(weights_mreg[[p]])) w.mid <- w_mreg[, p][control_indices]
      # update mreg[[p]]
      call_mreg[[p]]$weights <- w.mid
      call_mreg[[p]]$data <- data[control_indices, ]
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }

    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w_yreg
    if (is.null(weights_yreg)) weights_yreg <- w_yreg
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)

    rm(control_indices)

  } else {

    # not a case control design
    if (length(basec) != 0) {
      # update ereg
      call_ereg$weights <- weights_ereg[indices]
      call_ereg$data <- data
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
      wa <- wanom / wadenom
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)

    # calculate w_{m_p,i}=P(M_p=M_{p,i}|A=A_i)/P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
    wmdenom <- wmnom <- matrix(nrow = n, ncol = length(mediator))
    wmdenomreg <- wmnomreg <- list()
    for (p in 1:length(mediator)) {
      call_wmdenomreg <- call_wmreg[[p]]
      call_wmnomreg <- call_mreg[[p]]
      # update wmdenomreg[[p]] and wmnomreg[[p]]
      call_wmdenomreg$weights <- weights_wmreg[[p]][indices]
      call_wmdenomreg$data <- data
      call_wmnomreg$weights <- weights_mreg[[p]][indices]
      call_wmnomreg$data <- data
      wmdenomreg[[p]] <- eval.parent(call_wmdenomreg)
      wmnomreg[[p]] <- eval.parent(call_wmnomreg)
      m_lev <- levels(droplevels(as.factor(data_orig[, mediator[p]])))
      wm_data <- data[, mediator[p], drop = FALSE]
      wm_data[, mediator[p]] <- factor(wm_data[, mediator[p]], levels = m_lev)
      type <- ifelse(is_multinom_wmreg[p] | is_polr_wmreg[p], "probs", "response")
      wmdenom_prob <- as.matrix(predict(wmdenomreg[[p]], newdata = data, type = type))
      wmnom_prob <- as.matrix(predict(wmnomreg[[p]], newdata = data, type = type))
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
    rm(call_wmdenomreg, call_wmnomreg, m_lev, wm_data, type, wmdenom_prob, wmnom_prob,
       category, wmdenom, wmnom)
    # weights for each mediator regression
    w_mreg <- replicate(length(mediator), wa)
    # weights for outcome regression
    w_yreg <- wa
    for (p in 1:ncol(wm)) w_yreg <- w_yreg * wm[, p]

    for (p in 1:length(mediator)) {
      # weights for mreg[[p]]
      if (!is.null(weights_mreg[[p]])) w.mid <- (weights_mreg[[p]][indices] * w_mreg[, p])
      if (is.null(weights_mreg[[p]])) w.mid <- w_mreg[, p]
      # update mreg[[p]]
      call_mreg[[p]]$weights <- w.mid
      call_mreg[[p]]$data <- data
      mreg[[p]] <- eval.parent(call_mreg[[p]])
    }

    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w_yreg
    if (is.null(weights_yreg)) weights_yreg <- w_yreg
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)

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

  # design matrices for simulating mediator[1]
  mdesign_a <- data.frame(a_sim)
  mdesign_astar <- data.frame(astar_sim)
  colnames(mdesign_a) <- colnames(mdesign_astar) <- exposure

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
      mid_a <- rinvGauss(n, nu = mpred_a, lambda = lambda)
      mid_astar <- rinvGauss(n, nu = mpred_astar, lambda = lambda)

      rm(lambda)

      # poisson M
    } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && family_mreg[[p]]$family == "poisson") {

      mid_a <- rpois(n, lambda = mpred_a)
      mid_astar <- rpois(n, lambda = mpred_astar)

      # negative binomial M
    } else if ((is_lm_mreg[p] | is_glm_mreg[p]) && startsWith(family_reg[[p]]$family, "Negative Binomial")) {

      theta <- summary(mreg[[p]])$theta
      mid_a <- rnegbin(n, mu = mpred_a, theta = theta)
      mid_astar <- rnegbin(n, mu = mpred_astar, theta = theta)

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
      int <- (intref + intmed)/te
      pe <- (intref + intmed + pnie)/te
      est <- c(cde, pnde, tnde, pnie, tnie, te, intref, intmed, cde_prop, intref_prop, 
               intmed_prop, pnie_prop, pm, int, pe)

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
      int <- (ERRintref + ERRintmed)/ERRte
      pe <- (ERRintref + ERRintmed + ERRpnie)/ERRte
      est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte, 
               ERRcde, ERRintref, ERRintmed, ERRpnie,
               ERRcde_prop, ERRintref_prop, ERRintmed_prop, ERRpnie_prop,
               pm, int, pe)
    } else est <- c(logRRcde, logRRpnde, logRRtnde, logRRpnie, logRRtnie, logRRte)

  } else stop("Unsupported yreg")

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

