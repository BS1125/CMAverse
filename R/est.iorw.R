est.iorw <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {
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
    prob1 <- mean(data[, outcome] == y_case, na.rm = TRUE)
    w4casecon <- as.vector(ifelse(data[, outcome] == y_case, yprevalence / prob1, (1 - yprevalence) / (1 - prob1)))
    
    # weights for ereg
    if (!is.null(weights_ereg)) weights_ereg <- weights_ereg[indices] * w4casecon
    if (is.null(weights_ereg)) weights_ereg <- w4casecon
    # update ereg
    call_ereg$weights <- weights_ereg
    call_ereg$data <- data
    ereg <- eval.parent(call_ereg)
    # calculate w_{a,i}=P(A=0|M_i,C_i)/P(A=A_i|M_i,C_i)
    wadenom_prob <- as.matrix(predict(ereg, newdata = data,
                                      type = ifelse(is_multinom_ereg | is_polr_ereg, "probs", "response")))
    a_lev <- levels(droplevels(as.factor(data[, exposure])))
    wa_data <- data[, exposure, drop = FALSE]
    wa_data[, exposure] <- factor(wa_data[, exposure], levels = a_lev)
    if (dim(wadenom_prob)[2] == 1) {
      category <- as.numeric(wa_data[, 1]) - 1
      wadenom <-  wadenom_prob[, 1] ^ category * (1 - wadenom_prob[, 1]) ^ (1 - category)
      wanom <- 1 - wadenom_prob
    } else {
      category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)
      wadenom <- rowSums(category * wadenom_prob)
      wanom <- wadenom_prob[, 1]
    }
    wa <- as.vector(wanom / wadenom)
    rm(weights_ereg, wadenom_prob, a_lev, wa_data, category, wanom, wadenom)

    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon
    if (is.null(weights_yreg)) weights_yreg <- w4casecon
    # update yreg
    call_yreg_tot <- call_yreg
    call_yreg_tot$weights <- weights_yreg
    call_yreg_tot$data <- data
    yreg_tot <- eval.parent(call_yreg_tot)
    call_yreg_dir <- call_yreg
    call_yreg_dir$weights <- as.vector(weights_yreg * wa)
    call_yreg_dir$data <- data
    yreg_dir <- eval.parent(call_yreg_dir)
    rm(prob1, w4casecon, weights_ereg, weights_yreg)

  } else if (casecontrol && yrare) {
    # method 2 for a case control design
    # data from controls
    control_indices <- which(data[, outcome] == y_control)
    # update ereg
    call_ereg$weights <- as.vector(weights_ereg[indices][control_indices])
    call_ereg$data <- data[control_indices, ]
    ereg <- eval.parent(call_ereg)
    # calculate w_{a,i}=P(A=0|M_i,C_i)/P(A=A_i|M_i,C_i)
    wadenom_prob <- as.matrix(predict(ereg, newdata = data, 
                                      type = ifelse(is_multinom_ereg | is_polr_ereg, "probs", "response")))
    a_lev <- levels(droplevels(as.factor(data[, exposure])))
    wa_data <- data[, exposure, drop = FALSE]
    wa_data[, exposure] <- factor(wa_data[, exposure], levels = a_lev)
    if (dim(wadenom_prob)[2] == 1) {
      category <- as.numeric(wa_data[, 1]) - 1
      wadenom <-  wadenom_prob[, 1] ^ category * (1 - wadenom_prob[, 1]) ^ (1 - category)
      wanom <- 1 - wadenom_prob
    } else {
      category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)
      wadenom <- rowSums(category * wadenom_prob)
      wanom <- wadenom_prob[, 1]
    }
    wa <- as.vector(wanom / wadenom)
    rm(wadenom_prob, a_lev, wa_data, category, wanom, wadenom)

    # update yreg
    call_yreg_tot <- call_yreg
    call_yreg_tot$weights <- weights_yreg[indices]
    call_yreg_tot$data <- data
    yreg_tot <- eval.parent(call_yreg_tot)
    call_yreg_dir <- call_yreg
    if (!is.null(weights_yreg)) call_yreg_dir$weights <- weights_yreg[indices] * wa
    if (is.null(weights_yreg)) call_yreg_dir$weights <- wa
    call_yreg_dir$data <- data
    yreg_dir <- eval.parent(call_yreg_dir)
    rm(control_indices)

  } else {
    # not a case control design
    # update ereg
    call_ereg$weights <- weights_ereg[indices]
    call_ereg$data <- data
    ereg <- eval.parent(call_ereg)
    # calculate w_{a,i}=P(A=0|M_i,C_i)/P(A=A_i|M_i,C_i)
    wadenom_prob <- as.matrix(predict(ereg, newdata = data,
                                      type = ifelse(is_multinom_ereg | is_polr_ereg, "probs", "response")))
    a_lev <- levels(droplevels(as.factor(data[, exposure])))
    wa_data <- data[, exposure, drop = FALSE]
    wa_data[, exposure] <- factor(wa_data[, exposure], levels = a_lev)
    if (dim(wadenom_prob)[2] == 1) {
      category <- as.numeric(wa_data[, 1]) - 1
      wadenom <-  wadenom_prob[, 1] ^ category * (1 - wadenom_prob[, 1]) ^ (1 - category)
      wanom <- 1 - wadenom_prob
    } else {
      category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)
      wadenom <- rowSums(category * wadenom_prob)
      wanom <- wadenom_prob[, 1]
    }
    wa <- as.vector(wanom / wadenom)
    rm(wadenom_prob, a_lev, wa_data, category, wanom, wadenom)

    # update yreg
    call_yreg_tot <- call_yreg
    call_yreg_tot$weights <- weights_yreg[indices]
    call_yreg_tot$data <- data
    yreg_tot <- eval.parent(call_yreg_tot)
    call_yreg_dir <- call_yreg
    if (!is.null(weights_yreg)) call_yreg_dir$weights <- weights_yreg[indices] * wa
    if (is.null(weights_yreg)) call_yreg_dir$weights <- wa
    call_yreg_dir$data <- data
    yreg_dir <- eval.parent(call_yreg_dir)
  }

  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yregTot <- yreg_tot
    out$reg.output$yregDir <- yreg_dir
    out$reg.output$ereg <- ereg
  }

  # the index of the reference level for a categorical outcome
  if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                       startsWith(family_yreg$family, "Ordered Categorical"))) |
      is_multinom_yreg | is_polr_yreg) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    yref_index <- switch((yref %in% y_lev) + 1, "1" = NULL, "2" = which(y_lev == yref))
    rm(y_lev)
  }

  a_lev <- levels(droplevels(as.factor(data[, exposure])))
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

  # design matrices for outcome simulation
  totdesign0 <- dirdesign0 <- data.frame(astar_sim, basec_sim)
  totdesign1 <- dirdesign1 <- data.frame(a_sim, basec_sim)
  rm(a_sim, astar_sim, basec_sim)
  colnames(totdesign0) <- colnames(totdesign1) <-
    colnames(dirdesign0) <- colnames(dirdesign1) <- c(exposure, basec)

  # predict Y
  type <- ifelse(is_coxph_yreg, "risk", ifelse(is_multinom_yreg | is_polr_yreg, "probs", "response"))
  tot0_pred <- as.matrix(predict(yreg_tot, newdata = totdesign0, type = type))
  tot1_pred <- as.matrix(predict(yreg_tot, newdata = totdesign1, type = type))
  dir0_pred <- as.matrix(predict(yreg_dir, newdata = dirdesign0, type = type))
  dir1_pred <- as.matrix(predict(yreg_dir, newdata = dirdesign1, type = type))
  rm(type, totdesign0, totdesign1, dirdesign0, dirdesign1)

  # weights for calculating counterfactuals
  weightsEY_tot <- as.vector(model.frame(yreg_tot)$'(weights)')
  if (is.null(weightsEY_tot)) weightsEY_tot <- rep(1, n)
  weightsEY_dir <- as.vector(model.frame(yreg_dir)$'(weights)')
  if (is.null(weightsEY_dir)) weightsEY_dir <- rep(1, n)

  # categorical Y
  if ((is_glm_yreg && ((family_yreg$family %in% c("binomial", "quasibinomial", "multinom")) |
                       startsWith(family_yreg$family, "Ordered Categorical")))|
      is_multinom_yreg | is_polr_yreg) {
    if (!is.null(yref_index)) {
      if (dim(tot0_pred)[2] == 1) {
        EYtot0 <- weighted.mean(cbind(1 - tot0_pred, tot0_pred)[, yref_index], na.rm = TRUE, w = weightsEY_tot)
        EYtot1 <- weighted.mean(cbind(1 - tot1_pred, tot1_pred)[, yref_index], na.rm = TRUE, w = weightsEY_tot)
        EYdir0 <- weighted.mean(cbind(1 - dir0_pred, dir0_pred)[, yref_index], na.rm = TRUE, w = weightsEY_dir)
        EYdir1 <- weighted.mean(cbind(1 - dir1_pred, dir1_pred)[, yref_index], na.rm = TRUE, w = weightsEY_dir)
      } else {
        EYtot0 <- weighted.mean(tot0_pred[, yref_index], na.rm = TRUE, w = weightsEY_tot)
        EYtot1 <- weighted.mean(tot1_pred[, yref_index], na.rm = TRUE, w = weightsEY_tot)
        EYdir0 <- weighted.mean(dir0_pred[, yref_index], na.rm = TRUE, w = weightsEY_dir)
        EYdir1 <- weighted.mean(dir1_pred[, yref_index], na.rm = TRUE, w = weightsEY_dir)
      }
    } else EYtot0 <- EYtot1 <- EYdir0 <- EYdir1 <- 0
  } else {
    # non-categorical Y
    EYtot0 <- weighted.mean(tot0_pred, na.rm = TRUE, w = weightsEY_tot)
    EYtot1 <- weighted.mean(tot1_pred, na.rm = TRUE, w = weightsEY_tot)
    EYdir0 <- weighted.mean(dir0_pred, na.rm = TRUE, w = weightsEY_dir)
    EYdir1 <- weighted.mean(dir1_pred, na.rm = TRUE, w = weightsEY_dir)
  }
  rm(weightsEY_tot, weightsEY_dir, tot0_pred, tot1_pred, dir0_pred, dir1_pred)

  # output causal effects on the difference scale for continuous Y
  if ((is_lm_yreg | is_glm_yreg) &&
      (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
    tot <- EYtot1 - EYtot0
    dir <- EYdir1 - EYdir0
    ind <- tot - dir
    if (full) {
    pm <- ind / tot
    est <- c(tot, dir, ind, pm)
    } else est <- c(tot, dir, ind)
  } else if (((is_lm_yreg | is_glm_yreg) &&
              (family_yreg$family %in% c("binomial", "quasibinomial", "multinom", "poisson", "quasipoisson", "ziplss") |
               startsWith(family_yreg$family, "Negative Binomial") |
               startsWith(family_yreg$family, "Zero inflated Poisson") |
               startsWith(family_yreg$family, "Ordered Categorical"))) |
             is_multinom_yreg | is_polr_yreg | is_survreg_yreg | is_coxph_yreg){
    # output causal effects in ratio scale for non-continuous Y
    logRRtot <- log(EYtot1) - log(EYtot0)
    logRRdir <- log(EYdir1) - log(EYdir0)
    logRRind <- logRRtot - logRRdir
    if (full) {
    pm <- (exp(logRRdir) * (exp(logRRind) - 1)) / (exp(logRRtot) - 1)
    est <- c(logRRtot, logRRdir, logRRind, pm)
} else est <- c(logRRtot, logRRdir, logRRind)
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

