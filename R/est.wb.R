est.wb <- function(data = NULL, indices = NULL, outReg = FALSE, full = TRUE) {
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
      rm(weights_ereg, a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)

    # weights for yreg
    if (!is.null(weights_yreg)) weights_yreg <- weights_yreg[indices] * w4casecon
    if (is.null(weights_yreg)) weights_yreg <- w4casecon
    # update yreg
    call_yreg$weights <- weights_yreg
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)
    rm(prob1, w4casecon, weights_ereg, weights_yreg)

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
      wa <- as.vector(wanom / wadenom)
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)

    # update yreg
    call_yreg$weights <- weights_yreg[indices]
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
      wa <- as.vector(wanom / wadenom)
      rm(a_lev, wa_data, wanom, wadenom_prob, category, wadenom)
    } else wa <- rep(1, n)

    # update yreg
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    yreg <- eval.parent(call_yreg)
  }

  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yreg <- yreg
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

  # subjects with A = a and subjects with A = astar
  subj_astar <- which(data[, exposure] == astar)
  subj_a <- which(data[, exposure] == a)

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

  # simulate mstar for cde
  mstar_sim <- do.call(cbind, lapply(1:length(mediator), function(x)
    if (is.factor(data[, mediator[x]])) {
      data.frame(factor(rep(mval[[x]], n), levels = levels(data[, mediator[x]])))
    } else data.frame(rep(mval[[x]], n))))

  # design matrices for outcome simulation
  ydesign0m <- data.frame(astar_sim, mstar_sim, basec_sim)
  ydesign1m <- data.frame(a_sim, mstar_sim, basec_sim)
  ydesign01 <- data.frame(astar_sim, data[, mediator], basec_sim)[subj_a, ]
  ydesign10 <- data.frame(a_sim, data[, mediator], basec_sim)[subj_astar, ]
  rm(a_sim, astar_sim, mstar_sim, basec_sim)
  colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign01) <-
    colnames(ydesign10) <- c(exposure, mediator, basec)

  # predict Y
  type <- ifelse(is_coxph_yreg, "risk", ifelse(is_multinom_yreg | is_polr_yreg, "probs", "response"))
  EY0m_pred <- as.matrix(predict(yreg, newdata =  ydesign0m, type = type))
  EY1m_pred <- as.matrix(predict(yreg, newdata =  ydesign1m, type = type))
  EY01_pred <- as.matrix(predict(yreg, newdata =  ydesign01, type = type))
  EY10_pred <- as.matrix(predict(yreg, newdata =  ydesign10, type = type))
  rm(type, ydesign0m, ydesign1m, ydesign01, ydesign10)

  # weights for calculating counterfactuals
  weightsEY_cde <- as.vector(model.frame(yreg)$'(weights)')
  if (is.null(weightsEY_cde)) weightsEY_cde <- rep(1, n)
  weightsEY <- as.vector(model.frame(yreg)$'(weights)')
  if (!is.null(weightsEY)) weightsEY <- weightsEY * wa
  if (is.null(weightsEY)) weightsEY <- wa

  y_obs <- data[, outcome]
  # categorical Y
  if ((is_glm_yreg && ((family_yreg$family %in% c("binomial", "quasibinomial", "multinom")) |
                       startsWith(family_yreg$family, "Ordered Categorical")))|
      is_multinom_yreg | is_polr_yreg) {
    if (!is.null(yref_index)) {
      if (dim(EY0m_pred)[2] == 1) {
        EY0m <- weighted.mean(cbind(1 - EY0m_pred, EY0m_pred)[, yref_index], na.rm = TRUE, w = weightsEY_cde)
        EY1m <- weighted.mean(cbind(1 - EY1m_pred, EY1m_pred)[, yref_index], na.rm = TRUE, w = weightsEY_cde)
        EY00 <- weighted.mean(as.numeric(y_obs[subj_astar] == yref), na.rm = TRUE, w = weightsEY[subj_astar])
        EY01 <- weighted.mean(cbind(1 - EY01_pred, EY01_pred)[, yref_index], na.rm = TRUE, w = weightsEY[subj_a])
        EY10 <- weighted.mean(cbind(1 - EY10_pred, EY10_pred)[, yref_index], na.rm = TRUE, w = weightsEY[subj_astar])
        EY11 <- weighted.mean(as.numeric(y_obs[subj_a] == yref), na.rm = TRUE, w = weightsEY[subj_a])
      } else {
        EY0m <- weighted.mean(EY0m_pred[, yref_index], na.rm = TRUE, w = weightsEY_cde)
        EY1m <- weighted.mean(EY1m_pred[, yref_index], na.rm = TRUE, w = weightsEY_cde)
        EY00 <- weighted.mean(as.numeric(y_obs[subj_astar] == yref), na.rm = TRUE, w = weightsEY[subj_astar])
        EY01 <- weighted.mean(EY01_pred[, yref_index], na.rm = TRUE, w = weightsEY[subj_a])
        EY10 <- weighted.mean(EY10_pred[, yref_index], na.rm = TRUE, w = weightsEY[subj_astar])
        EY11 <- weighted.mean(as.numeric(y_obs[subj_a] == yref), na.rm = TRUE, w = weightsEY[subj_a])
      }
    } else EY0m <- EY1m <- EY00 <- EY01 <- EY10 <- EY11 <- 0
  } else {
    # non-categorical Y
    EY0m <- weighted.mean(EY0m_pred, na.rm = TRUE, w = weightsEY_cde)
    EY1m <- weighted.mean(EY1m_pred, na.rm = TRUE, w = weightsEY_cde)
    EY00 <- weighted.mean(y_obs[subj_astar], na.rm = TRUE, w = weightsEY[subj_astar])
    EY01 <- weighted.mean(EY01_pred, na.rm = TRUE, w = weightsEY[subj_a])
    EY10 <- weighted.mean(EY10_pred, na.rm = TRUE, w = weightsEY[subj_astar])
    EY11 <- weighted.mean(y_obs[subj_a], na.rm = TRUE, w = weightsEY[subj_a])
  }
  rm(weightsEY_cde, weightsEY, EY0m_pred, EY1m_pred, EY01_pred, EY10_pred, y_obs)

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
  } else if ((is_glm_yreg &&
              (family_yreg$family %in% c("binomial", "quasibinomial", "multinom", "poisson", "quasipoisson", "ziplss") |
               startsWith(family_yreg$family, "Negative Binomial") |
               startsWith(family_yreg$family, "Zero inflated Poisson") |
               startsWith(family_yreg$family, "Ordered Categorical"))) |
             is_multinom_yreg | is_polr_yreg | is_survreg_yreg | is_coxph_yreg){
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

