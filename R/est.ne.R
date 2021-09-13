est.ne <- function(data = NULL, indices = NULL, outReg = FALSE, full = NULL) {
  if (is.null(indices)) indices <- 1:n
  # resample data
  data <- data[indices, ]
  
  if (casecontrol && !is.null(yprevalence)) {
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
    rm(prob1, w4casecon)
  } else {
    call_yreg$weights <- weights_yreg[indices]
    call_yreg$data <- data
    if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_yreg$variance <- TRUE
    yreg <- eval.parent(call_yreg)
  }
  
  # the expanded dataset
  expdata <- eval(bquote(neImpute(yreg, nMed = length(mediator), nRep = .(nRep))))
  # natural effect model formula
  ne_a <- paste0(exposure, "0")
  ne_m <- paste0(exposure, "1")
  formula_yreg <- paste0(formula(yreg)[2], "~", formula(yreg)[3])
  ne_formula <- gsub(exposure, ne_a, formula_yreg)
  for (p in 1:length(mediator)) ne_formula <- gsub(mediator[p], ne_m, ne_formula)
  
  # natural effect model
  call_nereg <- call_yreg
  weights_yreg_new <- as.vector(model.frame(yreg)$'(weights)')
  if (!is.null(weights_yreg_new)) call_nereg$weights <- rep(weights_yreg_new, rep(nRep, length(weights_yreg_new)))
  call_nereg$formula <- as.formula(ne_formula)
  call_nereg$data <- expdata
  if (outReg && (inherits(yreg, "rcreg") | inherits(yreg, "simexreg"))) call_nereg$variance <- TRUE
  nereg <- eval.parent(call_nereg)
  rm(expdata, formula_yreg, ne_formula)
  
  # output list
  out <- list()
  if (outReg) {
    out$reg.output$yreg <- yreg
    out$reg.output$nereg <- nereg
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
  
  # simulate mstar for cde
  mstar_sim <- do.call(cbind, lapply(1:length(mediator), function(x)
    if (is.factor(data[, mediator[x]])) {
      data.frame(factor(rep(mval[[x]], n), levels = levels(data[, mediator[x]])))
    } else data.frame(rep(mval[[x]], n))))
  
  # design matrices for outcome simulation
  ydesign0m <- data.frame(astar_sim, mstar_sim, basec_sim)
  ydesign1m <- data.frame(a_sim, mstar_sim, basec_sim)
  ydesign00 <- data.frame(astar_sim, astar_sim, basec_sim)
  ydesign01 <- data.frame(astar_sim, a_sim, basec_sim)
  ydesign10 <- data.frame(a_sim, astar_sim, basec_sim)
  ydesign11 <- data.frame(a_sim, a_sim, basec_sim)
  rm(a_sim, astar_sim, mstar_sim, basec_sim)
  colnames(ydesign0m) <- colnames(ydesign1m) <- c(exposure, mediator, basec)
  colnames(ydesign00) <- colnames(ydesign01) <- colnames(ydesign10) <- colnames(ydesign11) <- c(ne_a, ne_m, basec)
  
  # predict Y
  EY0m_pred <- as.matrix(predict(yreg, newdata = ydesign0m, type = "response"))
  EY1m_pred <- as.matrix(predict(yreg, newdata = ydesign1m, type = "response"))
  EY00_pred <- as.matrix(predict(nereg, newdata = ydesign00, type = "response"))
  EY01_pred <- as.matrix(predict(nereg, newdata = ydesign01, type = "response"))
  EY10_pred <- as.matrix(predict(nereg, newdata = ydesign10, type = "response"))
  EY11_pred <- as.matrix(predict(nereg, newdata = ydesign11, type = "response"))
  rm(ydesign0m, ydesign1m, ydesign00, ydesign01, ydesign10, ydesign11)
  
  # weights of yreg
  weightsEY <- weights_yreg_new
  if (is.null(weightsEY)) weightsEY <- rep(1, n)
  
  # categorical Y
  if (family_yreg$family %in% c("binomial", "quasibinomial")) {
    if (!is.null(yval_index)) {
      EY0m <- weighted.mean(cbind(1 - EY0m_pred, EY0m_pred)[, yval_index], na.rm = TRUE, w = weightsEY)
      EY1m <- weighted.mean(cbind(1 - EY1m_pred, EY1m_pred)[, yval_index], na.rm = TRUE, w = weightsEY)
      EY00 <- weighted.mean(cbind(1 - EY00_pred, EY00_pred)[, yval_index], na.rm = TRUE, w = weightsEY)
      EY01 <- weighted.mean(cbind(1 - EY01_pred, EY01_pred)[, yval_index], na.rm = TRUE, w = weightsEY)
      EY10 <- weighted.mean(cbind(1 - EY10_pred, EY10_pred)[, yval_index], na.rm = TRUE, w = weightsEY)
      EY11 <- weighted.mean(cbind(1 - EY11_pred, EY11_pred)[, yval_index], na.rm = TRUE, w = weightsEY)
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
  if (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi")) {
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
  } else if (family_yreg$family %in% c("binomial", "quasibinomial", "poisson", "quasipoisson")) {
    # output causal effects on the ratio scale for non-continuous Y
    
    ## output effects on the odds ratio scale for logistic regressions
    if (is_glm_yreg && family_yreg$family %in% c("binomial", "quasibinomial") &&
        yreg$family$link == "logit") {
      logRRcde <- log(EY1m/(1-EY1m)) - log(EY0m/(1-EY0m))
      logRRpnde <- log(EY10/(1-EY10)) - log(EY00/(1-EY00))
      logRRtnde <- log(EY11/(1-EY11)) - log(EY01/(1-EY01))
      logRRpnie <- log(EY01/(1-EY01)) - log(EY00/(1-EY00))
      logRRtnie <- log(EY11/(1-EY11)) - log(EY10/(1-EY10))
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
      
      ## otherwise on the risk ratio scale
    } else {
      logRRcde <- log(EY1m) - log(EY0m)
      logRRpnde <- log(EY10) - log(EY00)
      logRRtnde <- log(EY11) - log(EY01)
      logRRpnie <- log(EY01) - log(EY00)
      logRRtnie <- log(EY11) - log(EY10)
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
