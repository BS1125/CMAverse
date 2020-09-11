estinf <- function() {
  # restrict data types of variables
  allvar <- c(outcome, exposure, mediator, postc, basec)
  for (i in 1:length(allvar))
    if (!(is.numeric(data[, allvar[i]]) | is.factor(data[, allvar[i]]) |
          is.character(data[, allvar[i]]))) stop(paste0("The variable ", allvar[i], " should be numeric, factor or character"))
  # output list
  out <- list()
  out$multimp <- list(multimp = multimp)
  # obtain calls, weights, classes, and families of regs
  for (reg_name in c("yreg", "ereg", "mreg", "wmreg", "postcreg")) {
    reg <- get(reg_name)
    if (!is.null(reg)) {
      if (reg_name %in% c("yreg", "ereg")) {
        assign(paste0("call_", reg_name), getCall(reg))
        assign("reg_mid", switch(("rcreg" %in% class(reg) | "simexreg" %in% class(reg)) + 1, "1" = reg, "2" = reg$NAIVEreg))
        assign(paste0("is_lm_", reg_name), inherits(reg_mid, "lm"))
        assign(paste0("is_glm_", reg_name), inherits(reg_mid, "glm"))
        assign(paste0("is_svyglm_", reg_name), inherits(reg_mid, "svyglm"))
        assign(paste0("is_gam_", reg_name), inherits(reg_mid, "gam"))
        if (get(paste0("is_lm_", reg_name)) | get(paste0("is_glm_", reg_name))) assign(paste0("family_", reg_name), family(reg_mid))
        assign(paste0("is_multinom_", reg_name), inherits(reg_mid, "multinom"))
        assign(paste0("is_svymultinom_", reg_name), inherits(reg_mid, "svymultinom"))
        assign(paste0("is_polr_", reg_name), inherits(reg_mid, "polr"))
        if (reg_name == "yreg") {
          assign(paste0("is_survreg_", reg_name), inherits(reg_mid, "survreg"))
          assign(paste0("is_coxph_", reg_name), inherits(reg_mid, "coxph"))
        }
        assign(paste0("weights_", reg_name), model.frame(get(reg_name))$'(weights)')   
      } else {
        assign(paste0("call_", reg_name), lapply(1:length(reg), function(x) getCall(reg[[x]])))
        assign("reg_mid", lapply(1:length(reg), function(x)
          switch(("rcreg" %in% class(reg[[x]]) | "simexreg" %in% class(reg[[x]])) + 1, "1" = reg[[x]], "2" = reg[[x]]$NAIVEreg)))
        assign(paste0("is_lm_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "lm")))
        assign(paste0("is_glm_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "glm")))
        assign(paste0("is_svyglm_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "svyglm")))
        assign(paste0("is_gam_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "gam")))
        assign(paste0("family_", reg_name), lapply(1:length(reg_mid), function(x)
          if (get(paste0("is_lm_", reg_name))[x] | get(paste0("is_glm_", reg_name))[x]) family(reg_mid[[x]])))
        assign(paste0("is_multinom_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "multinom")))
        assign(paste0("is_svymultinom_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "svymultinom")))
        assign(paste0("is_polr_", reg_name), sapply(1:length(reg_mid), function(x) inherits(reg_mid[[x]], "polr")))
        assign(paste0("weights_", reg_name), lapply(1:length(reg), function(x) model.frame(get(reg_name)[[x]])$'(weights)'))
      }
      rm(reg_mid)
    }
  }
  
  # restrict classes and families of regression objects
  if (!(((is_lm_yreg | is_glm_yreg) && 
         (family_yreg$family %in% 
          c("gaussian", "inverse.gaussian", "quasi", "poisson", "quasipoisson", 
            "Gamma", "binomial", "quasibinomial", "multinom", "ziplss") |
          startsWith(family_yreg$family, "Negative Binomial") |
          startsWith(family_yreg$family, "Zero inflated Poisson") |
          startsWith(family_yreg$family, "Ordered Categorical"))) |
        is_multinom_yreg | is_polr_yreg | is_survreg_yreg | is_coxph_yreg |
        inference == "delta")) stop("Unsupported yreg")
  if (!is.null(ereg) && 
      !(((is_lm_ereg | is_glm_ereg) && 
          (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
           startsWith(family_ereg$family, "Ordered Categorical"))) |
         is_multinom_ereg | is_polr_ereg)) stop("Unsupported ereg")
  if (!is.null(mreg) && inference == "bootstrap") {
    for (p in 1:length(mreg)) {
      if (!(((is_lm_mreg[[p]] | is_glm_mreg[[p]]) && 
              (family_mreg[[p]]$family %in% 
               c("gaussian", "inverse.gaussian", "poisson", "quasipoisson", 
                 "Gamma", "binomial", "multinom") |
               startsWith(family_mreg[[p]]$family, "Negative Binomial") |
               startsWith(family_mreg[[p]]$family, "Ordered Categorical"))) |
             is_multinom_mreg[[p]] | is_polr_mreg[[p]])) stop(paste0("Unsupported mreg[[", p, "]]"))
    }
  }
  if (!is.null(wmreg)) {
    for (p in 1:length(wmreg)) {
      if (!((((is_lm_wmreg[[p]] | is_glm_wmreg[[p]]) && 
              (family_wmreg[[p]]$family %in% 
               c("binomial", "quasibinomial", "multinom") |
               startsWith(family_wmreg[[p]]$family, "Ordered Categorical"))) |
             is_multinom_wmreg[[p]] | is_polr_wmreg[[p]]))) stop(paste0("Unsupported wmreg[[", p, "]]"))
    }
  }
  if (!is.null(postcreg)) {
    for (p in 1:length(postcreg)) {
      if (!(((is_lm_postcreg[[p]] | is_glm_postcreg[[p]]) && 
              (family_postcreg[[p]]$family %in% 
               c("gaussian", "inverse.gaussian", "poisson", "quasipoisson", 
                 "Gamma", "binomial", "multinom") |
               startsWith(family_postcreg[[p]]$family, "Negative Binomial") |
               startsWith(family_postcreg[[p]]$family, "Ordered Categorical"))) |
             is_multinom_postcreg[[p]] | is_polr_postcreg[[p]])) stop(paste0("Unsupported postcreg[[", p, "]]"))
    }
  }
  
  # reference values of the exposure
  if (is.factor(data[, exposure]) | is.character(data[, exposure])) {
    a_lev <- levels(droplevels(as.factor(data[, exposure])))
    if (!a %in% a_lev) {
      a <- a_lev[length(a_lev)]
      warning(paste0("a is not a value of the exposure; ", a, " is used"))
    }
    if (!astar %in% a_lev) {
      astar <- a_lev[1]
      warning(paste0("astar is not a value of the exposure; ", astar, " is used"))
    }
  }
  out$ref$a <- a
  out$ref$astar <- astar
  
  # yref: the reference level for a categorical outcome
  if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                       startsWith(family_yreg$family, "Ordered Categorical"))) |
      is_multinom_yreg | is_polr_yreg) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    # if yref is not provided or yref provided is not a value of the outcome, use the last level of the outcome
    if (is.null(yref)) {
      yref <- y_lev[length(y_lev)]
      warning(paste0("yref is not specified; ", yref, " is used"))
    }
    if (!yref %in% y_lev) {
      yref <- y_lev[length(y_lev)]
      warning(paste0("yref is not a value of the outcome; ", yref, " is used"))
    }
    out$ref$yref <- yref
  }
  
  # get the level of the case and the level of the control
  if (casecontrol) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    if (length(y_lev) != 2) stop("outcome with more than 2 levels")
    y_control <- y_lev[1]
    y_case <- y_lev[2]
  }

  if (model == "rb") {
    ###################################################################################################
    ############################################Regression-based Approach##############################
    ###################################################################################################
    if (is.null(yreg)) stop("yreg is required")
    # a regression is required for each mediator
    if (is.null(mreg)) stop("mreg is required for model = 'rb'")
    if (!is.list(mreg)) stop("mreg should be a list")
    if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
    for (p in 1:length(mediator)) if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
    out$ref$mval <- mval
    
    # closed-form parameter function estimation
    if (estimation == "paramfunc") {
      # create a list of covariate values to calculate conditional causal effects
      if (length(basec) != 0) {
        if (!is.null(basecval)) {
          if (!is.list(basecval)) stop("basecval should be a list")
          if (length(basecval) != length(basec)) stop("length(basecval) != length(basec)")
        }
        if (is.null(basecval)) basecval <- rep(list(NULL), length(basec))
        # if NULL, set basecval[[i]] to be the mean value of basec[i]
        for (i in 1:length(basec)) {
          if (is.factor(data[, basec[i]]) | is.character(data[, basec[i]])) {
            c_lev <- levels(droplevels(as.factor(data[, basec[i]])))
            if (is.null(basecval[[i]])) {
              c_data <- data[, basec[i], drop = FALSE]
              c_data[, basec[i]] <- factor(c_data[, basec[i]], levels = c_lev)
              # set basecval[[i]] to be the mean values of dummy variables
              basecval[[i]] <- unname(colMeans(as.matrix(model.matrix(as.formula(paste0("~", basec[i])),
                                                                      data = model.frame(~., data = c_data, 
                                                                                         na.action = na.pass))[, -1]), 
                                               na.rm = TRUE))
              rm(c_data)
              # dummy values of basecval[[i]]
            } else basecval[[i]] <- as.numeric(c_lev == basecval[[i]])[-1]
            rm(c_lev)
          } else if (is.numeric(data[, basec[i]])) {
            if (is.null(basecval[[i]])) {
              # set basecval[[i]] to be the mean value of basec[i]
              basecval[[i]] <- mean(data[, basec[i]], na.rm = TRUE)
            } else basecval[[i]] <- basecval[[i]]
          } 
        }
        out$ref$basecval <- basecval
      }
    }
    
    # estimation and inference
    environment(est.rb) <- environment()
    if (!multimp) {
      # point estimates of causal effects
      est <- est.rb(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      
      if (inference == "bootstrap") {
        # bootstrap results
        boots <- boot(data = data, statistic = est.rb, R = nboot, outReg = FALSE, full = full)
        # bootstrap percentile CIs
        effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm = TRUE))
        effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm = TRUE))
        # bootstrap p-values
        effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      } else if (inference == "delta") {
        if (estimation != "paramfunc") stop("When inference = 'delta', estimation can only be 'paramfunc'")
        yreg <- est$reg.output$yreg
        mreg <- est$reg.output$mreg[[1]]
        # standard errors by the delta method
        environment(inf.delta) <- environment()
        effect.se <- inf.delta(data = data, yreg = yreg, mreg = mreg)
        # critical value
        z0 <- qnorm(1 - (1 - 0.95) / 2)
        z <- effect.pe/effect.se
        # delta method CIs
        effect.ci.low <- effect.pe - z0 * effect.se
        effect.ci.high <- effect.pe + z0 * effect.se
        # delta method p-values
        effect.pval <- 2 * (1 - pnorm(abs(z)))
      } else stop("Select inference from 'delta', 'bootstrap'")
      
    } else {
      
      # arguments for mice()
      if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
      args_mice$data <- data
      out$multimp$args_mice <- args_mice
      
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.rb(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      if (inference == "bootstrap") {
        boot.step <- function(data = NULL, indices = NULL) {
          data_boot <- data[indices, ]
          args_mice$data <- data_boot
          data_imp <- complete(do.call(mice, args_mice), action = "all")
          curVal <- get("counter", envir = env)
          assign("counter", curVal + 1, envir = env)
          setTxtProgressBar(get("progbar", envir = env), curVal + 1)
          return(colMeans(do.call(rbind, lapply(1:m, function(x)
            est.rb(data = data_imp[[x]], outReg = FALSE, full = full)))))
        }
        environment(boot.step) <- environment()
        # bootstrap results
        boots <- boot(data = data, statistic = boot.step, R = nboot)
        # bootstrap percentile CIs
        effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm=TRUE))
        effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm=TRUE))
        # bootstrap p-values
        effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      } else if (inference == "delta") {
        if (estimation != "paramfunc") stop("When inference = 'delta', estimation can only be 'paramfunc'")
        environment(inf.delta) <- environment()
        # standard errors by the delta method
        se_imp <- do.call(rbind, lapply(1:m, function(x)
          inf.delta(data = data_imp[[x]], yreg = est_imp[[x]]$reg.output$yreg,
                    mreg = est_imp[[x]]$reg.output$mreg[[1]])))
        # pool the results by Rubin's rule
        var_within <- colMeans(se_imp ^ 2)
        var_between <- colSums((est_imp_df - t(replicate(m, effect.pe)))^2)/(m - 1)
        effect.se <- sqrt(var_within + var_between * (m + 1) / m)
        z0 <- qnorm(1 - (1 - 0.95) / 2)
        z <- effect.pe/effect.se
        effect.ci.low <- effect.pe - z0 * effect.se
        effect.ci.high <- effect.pe + z0 * effect.se
        effect.pval <- 2 * (1 - pnorm(abs(z)))
      } else stop("Select inference from 'delta', 'bootstrap'")
      
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi", "gaulss", "gevlss") |
         startsWith(family_yreg$family, "Tweedie") |
         startsWith(family_yreg$family, "Beta regression") |
         startsWith(family_yreg$family, "Scaled t"))) {
      # standard errors by bootstrapping
      if (inference == "bootstrap") effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                 "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                 "pm", "int", "pe")
      if (!full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
    } else {
      # transform standard errors of effects in log scale
      if (inference == "bootstrap") effect.se <- sapply(1:n_effect, function(x)
        ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x]))) #
      if (inference == "delta") effect.se[1:6] <- sapply(1:6, function(x)
        deltamethod(as.formula("~exp(x1)"), effect.pe[x], effect.se[x]^2))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                 "ERcde", "ERintref", "ERintmed", "ERpnie",
                                 "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                 "pm", "int", "pe")
      if (!full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "gformula") {
    
    ###################################################################################################
    #############################################G-formula Approach####################################
    ###################################################################################################
    
    if (is.null(yreg)) stop("yreg is required")
    # a regression is required for each mediator
    if (is.null(mreg)) stop("mreg is required for model = 'gformula'")
    if (!is.list(mreg)) stop("mreg should be a list")
    if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
    for (p in 1:length(mediator)) if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
    # a regression is required for each post-exposure confounder
    if (length(postc) != 0) {
      if (is.null(postcreg)) stop("postcreg is required for model = 'gformula' when length(postc) != 0")
      if (!is.list(postcreg)) stop("postcreg should be a list")
      if (length(postcreg) != length(postc)) stop("length(postcreg) != length(postc)")
      for (p in 1:length(postc)) if (is.null(postcreg[[p]])) stop(paste0("Unspecified postcreg[[", p, "]]"))
    } else if (!is.null(postcreg)) warning("postcreg is ignored when length(postc) = 0")
    
    out$ref$mval <- mval
    
    environment(est.gformula) <- environment()
    if (!multimp) {
      
      # point estimates of causal effects
      est <- est.gformula(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      
      # bootstrap results
      boots <- boot(data = data, statistic = est.gformula, R = nboot, outReg = FALSE, full = full)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm = TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm = TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    } else {
      
      # arguments for mice()
      if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
      args_mice$data <- data
      out$multimp$args_mice <- args_mice
      
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.gformula(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.gformula(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm=TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm=TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi", "gaulss", "gevlss") |
         startsWith(family_yreg$family, "Tweedie") |
         startsWith(family_yreg$family, "Beta regression") |
         startsWith(family_yreg$family, "Scaled t"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (length(postc) == 0 && full) effect_name <-
          c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
            "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
            "pm", "int", "pe")
      if (length(postc) == 0 && !full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
      if (length(postc) != 0 && full) effect_name <-
          c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te", 
            "rintref", "rintmed", "cde(prop)", "rintref(prop)", "rintmed(prop)", "rpnie(prop)",
            "rpm", "rint", "rpe")
      if (length(postc) != 0 && !full) effect_name <- c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te")
    } else {
      # transform standard errors of effects in log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (length(postc) == 0 && full) effect_name <-
        c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
          "ERcde", "ERintref", "ERintmed", "ERpnie",
          "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
          "pm", "int", "pe")
      if (length(postc) == 0 && !full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
      if (length(postc) != 0 && full) effect_name <-
        c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte", 
          "ERcde", "rERintref", "rERintmed", "rERpnie",
          "ERcde(prop)", "rERintref(prop)", "rERintmed(prop)", "rERpnie(prop)",
          "rpm", "rint", "rpe")
      if (length(postc) != 0 && !full) effect_name <- c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte")
      
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "wb") {
    
    ###################################################################################################
    ############################################Weighting-based Approach##############################
    ###################################################################################################
    
    if (is.null(yreg)) stop("yreg is required")
    if (length(basec) != 0 && is.null(ereg)) stop("ereg is required for model = 'wb' when length(basec) != 0")
    if (length(basec) != 0 && (!((is_glm_ereg && (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
                                                  startsWith(family_ereg$family, "Ordered Categorical"))) |
                                 is_multinom_ereg | is_polr_ereg))) stop(
                                   "model = 'wb' only supports categorical exposure when length(basec) != 0")
    if (is_survreg_yreg | is_coxph_yreg) stop("model = 'wb' doesn't support survival outcomes")
    out$ref$mval <- mval
    
    environment(est.wb) <- environment()
    if (!multimp) {
      
      # point estimates of causal effects
      est <- est.wb(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      
      # bootstrap results
      boots <- boot(data = data, statistic = est.wb, R = nboot, outReg = FALSE, full = full)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm = TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm = TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    } else {
      
      # arguments for mice()
      if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
      args_mice$data <- data
      out$multimp$args_mice <- args_mice
      
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.wb(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.wb(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm=TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm=TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi", "gaulss", "gevlss") |
         startsWith(family_yreg$family, "Tweedie") |
         startsWith(family_yreg$family, "Beta regression") |
         startsWith(family_yreg$family, "Scaled t"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                 "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                 "pm", "int", "pe")
      if (!full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
    } else {
      # transform standard errors of effects in log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                 "ERcde", "ERintref", "ERintmed", "ERpnie",
                                 "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                 "pm", "int", "pe")
      if (!full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "iorw") {
    
    ###################################################################################################
    ###################################Inverse Odds Ratio Weighting Approach###########################
    ###################################################################################################
    
    if (is.null(yreg)) stop("yreg is required")
    if (is.null(ereg)) stop("ereg is required for model = 'iorw'")
    if (!((is_glm_ereg && (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
                           startsWith(family_ereg$family, "Ordered Categorical"))) |
          is_multinom_ereg | is_polr_ereg)) stop("model = 'iorw' only supports categorical exposure")
    
    environment(est.iorw) <- environment()
    if (!multimp) {
      
      # point estimates of causal effects
      est <- est.iorw(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      
      # bootstrap results
      boots <- boot(data = data, statistic = est.iorw, R = nboot, outReg = FALSE,
                    full = full)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm = TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm = TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    } else {
      
      # arguments for mice()
      if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
      args_mice$data <- data
      out$multimp$args_mice <- args_mice
      
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.iorw(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.iorw(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm=TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm=TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi", "gaulss", "gevlss") |
         startsWith(family_yreg$family, "Tweedie") |
         startsWith(family_yreg$family, "Beta regression") |
         startsWith(family_yreg$family, "Scaled t"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) effect_name <- c("te", "pnde", "tnie", "pm")
      if (!full) effect_name <- c("te", "pnde", "tnie")
    } else {
      # transform standard errors of effects in log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 3, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:3] <- exp(effect.pe[1:3])
      effect.ci.low[1:3] <- exp(effect.ci.low[1:3])
      effect.ci.high[1:3] <- exp(effect.ci.high[1:3])
      # effect names
      if (full) effect_name <- c("Rte", "Rpnde", "Rtnie", "pm")
      if (!full) effect_name <- c("Rte", "Rpnde", "Rtnie")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "msm") {
    
    ###################################################################################################
    #########################################Marginal Structural Model#################################
    ###################################################################################################
    
    if (is.null(yreg)) stop("yreg is required")
    if (length(basec) != 0 && is.null(ereg)) stop("ereg is required for model = 'msm' when length(basec) != 0")
    if (length(basec) != 0 && (!((is_glm_ereg && (family_ereg$family %in% c("binomial", "quasibinomial", "multinom") |
                                                  startsWith(family_ereg$family, "Ordered Categorical"))) |
                                 is_multinom_ereg | is_polr_ereg))) stop(
                                   "model = 'msm' only supports categorical exposure when length(basec) != 0")
    # a regression is required for each mediator
    if (is.null(mreg)) stop("mreg is required for model = 'msm'")
    if (!is.list(mreg)) stop("mreg should be a list")
    if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
    for (p in 1:length(mediator)) {
      if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
      if (!((is_glm_mreg[p] && (family_mreg[[p]]$family %in% c("binomial", "multinom") |
                                 startsWith(family_mreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_mreg[p] | is_polr_mreg[p])) stop(
              "model = 'msm' only supports categorical mediators")
    }
    # a regression for calculating weights is required for each mediator
    if (is.null(wmreg)) stop("wmreg is required for model = 'msm'")
    if (!is.list(wmreg)) stop("wmreg should be a list")
    if (length(wmreg) != length(mediator)) stop("length(wmreg) != length(mediator)")
    for (p in 1:length(mediator)) {
      if (is.null(wmreg[[p]])) stop(paste0("Unspecified wmreg[[", p, "]]"))
      if (!((is_glm_wmreg[p] && (family_wmreg[[p]]$family %in% c("binomial", "quasibinomial", "multinom") |
                                 startsWith(family_wmreg[[p]]$family, "Ordered Categorical"))) |
            is_multinom_wmreg[p] | is_polr_wmreg[p])) stop(
              "model = 'msm' only supports categorical mediators")
    }
    
    environment(est.msm) <- environment()
    if (!multimp) {
      
      # point estimates of causal effects
      est <- est.msm(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      
      # bootstrap results
      boots <- boot(data = data, statistic = est.msm, R = nboot, outReg = FALSE, full = full)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm = TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm = TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    } else {
      
      # arguments for mice()
      if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
      args_mice$data <- data
      out$multimp$args_mice <- args_mice
      
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.msm(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.msm(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm=TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm=TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    }
    
    if ((is_lm_yreg | is_glm_yreg) &&
        (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi", "gaulss", "gevlss") |
         startsWith(family_yreg$family, "Tweedie") |
         startsWith(family_yreg$family, "Beta regression") |
         startsWith(family_yreg$family, "Scaled t"))) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (length(postc) == 0 && full) effect_name <-
          c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
            "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
            "pm", "int", "pe")
      if (length(postc) == 0 && !full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
      if (length(postc) != 0 && full) effect_name <-
          c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te", 
            "rintref", "rintmed", "cde(prop)", "rintref(prop)", "rintmed(prop)", "rpnie(prop)",
            "rpm", "rint", "rpe")
      if (length(postc) != 0 && !full) effect_name <- c("cde", "rpnde", "rtnde", "rpnie", "rtnie", "te")
    } else {
      # transform standard errors of effects in log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (length(postc) == 0 && full) effect_name <-
        c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
          "ERcde", "ERintref", "ERintmed", "ERpnie",
          "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
          "pm", "int", "pe")
      if (length(postc) == 0 && !full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
      if (length(postc) != 0 && full) effect_name <-
        c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte", 
          "ERcde", "rERintref", "rERintmed", "rERpnie",
          "ERcde(prop)", "rERintref(prop)", "rERintmed(prop)", "rERpnie(prop)",
          "rpm", "rint", "rpe")
      if (length(postc) != 0 && !full) effect_name <- c("Rcde", "rRpnde", "rRtnde", "rRpnie", "rRtnie", "Rte")
      
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  } else if (model == "ne") {
    
    ###################################################################################################
    #########################################Natural Effect Model######################################
    ###################################################################################################
    
    if (is.null(yreg)) stop("yreg is required")
    if (!identical(class(yreg), c("glm", "lm"))) stop("model = 'ne' only supports yreg fitted via glm()")
    if (is.character(data[, exposure])) data[, exposure] <- as.factor(data[, exposure])
    
    out$ref$mval <- mval
    
    environment(est.ne) <- environment()
    if (!multimp) {
      
      # point estimates of causal effects
      est <- est.ne(data = data, indices = NULL, outReg = TRUE, full = full)
      effect.pe <- est$est
      n_effect <- length(effect.pe)
      out$reg.output <- est$reg.output
      
      # bootstrap results
      boots <- boot(data = data, statistic = est.ne, R = nboot, outReg = FALSE, full = full)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm = TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm = TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    } else {
      
      # arguments for mice()
      if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
      args_mice$data <- data
      out$multimp$args_mice <- args_mice
      
      # imputed data sets
      data_imp <- complete(do.call(mice, args_mice), action = "all")
      m <- length(data_imp)
      # estimate causal effects for each imputed data set
      est_imp <- lapply(1:m, function(x)
        est.ne(data = data_imp[[x]], indices = NULL, outReg = TRUE, full = full))
      est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
      effect.pe <- colMeans(est_imp_df)
      n_effect <- length(effect.pe)
      out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
      
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est.ne(data = data_imp[[x]], outReg = FALSE, full = full)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap percentile CIs
      effect.ci.low <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = (1 - 0.95)/2, na.rm=TRUE))
      effect.ci.high <- sapply(1:n_effect, function(x) quantile(x = boots$t[, x], probs = 1 - (1 - 0.95)/2, na.rm=TRUE))
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
      
    }
    
    if (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi")) {
      # standard errors by bootstrapping
      effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x]))
      # effect names
      if (full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                 "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                 "pm", "int", "pe")
      if (!full) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te")
    } else {
      # transform standard errors of effects in log scale
      effect.se <- sapply(1:n_effect, function(x) ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x])))
      # transform effects in log ratio scale into effects in ratio scale
      effect.pe[1:6] <- exp(effect.pe[1:6])
      effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
      effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
      # effect names
      if (full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                 "ERcde", "ERintref", "ERintmed", "ERpnie",
                                 "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                 "pm", "int", "pe")
      if (!full) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte")
    }
    
    names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
      names(effect.pval) <- effect_name
    out$effect.pe <- effect.pe
    out$effect.se <- effect.se
    out$effect.ci.low <- effect.ci.low
    out$effect.ci.high <- effect.ci.high
    out$effect.pval <- effect.pval
    
  }
  
  return(out)
}


boot.pval <- function(boots, pe){
  if (pe == 0) out <- 1
  if (pe != 0) out <- 2 * min(sum(boots > 0), sum(boots < 0)) / length(boots)
  return(out)
}


rqpois = function(n, lambda, phi) {
  r = stats::rnbinom(n, mu = lambda, size = lambda/(phi-1))
  return(r)
}

