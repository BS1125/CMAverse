allreg <- function(model = NULL, data = NULL,
                   yreg = NULL, mreg = NULL, ereg = NULL, postcreg = NULL,
                   wmnomreg = NULL, wmdenomreg = NULL,
                   outcome = NULL, event = NULL,
                   exposure = NULL, mediator = NULL, EMint = FALSE,
                   prec = NULL, postc = NULL) {

  if (is.null(data)) stop("Unspecified data")

  if (is.null(model)) {

    stop("Unspecified model")

  } else if (!model %in% c("rb", "wb", "iorw", "msm", "gformula", "ne")) {

    stop("Unsupported model")

  } else if (!model %in% c("msm", "gformula") && !is.null(postc)) {

    stop("The selected model doesn't support post-exposure confounding")

  }

  n <- nrow(data)


  ##################################################################################################
  #######################################Classes and Families of regs###############################
  ##################################################################################################

  # classes of regs; mreg_class, wmnomreg_class, wmdenomreg_class and postcreg_class are lists
  yreg_class <- switch((class(yreg) %in% c("rcreg", "simexreg")) + 1,
                       "1" = class(yreg), "2" = class(yreg$NAIVEreg))

  ereg_class <- switch((class(ereg) %in% c("rcreg", "simexreg")) + 1,
                       "1" = class(ereg), "2" = class(ereg$NAIVEreg))

  mreg_class <- switch((length(mreg) != 0) + 1, "1" = "NULL",
                       "2" = lapply(1:length(mreg), function(x)
                         switch((class(mreg[[x]]) %in% c("rcreg", "simexreg")) + 1,
                                "1" = class(mreg[[x]]), "2" = class(mreg[[x]]$NAIVEreg))))

  wmnomreg_class <- switch((length(wmnomreg) != 0) + 1, "1" = "NULL",
                           "2" = lapply(1:length(wmnomreg), function(x)
                             switch((class(wmnomreg[[x]]) %in% c("rcreg", "simexreg")) + 1,
                                    "1" = class(wmnomreg[[x]]), "2" = class(wmnomreg[[x]]$NAIVEreg))))

  wmdenomreg_class <- switch((length(wmdenomreg) != 0) + 1, "1" = "NULL",
                             "2" = lapply(1:length(wmdenomreg), function(x)
                               switch((class(wmdenomreg[[x]]) %in% c("rcreg", "simexreg")) + 1,
                                      "1" = class(wmdenomreg[[x]]), "2" = class(wmdenomreg[[x]]$NAIVEreg))))

  postcreg_class <- switch((length(postcreg) != 0) + 1, "1" = "NULL",
                           "2" = lapply(1:length(postcreg), function(x)
                             switch((class(postcreg[[x]]) %in% c("rcreg", "simexreg")) + 1,
                                    "1" = class(postcreg[[x]]), "2" = class(postcreg[[x]]$NAIVEreg))))

  # indicators of classes of regs and families of regs
  for (reg in c("yreg", "mreg", "wmnomreg", "wmdenomreg", "ereg", "postcreg")) {

    reg_class <- get(paste0(reg, "_class"))

    if (!is.list(reg_class)) {

      assign(paste0("is_lm_", reg), ifelse(identical(reg_class, "lm"), TRUE, FALSE))

      assign(paste0("is_glm_", reg), ifelse(identical(reg_class, c("glm", "lm")), TRUE, FALSE))

      assign(paste0("is_gam_", reg), ifelse(identical(reg_class, c("gam", "glm", "lm")), TRUE, FALSE))

      if (get(paste0("is_glm_", reg)) | get(paste0("is_gam_", reg))) {

        assign(paste0("family_", reg), switch((class(get(reg)) %in% c("rcreg", "simexreg")) + 1,
                                              "1" = family(get(reg))$family, "2" = family(get(reg)$NAIVEreg)$family))

      } else assign(paste0("family_", reg), "NULL")

      assign(paste0("is_multinom_", reg), ifelse(identical(reg_class, c("multinom", "nnet")), TRUE, FALSE))

      assign(paste0("is_polr_", reg), ifelse(identical(reg_class, "polr"), TRUE, FALSE))

      assign(paste0("is_glmnb_", reg), ifelse(identical(reg_class, c("negbin", "glm", "lm")), TRUE, FALSE))

      assign(paste0("is_survreg_", reg), ifelse(identical(reg_class, "survreg"), TRUE, FALSE))

      assign(paste0("is_coxph_", reg), ifelse(identical(reg_class, "coxph"), TRUE, FALSE))

      assign(paste0("is_chara_", reg), ifelse(identical(reg_class, "character"), TRUE, FALSE))

    } else {

      assign(paste0("is_lm_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], "lm"), TRUE, FALSE)))

      assign(paste0("is_glm_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], c("glm", "lm")), TRUE, FALSE)))

      assign(paste0("is_gam_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], c("gam", "glm", "lm")), TRUE, FALSE)))

      assign(paste0("family_", reg), lapply(1:length(reg_class), function(x)
        if (get(paste0("is_glm_", reg))[x] | get(paste0("is_gam_", reg))[x]) {
          if (class(get(reg)[[i]]) %in% c("rcreg", "simexreg")) {
            family(get(reg)[[i]]$NAIVEreg)$family
          } else family(get(reg)[[i]])$family
        } else "NULL"))

      assign(paste0("is_multinom_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], c("multinom", "nnet")), TRUE, FALSE)))

      assign(paste0("is_polr_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], "polr"), TRUE, FALSE)))

      assign(paste0("is_glmnb_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], c("negbin", "glm", "lm")), TRUE, FALSE)))

      assign(paste0("is_survreg_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], "survreg"), TRUE, FALSE)))

      assign(paste0("is_coxph_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], "coxph"), TRUE, FALSE)))

      assign(paste0("is_chara_", reg), sapply(1:length(reg_class), function(x)
        ifelse(identical(reg_class[[x]], "character"), TRUE, FALSE)))

    }

  }


  ##################################################################################################
  ########################################Exposure Regression#######################################
  ##################################################################################################

  # for wb and msm, the exposure regression is needed for calculating weights if prec is not empty;
  # for iorw, the exposure regression is needed for calculating weights
  if ((model %in% c("wb", "msm") && !is.null(prec)) | model == "iorw") {

    if (is.null(ereg)) stop("Unspecified ereg")

    if (is_chara_ereg) {

      if (is.null(exposure)) stop("Unspecified exposure")

      if (is.null(mediator) && model == "iorw") stop("Unspecified mediator")

      # create exposure regression formula
      switch(model,
             iorw = exposure_formula <- paste0(exposure, "~", paste0(c(mediator, prec), collapse = "+")),
             wb = exposure_formula <- paste0(exposure, "~", paste0(prec, collapse = "+")),
             msm = exposure_formula <- paste0(exposure, "~", paste0(prec, collapse = "+")))

      # run exposure regression; exposure weights are only applicable when the exposure is categorical
      if (!ereg %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop("Unsupported ereg")

      switch(ereg,
             logistic = exposure_regression <- glm(exposure_formula, family = binomial(),
                                                   data = data),
             loglinear = exposure_regression <- glm(exposure_formula, family = binomial("log"), data = data),
             multinomial = exposure_regression <- nnet::multinom(exposure_formula, data = data, trace = FALSE),
             ordinal = exposure_regression <- MASS::polr(exposure_formula, data = data))

    } else if (((is_glm_ereg | is_gam_ereg) &&
                (family_ereg %in% c("multinom", "binomial", "quasibinomial") |
                 startsWith(family_ereg,"Ordered Categorical")))|
               is_multinom_ereg | is_polr_ereg) {

      exposure_regression <- update(ereg, data = data,
                                    formula. = formula(ereg))

    } else stop("Unsupported ereg")

  } else exposure_regression <- NULL


  ##################################################################################################
  ###########################################Mediator Regression####################################
  ##################################################################################################

  # for rb, msm and gformula, a mediator regression is needed for each mediator
  if (model %in% c("rb", "msm", "gformula")) {

    if (is.null(mreg)) stop("Unspecified mreg")

    if (is.null(mediator)|is.null(exposure)) stop("Unspecified mediator or exposure")

    if (length(mediator) != length(mreg)) stop("length(mediator) != length(mreg)")

    mediator_regression <- list()

    for (i in 1:length(mreg)) {

      if (is.null(mreg[[i]])) stop("Unspecified mreg[[i]]")

      if (is_chara_mreg[i]) {

        if (!mreg[[i]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                              "negbin", "multinomial", "ordinal")) stop("Unspecified mreg[[i]]")

        if ((mreg[[i]] %in% c("linear", "poisson", "quasipoisson", "negbin")) &&
            model == "msm") stop("Unsupported mreg[[i]] for msm")

        switch(model,
               rb = mediator_formula <- paste0(mediator[i], "~",
                                               paste(c(exposure, mediator[0:(i-1)], prec), collapse = "+")),
               msm = mediator_formula <- paste0(mediator[i], "~", paste(c(exposure, mediator[0:(i-1)]),
                                                                        collapse = "+")),
               gformula = mediator_formula <- paste0(mediator[i], "~", paste(c(exposure, mediator[0:(i-1)],
                                                                               prec, postc), collapse = "+")))

        switch(mreg[[i]],
               linear = mediator_regression[[i]] <- lm(mediator_formula, data = data),
               logistic = mediator_regression[[i]] <- glm(mediator_formula, family = binomial(),
                                                          data = data),
               loglinear = mediator_regression[[i]] <- glm(mediator_formula,
                                                           family = binomial("log"), data = data),
               poisson = mediator_regression[[i]]  <- glm(mediator_formula, family = poisson(),
                                                          data = data),
               quasipoisson = mediator_regression[[i]] <- glm(mediator_formula,
                                                              family = quasipoisson(), data = data),
               negbin = mediator_regression[[i]] <- MASS::glm.nb(mediator_formula, data = data),
               multinomial = mediator_regression[[i]] <- nnet::multinom(mediator_formula, data = data,
                                                                        trace = FALSE),
               ordinal = mediator_regression[[i]] <- MASS::polr(mediator_formula, data = data))

      } else if ((is_lm_mreg[i] | is_glm_mreg[i] | is_gam_mreg[i] | is_multinom_mreg[i]|
                  is_polr_mreg[i] | is_glmnb_mreg[i])) {

        if ((!((((is_glm_mreg[i] | is_gam_mreg[i]) &&
                 (family_mreg[[i]] %in% c("multinom", "binomial", "quasibinomial") |
                  startsWith(family_mreg[[i]],"Ordered Categorical")))|
                is_multinom_mreg[i] | is_polr_mreg[i]))) && model == "msm") stop("Unsupported mreg[[i]] for msm")

        mediator_regression[[i]] <- update(mreg[[i]], data = data,
                                           formula. = formula(mreg[[i]]))

      } else stop("Unsupported mreg[[i]]")

    }

  } else mediator_regression <- NULL


  ##################################################################################################
  ####################################Mediator Regression For Weights###############################
  ##################################################################################################

  if (model == "msm") {

    if (length(wmnomreg) != length(wmdenomreg)) {

      stop("length(wmnomreg) != length(wmdenomreg)")

    } else if (length(wmdenomreg) != length(mediator)) {

      stop("length(wmdenomreg) != length(mediator)")

    }

    # regressions for P(Mp|A)
    wmnom_regression <- list()

    for (i in 1:length(wmnomreg)) {

      if (is_chara_wmnomreg[i]) {

        if (wmnomreg[[i]] %in% c("linear", "poisson", "quasipoisson", "negbin")) stop("Unsupported wmnomreg[[i]] for msm")

        wmnom_formula <- paste0(mediator[i], "~", exposure)

        switch(wmnomreg[[i]],
               logistic = wmnom_regression[[i]] <- glm(wmnom_formula, family = binomial(),
                                                       data = data),
               loglinear = wmnom_regression[[i]] <- glm(wmnom_formula,
                                                        family = binomial("log"), data = data),
               multinomial = wmnom_regression[[i]] <- nnet::multinom(wmnom_formula, data = data,
                                                                     trace = FALSE),
               ordinal = wmnom_regression[[i]] <- MASS::polr(wmnom_formula, data = data))

      } else if (((is_glm_wmnomreg[i] | is_gam_wmnomreg[i]) &&
                  (family_wmnomreg[[i]] %in% c("multinom", "binomial", "quasibinomial") |
                   startsWith(family_wmnomreg[[i]],"Ordered Categorical")))|
                 is_multinom_wmnomreg[i] | is_polr_wmnomreg[i]) {

        wmnom_regression[[i]] <- update(wmnomreg[[i]], data = data,
                                        formula. = formula(wmnomreg[[i]]))

      } else if (is.null(wmnomreg[[i]])) {

        wmnom_regression[[i]] <- NULL

      } else stop("Unsupported wmnomreg[[i]]")

    }

    # regressions for P(Mp|A,M1,...,Mp-1,C,L)
    wmdenom_regression <- list()

    for (i in 1:length(wmdenomreg)) {

      if (is_chara_wmdenomreg[i]) {

        if (wmdenomreg[[i]] %in% c("linear", "poisson", "quasipoisson", "negbin")) stop("Unsupported wmdenomreg[[i]] for msm")

        wmdenom_formula <- paste0(mediator[i], "~",
                                  paste(c(exposure, mediator[0:(i-1)], prec, postc),
                                        collapse = "+"))

        switch(wmdenomreg[[i]],
               logistic = wmdenom_regression[[i]] <- glm(wmdenom_formula, family = binomial(),
                                                         data = data),
               loglinear = wmdenom_regression[[i]] <- glm(wmdenom_formula,
                                                          family = binomial("log"), data = data),
               multinomial = wmdenom_regression[[i]] <- nnet::multinom(wmdenom_formula, data = data,
                                                                       trace = FALSE),
               ordinal = wmdenom_regression[[i]] <- MASS::polr(wmdenom_formula, data = data))

      } else if (((is_glm_wmdenomreg[i] | is_gam_wmdenomreg[i]) &&
                  (family_wmdenomreg[[i]] %in% c("multinom", "binomial", "quasibinomial") |
                   startsWith(family_wmdenomreg[[i]],"Ordered Categorical")))|
                 is_multinom_wmdenomreg[i] | is_polr_wmdenomreg[i]) {

        wmdenom_regression[[i]] <- update(wmdenomreg[[i]], data = data,
                                          formula. = formula(wmdenomreg[[i]]))

      } else if (is.null(wmdenomreg[[i]])) {

        wmdenom_regression[[i]] <- NULL

      } else stop("Unsupported wmdenomreg[[i]]")

    }

  } else wmnom_regression <- wmdenom_regression <- NULL


  ##################################################################################################
  ####################################Post-exposure Confounder Regression###########################
  ##################################################################################################

  # regressions for post-exposure confounders for gformula
  if (model == "gformula" && !is.null(postcreg)) {

    if (is.null(postc)) stop("Unspecified postc")

    if (length(postc) != length(postcreg)) stop("length(postc) != length(postcreg)")

    postc_regression <- list()

    for (i in 1:length(postcreg)) {

      if (is.null(postcreg[[i]])) stop("Unspecified postcreg[[i]]")

      if (is_chara_postcreg[i]) {

        if (!postcreg[[i]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                                  "negbin", "multinomial", "ordinal")) stop("Unspecified postcreg[[i]]")

        postc_formula <- paste0(postc[i], "~",
                                paste(c(exposure, prec, postc[0:(i-1)]),
                                      collapse = "+"))

        switch(postcreg[[i]],
               linear = postc_regression[[i]] <- lm(postc_formula, data = data),
               logistic = postc_regression[[i]] <- glm(postc_formula, family = binomial(),
                                                       data = data),
               loglinear = postc_regression[[i]] <- glm(postc_formula,
                                                        family = binomial("log"), data = data),
               poisson = postc_regression[[i]]  <- glm(postc_formula, family = poisson(),
                                                       data = data),
               quasipoisson = postc_regression[[i]] <- glm(postc_formula,
                                                           family = quasipoisson(), data = data),
               negbin = postc_regression[[i]] <- MASS::glm.nb(postc_formula, data = data),
               multinomial = postc_regression[[i]] <- nnet::multinom(postc_formula, data = data,
                                                                     trace = FALSE),
               ordinal = postc_regression[[i]] <- MASS::polr(postc_formula, data = data))

      } else if ((is_lm_postcreg[i] | is_glm_postcreg[i] | is_gam_postcreg[i] | is_multinom_postcreg[i]|
                  is_polr_postcreg[i] | is_glmnb_postcreg[i])) {

        postc_regression[[i]] <- update(postcreg[[i]], data = data,
                                        formula. = formula(postcreg[[i]]))

      } else stop("Unsupported postcreg[[i]]")

    }

  } else postc_regression <- NULL


  ##################################################################################################
  ###########################################Outcome Regression#####################################
  ##################################################################################################

  if (is.null(yreg)) stop("Unspecified yreg")

  if (is_chara_yreg) {

    if (is.null(outcome)|is.null(exposure)) stop("Unspecified outcome or exposure")

    if (is.null(mediator) && model != "iorw") stop("Unspecified mediator")

    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson", "negbin",
                     "multinomial", "ordinal", "coxph", "aft_exp", "aft_weibull")) stop("Unsupported yreg")

    int.terms <- switch((model != "iorw" && EMint) + 1, "1" = NULL, "2" = paste(exposure, mediator, sep = "*"))

    # create outcome regression formula
    switch(model,
           rb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator,
                                                                int.terms, prec), collapse = "+")),
           wb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator,
                                                                int.terms, prec), collapse = "+")),
           ne = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator,
                                                                int.terms, prec), collapse = "+")),
           iorw = exposure_formula <- paste0(outcome, "~", paste(c(exposure, prec), collapse = "+")),
           msm = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms),
                                                               collapse = "+")),
           gformula = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms,
                                                                      prec, postc), collapse = "+")))

    if (yreg %in% c("coxph","aft_exp","aft_weibull")) {

      if (is.null(event)) stop("Unspecified event")

      outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                               strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")

    }

    # run outcome regression
    switch(yreg,
           linear = outcome_regression <- lm(outcome_formula, data = data),
           logistic = outcome_regression <- glm(outcome_formula, family = binomial(),
                                                data = data),
           loglinear = outcome_regression <- glm(outcome_formula,
                                                 family = binomial("log"), data = data),
           poisson = outcome_regression  <- glm(outcome_formula, family = poisson(),
                                                data = data),
           quasipoisson = outcome_regression <- glm(outcome_formula,
                                                    family = quasipoisson(), data = data),
           negbin = outcome_regression <- MASS::glm.nb(outcome_formula, data = data),
           multinomial = outcome_regression <- nnet::multinom(outcome_formula, data = data,
                                                              trace = FALSE),
           ordinal = outcome_regression <- MASS::polr(outcome_formula, data = data),
           coxph = outcome_regression <- survival::coxph(outcome_formula, data = data),
           aft_exp = outcome_regression <- survival::survreg(outcome_formula,
                                                             dist = "exponential", data = data),
           aft_weibull = outcome_regression <- survival::survreg(outcome_formula,
                                                                 dist = "weibull", data = data))

  } else if ((is_lm_yreg | is_glm_yreg | is_gam_yreg | is_multinom_yreg|
              is_polr_yreg | is_glmnb_yreg | is_survreg_yreg | is_coxph_yreg)) {

    outcome_regression <- update(yreg, data = data, formula. = formula(yreg))

  } else stop("Unsupported yreg")


  ##################################################################################################
  ###########################################Return Regression List#################################
  ##################################################################################################

  if (model == "rb") {

    outcome_regression <- update(outcome_regression, formula. = formula(outcome_regression))

    for (i in 1:length(mediator_regression)) mediator_regression[[i]] <-
        update(mediator_regression[[i]], formula. = formula(mediator_regression[[i]]))

    regressions <- list(mediator_regression = mediator_regression, outcome_regression =  outcome_regression)

  } else if (model == "msm") {

    # calculate P(A=ai)/P(A=ai|C=ci)
    if (is.null(exposure_regression)) {

      wa <- NULL

    } else {

      wanom <- left_join(select(data, exposure),
                         count(data, !!as.name(exposure)),
                         by = exposure)[, "n"]/n

      if (ereg %in% c("logistic", "loglinear") |
          (is_glm_ereg && family_ereg %in% c("binomial", "quasibinomial")) |
          ((is_multinom_ereg | ereg == "multinomial") &&
           length(unique(data[, exposure])) == 2)) {

        category <- as.numeric(as.factor(data[, exposure])) - 1

        wadenom_prob <- predict(exposure_regression, newdata = data,
                                type = ifelse(is_multinom_ereg | ereg == "multinomial", "probs", "response"))

        wadenom <-  wadenom_prob ^ category * (1 - wadenom_prob)^(1 - category)

      } else {

        wa_data <- data[, exposure, drop = FALSE]

        wa_data[, exposure] <- as.factor(wa_data[, exposure])

        category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)

        wadenom_prob <- predict(exposure_regression, newdata = data,
                                type = ifelse(is_multinom_ereg | ereg == "multinomial" |
                                                is_polr_ereg | ereg == "ordinal", "probs", "response"))

        wadenom <- rowSums(category * wadenom_prob)

      }

      wa <- wanom/wadenom

    }

    # calculate P(Mp=Mpi|A=ai)/P(Mp=Mpi|A=ai,M1=M1i,...,Mp-1=Mp-1i,C=ci,L=Li)
    if (length(wmdenom_regression) == 0) {

      wm <- NULL

    } else {

      wmnom <- matrix(nrow = n, ncol = length(wmnom_regression))

      wmdenom <- matrix(nrow = n, ncol = length(wmdenom_regression))

      for (i in 1:length(wmdenom_regression)) {

        if (is.null(wmdenom_regression[[i]])) {

          wmdenom[, i] <- rep(1, n)

        } else if (wmdenomreg[[i]] %in% c("logistic", "loglinear") |
                   (is_glm_wmdenomreg[i] && family_wmdenomreg[[i]] %in% c("binomial", "quasibinomial")) |
                   ((is_multinom_wmdenomreg[i] | wmdenomreg[[i]] == "multinomial") &&
                    length(unique(data[, mediator[i]])) == 2)) {

          category <- as.numeric(as.factor(data[, mediator[i]])) - 1

          wmdenom_prob <- predict(wmdenom_regression[[i]], newdata = data,
                                  type = ifelse(is_multinom_wmdenomreg[i] |
                                                  wmdenomreg[[i]] == "multinomial",
                                                "probs", "response"))

          wmdenom[, i] <-  wmdenom_prob ^ category * (1 - wmdenom_prob)^(1 - category)

        } else {

          wm_data <- data[, mediator[i], drop = FALSE]

          wm_data[, mediator[i]] <- as.factor(wm_data[, mediator[i]])

          category <- model.matrix(as.formula(paste("~0+", mediator[i], sep = "")), data = wm_data)

          wmdenom_prob <- predict(wmdenom_regression[[i]], newdata = data,
                                  type = ifelse(is_multinom_wmdenomreg[i] | wmdenomreg[[i]] == "multinomial" |
                                                  is_polr_wmdenomreg[i] | wmdenomreg[[i]] == "ordinal",
                                                "probs", "response"))

          wmdenom[, i] <- rowSums(category * wmdenom_prob)

        }

        if (is.null(wmnom_regression[[i]])) {

          wmnom[, i] <- rep(1, n)

        } else if (wmnomreg[[i]] %in% c("logistic", "loglinear") |
                   (is_glm_wmnomreg[i] && family_wmnomreg[[i]] %in% c("binomial", "quasibinomial")) |
                   ((is_multinom_wmnomreg[i] | wmnomreg[[i]] == "multinomial") &&
                    length(unique(data[, mediator])) == 2)) {

          category <- as.numeric(as.factor(data[, mediator[i]])) - 1

          wmnom_prob <- predict(wmnom_regression[[i]], newdata = data,
                                type = ifelse(is_multinom_wmnomreg[i] |
                                                wmnomreg[[i]] == "multinomial",
                                              "probs", "response"))

          wmnom[, i] <-  wmnom_prob ^ category * (1 - wmnom_prob)^(1 - category)

        } else {

          wm_data <- data[, mediator[i], drop = FALSE]

          wm_data[, mediator[i]] <- as.factor(wm_data[, mediator[i]])

          category <- model.matrix(as.formula(paste("~0+", mediator[i], sep = "")), data = wm_data)

          wmnom_prob <- predict(wmnom_regression[[i]], newdata = data,
                                type = ifelse(is_multinom_wmnomreg[i] | wmnomreg[[i]] == "multinomial" |
                                                is_polr_wmnomreg[i] | wmnomreg[[i]] == "ordinal",
                                              "probs", "response"))

          wmnom[, i] <- rowSums(category * wmnom_prob)

        }

      }

      wm <- wmnom/wmdenom

    }

    # calculate weights for each mediator regression and outcome regression
    if (is.null(wa) && is.null(wm)) {

      w_yreg <- w_mreg <- NULL

    } else if (!is.null(wa) && is.null(wm)) {

      w_yreg <- wa

      w_mreg <- replicate(length(mediator), wa)

    } else if (!is.null(wm)) {

      w_yreg <- switch((is.null(wa))+1, "1" = wa, "2" = rep(1, n))

      for (i in 1:ncol(wm)) w_yreg <- w_yreg * wm[, i]

      w_mreg <- switch((is.null(wa))+1, "1" = replicate(length(mediator), wa),
                       "2" = replicate(length(mediator), rep(1, n)))

      for (i in 2:length(mediator)) {

        for (j in 1:(i - 1)) w_mreg[, i] <- w_mreg[, i] * wm[, j]

      }

    }

    for (i in 1:length(mediator_regression)) mediator_regression[[i]] <-
        update(mediator_regression[[i]], weights = w_mreg[, i], formula. = formula(mediator_regression[[i]]))

    outcome_regression <- update(outcome_regression, weights = w_yreg, formula. = formula(outcome_regression))

    regressions <- list(mediator_regression = mediator_regression, outcome_regression =  outcome_regression)

  } else if (model == "wb") {

    outcome_regression <- update(outcome_regression, formula. = formula(outcome_regression))

    exposure_regression <- update(exposure_regression, formula. = formula(exposure_regression))

    regressions <- list(exposure_regression = exposure_regression, outcome_regression =  outcome_regression)


  } else if (model == "gformula") {

    outcome_regression <- update(outcome_regression, formula. = formula(outcome_regression))

    for (i in 1:length(mediator_regression)) mediator_regression[[i]] <-
        update(mediator_regression[[i]], formula. = formula(mediator_regression[[i]]))

    if (!is.null(postc_regression)) {

      for (i in 1:length(postc_regression)) postc_regression[[i]] <-
          update(postc_regression[[i]], formula. = formula(postc_regression[[i]]))

    }

    regressions <- list(postc_regression = postc_regression,
                        mediator_regression = mediator_regression,
                        outcome_regression = outcome_regression)

  } else if (model == "iorw") {

    if (ereg %in% c("logistic", "loglinear") |
        (is_glm_ereg && family_ereg %in% c("binomial", "quasibinomial")) |
        ((is_multinom_ereg | ereg == "multinomial") &&
         length(unique(data[, exposure])) == 2)) {

      category <- as.numeric(as.factor(data[, exposure])) - 1

      wadenom.prob <- predict(exposure_regression, newdata = data,
                              type = ifelse(is_multinom_ereg | ereg == "multinomial", "probs", "response"))

      wanom <- 1 - wadenom.prob

      wadenom <-  wadenom.prob ^ category * (1 - wadenom.prob)^(1 - category)

    } else {

      wa_data <- data[, exposure, drop = FALSE]

      wa_data[, exposure] <- as.factor(wa_data[, exposure])

      category <- model.matrix(as.formula(paste("~0+", exposure, sep = "")), data = wa_data)

      wadenom.prob <- predict(exposure_regression, newdata = data,
                              type = ifelse(is_multinom_ereg | ereg == "multinomial" |
                                              is_polr_ereg | ereg == "ordinal", "probs", "response"))

      wanom <- wadenom.prob[, 1]

      wadenom <- rowSums(category * wadenom.prob)

    }

    wa <- wanom/wadenom

    tot_regression <- update(outcome_regression, formula. = formula(outcome_regression))

    dir_regression <- update(outcome_regression, weights = wa, formula. = formula(outcome_regression))

    regressions <- list(tot_regression = tot_regression, dir_regression = dir_regression)

  } else if (model == "ne") {

    outcome_regression <- update(outcome_regression, formula. = formula(outcome_regression))

    regressions <- list(outcome_regression)

  }

  return(regressions)

}
