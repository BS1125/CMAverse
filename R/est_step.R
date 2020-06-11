est_step <- function(data = NULL, indices = NULL, model = NULL,
                     outcome = NULL, event, exposure = NULL, mediator = NULL, EMint = NULL,
                     prec = NULL, postc = NULL,
                     yreg = NULL, mreg = NULL, ereg = NULL, postcreg = NULL,
                     wmnomreg = NULL, wmdenomreg = NULL,
                     astar = 0, a = 1, mval = NULL, yref = NULL, vecc = NULL,
                     estimation = "imputation") {

  if (is.null(data)) stop("Unspecified data")

  if (is.null(indices)) indices <- 1:nrow(data)

  data_boot <- data[indices, ]

  n <- nrow(data_boot)

  if (is.null(model)) stop("Unspecified model")

  if (!model %in% c("rb", "wb", "iorw", "msm", "g-formula")) stop("Unsupported model")

  # run all regressions needed
  formulas <- allformula(model = model,
                         yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                         wmnomreg = wmnomreg, wmdenomreg = wmdenomreg,
                         outcome = outcome, event = event,
                         exposure = exposure, mediator = mediator, EMint = EMint,
                         prec = prec, postc = postc)

  regressions <- allreg(formulas = formulas, data = data_boot, model = model,
                        yreg = yreg, mreg = mreg, ereg = ereg, postcreg = postcreg,
                        wmnomreg = wmnomreg, wmdenomreg = wmdenomreg)

  outcome_regression <- regressions$outcome_regression

  tot_regression <- regressions$tot_regression

  dir_regression <- regressions$dir_regression

  mediator_regression <- regressions$mediator_regression

  exposure_regression <- regressions$exposure_regression

  postc_regression <- regressions$postc_regression

  # retrieve classes and families of regressions and create class indicators for them
  if (class(outcome_regression) %in% c("rcreg", "simexreg")) {

    yreg_class <- class(outcome_regression$NAIVEreg)

  } else yreg_class <- class(outcome_regression)

  is_lm_yreg <- ifelse(identical(yreg_class, "lm"), TRUE, FALSE)

  is_glm_yreg <- ifelse(identical(yreg_class, c("glm", "lm")), TRUE, FALSE)

  is_gam_yreg <- ifelse(identical(yreg_class, c("gam", "glm", "lm")), TRUE, FALSE)

  if (is_glm_yreg|is_gam_yreg) {

    if (class(outcome_regression) %in% c("rcreg", "simexreg")) {

      family_yreg <- family(outcome_regression$NAIVEreg)$family

    } else family_yreg <- family(outcome_regression)$family

  }

  is_multinom_yreg <- ifelse(identical(yreg_class, c("multinom", "nnet")), TRUE, FALSE)

  is_polr_yreg <- ifelse(identical(yreg_class, "polr"), TRUE, FALSE)

  is_glmnb_yreg <- ifelse(identical(yreg_class, c("negbin", "glm", "lm")), TRUE, FALSE)

  is_survreg_yreg <- ifelse(identical(yreg_class, "survreg"), TRUE, FALSE)

  is_coxph_yreg <- ifelse(identical(yreg_class, "coxph"), TRUE, FALSE)


  if (class(tot_regression) %in% c("rcreg", "simexreg")) {

    totreg_class <- class(tot_regression$NAIVEreg)

  } else totreg_class <- class(tot_regression)

  is_lm_totreg <- ifelse(identical(totreg_class, "lm"), TRUE, FALSE)

  is_glm_totreg <- ifelse(identical(totreg_class, c("glm", "lm")), TRUE, FALSE)

  is_gam_totreg <- ifelse(identical(totreg_class, c("gam", "glm", "lm")), TRUE, FALSE)

  if (is_glm_totreg|is_gam_totreg) {

    if (class(tot_regression) %in% c("rcreg", "simexreg")) {

      family_totreg <- family(tot_regression$NAIVEreg)$family

    } else family_totreg <- family(tot_regression)$family

  } else family_totreg <- NULL

  is_multinom_totreg <- ifelse(identical(totreg_class, c("multinom", "nnet")), TRUE, FALSE)

  is_polr_totreg <- ifelse(identical(totreg_class, "polr"), TRUE, FALSE)

  is_glmnb_totreg <- ifelse(identical(totreg_class, c("negbin", "glm", "lm")), TRUE, FALSE)

  is_survreg_totreg <- ifelse(identical(totreg_class, "survreg"), TRUE, FALSE)

  is_coxph_totreg <- ifelse(identical(totreg_class, "coxph"), TRUE, FALSE)


  if (class(dir_regression) %in% c("rcreg", "simexreg")) {

    dirreg_class <- class(dir_regression$NAIVEreg)

  } else dirreg_class <- class(dir_regression)

  is_lm_dirreg <- ifelse(identical(dirreg_class, "lm"), TRUE, FALSE)

  is_glm_dirreg <- ifelse(identical(dirreg_class, c("glm", "lm")), TRUE, FALSE)

  is_gam_dirreg <- ifelse(identical(dirreg_class, c("gam", "glm", "lm")), TRUE, FALSE)

  if (is_glm_dirreg|is_gam_dirreg) {

    if (class(dir_regression) %in% c("rcreg", "simexreg")) {

      family_dirreg <- family(dir_regression$NAIVEreg)$family

    } else family_dirreg <- family(dir_regression)$family

  } else family_dirreg <- NULL

  is_multinom_dirreg <- ifelse(identical(dirreg_class, c("multinom", "nnet")), TRUE, FALSE)

  is_polr_dirreg <- ifelse(identical(dirreg_class, "polr"), TRUE, FALSE)

  is_glmnb_dirreg <- ifelse(identical(dirreg_class, c("negbin", "glm", "lm")), TRUE, FALSE)

  is_survreg_dirreg <- ifelse(identical(dirreg_class, "survreg"), TRUE, FALSE)

  is_coxph_dirreg <- ifelse(identical(dirreg_class, "coxph"), TRUE, FALSE)

  if (!identical(c(is_lm_totreg, is_glm_totreg, is_gam_totreg, is_multinom_totreg, is_polr_totreg,
                   is_glmnb_totreg, is_survreg_totreg, is_coxph_totreg),
                 c(is_lm_dirreg, is_glm_dirreg, is_gam_dirreg, is_multinom_dirreg, is_polr_dirreg,
                   is_glmnb_dirreg, is_survreg_dirreg, is_coxph_dirreg))) {

    stop("different classes for tot_regression and dir_regression")

  } else if (family_dirreg != family_totreg) stop("different families for tot_regression and dir_regression")


  if (length(mediator_regression) != 0) {

    mreg_class <- is_lm_mreg <- is_glm_mreg <- is_gam_mreg <- is_multinom_mreg <- is_polr_mreg <-
      is_glmnb_mreg <- is_survreg_mreg <- is_coxph_mreg <- family_mreg <- c()

    for (i in 1:length(mediator_regression)) {

      if (class(mediator_regression[[i]]) %in% c("rcreg", "simexreg")) {

        mreg_class <- c(mreg_class, class(mediator_regression[[i]]$NAIVEreg))

      } else mreg_class <- c(mreg_class, class(mediator_regression[[i]]))

      is_lm_mreg <- c(is_lm_mreg, ifelse(identical(mreg_class[i], "lm"), TRUE, FALSE))

      is_glm_mreg <- c(is_glm_mreg, ifelse(identical(mreg_class[i], c("glm", "lm")), TRUE, FALSE))

      is_gam_mreg <- c(is_gam_mreg, ifelse(identical(mreg_class[i], c("gam", "glm", "lm")), TRUE, FALSE))

      if (is_glm_mreg[i]|is_gam_mreg[i]) {

        if (class(mediator_regression[[i]]) %in% c("rcreg", "simexreg")) {

          family_mreg <- c(family_mreg, family(mediator_regression[[i]]$NAIVEreg)$family)

        } else family_mreg <- c(family_mreg, family(mediator_regression[[i]])$family)

      }

      is_multinom_mreg <- c(is_multinom_mreg, ifelse(identical(mreg_class[i], c("multinom", "nnet")), TRUE, FALSE))

      is_polr_mreg <- c(is_polr_mreg, ifelse(identical(mreg_class[i], "polr"), TRUE, FALSE))

      is_glmnb_mreg <- c(is_glmnb_mreg, ifelse(identical(mreg_class[i], c("negbin", "glm", "lm")), TRUE, FALSE))

      is_survreg_mreg <- c(is_survreg_mreg, ifelse(identical(mreg_class[i], "survreg"), TRUE, FALSE))

      is_coxph_mreg <- c(is_coxph_mreg, ifelse(identical(mreg_class[i], "coxph"), TRUE, FALSE))

    }

  } else mreg_class <- is_lm_mreg <- is_glm_mreg <- is_gam_mreg <- is_multinom_mreg <-
    is_polr_mreg <- is_glmnb_mreg <- is_survreg_mreg <- is_coxph_mreg <-NULL


  if (class(exposure_regression) %in% c("rcreg", "simexreg")) {

    ereg_class <- class(exposure_regression$NAIVEreg)

  } else ereg_class <- class(exposure_regression)

  is_lm_ereg <- ifelse(identical(ereg_class, "lm"), TRUE, FALSE)

  is_glm_ereg <- ifelse(identical(ereg_class, c("glm", "lm")), TRUE, FALSE)

  is_gam_ereg <- ifelse(identical(ereg_class, c("gam", "glm", "lm")), TRUE, FALSE)

  if (is_glm_ereg|is_gam_ereg) {

    if (class(exposure_regression) %in% c("rcreg", "simexreg")) {

      family_ereg <- family(exposure_regression$NAIVEreg)$family

    } else family_ereg <- family(exposure_regression)$family

  }

  is_multinom_ereg <- ifelse(identical(ereg_class, c("multinom", "nnet")), TRUE, FALSE)

  is_polr_ereg <- ifelse(identical(ereg_class, "polr"), TRUE, FALSE)

  is_glmnb_ereg <- ifelse(identical(ereg_class, c("negbin", "glm", "lm")), TRUE, FALSE)

  is_survreg_ereg <- ifelse(identical(ereg_class, "survreg"), TRUE, FALSE)

  is_coxph_ereg <- ifelse(identical(ereg_class, "coxph"), TRUE, FALSE)


  if (length(postc_regression) != 0) {

    postcreg_class <- c()

    is_lm_postcreg <- is_glm_postcreg <- is_gam_postcreg <- is_multinom_postcreg <- is_polr_postcreg <-
      is_glmnb_postcreg <- is_survreg_postcreg <- is_coxph_postcreg <- family_postcreg <- c()

    for (i in 1:length(postc_regression)) {

      if (class(postc_regression[[i]]) %in% c("rcreg", "simexreg")) {

        postcreg_class <- c(postcreg_class, class(postc_regression[[i]]$NAIVEreg))

      } else postcreg_class <- c(postcreg_class, class(postc_regression[[i]]))

      is_lm_postcreg <- c(is_lm_postcreg, ifelse(identical(postcreg_class[i], "lm"), TRUE, FALSE))

      is_glm_postcreg <- c(is_glm_postcreg, ifelse(identical(postcreg_class[i], c("glm", "lm")), TRUE, FALSE))

      is_gam_postcreg <- c(is_gam_postcreg, ifelse(identical(postcreg_class[i], c("gam", "glm", "lm")), TRUE, FALSE))

      if (is_glm_postcreg[i]|is_gam_postcreg[i]) {

        if (class(postc_regression[[i]]) %in% c("rcreg", "simexreg")) {

          family_postcreg <- c(family_postcreg, family(postc_regression[[i]]$NAIVEreg)$family)

        } else family_postcreg <- c(family_postcreg, family(postc_regression[[i]])$family)

      }

      is_multinom_postcreg <- c(is_multinom_postcreg, ifelse(identical(postcreg_class[i], c("multinom", "nnet")), TRUE, FALSE))

      is_polr_postcreg <- c(is_polr_postcreg, ifelse(identical(postcreg_class[i], "polr"), TRUE, FALSE))

      is_glmnb_postcreg <- c(is_glmnb_postcreg, ifelse(identical(postcreg_class[i], c("negbin", "glm", "lm")), TRUE, FALSE))

      is_survreg_postcreg <- c(is_survreg_postcreg, ifelse(identical(postcreg_class[i], "survreg"), TRUE, FALSE))

      is_coxph_postcreg <- c(is_coxph_postcreg, ifelse(identical(postcreg_class[i], "coxph"), TRUE, FALSE))

    }

  } else postcreg_class <- is_glm_postcreg <- is_gam_postcreg <- is_multinom_postcreg <-
    is_polr_postcreg <- is_glmnb_postcreg <- is_survreg_postcreg <- is_coxph_postcreg <- NULL


  # the index of the reference level for a categorical outcome
  if ((is_gam_yreg && (family_yreg == "multinom" | startsWith(family_yreg, "Ordered Categorical")))|
      (is_multinom_yreg && length(unique(data_boot[, outcome])) > 2)| is_polr_yreg) {

    if (is.null(yref)) {

      warnings("Unspecified yref, 1 is used")

      yref = 1

    } else {

      if (is.null(outcome)) stop("Unspecified outcome")

      yref = which(levels(as.factor(data_boot[, outcome])) == yref)

    }

  }


  ###################################################################################################
  #################################Closed-form Parameter Function Estimation#########################
  ###################################################################################################

  if (estimation == "paramfunc") {

    if (is.null(exposure) | is.null(mediator)) stop("Unspecified exposure or mediator")

    if (!(model %in% c("rb", "iorw"))) {

      stop("Closed-form parameter function estimation doesn't support the selected model")

    } else if (length(mediator) > 1 && model == "rb") {

      stop("Closed-form parameter function estimation only supports a single mediator for rb")

    } else if (!(yreg %in% c("linear", "logistic", "loglinear", "poisson",
                             "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")|
                 !(mreg[[1]] %in% c("linear", "logistic", "multinominal"))) && model == "rb") {

      stop("Closed-form parameter function estimation doesn't support the selected yreg or mreg for rb")

    } else if (!(yreg %in% c("linear", "logistic", "loglinear", "poisson",
                             "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) && model == "iorw") {

      stop("Closed-form parameter function estimation doesn't support the selected yreg for iorw")

    }

    # create indicator vector a and astar for categorical exposure
    if (is.factor(data_boot[, exposure]) | is.character(data_boot[, exposure])) {

      a_lev <- levels(as.factor(data_boot[, exposure]))

      a_lev <- a_lev[which(a_lev %in% unique(data_boot[, exposure]))]

      a <- as.numeric(a_lev == a)[-1]

      astar <- as.numeric(a_lev == astar)[-1]

      elevel <- length(a_lev)

    } else elevel <- 2

    if (model == "rb"){

      if(length(prec) != 0) {

        if (length(vecc) == 0) stop("Unspecified vecc")

        # after resampling, some levels of prec might be lost; recreate vecc
        vecc_new <- c()

        for (c in 1:length(prec)) {

          if (is.factor(data[, prec[c]]) | is.character(data[, prec[c]])) {

            c_lev <- levels(as.factor(data[, prec[c]]))

            c_lev_new_index <- which(c_lev %in% unique(data_boot[, prec[c]]))

            vecc_new <- c(vecc_new, vecc[1:(length(c_lev) - 1)][c_lev_new_index])

            vecc <- vecc[-c(1:(length(c_lev) - 1))]

          } else {

            vecc_new <- c(vecc_new, vecc[1])

            vecc <- vecc[-1]

          }

        }

        vecc <- vecc_new

      } else vecc <- null

      # create indicator vector mstar for categorical mediator
      if (is.factor(data_boot[, mediator]) | is.character(data_boot[, mediator])) {

        m_lev <- levels(as.factor(data_boot[, mediator]))

        m_lev <- m_lev[which(m_lev %in% unique(data_boot[, mediator]))]

        mstar <- as.numeric(m_lev == mval[[1]])[-1]

        mlevel <- length(m_lev)

      } else {

        mstar <- mval[[1]]

        mlevel <- 2

      }

      mediator_regression <- mediator_regression[[1]]

      # coefficients for the outcome regression
      thetas <- coef(outcome_regression)

      # coefficients for the mediator regression
      betas  <- as.vector(t(coef(mediator_regression)))

      # sigma hat
      if (mreg == "linear") { variance <- sigma(mediator_regression)^2
      } else variance <- NULL

      # intercept coefficient for the outcome regression
      theta0 <- thetas[1]

      # exposure coefficient for the outcome regression
      theta1 <- thetas[2:elevel]

      # mediator coefficient for the outcome regression
      theta2 <- thetas[(elevel+1):(elevel + mlevel - 1)]

      # exposure-mediator interaction coefficient for the outcome regression
      if (EMint == TRUE) {

        theta3 <- t(matrix(thetas[length(thetas)-(((elevel-1)*(mlevel-1)-1):0)],
                           ncol = mlevel - 1))

      } else theta3 <- t(matrix(rep(0, (elevel-1)*(mlevel-1)),
                                ncol = mlevel - 1))

      # intercept coefficient for the mediator regression
      beta0 <- betas[1+(0:(mlevel-2))*length(betas)/(mlevel-1)]

      # exposure coefficient for the mediator regression
      beta1 <- t(matrix(betas[rowSums(expand.grid(2:elevel, (0:(mlevel-2)) *
                                                    length(betas)/(mlevel-1)))], ncol = mlevel - 1))

      # t(vecc)%*%beta'2 for the mediator regression
      covariatesTerm <- sapply(0:(mlevel-2), function(x)
        ifelse(length(prec) == 0, 0, sum(betas[elevel + 1:length(vecc) +
                                                 x * length(betas)/(mlevel-1)] * vecc)))

      # closed-form parameter function estimation
      if (yreg == "linear") {

        if (mreg == "linear") {

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

        } else if (mreg %in% c("logistic", "multinomial")) {

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

        pm <- tnie / te

        intref <- pnde - cde

        intmed <- tnie - pnie

        pie <- pnie

        cde_prop <- cde/te

        intref_prop <- intref/te

        intmed_prop <- intmed/te

        pie_prop <- pie/te

        overall_pm <- (pie+intmed)/te

        overall_int <- (intref+intmed)/te

        overall_pe <- (intref+intmed+pie)/te

        out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie, te = te, pm = pm,
                 intref = intref, intmed = intmed, pie = pie,
                 cde_prop = cde_prop, intref_prop = intref_prop,
                 intmed_prop = intmed_prop, pie_prop = pie_prop,
                 overall_pm = overall_pm, overall_int = overall_int, overall_pe = overall_pe)

      } else if (!yreg == "linear") {

        if (mreg == "linear") {

          RRcde <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                (sum(theta3 * a) - sum(theta3 * astar)) * mstar))

          RRpnde <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                 (sum(theta3 * a) - sum(theta3 * astar)) *
                                 (beta0 + sum(beta1 * astar) +
                                    covariatesTerm + theta2  * variance) +
                                 0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar))))

          RRtnde <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                 (sum(theta3 * a) - sum(theta3 * astar)) *
                                 (beta0 + sum(beta1 * a) +
                                    covariatesTerm + theta2  * variance) +
                                 0.5 * variance * (sum(theta3 ^ 2 * a) - sum(theta3 ^ 2 * astar))))

          RRpnie <- unname(exp(theta2 * (sum(beta1 * a) - sum(beta1 * astar)) +
                                 sum(theta3  * astar) * (sum(beta1 * a) - sum(beta1 * astar))))

          RRtnie <- unname(exp(theta2 * (sum(beta1 * a) - sum(beta1 * astar)) +
                                 sum(theta3  * a) * (sum(beta1 * a) - sum(beta1 * astar))))

          ERRcde <- unname((exp(sum(theta1 * a) - sum(theta1 * astar) +
                                  sum(theta3  * a) * mstar) -
                              exp(sum(theta3 * astar) * mstar)) *
                             exp(theta2 * mstar - (theta2 + sum(theta3 * astar)) *
                                   (beta0 + sum(beta1 * astar) + covariatesTerm) -
                                   0.5 * (theta2 + sum(theta3 * astar)) ^ 2 * variance))

        } else if (mreg %in% c("logistic", "multinomial")) {

          RRcde <- unname(exp(sum(theta1 * a) - sum(theta1 * astar) +
                                ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a) -
                                ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar)))

          RRpnde <- unname((exp(sum(theta1 * a) - sum(theta1 * astar)) *
                              (1 + sum(exp(theta2 + theta3 %*% a +
                                             beta0 + beta1 %*% astar + covariatesTerm)))) /
                             (1 + sum(exp(theta2 + theta3 %*% astar +
                                            beta0 + beta1 %*% astar + covariatesTerm))))

          RRtnde <- unname((exp(sum(theta1 * a) - sum(theta1 * astar)) *
                              (1 + sum(exp(theta2 + theta3 %*% a +
                                             beta0 + beta1 %*% a + covariatesTerm)))) /
                             (1 + sum(exp(theta2 + theta3 %*% astar +
                                            beta0 + beta1 %*% a + covariatesTerm))))

          RRpnie <- unname(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                              (1 + sum(exp(theta2 + theta3 %*% astar +
                                             beta0 + beta1 %*% a + covariatesTerm)))) /
                             ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                                (1 + sum(exp(theta2 + theta3 %*% astar +
                                               beta0 + beta1 %*% astar + covariatesTerm)))))

          RRtnie <- unname(((1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) *
                              (1 + sum(exp(theta2 + theta3 %*% a +
                                             beta0 + beta1 %*% a + covariatesTerm)))) /
                             ((1 + sum(exp(beta0 + beta1 %*% a + covariatesTerm))) *
                                (1 + sum(exp(theta2 + theta3 %*% a +
                                               beta0 + beta1 %*% astar + covariatesTerm)))))

          ERRcde <- unname(exp(sum(theta2*mstar)) *
                             (exp(sum(theta1 * a) - sum(theta1 * astar) +
                                    ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% a)) -
                                exp(ifelse(sum(mstar) == 0, 0, theta3[which(mstar == 1),] %*% astar))) *
                             (1 + sum(exp(beta0 + beta1 %*% astar + covariatesTerm))) /
                             (1+ sum(exp(theta2 + theta3 %*% astar + beta0 +
                                           beta1 %*% astar + covariatesTerm))))

        }

        ERRintref <- RRpnde - 1 - ERRcde

        ERRintmed <- RRtnie * RRpnde - RRpnde - RRpnie + 1

        ERRpie <- RRpnie - 1

        RRte <- RRtnie * RRpnde

        pm <- (RRpnde * (RRtnie - 1)) / (RRte - 1)

        ERRte <- RRte - 1

        ERRcde_prop <- ERRcde/ERRte

        ERRintmed_prop <- ERRintmed/ERRte

        ERRintref_prop <- ERRintref/ERRte

        ERRpie_prop <- ERRpie/ERRte

        overall_pm <- (ERRpie+ERRintmed)/ERRte

        overall_int <- (ERRintref+ERRintmed)/ERRte

        overall_pe <- (ERRintref+ERRintmed+ERRpie)/ERRte

        out <- c(RRcde = RRcde, RRpnde = RRpnde, RRtnde = RRtnde, RRpnie = RRpnie,
                 RRtnie = RRtnie, RRte = RRte, pm = pm,
                 ERRcde = ERRcde, ERRintref = ERRintref,
                 ERRintmed = ERRintmed, ERRpie = ERRpie,
                 ERRcde_prop = ERRcde_prop, ERRintref_prop = ERRintref_prop,
                 ERRintmed_prop = ERRintmed_prop, ERRpie_prop = ERRpie_prop,
                 overall_pm = overall_pm, overall_int = overall_int,
                 overall_pe = overall_pe)

      }

    } else if (model == "iorw") {

      dir_coef <- coef(dir_regression)

      tot_coef <- coef(tot_regression)

      if (is.factor(data_boot[, exposure]) | is.character(data_boot[, exposure])) {

        dir_coef <- unname(dir_coef[2 + (0:(length(a_lev) - 2))])

        tot_coef <- unname(tot_coef[2 + (0:(length(a_lev) - 2))])

      } else {

        dir_coef <- unname(coef(dir_regression)[2])

        tot_coef <- unname(coef(tot_regression)[2])

      }

      if (yreg == "linear") {

        dir <- unname(sum(dir_coef %*% a) - sum(dir_coef %*% astar))

        dir <- unname(sum(tot_coef %*% a) - sum(tot_coef %*% astar))

        ind <- tot - dir

        out <- c(tot = tot, dir = dir, ind = ind)

      } else {

        RRdir <- unname(exp(sum(dir_coef %*% a) - sum(dir_coef %*% astar)))

        RRtot <- unname(exp(sum(tot_coef %*% a) - sum(tot_coef %*% astar)))

        RRind <- RRtot / RRdir

        out <- c(RRtot = RRtot, RRdir = RRdir, RRind = RRind)

      }

    }


    ###################################################################################################
    ###############################Direct Counterfactual Imputation Estimation#########################
    ###################################################################################################

  } else if (estimation == "imputation") {

    if (model %in% c("rb", "msm", "g-formula", "wb")) {

      if (model %in% c("rb", "msm", "g-formula")) {

        # simulate postc
        if (model == "g-formula" && !is.null(postc)) {

          # design matrices for simulating the first postc
          postcdesign_a <- data.frame(c(rep(a, n)),
                                      data_boot[, prec])

          postcdesign_astar <- data.frame(c(rep(astar, n)),
                                          data_boot[, prec])

          colnames(postcdesign_a) <- colnames(postcdesign_astar) <- c(exposure, prec)

          if (is.factor(data_boot[, exposure])) {

            postcdesign_a[, exposure] <- factor(postcdesign_a[, exposure], levels = levels(data_boot[, exposure]))

            postcdesign_astar[, exposure] <- factor(postcdesign_astar[, exposure], levels = levels(data_boot[, exposure]))

          }

          postc_a <- data.frame(matrix(nrow = n, ncol = length(postc)))

          postc_astar <- data.frame(matrix(nrow = n, ncol = length(postc)))

          colnames(postc_a) <- colnames(postc_astar) <- postc

          for (i in 1:length(postc)) {

            # design matrices for postc[i] simulation
            postcdesign_a <- cbind(postcdesign_a, postc_a[, i-1, drop = FALSE])

            postcdesign_astar <- cbind(postcdesign_astar, postc_astar[, i-1, drop = FALSE])

            if (((is_glm_postcreg[i] | is_gam_postcreg[i]) &&
                 (family_postcreg[i] %in% c("binomial", "quasibinomial")))|
                (is_multinom_postcreg[i] && length(unique(data_boot[, postc[i]])) == 2)) {

              prob_a <- predict(postc_regression[[i]], newdata = postcdesign_a,
                                type = ifelse(is_multinom_postcreg[i], "probs", "response"))

              # simulated postc[i] for exposure=a; indexed by 0 or 1
              pred_a <- rbinom(n, size = 1, prob = prob_a)

              prob_astar <- predict(postc_regression[[i]], newdata = postcdesign_astar,
                                    type = ifelse(is_multinom_postcreg[i], "probs", "response"))

              # simulated postc[i] for exposure=astar; indexed by 0 or 1
              pred_astar <- rbinom(n, size = 1, prob = prob_astar)

              postc_lev <- levels(as.factor(data_boot[, postc[i]]))

              postc_lev_new <- postc_lev[which(postc_lev %in% unique(data_boot[, postc[i]]))]

              # mid_a: simulated postc[i] for exposure=a; mid_astar: simulated postc[i] for exposure=astar
              if (is.numeric(data_boot[, postc[i]])) {

                mid_a <- as.numeric(postc_lev_new[pred_a + 1])

                mid_astar <- as.numeric(postc_lev_new[pred_astar + 1])

              } else if (is.factor(data_boot[, postc[i]])) {

                mid_a <- factor(postc_lev_new[pred_a + 1], levels = postc_lev)

                mid_astar <- factor(postc_lev_new[pred_astar + 1], levels = postc_lev)

              } else if (is.character(data_boot[, postc[i]])) {

                mid_a <- postc_lev_new[pred_a + 1]

                mid_astar <- postc_lev_new[pred_astar + 1]

              }

            } else if ((is_gam_postcreg[i] && (family_postcreg[i] == "multinom" |
                                               startsWith(family_postcreg[i],"Ordered Categorical")))|
                       (is_multinom_postcreg[i] && length(unique(data_boot[, postc[i]])) > 2)|
                       is_polr_postcreg[i]) {

              prob_a <- predict(postc_regression[[i]], newdata = postcdesign_a,
                                type = ifelse(is_gam_postcreg[i], "response", "probs"))

              # indices of the simulated postc[i] for exposure=a
              pred_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

              prob_astar <- predict(postc_regression[[i]], newdata = postcdesign_astar,
                                    type = ifelse(is_gam_postcreg[i], "response", "probs"))

              # indices of the simulated postc[i] for exposure=astar
              pred_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

              postc_lev <- levels(as.factor(data_boot[, postc[i]]))

              postc_lev_new <- postc_lev[which(postc_lev %in% unique(data_boot[, postc[i]]))]

              # mid_a: simulated postc[i] for exposure=a; mid_astar: simulated postc[i] for exposure=astar
              if (is.numeric(data_boot[, postc[i]])) {

                mid_a <- as.numeric(postc_lev_new[pred_a])

                mid_astar <- as.numeric(postc_lev_new[pred_astar])

              } else if (is.factor(data_boot[, postc[i]])) {

                mid_a <- factor(postc_lev_new[pred_a],
                                levels = postc_lev)

                mid_astar <- factor(postc_lev_new[pred_astar],
                                    levels = postc_lev)

              } else if (is.character(data_boot[, postc[i]])) {

                mid_a <- postc_lev_new[pred_a]

                mid_astar <- postc_lev_new[pred_astar]

              }

            } else {

              # mid_a: simulated postc[i] for exposure=a; mid_astar: simulated postc[i] for exposure=astar
              mid_a <- predict(postc_regression[[i]], newdata = postcdesign_a, type = "response")

              mid_astar <- predict(postc_regression[[i]], newdata = postcdesign_astar, type = "response")

            }

            postc_a[, i] <- mid_a

            postc_astar[, i] <- mid_astar

          }

        } else {

          postc_a <- postc_astar <- data.frame()[1:n, ]

        }

        # simulate mediator
        # design matrices for simulating the first mediator
        mdesign_a <- data.frame(c(rep(a, n)), data_boot[, prec], postc_a)

        mdesign_astar <- data.frame(c(rep(astar, n)), data_boot[,prec], postc_astar)

        colnames(mdesign_a) <- colnames(mdesign_astar) <- c(exposure, prec, postc)

        if (is.factor(data_boot[, exposure])) {

          mdesign_a[, exposure] <- factor(mdesign_a[, exposure], levels = levels(data_boot[, exposure]))

          mdesign_astar[, exposure] <- factor(mdesign_astar[, exposure], levels = levels(data_boot[, exposure]))

        }

        m_a <- data.frame(matrix(nrow = n, ncol = length(mediator)))

        m_astar <- data.frame(matrix(nrow = n, ncol = length(mediator)))

        colnames(m_a) <- colnames(m_astar) <- mediator

        for (i in 1:length(mediator)) {

          # design matrices for mediator[i] simulation
          mdesign_a <- cbind(mdesign_a, m_a[, i-1, drop = FALSE])

          mdesign_astar <- cbind(mdesign_astar, m_astar[, i-1, drop = FALSE])

          if (((is_glm_mreg[i] | is_gam_mreg[i]) && (family_mreg[i] %in% c("binomial", "quasibinomial")))|
              (is_multinom_mreg[i] && length(unique(data_boot[, mediator[i]])) == 2)) {

            prob_a <- predict(mediator_regression[[i]], newdata = mdesign_a,
                              type = ifelse(is_multinom_mreg[i], "probs", "response"))

            # simulated mediator[i] for exposure=a; indexed by 0 or 1
            pred_a <- rbinom(n, size = 1, prob = prob_a)

            prob_astar <- predict(mediator_regression[[i]], newdata = mdesign_astar,
                                  type = ifelse(is_multinom_mreg[i], "probs", "response"))

            # simulated mediator[i] for exposure=astar; indexed by 0 or 1
            pred_astar <- rbinom(n, size = 1, prob = prob_astar)

            m_lev <- levels(as.factor(data_boot[, mediator[i]]))

            m_lev_new <- m_lev[which(m_lev %in% unique(data_boot[, mediator[i]]))]

            # mid_a: simulated mediator[i] for exposure=a; mid_astar: simulated mediator[i] for exposure=astar
            if (is.numeric(data_boot[, mediator[i]])) {

              mid_a <- as.numeric(m_lev_new[pred_a + 1])

              mid_astar <- as.numeric(m_lev_new[pred_astar + 1])

            } else if (is.factor(data_boot[, mediator[i]])) {

              mid_a <- factor(m_lev_new[pred_a + 1], levels = m_lev)

              mid_astar <- factor(m_lev_new[pred_astar + 1], levels = m_lev)

            } else if (is.character(data_boot[, mediator[i]])) {

              mid_a <- m_lev_new[pred_a + 1]

              mid_astar <- m_lev_new[pred_astar + 1]

            }

          } else if ((is_gam_mreg[i] && (family_mreg[i] == "multinom" |
                                         startsWith(family_mreg[i], "Ordered Categorical"))) |
                     (is_multinom_mreg[i] && length(unique(data_boot[, mediator[i]])) > 2) |
                     is_polr_mreg[i]) {

            prob_a <- predict(mediator_regression[[i]], newdata = mdesign_a,
                              type = ifelse(is_gam_mreg[i], "response", "probs"))

            # indices of the simulated mediator[i] for exposure=a
            pred_a <- apply(prob_a, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

            prob_astar <- predict(mediator_regression[[i]], newdata = mdesign_astar,
                                  type = ifelse(is_gam_mreg[i], "response", "probs"))

            # indices of the simulated mediator[i] for exposure=astar
            pred_astar <- apply(prob_astar, 1, FUN = function(x) apply(t(rmultinom(1, 1, prob = x)), 1, which.max))

            m_lev <- levels(as.factor(data_boot[, mediator[i]]))

            m_lev_new <- m_lev[which(m_lev %in% unique(data_boot[, mediator[i]]))]

            # mid_a: simulated mediator[i] for exposure=a; mid_astar: simulated mediator[i] for exposure=astar
            if (is.numeric(data_boot[, mediator[i]])) {

              mid_a <- as.numeric(m_lev_new[pred_a])

              mid_astar <- as.numeric(m_lev_new[pred_astar])

            } else if (is.factor(data_boot[, mediator[i]])) {

              mid_a <- factor(m_lev_new[pred_a], levels = m_lev)

              mid_astar <- factor(m_lev_new[pred_astar], levels = m_lev)

            } else if (is.character(data_boot[, mediator[i]])) {

              mid_a <- m_lev_new[pred_a]

              mid_astar <- m_lev_new[pred_astar]

            }

          } else {

            # mid_a: simulated mediator[i] for exposure=a; mid_astar: simulated mediator[i] for exposure=astar
            mid_a <- predict(mediator_regression_pred[[i]], newdata = mdesign_a, type = "response")

            mid_astar <- predict(mediator_regression_pred[[i]], newdata = mdesign_astar, type = "response")

          }

          m_a[, i] <- mid_a

          m_astar[, i] <- mid_astar

        }

        # shuffle rows of the simulated mediator for the g formula approach
        if (model == "g-formula" && !is.null(prec)) {

          m_a <- m_a[sample(1:n, replace = FALSE), ]

          m_astar <- m_astar[sample(1:n, replace = FALSE), ]

        }

        # design matrices for outcome simulation
        ydesign0m <- data.frame(rep(astar, n), as.data.frame(mval)[rep(1, n),],
                                data_boot[,prec], postc_astar)

        ydesign1m <- data.frame(rep(a, n), as.data.frame(mval)[rep(1, n),],
                                data_boot[,prec], postc_a)

        ydesign00 <- data.frame(rep(astar, n), m_astar,
                                data_boot[,prec], postc_astar)

        ydesign01 <- data.frame(rep(astar, n), m_a,
                                data_boot[,prec], postc_astar)

        ydesign10 <- data.frame(rep(a, n), m_astar,
                                data_boot[,prec], postc_a)

        ydesign11 <- data.frame(rep(a, n), m_a,
                                data_boot[,prec], postc_a)

        colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
          colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, prec, postc)

        if (is.factor(data_boot[, exposure])) {

          ydesign0m[, exposure] <- factor(ydesign0m[, exposure], levels = levels(data_boot[, exposure]))

          ydesign1m[, exposure] <- factor(ydesign1m[, exposure], levels = levels(data_boot[, exposure]))

          ydesign00[, exposure] <- factor(ydesign00[, exposure], levels = levels(data_boot[, exposure]))

          ydesign10[, exposure] <- factor(ydesign10[, exposure], levels = levels(data_boot[, exposure]))

          ydesign01[, exposure] <- factor(ydesign01[, exposure], levels = levels(data_boot[, exposure]))

          ydesign11[, exposure] <- factor(ydesign11[, exposure], levels = levels(data_boot[, exposure]))

        }

        for (i in 1:length(mediator)) {

          if (is.factor(data_boot[, mediator[i]])) {

            ydesign0m[,1+i] <- factor(ydesign0m[, 1+i], levels = levels(data_boot[, mediator[i]]))

            ydesign1m[,1+i] <- factor(ydesign1m[, 1+i], levels = levels(data_boot[, mediator[i]]))

          }

        }

        if ((is_gam_yreg && (family_yreg == "multinom" |
                             startsWith(family_yreg, "Ordered Categorical")))|
            (is_multinom_yreg && length(unique(data_boot[, outcome])) > 2)| is_polr_yreg) {

          type = ifelse(inherits(outcome_regression, "gam"), "response", "probs")

          EY0m <- mean(predict(outcome_regression_pred, newdata =  ydesign0m, type = type)[, yref], na.rm = TRUE)

          EY1m <- mean(predict(outcome_regression_pred, newdata =  ydesign1m, type = type)[, yref], na.rm = TRUE)

          EY00 <- mean(predict(outcome_regression_pred, newdata =  ydesign00, type = type)[, yref], na.rm = TRUE)

          EY01 <- mean(predict(outcome_regression_pred, newdata =  ydesign01, type = type)[, yref], na.rm = TRUE)

          EY10 <- mean(predict(outcome_regression_pred, newdata =  ydesign10, type = type)[, yref], na.rm = TRUE)

          EY11 <- mean(predict(outcome_regression_pred, newdata =  ydesign11, type = type)[, yref], na.rm = TRUE)

        } else {

          type <- ifelse(is_coxph_yreg, "risk", ifelse(is_multinom_yreg, "probs", "response"))

          EY0m <- mean(predict(outcome_regression_pred, newdata =  ydesign0m, type = type), na.rm = TRUE)

          EY1m <- mean(predict(outcome_regression_pred, newdata =  ydesign1m, type = type), na.rm = TRUE)

          EY00 <- mean(predict(outcome_regression_pred, newdata =  ydesign00, type = type), na.rm = TRUE)

          EY01 <- mean(predict(outcome_regression_pred, newdata =  ydesign01, type = type), na.rm = TRUE)

          EY10 <- mean(predict(outcome_regression_pred, newdata =  ydesign10, type = type), na.rm = TRUE)

          EY11 <- mean(predict(outcome_regression_pred, newdata =  ydesign11, type = type), na.rm = TRUE)

        }

      } else if (model == "wb") {

        subj_astar <- which(data_boot[, exposure] == astar)

        subj_a <- which(data_boot[, exposure] == a)

        a_lev <- levels(as.factor(data_boot[, exposure]))

        a_lev <- a_lev[which(a_lev %in% unique(data_boot[, exposure]))]

        astar_index <- which(a_lev == astar)

        a_index <- which(a_lev == a)

        # P(A=a) and P(A=astar)
        wnom <- left_join(select(data_boot, exposure),
                          count(data_boot, !!as.name(exposure)),
                          by = exposure)[, "n"]/n

        # wdenom1: P(A=a|Ci) and wdenom0: P(A=astar|Ci)
        wdenom.prob <- predict(exposure_regression, newdata = data_boot,
                               type = ifelse(inherits(exposure_regression, "multinom")|
                                               inherits(exposure_regression, "polr"),
                                             "probs", "response"))

        if (((is_glm_ereg | is_gam_ereg) && family_ereg %in% c("binomial", "quasibinomial"))|
            (is_multinom_ereg && length(unique(data_boot[, exposure])) == 2)) {

          wdenom0 <- wdenom.prob^(astar_index - 1)*(1 - wdenom.prob)^(2 - astar_index)

          wdenom1 <- wdenom.prob^(a_index - 1)*(1 - wdenom.prob)^(2 - a_index)

        } else {

          wdenom0 <- wdenom.prob[, astar_lev]

          wdenom1 <- wdenom.prob[, a_lev]

        }

        #design matrices for outcome simulation
        ydesign0m <- data.frame(rep(astar, length(index_astar)), as.data.frame(mval)[rep(1,length(index_astar)),],
                                data_boot[index_astar,prec])

        ydesign1m <- data.frame(rep(a, length(index_a)), as.data.frame(mval)[rep(1,length(index_a)),],
                                data_boot[index_a,prec])

        ydesign01 <- data.frame(rep(astar, length(index_a)), data_boot[index_a, mediator],
                                data_boot[index_a, prec])

        ydesign10 <- data.frame(rep(a, length(index_astar)), data_boot[index_astar, mediator],
                                data_boot[index_astar, prec])

        colnames(ydesign0m) <- colnames(ydesign1m) <-
          colnames(ydesign10) <- colnames(ydesign01) <- c(exposure, mediator, prec)

        if (is.factor(data_boot[, exposure])) {

          ydesign0m[, exposure] <- factor(ydesign0m[, exposure], levels = levels(data_boot[, exposure]))

          ydesign1m[, exposure] <- factor(ydesign1m[, exposure], levels = levels(data_boot[, exposure]))

          ydesign10[, exposure] <- factor(ydesign10[, exposure], levels = levels(data_boot[, exposure]))

          ydesign01[, exposure] <- factor(ydesign01[, exposure], levels = levels(data_boot[, exposure]))

        }

        for (i in 1:length(mediator)) {

          if (is.factor(data_boot[, mediator[i]])) {

            ydesign0m[,1+i] <- factor(ydesign0m[, 1+i], levels = levels(data_boot[, mediator[i]]))

            ydesign1m[,1+i] <- factor(ydesign1m[, 1+i], levels = levels(data_boot[, mediator[i]]))

          }

        }

        if ((is_gam_yreg && (family_yreg == "multinom" |
                             startsWith(family_yreg, "Ordered Categorical")))|
            (is_multinom_yreg && length(unique(data_boot[, outcome])) > 2)|
            is_polr_yreg) {

          stop("Unsupported yreg for wb")

        } else {

          type <- ifelse(is_coxph_yreg, "risk", ifelse(is_multinom_yreg, "probs", "response"))

          EY0m <- weighted.mean(predict(outcome_regression, newdata = ydesign0m, type = type),
                                w = (wnom/wdenom0)[index_astar], na.rm = TRUE)

          EY1m <- weighted.mean(predict(outcome_regression, newdata = ydesign1m, type = type),
                                w = (wnom/wdenom1)[index_a], na.rm = TRUE)

          EY00 <- weighted.mean(data_boot[index_astar, outcome],
                                w = (wnom/wdenom0)[index_astar], na.rm = TRUE)

          EY11 <- weighted.mean(data_boot[index_a, outcome],
                                w = (wnom/wdenom1)[index_a], na.rm = TRUE)

          EY01 <- weighted.mean(predict(outcome_regression, newdata = ydesign01, type = type),
                                w = (wnom/wdenom1)[index_a], na.rm = TRUE)

          EY10 <- weighted.mean(predict(outcome_regression, newdata = ydesign10, type = type),
                                w = (wnom/wdenom0)[index_astar], na.rm = TRUE)

        }

      }

      if (is_lm_yreg | ((is_glm_yreg | is_gam_yreg) &&
                        (family_yreg %in% c("gaussian","Gamma","inverse.gaussian","quasi", "gaulss", "gevlss") |
                         startsWith(family_yreg, "Tweedie") | startsWith(family_yreg, "Beta regression") |
                         startsWith(family_yreg, "Scaled t")))) {

        cde <- EY1m - EY0m

        pnde <- EY10 - EY00

        tnde <- EY11 - EY01

        pnie <- EY01 - EY00

        tnie <- EY11 - EY10

        te <- tnie + pnde

        pm <- tnie / te

        intref <- pnde - cde

        intmed <- tnie - pnie

        pie <- pnie

        cde_prop <- cde/te

        intref_prop <- intref/te

        intmed_prop <- intmed/te

        pie_prop <- pie/te

        overall_pm <- (pie + intmed)/te

        overall_int <- (intref + intmed)/te

        overall_pe <- (intref + intmed + pie)/te

        out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie,
                 te = te, pm = pm,
                 intref = intref, intmed = intmed, pie = pie,
                 cde_prop = cde_prop, intref_prop = intref_prop,
                 intmed_prop = intmed_prop, pie_prop = pie_prop,
                 overall_pm = overall_pm, overall_int = overall_int, overall_pe = overall_pe)

        if (!is.null(postc) && model %in% c("msm", "g-formula")) {

          names(out) <- c("cde","rpnde","rtnde","rpnie","rtnie","te","pm",
                          "rintref","rintmed","rpie",
                          "cde_prop","intref_prop","intmed_prop","pie_prop",
                          "overall_pm", "overall_int", "overall_pe")

        }

      } else {

        RRcde <- EY1m/EY0m

        RRpnde <- EY10/EY00

        RRtnde <- EY11/EY01

        RRpnie <- EY01/EY00

        RRtnie <- EY11/EY10

        RRte <- RRtnie * RRpnde

        pm <- (RRpnde * (RRtnie - 1)) / (RRte - 1)

        ERRcde <- (EY1m-EY0m)/EY00

        ERRintref <- RRpnde - 1 - ERRcde

        ERRintmed <- RRtnie * RRpnde - RRpnde - RRpnie + 1

        ERRpie <- RRpnie - 1

        ERRte <- RRte - 1

        ERRcde_prop <- ERRcde/ERRte

        ERRintmed_prop <- ERRintmed/ERRte

        ERRintref_prop <- ERRintref/ERRte

        ERRpie_prop <- ERRpie/ERRte

        overall_pm <- (ERRpie + ERRintmed)/ERRte

        overall_int <- (ERRintref + ERRintmed)/ERRte

        overall_pe <- (ERRintref + ERRintmed + ERRpie)/ERRte

        out <- c(RRcde = RRcde, RRpnde = RRpnde, RRtnde = RRtnde, RRpnie = RRpnie,
                 RRtnie = RRtnie, RRte = RRte, pm = pm,
                 ERRcde = ERRcde, ERRintref = ERRintref,
                 ERRintmed = ERRintmed, ERRpie = ERRpie,
                 ERRcde_prop = ERRcde_prop, ERRintref_prop = ERRintref_prop,
                 ERRintmed_prop = ERRintmed_prop, ERRpie_prop = ERRpie_prop,
                 overall_pm = overall_pm, overall_int = overall_int,
                 overall_pe = overall_pe)

        if (!is.null(postc) && model %in% c("msm", "g-formula")) {

          names(out) <- c("RRcde","RRrpnde","RRrtnde","RRrpnie","RRrtnie","RRte","pm",
                          "ERRcde","ERRrintref","ERRrintmed","ERRrpie",
                          "ERRcde_prop","ERRrintref_prop","ERRrintmed_prop","ERRrpie_prop",
                          "overall_pm", "overall_int", "overall_pe")

        }

      }

    } else if (model == "iorw") {

      #design matrices for outcome simulation
      totdesign0 <- dirdesign0 <- data.frame(rep(astar, n), data_boot[, prec])

      totdesign1 <- dirdesign1 <- data.frame(rep(a, n), data_boot[, prec])

      colnames(totdesign0) <- colnames(totdesign1) <-
        colnames(dirdesign0) <- colnames(dirdesign1) <- c(exposure, prec)

      if (is.factor(data_boot[, exposure])) {

        totdesign0[, exposure] <- factor(totdesign0[, exposure], levels = levels(data_boot[, exposure]))

        totdesign1[, exposure] <- factor(totdesign1[, exposure], levels = levels(data_boot[, exposure]))

        dirdesign0[, exposure] <- factor(dirdesign0[, exposure], levels = levels(data_boot[, exposure]))

        dirdesign1[, exposure] <- factor(dirdesign1[, exposure], levels = levels(data_boot[, exposure]))

      }

      if ((is_gam_totreg && (family_totreg == "multinom" |
                             startsWith(family_totreg, "Ordered Categorical")))|
          (is_multinom_totreg && length(unique(data_boot[, outcome])) > 2)|
          is_polr_totreg) {

        type = ifelse(is_gam_totreg, "response", "probs")

        EYtot0 <- mean(predict(tot_regression, newdata = totdesign0, type = type)[, yref], na.rm = TRUE)

        EYtot1 <- mean(predict(tot_regression, newdata = totdesign1, type = type)[, yref], na.rm = TRUE)

        EYdir0 <- mean(predict(dir_regression, newdata = dirdesign0, type = type)[, yref], na.rm = TRUE)

        EYdir1 <- mean(predict(dir_regression, newdata = dirdesign1, type = type)[, yref], na.rm = TRUE)

      } else {

        type <- ifelse(is_coxph_totreg, "risk",
                       ifelse(is_multinom_totreg, "probs", "response"))

        EYtot0 <- mean(predict(tot_regression, newdata = totdesign0, type = type), na.rm = TRUE)

        EYtot1 <- mean(predict(tot_regression, newdata = totdesign1, type = type), na.rm = TRUE)

        EYdir0 <- mean(predict(dir_regression, newdata = dirdesign0, type = type), na.rm = TRUE)

        EYdir1 <- mean(predict(dir_regression, newdata = dirdesign1, type = type), na.rm = TRUE)

      }

      if (is_lm_totreg |
          ((is_glm_totreg|is_gam_totreg) &&
           (family_totreg %in% c("gaussian","Gamma","inverse.gaussian","quasi", "gaulss", "gevlss")|
            startsWith(family_totreg, "Tweedie") | startsWith(family_totreg, "Beta regression") |
            startsWith(family_totreg, "Scaled t")))) {

        tot <- EYtot1 - EYtot0

        dir <- EYdir1 - EYdir0

        ind <- tot - dir

        out <- c(tot = tot, dir = dir, ind = ind)

      } else {

        RRtot <- EYtot1 / EYtot0

        RRdir <- EYdir1 / EYdir0

        RRind <- tot / dir

        out <- c(RRtot = RRtot, RRdir = RRdir, RRind = RRind)

      }

    }

  }

  return(out)

}
