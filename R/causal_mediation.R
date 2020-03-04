causal_mediation <- function(data = NULL, outcome = NULL, event = NULL,
                             exposure = NULL, exposure.type = "binary",
                             mediator = NULL, covariates.pre = NULL, covariates.post = NULL,
                             covariates.post.type = NULL, cval = NULL,
                             EMint = FALSE, MMint = FALSE, EMMint = FALSE,
                             EMint.terms = NULL, MMint.terms = NULL, EMMint.terms = NULL,
                             mreg = "linear", yreg = "linear",
                             mval = NULL, a_star = 0, a = 1,
                             model = "rb", est.method = "paramfunc", inf.method = "delta",
                             nboot = 200, nrep = 5) {

  require(dplyr)

  if (is.null(exposure) | is.null(outcome) | is.null(mediator)) {
    stop("Unspecified exposure, mediator, or outcome")
  }

  if (model %in% c("rb", "msm", "g-formula") &&
      ((length(mreg) < length(mediator)))) {
    stop("Unspecified mediator model")
  }

  if (is.null(yreg)) {
    stop("Unspecified outcome model")
  }

  if (!(model %in% c("rb", "wb", "ne", "msm", "iorw", "g-formula"))) {
    stop("Unsupported causal mediation model")
  }

  if (is.null(covariates.post) == FALSE && !(model %in% c("msm", "g-formula"))) {
    stop("Unsupported causal mediation model when post-exposure covariates exist.")
  }

  ##########################################reorganize the dataset################################

  if (exposure.type == "binary") {

    data[, exposure] <- as.factor(data[, exposure])

    a <- which(levels(data[,exposure]) == a) - 1

    a_star <- which(levels(data[,exposure]) == a_star) - 1

    levels(data[,exposure]) <- 0:1

    data[,exposure] <- as.numeric(data[,exposure]) - 1

  }

  if (!is.null(mval) && (length(mval) != length(mediator))) {
    stop("length(mval) != length(mediator)")
  } else if (is.null(mval)) mval <- list(rep(0, length(mediator)))

  m_star <- c()

  for (i in 1:length(mreg)) {

    if (mreg[[i]] %in% c("logistic", "multinomial")|inherits(mreg[[i]], "multinom")|
        (!is.character(mreg[[i]]) && (family(mreg[[i]])$family == "binomial"))) {

      data[,mediator[i]] <- as.factor(data[,mediator[i]])

      if (model == "wb") {
        m_star <- c(m_star, which(levels(data[,mediator[i]]) == mval[[i]])-1)
      } else {m_star <- c(m_star, as.numeric(levels(data[,mediator[i]]) == mval[[i]])[-1])}

      levels(data[,mediator[i]]) <- 0:(length(levels(data[,mediator[i]])) - 1)

      if (mreg[[i]] == "logistic"|
          (!is.character(mreg[[i]]) && (family(mreg[[i]])$family == "binomial"))) {
        data[,mediator[i]] <- as.numeric(data[,mediator[i]]) - 1
      } else {
        data <- data.frame(data, model.matrix(as.formula(paste("~", mediator[i])),
                                              data = data)[,-1])
      }

    } else {

      m_star <- c(m_star, mval[[i]])

    }
  }

  if (!is.null(covariates.post)) {

    for (i in 1:length(covariates.post)) {

      if (covariates.post.type[i] %in% c("binary", "categorical")) {

        data[,covariates.post[i]] <- as.factor(data[,covariates.post[i]])

        levels(data[,covariates.post[i]]) <- 0:(length(levels(data[,covariates.post[i]])) - 1)

        if (covariates.post.type[i] == "binary") {
          data[,covariates.post[i]] <- as.numeric(data[,covariates.post[i]]) - 1
        } else {
          data <- data.frame(data, model.matrix(as.formula(paste("~", covariates.post[i])),
                                                data = data)[,-1])
        }
      }
    }
  }

  if (model != "ne") {

    if (model == "rb" && est.method == "paramfunc") {

      if (is.null(covariates.pre)) {
        vecc <- NULL
      } else if (!is.null(covariates.pre)) {

        vecc <- c()

        if (is.null(cval)) cval <- as.list(rep(NA, length(covariates.pre)))

        for (i in 1:length(covariates.pre)) {

          if (is.numeric(data[, covariates.pre[i]])) {

            if (!is.null(cval[[i]])&&!is.na(cval[[i]])){vecc <- c(vecc, cval[[i]])
            } else  {vecc <- c(vecc, mean(data[, covariates.pre[i]],na.rm=TRUE))}

          } else if (!is.numeric(data[, covariates.pre[i]])) {
            data[,covariates.pre[i]] <- as.factor(data[,covariates.pre[i]])
            mid <- which(levels(data[,covariates.pre[i]]) == cval[[i]])
            levels(data[,covariates.pre[i]]) <- 0:(length(levels(data[,covariates.pre[i]]))-1)

            if (!is.null(cval[[i]])&&!is.na(cval[[i]])) {
              vecc <- c(vecc, as.numeric(levels(data[,covariates.pre[i]]) == mid - 1)[-1])
            } else  {vecc <- c(vecc, unname(colMeans(as.matrix(model.matrix(as.formula(paste0("~", covariates.pre[i])), data=data)[,-1]),na.rm = TRUE)))}

          }
        }
      }
    }


    n <- nrow(data)

    effect_estimate <- est_step(data = data, indices = 1:n,
                                outcome = outcome, exposure = exposure,
                                exposure.type = exposure.type, mediator = mediator,
                                covariates.pre = covariates.pre, covariates.post = covariates.post,
                                covariates.post.type = covariates.post.type, vecc = vecc,
                                EMint = EMint, MMint = MMint, EMMint = EMMint, EMint.terms = EMint.terms,
                                MMint.terms = MMint.terms, EMMint.terms = EMMint.terms,
                                event = event, mreg = mreg, yreg = yreg, m_star = m_star,
                                a_star = a_star, a = a, est.method = est.method, model = model)


    effect_se <- inf_step(data = data, nboot = nboot, outcome = outcome, exposure = exposure,
                          exposure.type = exposure.type, mediator = mediator,
                          covariates.pre = covariates.pre, covariates.post = covariates.post,
                          covariates.post.type = covariates.post.type, vecc = vecc,
                          EMint = EMint, MMint = MMint, EMMint = EMMint, EMint.terms = EMint.terms,
                          MMint.terms = MMint.terms, EMMint.terms = EMMint.terms,
                          event = event, mreg = mreg, yreg = yreg, m_star = m_star,
                          a_star = a_star, a = a,
                          est.method = est.method, inf.method = inf.method, model = model)

    out <- list(effect_estimate = effect_estimate, effect_se = effect_se)


  } else if (model == "ne") {

    formulas <- create_formulas(outcome=outcome, exposure=exposure, mediator=mediator,
                                covariates.pre=covariates.pre, covariates.post=covariates.post,
                                EMint=EMint, MMint=MMint, EMMint=EMMint,
                                EMint.terms=EMint.terms, MMint.terms=MMint.terms, EMMint.terms=EMMint.terms,
                                event=event, mreg=mreg, yreg=yreg, model=model, data=data)

    regressions <- run_regressions(model = model, formulas = formulas,
                                   exposure = exposure, exposure.type = exposure.type,
                                   mediator = mediator, covariates.post = covariates.post,
                                   covariates.post.type = covariates.post.type,
                                   mreg = mreg, yreg = yreg, data = data)

    expData <- medflex::neImpute(regressions$outcome_regression, data = data,
                                 nMed = length(mediator), nrep = nrep)

    if (EMint == TRUE) {

      medflex_formula <- paste(c(paste(outcome, paste0(exposure, "0"), sep = "~"), paste0(exposure, "1"),
                                 paste(paste0(exposure, "0"), paste0(exposure, "1"), sep = "*"), covariates.pre),
                               collapse = "+")

    } else if (EMint == FALSE) {

      medflex_formula <- paste(c(paste(outcome, paste0(exposure, "0"), sep = "~"), paste0(exposure, "1"),
                                 covariates.pre),
                               collapse = "+")

    }

    assign("medflex_formula", medflex_formula, envir = globalenv())

    family <- regressions$outcome_regression$family

    out <- summary(medflex::neEffdecomp(medflex::neModel(as.formula(medflex_formula),
                                                         family = family,
                                                         expData = expData, se = "robust")))$coefficients[, c("Estimate", "Std. Error")]

    if (yreg != "linear") {

      for (i in 1:length(out[, "Std. Error"]))
        out[, "Std. Error"][i] <- msm::deltamethod(as.formula("~exp(x1)"), out[, "Estimate"][i],
                                                   (out[, "Std. Error"][i])^2)

      out[, "Estimate"] <- exp(out[, "Estimate"])

    }

    rm(medflex_formula, pos = ".GlobalEnv")

  }

  return(out)

}

