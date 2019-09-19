create_formulas <- function(outcome, exposure, mediator, covariates,
                            EMint, MMint, EMMint, EMint.terms, MMint.terms, EMMint.terms,
                            event, mreg, yreg) {

  mediator_formula <- paste(mediator, exposure, sep = " ~ ")

  if (length(mediator) == 1) {

    outcome_formula  <- paste(paste(outcome, exposure, sep = " ~ "),
                                     mediator,
                                     sep = " + ")

    if (EMint == TRUE) outcome_formula  <- paste(outcome_formula,
                                                       paste(exposure, mediator, sep = " * "),
                                                       sep = " + ")

  } else if (length(mediator) > 1) {

    if (EMint == TRUE&&length(EMint.terms) == 0) {
      stop("Exposure-mediator interaction terms need to be specified")
    } else if (EMint == FALSE) EMint.terms <- NULL

    if (MMint == TRUE&&length(MMint.terms) == 0) {
      stop("Mediator-mediator interaction terms need to be specified")
    } else if (MMint == FALSE) MMint.terms <- NULL

    if (EMMint == TRUE&&length(EMMint.terms) == 0) {
      stop("Exposure-mediator-mediator interaction terms need to be specified")
    } else if (EMMint == FALSE) EMMint.terms <- NULL

    if (!is.null(c(EMint.terms, MMint.terms, EMMint.terms))) {
    if(!all(unlist(strsplit(c(EMint.terms, MMint.terms, EMMint.terms),
                            split = "[*]")) %in% c(exposure, mediator))){
      stop("Variables other than exposure and mediators in interaction terms")
      }
    }

    if (!is.null(EMMint.terms)) {

      int.main.effect <- unlist(lapply(1:length(EMMint.terms),
                 FUN = function(i)
                   sapply(2:(length( unlist(strsplit(EMMint.terms[i], split = "[*]")))-1),
                          FUN = function(x)
            sapply(1:length(combn(unlist(strsplit(EMMint.terms[i], split = "[*]")),
                                  m = x,simplify = FALSE)),
                   FUN = function(y) paste(combn(unlist(strsplit(EMMint.terms[i], split = "[*]")),
                                                 m = x,simplify = FALSE)[[y]],collapse = "*")))))

      int.terms <- unique(c(EMint.terms, MMint.terms, EMMint.terms, int.main.effect))

      } else int.terms <- unique(c(EMint.terms, MMint.terms))

     outcome_formula  <- paste(c(paste(outcome, exposure, sep = " ~ "), mediator,
                                  int.terms), collapse = " + ")

  }

  if (length(covariates) != 0) {

    mediator_formula <- paste(mediator_formula,
                              paste(covariates, collapse = " + "),
                              sep = " + ")

    outcome_formula  <- paste(outcome_formula,
                              paste(covariates, collapse = " + "),
                              sep = " + ")

  }

  if (yreg %in% c("coxph","aft_exp","aft_weibull")) {

    l <- strsplit(outcome_formula, split = " ~ ")

    l[[1]][1] <- paste0("Surv(", outcome, ", ", event, ")")

    outcome_formula <- paste(l[[1]][1], l[[1]][2], sep = " ~ ")
  }

  formulas <- list(outcome_formula = outcome_formula,
                   mediator_formula = mediator_formula)

  return(formulas)

}


run_regressions <- function(formulas = NULL, model, yreg, mreg, data = NULL) {

  outcome_formula <- formulas$outcome_formula

  mediator_formula <- formulas$mediator_formula

  if (yreg == "linear") {
    outcome_regression  <- glm(outcome_formula, family = gaussian(), data = data)
  }

  if (yreg == "logistic") {
    outcome_regression  <- glm(outcome_formula, family = binomial(), data = data)
  }

  if (yreg == "loglinear") {
    outcome_regression  <- glm(outcome_formula, family = binomial("log"), data = data)
  }

  if (yreg == "poisson") {
    outcome_regression  <- glm(outcome_formula, family = poisson(), data = data)
  }

  if (yreg == "quasipoisson") {
    outcome_regression  <- glm(outcome_formula, family = quasipoisson(), data = data)
  }

  if (yreg == "negbin") {
    outcome_regression  <- glm.nb(outcome_formula, data = data)
  }

  if (yreg == "coxph") {
    outcome_regression <- survival::coxph(as.formula(outcome_formula), data = data)
  }

  if (yreg == "aft_exp") {
    outcome_regression <- survreg(as.formula(outcome_formula), dist = "exponential", data = data)
  }

  if (yreg == "aft_weibull") {
    outcome_regression <- survreg(as.formula(outcome_formula), dist = "weibull", data = data)
  }

  if (!(model == "ne" | model == "wb")) {
  mediator_regression <- lapply(1:length(mediator_formula), FUN = function(x) {
    if (mreg[x] == "linear") {
      glm(mediator_formula[x],
          family = gaussian(),data = data)
    } else if (mreg[x] == "logistic") {
      glm(mediator_formula[x],
          family = binomial("logit"), data = data)}})
   } else mediator_regression <- NULL

  regressions <- list(mediator_regression = mediator_regression,
                      outcome_regression = outcome_regression)

  return(regressions)

}


get_coef <- function(regressions = NULL, mreg, yreg) {

  mediator_regression <- regressions$mediator_regression

  outcome_regression <- regressions$outcome_regression

  ## Store coefficients from regressions
  betas  <- lapply(1:length(mediator_regression), FUN = function(x) coefficients(mediator_regression[[x]]))

  thetas <- coefficients(outcome_regression)

  ## Store covariances from regressions
  vcov_betas <- lapply(1:length(mediator_regression), FUN = function(x) vcov(mediator_regression[[x]]))

  vcov_thetas <- vcov(outcome_regression)

  variance <- lapply(1:length(mediator_regression), FUN = function(x) {
      if (mreg[x] == "linear") { sigma(mediator_regression[[x]])^2
      } else if (mreg[x] == "logistic") NULL })

  if (yreg == "aft_weibull")
    vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

  ## Build block diagonal matrix
  vcov_block <- lapply(1:length(mediator_regression), FUN = function(x) Matrix::bdiag(vcov_thetas, vcov_betas[[x]]))

  coef <- list(betas = betas, thetas = thetas,
               vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
               variance = variance)

  return(coef)
}


bootstrap_step_singleM <- function(data, indices, outcome, exposure, exposure.type,
                           mediator, covariates, vecc, EMint, model,
                           event, mreg, yreg, m_star, a_star, a) {

  data_boot <- data[indices, ]

  formulas <- create_formulas(outcome = outcome, exposure = exposure,
                               mediator = mediator, covariates = covariates,
                               EMint = EMint, MMint = F, EMMint = F, EMint.terms = NULL,
                               MMint.terms = NULL, EMMint.terms = NULL,
                               event = event, mreg = mreg, yreg = yreg)

  regressions <- run_regressions(formulas = formulas, model = model,
                                 mreg = mreg, yreg = yreg, data = data_boot)

  coef <- get_coef(regressions = regressions, mreg = mreg, yreg = yreg)

  thetas <- coef$thetas

  betas <- coef$betas

  variance <- coef$variance

  covariatesTerm <- ifelse(!is.null(vecc), sum(betas[[1]][3:(2+length(vecc))]*t(vecc)), 0)

  EMintTerm <- ifelse(EMint, thetas[paste(exposure, mediator, sep = ":")], 0)

  if (yreg == "linear") {

    if (mreg == "linear") {

      cde <- unname((thetas[exposure] + EMintTerm * m_star) * (a - a_star))

      pnde <- unname((thetas[exposure] + EMintTerm * (betas[[1]]["(Intercept)"] +
                                  betas[[1]][exposure] * a_star + covariatesTerm)) * (a - a_star))

      tnde <- unname((thetas[exposure] + EMintTerm *
                        (betas[[1]]["(Intercept)"] + betas[[1]][exposure]*a +
                           covariatesTerm)) * (a - a_star))

      pnie <- unname((betas[[1]][exposure] * thetas[mediator] +
                        EMintTerm * betas[[1]][exposure] * a_star)*(a - a_star))

      tnie <- unname((betas[[1]][exposure] * thetas[mediator] +
                        EMintTerm * betas[[1]][exposure] * a)*(a - a_star))

    }

    if (mreg == "logistic") {

      cde <- unname((thetas[exposure] + EMintTerm * m_star * (a - a_star)))

      pnde <- unname(thetas[exposure] * (a - a_star) + EMintTerm*(a - a_star) *
                       (exp(betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm) /
                          (1 + exp(betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm))))

      tnde <- unname(thetas[exposure] * (a - a_star) + EMintTerm*(a - a_star) *
                       (exp(betas[[1]][1] + betas[[1]][exposure] * a + covariatesTerm) /
                          (1 + exp(betas[[1]][1] + betas[[1]][exposure] * a + covariatesTerm))))

      pnie <- unname((thetas[mediator]+EMintTerm*a_star) * (exp(betas[[1]][1] +
                                                                        betas[[1]][exposure] * a + covariatesTerm) / (1 + exp(betas[[1]][1] +
                                                                                                                            betas[[1]][exposure] * a + covariatesTerm)) - exp(betas[[1]][1] + betas[[1]][exposure] * a_star +
                                                                                                                                                                            covariatesTerm) / (1 + exp(betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm))))

      tnie <- unname((thetas[mediator]+EMintTerm*a) * (exp(betas[[1]][1] +
                                                                   betas[[1]][exposure] * a + covariatesTerm) / (1 + exp(betas[[1]][1] +
                                                                                                                       betas[[1]][exposure] * a + covariatesTerm)) - exp(betas[[1]][1] + betas[[1]][exposure] * a_star +
                                                                                                                                                                       covariatesTerm) / (1 + exp(betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm))))

      }

    intref <- pnde - cde

    intmed <- tnie - pnie

    pie <- pnie

    te <- tnie + pnde

    pm <- tnie / (pnde + te)

    cde_prop <- cde/te

    intref_prop <- intref/te

    intmed_prop <- intmed/te

    pie_prop <- pie/te

    overall_pm <- (pie+intmed)/te

    overall_int <- (intref+intmed)/te

    overall_pe <- (intref+intmed+pie)/te

    out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie,
             intref = intref, intmed = intmed, pie = pie,
             te = te, pm = pm, cde_prop = cde_prop, intref_prop = intref_prop,
             intmed_prop = intmed_prop, pie_prop = pie_prop,
             overall_pm = overall_pm, overall_int = overall_int, overall_pe = overall_pe)

  }

  if (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                  "negbin", "coxph", "aft_exp", "aft_weibull")) {

    if (mreg == "linear") {

      cde_rr <- unname(exp((thetas[exposure] + EMintTerm * m_star) * (a - a_star)))

      pnde_rr <- unname(exp((thetas[exposure] + EMintTerm * (betas[[1]][1] + betas[[1]][exposure] * a_star +
                                                                   covariatesTerm + thetas[mediator]  * variance[[1]])) * (a - a_star) +
                           0.5 * EMintTerm ^ 2 * variance[[1]] * (a^2 - a_star ^ 2)))

      tnde_rr <- unname(exp((thetas[exposure] + EMintTerm * (betas[[1]][1] + betas[[1]][exposure] * a +
                                                                   covariatesTerm + thetas[mediator]  * variance[[1]])) * (a - a_star) +
                           0.5 * EMintTerm ^ 2 * variance[[1]] * (a^2 - a_star ^ 2)))

      pnie_rr <- unname(exp((thetas[mediator] * betas[[1]][exposure] +
                            EMintTerm * betas[[1]][exposure] * a_star) * (a - a_star)))

      tnie_rr <- unname(exp((thetas[mediator] * betas[[1]][exposure] +
                            EMintTerm * betas[[1]][exposure] * a) * (a - a_star)))

      cde_err <- unname(exp(thetas[exposure]*(a-a_star)+thetas[mediator]*m_star+
                              EMintTerm*a*m_star- (thetas[mediator]+EMintTerm*a_star)*
                               (betas[[1]][1]+betas[[1]][exposure]*a_star+
                                  covariatesTerm)-
                               0.5*(thetas[mediator]+EMintTerm*a_star)^2*variance[[1]])-
                           exp(thetas[mediator]*m_star+EMintTerm*a_star*m_star-
                                 (thetas[mediator]+EMintTerm*a_star)*(betas[[1]][1]+
                                 betas[[1]][exposure]*a_star+covariatesTerm)-
                                 0.5*(thetas[mediator]+EMintTerm*a_star)^2*variance[[1]]))

    }

    if (mreg == "logistic") {

      cde_rr <- unname(exp((thetas[exposure] + EMintTerm*m_star) * (a - a_star)))

      pnde_rr <- unname((exp(thetas[exposure] * (a - a_star)) * (1 + exp(thetas[mediator] +
                                                                         EMintTerm * a + betas[[1]][1] + betas[[1]][exposure] * a_star +
                                                                         covariatesTerm))) /(1 + exp(thetas[mediator] + EMintTerm * a_star +
                                                                                                       betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm)))

      tnde_rr <- unname((exp(thetas[exposure] * (a - a_star)) * (1 + exp(thetas[mediator] +
                                                                         EMintTerm * a + betas[[1]][1] + betas[[1]][exposure] * a + covariatesTerm))) /
                       (1 + exp(thetas[mediator] + EMintTerm * a_star +  betas[[1]][1] +
                                  betas[[1]][exposure] * a + covariatesTerm)))

      pnie_rr <- unname(((1 + exp(betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm)) *
                        (1 + exp(thetas[mediator] + EMintTerm * a_star + betas[[1]][1] +
                                   betas[[1]][exposure] * a + covariatesTerm))) /
                       ((1 + exp(betas[[1]][1] + betas[[1]][exposure] * a + covariatesTerm)) *
                          (1+ exp(thetas[mediator] + EMintTerm * a_star + betas[[1]][1] +
                                    betas[[1]][exposure] * a_star + covariatesTerm))))

      tnie_rr <- unname(((1 + exp(betas[[1]][1] + betas[[1]][exposure] * a_star + covariatesTerm)) *
                        (1 + exp(thetas[mediator] + EMintTerm * a + betas[[1]][1] + betas[[1]][exposure] * a +
                                   covariatesTerm))) / ((1 + exp(betas[[1]][1] + betas[[1]][exposure] * a + covariatesTerm)) *
                                                          (1 + exp(thetas[mediator] + EMintTerm * a + betas[[1]][1] + betas[[1]][exposure] * a_star +
                                                                     covariatesTerm))))

      cde_err <- unname((exp(thetas[exposure]*(a-a_star)+thetas[mediator]*m_star+
                               EMintTerm*a*m_star)*(1+exp(betas[[1]]["(Intercept)"]+
                  betas[[1]][exposure]*a_star+covariatesTerm))/(1+exp(betas[[1]]["(Intercept)"]+
                  betas[[1]][exposure]*a_star+covariatesTerm+thetas[mediator]+
                    EMintTerm*a_star))-exp(thetas[mediator]*m_star+
                                                   EMintTerm*a_star*m_star)*(1+exp(betas[[1]]["(Intercept)"]+
                  betas[[1]][exposure]*a_star+covariatesTerm))/(1+exp(betas[[1]]["(Intercept)"]+
                  betas[[1]][exposure]*a_star+covariatesTerm+thetas[mediator]+
                    EMintTerm*a_star))))

    }

    intref_err <- pnde_rr - 1 - cde_err

    intmed_err <- tnie_rr * pnde_rr - pnde_rr - pnie_rr + 1

    pie_err <- pnie_rr - 1

    te_rr <- tnie_rr * pnde_rr

    pm <- (pnde_rr * (tnie_rr - 1)) / (pnde_rr * tnie_rr - 1)

    total_err <- te_rr - 1

    cde_err_prop <- cde_err/total_err

    intmed_err_prop <- intmed_err/total_err

    intref_err_prop <- intref_err/total_err

    pie_err_prop <- pie_err/total_err

    overall_pm <- (pie_err+intmed_err)/total_err

    overall_int <- (intref_err+intmed_err)/total_err

    overall_pe <- (intref_err+intmed_err+pie_err)/total_err

    out <- c(cde_rr = cde_rr, pnde_rr = pnde_rr, tnde_rr = tnde_rr, pnie_rr = pnie_rr,
             tnie_rr = tnie_rr, te_rr = te_rr, pm = pm,
             cde_err = cde_err, intref_err = intref_err,
             intmed_err = intmed_err, pie_err = pie_err, total_err = total_err,
             cde_err_prop = cde_err_prop, intref_err_prop = intref_err_prop,
             intmed_err_prop = intmed_err_prop, pie_err_prop = pie_err_prop,
             overall_pm = overall_pm, overall_int = overall_int,
             overall_pe = overall_pe)


  }

  return(out)

}


bootstrap_step_multipleM <- function(data, indices, outcome, exposure,
                                     mediator, covariates, vecc, model,
                                     EMint, MMint, EMMint, EMint.terms, MMint.terms, EMMint.terms,
                                     event, mreg, yreg, m_star, a_star, a) {

  data_boot <- data[indices, ]

  if (length(mediator) == 1&&EMint == TRUE) EMint.terms <- c(paste(exposure, mediator, sep = "*"))

  if (model == "rb") {

  formulas <- create_formulas(outcome = outcome, exposure = exposure,
                              mediator = mediator, covariates = covariates,
                              EMint = EMint, MMint = F, EMMint = F,
                              EMint.terms = EMint.terms, MMint.terms = NULL, EMMint.terms = NULL,
                              event = event, mreg = mreg, yreg = yreg)

  if (!is.null(EMint.terms)) {

  EMint.index <- sapply(1:length(EMint.terms),
                        FUN = function(x) which(mediator %in% strsplit(EMint.terms, split = "[*]")[[x]]))

  }

  regressions <- run_regressions(formulas = formulas, model = model, mreg = mreg, yreg = yreg, data = data_boot)

  coef <- get_coef(regressions = regressions, mreg = mreg, yreg = yreg)

  thetas <- coef$thetas

  betas <- coef$betas

  variance <- coef$variance

  covariatesTerm <- lapply(1:length(mediator), FUN = function (x) ifelse(!is.null(vecc),
                                                                         sum(betas[[x]][3:(2+length(vecc))]*t(vecc)), 0))

  if (EMint) { EMintTerm <- thetas[stringr::str_replace_all(EMint.terms,
                                                            pattern = "[*]", replacement = ":")]
  } else if (EMint == FALSE) EMintTerm <- 0

  if (yreg == "linear") {

    cde <- tnde <- pnde <- unname(thetas[exposure] * (a - a_star))

    tnie <- pnie <-  unname(sum(sapply(1:length(mediator), FUN = function(x) if (mreg[x] == "linear") {
      betas[[x]][exposure] * thetas[mediator[x]] * (a - a_star)
    } else if (mreg[x] == "logistic") {
      (exp(betas[[x]][1] + betas[[x]][exposure] * a + covariatesTerm[[x]]) / (1 + exp(betas[[x]][1] +
                                                                                        betas[[x]][exposure] * a + covariatesTerm[[x]])) -
         exp(betas[[x]][1] + betas[[x]][exposure] * a_star + covariatesTerm[[x]]) /
         (1 + exp(betas[[x]][1] + betas[[x]][exposure] * a_star + covariatesTerm[[x]]))) *
        thetas[mediator[x]]})))

    if (EMint == TRUE) {

      cde <- cde + sum(EMintTerm * m_star[EMint.index] * (a - a_star))

      pnde <- pnde + unname(sum(sapply(1:length(EMint.terms), FUN = function(x)
        if(mreg[EMint.index[x]] == "linear") {
          EMintTerm[x] * (betas[[EMint.index[x]]]["(Intercept)"] +
                            betas[[EMint.index[x]]][exposure] * a_star +
                            covariatesTerm[[EMint.index[x]]]) * (a - a_star)
        } else if (mreg[EMint.index[x]] == "logistic") {
          EMintTerm[x] * (exp(betas[[EMint.index[x]]][1] + betas[[EMint.index[x]]][exposure] * a_star +
                                covariatesTerm[[EMint.index[x]]]) / (1 + exp(betas[[EMint.index[x]]][1] +
                                                                               betas[[EMint.index[x]]][exposure] * a_star +
                                                                               covariatesTerm[[EMint.index[x]]]))) * (a - a_star)})))

      tnde <- tnde + unname(sum(sapply(1:length(EMint.terms), FUN = function(x)
        if(mreg[EMint.index[x]] == "linear") {
          EMintTerm[x] * (betas[[EMint.index[x]]][1] +
                            betas[[EMint.index[x]]][exposure] * a +
                            covariatesTerm[[EMint.index[x]]]) * (a - a_star)
        } else if (mreg[EMint.index[x]] == "logistic") {
          EMintTerm[x] * (exp(betas[[EMint.index[x]]][1] + betas[[EMint.index[x]]][exposure] * a +
                                covariatesTerm[[EMint.index[x]]]) / (1 + exp(betas[[EMint.index[x]]][1] +
                                                                               betas[[EMint.index[x]]][exposure] * a +
                                                                               covariatesTerm[[EMint.index[x]]]))) * (a - a_star)})))

      tnie <- tnie + unname(sum(sapply(1:length(EMint.terms), FUN = function(x)
        if(mreg[EMint.index[x]] == "linear") {
          EMintTerm[x]*betas[[EMint.index[x]]][exposure] * a * (a - a_star)
        } else if (mreg[EMint.index[x]] == "logistic") {
          EMintTerm[x] * (exp(betas[[EMint.index[x]]][1] + betas[[EMint.index[x]]][exposure] * a +
                                covariatesTerm[[EMint.index[x]]]) / (1 + exp(betas[[EMint.index[x]]][1] +
                                                                               betas[[EMint.index[x]]][exposure] * a + covariatesTerm[[EMint.index[x]]])) -
                            exp(betas[[EMint.index[x]]][1] + betas[[EMint.index[x]]][exposure] * a_star +
                                  covariatesTerm[[EMint.index[x]]]) / (1 + exp(betas[[EMint.index[x]]][1] +
                                                                                 betas[[EMint.index[x]]][exposure] * a_star +
                                                                                 covariatesTerm[[EMint.index[x]]]))) * a})))

      pnie <- pnie + unname(sum(sapply(1:length(EMint.terms), FUN = function(x)
        if(mreg[EMint.index[x]] == "linear") {
          EMintTerm[x]*betas[[EMint.index[x]]][exposure] * a_star * (a - a_star)
        } else if (mreg[EMint.index[x]] == "logistic") {
          EMintTerm[x] * (exp(betas[[EMint.index[x]]][1] + betas[[EMint.index[x]]][exposure] * a +
                                covariatesTerm[[EMint.index[x]]]) / (1 + exp(betas[[EMint.index[x]]][1] +
                                                                               betas[[EMint.index[x]]][exposure] * a + covariatesTerm[[EMint.index[x]]])) -
                            exp(betas[[EMint.index[x]]][1] + betas[[EMint.index[x]]][exposure] * a_star +
                                  covariatesTerm[[EMint.index[x]]]) / (1 + exp(betas[[EMint.index[x]]][1] +
                                                                                 betas[[EMint.index[x]]][exposure] * a_star +
                                                                                 covariatesTerm[[EMint.index[x]]]))) * a_star})))
    }

    te <- pnde + tnie

    pm <- tnie / (pnde + te)

    out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie, te = te, pm = pm)

  }


  if (yreg %in% c("logistic", "loglinear", "poisson", "quasipoisson",
                  "negbin", "coxph", "aft_exp", "aft_weibull")) {

    RRcde <- RRtnde <- RRpnde <- unname(exp(thetas[exposure] * (a - a_star)))

    RRtnie <- RRpnie <-  unname(exp(sum(sapply(1:length(mediator), FUN = function(x)
      betas[[x]][exposure] * thetas[mediator[x]] * (a - a_star)))))

    if (EMint == TRUE) {

      RRcde <- RRcde * exp(sum(EMintTerm * m_star[EMint.index] * (a - a_star)))

      if (length(mediator)==1) {Mcovariance <- unlist(variance)
      } else  Mcovariance <- diag(sqrt(unlist(variance)))%*%cor(data_boot[, mediator])%*%diag(sqrt(unlist(variance)))

      RRpnde <- RRpnde * unname(exp(sum(sapply(1:length(EMint.terms), FUN = function(x)
        EMintTerm[x] * (betas[[EMint.index[x]]]["(Intercept)"] + betas[[EMint.index[x]]][exposure] * a_star +
                          covariatesTerm[[EMint.index[x]]] +
                          thetas[mediator[EMint.index[x]]] * variance[[EMint.index[x]]]) * (a - a_star))) +
          as.numeric(0.5 * t(sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                      ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a, 0))) %*% Mcovariance %*%
                       (sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                 ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a, 0)))) - as.numeric(0.5 *
                                                                                                                         t(sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                                                                                                                    ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a_star, 0))) %*% Mcovariance %*%
                                                                                                                         (sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                                                                                                                   ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a_star, 0))))))

      RRtnde <- RRtnde *unname(exp(sum(sapply(1:length(EMint.terms), FUN = function(x)
        EMintTerm[x] * (betas[[EMint.index[x]]]["(Intercept)"] + betas[[EMint.index[x]]][exposure] * a +
                          covariatesTerm[[EMint.index[x]]] +
                          thetas[mediator[EMint.index[x]]] * variance[[EMint.index[x]]]) * (a - a_star))) +
          as.numeric(0.5 * t(sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                      ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a, 0))) %*% Mcovariance %*%
                       (sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                 ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a, 0)))) -
          as.numeric(0.5 * t(sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                      ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a_star, 0))) %*% Mcovariance %*%
                       (sapply(1:length(mediator), FUN = function(x) thetas[mediator[x]] +
                                 ifelse(x %in% EMint.index, EMintTerm[which(EMint.index == x)] * a_star, 0))))))

      RRtnie <- RRtnie * unname(exp(sum(sapply(1:length(EMint.terms), FUN = function(x)
        EMintTerm[x]*betas[[EMint.index[x]]][exposure] * a * (a - a_star)))))


      RRpnie <- RRpnie * unname(exp(sum(sapply(1:length(EMint.terms), FUN = function(x)
        EMintTerm[x]*betas[[EMint.index[x]]][exposure] * a_star * (a - a_star)))))
    }

    RRte <- RRpnde * RRtnie

    pm <- (RRpnde * (RRtnie - 1)) / (RRpnde * RRtnie - 1)

    out <- c(RRcde = RRcde, RRpnde = RRpnde, RRtnde = RRtnde, RRpnie = RRpnie, RRtnie = RRtnie,
               RRte = RRte, pm = pm)

    }
  } else if (model == "wb") {

    exposure_regression <- glm(paste(exposure, paste(covariates, collapse = "+"), sep = "~"),
                               family = binomial("logit"), data = data_boot)

    formulas <- create_formulas(outcome = outcome, exposure = exposure,
                                mediator = mediator, covariates = covariates,
                                EMint = EMint, MMint = MMint, EMMint = EMMint,
                                EMint.terms = EMint.terms,
                                MMint.terms = MMint.terms, EMMint.terms = EMMint.terms,
                                event = event, mreg = mreg, yreg = yreg)

    regressions <- run_regressions(formulas = formulas, model = model,
                                   mreg = mreg, yreg = yreg, data = data_boot)

    prob = predict.glm(exposure_regression, type = "response",
                       newdata = setNames(as.data.frame(data_boot[, covariates]), nm = covariates))

    w <- 1/dbinom(x = data_boot[, exposure], size = 1, prob = prob)

    data_cde <- data_boot[1,]
    data_cde[, exposure] <- a
    data_cde[, mediator] <- m_star
    EY1m <- as.numeric(predict(regressions$outcome_regression, newdata = data_cde,
                               type = ifelse(yreg != "coxph", "response", "risk")))

    data_cde[, exposure] <- a_star
    EY0m <- as.numeric(predict(regressions$outcome_regression, newdata = data_cde,
                               type = ifelse(yreg != "coxph", "response", "risk")))

    Yfitted <- predict(regressions$outcome_regression, newdata = data_boot,
                       type = ifelse(yreg != "coxph", "response", "risk"))

    data_cf <- data_boot
    data_cf[, exposure] <- 1 - data_boot[, exposure]
    Ycf <- predict(regressions$outcome_regression, newdata = data_cf,
                   type = ifelse(yreg != "coxph", "response", "risk"))

    EY00 <- weighted.mean(x = Yfitted[which(data_boot[, exposure] == a_star)],
                          w = w[which(data_boot[, exposure] == a_star)])
    EY11 <- weighted.mean(x = Yfitted[which(data_boot[, exposure] == a)],
                          w = w[which(data_boot[, exposure] == a)])
    EY10 <- weighted.mean(x = Ycf[which(data_boot[, exposure] == a_star)],
                          w = w[which(data_boot[, exposure] == a_star)])
    EY01 <- weighted.mean(x = Ycf[which(data_boot[, exposure] == a)],
                          w = w[which(data_boot[, exposure] == a)])

    if (yreg == "linear") {

      cde <- EY1m - EY0m

      pnde <- EY10 - EY00

      tnie <- EY11 - EY10

      tnde <- EY11 - EY01

      pnie <- EY01 - EY00

      te <- tnie + pnde

      pm <- tnie / (pnde + te)

      out <- c(cde = cde, pnde = pnde, tnde = tnde, pnie = pnie, tnie = tnie, te = te, pm = pm)

    } else if (yreg != "linear") {

      RRcde <- EY1m / EY0m

      RRpnde <- EY10 / EY00

      RRtnie <- EY11 / EY10

      RRtnde <- EY11 / EY01

      RRpnie <- EY01 / EY00

      RRte <- RRtnie * RRpnde

      pm <- (RRpnde * (RRtnie - 1)) / (RRpnde * RRtnie - 1)

      out <- c(RRcde = RRcde, RRpnde = RRpnde, RRtnde = RRtnde, RRpnie = RRpnie, RRtnie = RRtnie,
               RRte = RRte, pm = pm)

    }
   }

  return(out)

}




z_p <- function(s, n, nway, yreg) {

  z <- NULL

  pval <- NULL

  if (nway == 3) {

    if (yreg == "linear") {

      z <- s$estimate / s$std.error

      pval <- 2 * pt(-abs(z), n - 1)

      } else {

      z[1:6] <- (s$estimate[1:6]-1) / s$std.error[1:6]

      z[7] <- s$estimate[7] / s$std.error[7]

      pval <- 2 * pt(-abs(z), n - 1)

      }
    } else if (nway == 4) {

        z <- s$estimate / s$std.error

        pval <- 2 * pt(-abs(z), n - 1)

        }

  return(data.frame(z = z, pval = pval))

}

format_row_boot <- function(boot.out, index = 1, conf = 0.95) {
  ci <- boot::boot.ci(boot.out, index = index, type = "perc", conf = conf)
  d <- cbind(as.data.frame(t(tail(ci$percent[1, ], 2))))
  return(d)
}


format_df_boot <- function(boot.out, conf = 0.95, n, nway, yreg) {

  d_all <- NULL

  for (i in 1:length(boot.out$t0)) {

    d <- format_row_boot(boot.out, i, conf = conf)

    d_all <- rbind(d_all, d)

  }

  estimate <- boot.out$t0

  bias <- apply(boot.out$t, 2, mean) - estimate

  std.error <- apply(boot.out$t, 2, sd)

  d_all <- cbind(estimate, bias, std.error, d_all)

  rownames(d_all) <- names(boot.out$t0)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "bias", "std.error", label_CI[1], label_CI[2])


  return(cbind(d_all, z_p(d_all, n = n, nway = nway, yreg = yreg)))

}

format_df_delta <- function(delta.out, conf = 0.95, n, nway, yreg) {

  d_all <- data.frame(matrix(NA, length(delta.out)/2, 4))

  alpha <- (1 - conf) / 2

  z <- qnorm(1 - alpha)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "std.error", label_CI[1], label_CI[2])

  rownames(d_all) <- names(delta.out)[2*(1:(length(delta.out)/2))-1]

  for (name in rownames(d_all)) {

    d_all[name, ] <- c(unname(delta.out[name]), unname(delta.out[stringr::str_c(name,"se_delta",sep="_")]),
                       unname(delta.out[name]-z*delta.out[stringr::str_c(name,"se_delta",sep="_")]),
                       unname(delta.out[name]+z*delta.out[stringr::str_c(name,"se_delta",sep="_")]))

  }

  return(cbind(d_all, z_p(d_all, n = n, nway = nway, yreg = yreg)))

}

print.delta_out <- function(delta_out, digits = 2) {

  printCoefmat(delta_out, digits = digits, has.Pvalue = TRUE)

}

print.boot_out <- function(boot_out, digits = 2) {

  printCoefmat(boot_out, digits = digits, has.Pvalue = TRUE)

}
