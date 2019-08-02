create_formulas <- function(outcome = NULL, treatment = NULL, mediator = NULL, covariates = c(),
                            interaction = NULL, event = NULL,
                            mreg = c("linear", "logistic"),
                            yreg = c("linear", "logistic", "loglinear", "poisson","quasipoisson",
                                     "negbin", "coxph", "aft_exp", "aft_weibull")) {


  mediator_formula_basic <- paste(mediator, treatment, sep = " ~ ")

  outcome_formula_basic  <- paste(paste(outcome, treatment, sep = " ~ "),
                                  mediator,
                                  sep = " + ")

  if (interaction) {

    outcome_formula_basic <- paste(outcome_formula_basic,
                                   paste(treatment, mediator, sep = " * "),
                                   sep = " + ")
  }

  if (length(covariates) == 0) {

    mediator_formula <- mediator_formula_basic

    outcome_formula  <- outcome_formula_basic

  } else {

    mediator_formula <- paste(mediator_formula_basic,
                              paste(covariates, collapse = " + "),
                              sep = " + ")

    outcome_formula  <- paste(outcome_formula_basic,
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


run_regressions <- function(formulas = NULL, outcome, treatment, mediator, covariates,
                            interaction, event, yreg, mreg, data = NULL) {

  outcome_formula <- formulas$outcome_formula

  mediator_formula<- formulas$mediator_formula

  if (yreg == "linear") {
    outcome_regression  <- lm(outcome_formula, data = data)
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

  if (mreg == "linear") {
    mediator_regression <- lm(mediator_formula, data = data)
  }

  if (mreg == "logistic") {
    mediator_regression <- glm(mediator_formula, family = binomial(), data = data)
  }

  regressions <- list(mediator_regression = mediator_regression,
                      outcome_regression = outcome_regression)

  return(regressions)
}


get_coef <- function(regressions = NULL, outcome, treatment, mediator, covariates,
                     interaction, event, mreg, yreg) {

  mediator_regression <- regressions$mediator_regression

  outcome_regression <- regressions$outcome_regression

  ## Store coefficients from regressions
  betas  <- coefficients(mediator_regression)

  thetas <- coefficients(outcome_regression)

  ## Store covariances from regressions
  vcov_betas <- vcov(mediator_regression)

  vcov_thetas <- vcov(outcome_regression)

  if (mreg == "linear")
    variance = summary(mediator_regression)$sigma^2
  else variance = NULL

  if (yreg == "aft_weibull")
    vcov_thetas <- vcov_thetas[-nrow(vcov_thetas), -ncol(vcov_thetas)]

  ## Build block diagonal matrix
  vcov_block <- Matrix::bdiag(vcov_thetas, vcov_betas)

  coef <- list(outcome = outcome, treatment = treatment, mediator = mediator,
               covariates = covariates, interaction = interaction,
               event = event, mediator_reg = mreg, outcome_reg = yreg,betas = betas, thetas=thetas,
               vcov_betas = vcov_betas, vcov_thetas = vcov_thetas, vcov_block = vcov_block,
               variance = variance)

  return(coef)
}


bootstrap_step <- function(data, indices, outcome, treatment, mediator, covariates, vecc,
                                interaction, event, mreg, yreg, m_star, a_star, a) {

  data_boot <- data[indices, ]

  formulas <- create_formulas(outcome = outcome, treatment = treatment,
                              mediator = mediator, covariates = covariates, interaction = interaction,
                              event = event, mreg = mreg, yreg = yreg)

  regressions <- run_regressions(formulas = formulas, outcome = outcome, treatment = treatment,
                                 mediator = mediator, covariates = covariates, interaction = interaction,
                                 event = event, mreg = mreg, yreg = yreg, data = data_boot)

  coef <- get_coef(regressions = regressions, outcome = outcome, treatment = treatment,
                   mediator = mediator, covariates = covariates, interaction = interaction,
                   event = event, mreg = mreg, yreg = yreg)

  thetas <- coef$thetas

  betas <- coef$betas

  variance <- coef$variance

  covariatesTerm <- ifelse(!is.null(vecc), sum(betas[covariates]*t(vecc)), 0)

  interactionTerm <- ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")], 0)

  if (yreg == "linear") {

    if (mreg == "linear") {

      cde <- unname((thetas[treatment] + interactionTerm * m_star) * (a - a_star))

      pnde <- unname((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star +
                                                               covariatesTerm))*(a - a_star))

      tnde <- unname((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a +
                                                               covariatesTerm))*(a - a_star))

      pnie <- unname((thetas[mediator] * betas[treatment] +
                        interactionTerm * betas[treatment] * a_star) * (a - a_star))

      tnie <- unname((thetas[mediator] * betas[treatment] +
                        interactionTerm * betas[treatment] * a) * (a - a_star))

    }

    if (mreg == "logistic") {

      cde <- unname((thetas[treatment] + interactionTerm * m_star * (a - a_star)))

      pnde <- unname(thetas[treatment] * (a - a_star) + interactionTerm*(a - a_star) *
                       (exp(betas[1] + betas[treatment] * a_star + covariatesTerm) /
                          (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))

      tnde <- unname(thetas[treatment] * (a - a_star) + interactionTerm*(a - a_star) *
                       (exp(betas[1] + betas[treatment] * a + covariatesTerm) /
                          (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm))))

      pnie <- unname((thetas[mediator]+interactionTerm*a_star) * (exp(betas[1] +
                                                                        betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] +
                                                                                                                            betas[treatment] * a + covariatesTerm)) - exp(betas[1] + betas[treatment] * a_star +
                                                                                                                                                                            covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))

      tnie <- unname((thetas[mediator]+interactionTerm*a) * (exp(betas[1] +
                                                                   betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] +
                                                                                                                       betas[treatment] * a + covariatesTerm)) - exp(betas[1] + betas[treatment] * a_star +
                                                                                                                                                                       covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))

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

      cde_rr <- unname(exp(thetas[treatment] + interactionTerm * m_star * (a - a_star)))

      pnde_rr <- unname(exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star +
                                                                   covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                           0.5 * interactionTerm ^ 2 * variance * (a^2 - a_star ^ 2)))

      tnde_rr <- unname(exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a +
                                                                   covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                           0.5 * interactionTerm ^ 2 * variance * (a^2 - a_star ^ 2)))

      pnie_rr <- unname(exp((thetas[mediator] * betas[treatment] +
                            interactionTerm * betas[treatment] * a_star) * (a - a_star)))

      tnie_rr <- unname(exp((thetas[mediator] * betas[treatment] +
                            interactionTerm * betas[treatment] * a) * (a - a_star)))

      cde_err <- unname(exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                              interactionTerm*a*m_star- (thetas[mediator]+interactionTerm*a_star)*
                               (betas['(Intercept)']+betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                               0.5*(thetas[mediator]+interactionTerm*a_star)^2*coef$variance)-
                           exp(thetas[mediator]*m_star+interactionTerm*a_star*m_star-
                                 (thetas[mediator]+interactionTerm*a_star)*(betas['(Intercept)']+
                                                                                betas[treatment]*a_star+sum(betas[covariates]*t(vecc)))-
                                 0.5*(thetas[mediator]+interactionTerm*a_star)^2*coef$variance))

    }

    if (mreg == "logistic") {

      cde_rr <- unname(exp(thetas[treatment] + interactionTerm*m_star * (a - a_star)))

      pnde_rr <- unname((exp(thetas[treatment] * (a - a_star)) * (1 + exp(thetas[mediator] +
                                                                         interactionTerm * a + betas[1] + betas[treatment] * a_star +
                                                                         covariatesTerm))) /(1 + exp(thetas[mediator] + interactionTerm * a_star +
                                                                                                       betas[1] + betas[treatment] * a_star + covariatesTerm)))

      tnde_rr <- unname((exp(thetas[treatment] * (a - a_star)) * (1 + exp(thetas[mediator] +
                                                                         interactionTerm * a + betas[1] + betas[treatment] * a + covariatesTerm))) /
                       (1 + exp(thetas[mediator] + interactionTerm * a_star +  betas[1] +
                                  betas[treatment] * a + covariatesTerm)))

      pnie_rr <- unname(((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
                        (1 + exp(thetas[mediator] + interactionTerm * a_star + betas[1] +
                                   betas[treatment] * a + covariatesTerm))) /
                       ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) *
                          (1+ exp(thetas[mediator] + interactionTerm * a_star + betas[1] +
                                    betas[treatment] * a_star + covariatesTerm))))

      tnie_rr <- unname(((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
                        (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a +
                                   covariatesTerm))) / ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) *
                                                          (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a_star +
                                                                     covariatesTerm))))

      cde_err <- unname((exp(thetas[treatment]*(a-a_star)+thetas[mediator]*m_star+
                               interactionTerm*a*m_star)*(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    interactionTerm*a_star))-exp(thetas[mediator]*m_star+
                                                   interactionTerm*a_star*m_star)*(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))))/(1+exp(betas['(Intercept)']+
                  betas[treatment]*a_star+sum(betas[covariates]*t(vecc))+thetas[mediator]+
                    interactionTerm*a_star))))

    }

    intref_err <- pnde_rr - 1 - cde_err

    intmed_err <- tnie_rr * pnde_rr - pnde_rr - pnie_rr + 1

    pie_err <- pnie_rr - 1

    te_rr <- tnie_rr * pnde_rr

    pm <-  (pnde_rr * (tnie_rr - 1)) / (pnde_rr * tnie_rr - 1)

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
