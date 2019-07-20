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
    outcome_regression <- coxph(as.formula(outcome_formula), data = data)
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


z_p <- function(s, n, yreg) {

  if (yreg == "linear") {

  z <- s$estimate / s$std.error

  pval <- 2 * pt(-abs(z), n - 1)

  } else {

    z <- (s$estimate-1) / s$std.error

    pval <- 2 * pt(-abs(z), n - 1)

  }

  return(data.frame(z = z, pval = pval))

}

format_row_boot <- function(boot.out, index = 1, conf = 0.95) {
  ci <- boot::boot.ci(boot.out, index = index, type = "perc", conf = conf)
  d <- cbind(as.data.frame(t(tail(ci$percent[1, ], 2))))
  return(d)
}


format_df_boot <- function(boot.out, conf = 0.95, n, yreg) {

  d_all <- NULL

  for (i in 1:7) {

    d <- format_row_boot(boot.out, i, conf = conf)

    d_all <- rbind(d_all, d)

  }

  estimate <- boot.out$t0

  bias <- apply(boot.out$t, 2, mean) - estimate

  std.error <- apply(boot.out$t, 2, sd)

  d_all <- cbind(estimate, bias, std.error, d_all)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "bias", "std.error", label_CI[1], label_CI[2])

  if (yreg == "linear") { rownames(d_all) <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
  } else { rownames(d_all) <- c("ORcde", "ORpnde", "ORtnde", "ORpnie", "ORtnie", "ORte", "ORpm") }

  return(cbind(d_all, z_p(d_all, n = n, yreg = yreg)))

}

format_df_delta <- function(delta.out, conf = 0.95, n, yreg) {

  d_all <- data.frame(matrix(NA, 7, 4))

  alpha <- (1 - conf) / 2

  z <- qnorm(1 - alpha)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "std.error", label_CI[1], label_CI[2])

  rownames(d_all) <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")

  d_all["cde", ] <- c(delta.out$cde, delta.out$se_cde,
                      delta.out$cde - z * delta.out$se_cde,
                      delta.out$cde + z * delta.out$se_cde)

  d_all["pnde", ] <- c(delta.out$pnde, delta.out$se_pnde,
                       delta.out$pnde - z * delta.out$se_pnde,
                       delta.out$pnde + z * delta.out$se_pnde)

  d_all["tnde", ] <- c(delta.out$tnde, delta.out$se_tnde,
                       delta.out$tnde - z * delta.out$se_tnde,
                       delta.out$tnde + z * delta.out$se_tnde)

  d_all["pnie", ] <- c(delta.out$pnie, delta.out$se_pnie,
                       delta.out$pnie - z * delta.out$se_pnie,
                       delta.out$pnie + z * delta.out$se_pnie)

  d_all["tnie", ] <- c(delta.out$tnie, delta.out$se_tnie,
                       delta.out$tnie - z * delta.out$se_tnie,
                       delta.out$tnie + z * delta.out$se_tnie)

  d_all["te", ] <- c(delta.out$te, delta.out$se_te,
                     delta.out$te - z * delta.out$se_te,
                     delta.out$te + z * delta.out$se_te)

  d_all["pm", ] <- c(delta.out$pm, delta.out$se_pm,
                     delta.out$pm - z * delta.out$se_pm,
                     delta.out$pm + z * delta.out$se_pm)

  if (yreg == "linear") { rownames(d_all) <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
   } else { rownames(d_all) <- c("ORcde", "ORpnde", "ORtnde", "ORpnie", "ORtnie", "ORte", "ORpm") }

  return(cbind(d_all, z_p(d_all, n = n, yreg = yreg)))

}

print.delta_out <- function(delta_out, digits = 2) {

  printCoefmat(delta_out, digits = digits, has.Pvalue = TRUE)

}

print.boot_out <- function(boot_out, digits = 2) {

  printCoefmat(boot_out, digits = digits, has.Pvalue = TRUE)

}
