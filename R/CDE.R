CDE_bin <- function(thetas, treatment, mediator, m, a_star, a, interaction) {

  ORcde <- exp(thetas[treatment] + ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")] * m, 0) * (a - a_star))

  unname(ORcde)
}

CDE_cont <- function(thetas, treatment, mediator, m, a_star, a, interaction) {

  cde <- (thetas[treatment] + ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")] * m, 0)) * (a - a_star)

  unname(cde)
}

CDE_bin_delta <- function(thetas, treatment, mediator, m, a_star, a, interaction) {

  s <- ifelse(interaction,
              paste0("~ exp((x2 + x", length(thetas), " * m) * (a - a_star))"),
              paste0(" ~exp(x2 * (a - a_star))"))

  s <- stringr::str_replace_all(s, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a),
                                               "\\bm\\b" = as.character(m)))

  return(as.formula(s))
}

CDE_cont_delta <- function(thetas, treatment, mediator, m, a_star, a, interaction) {

  s <- ifelse(interaction,
              paste0("~ (x2 + x", length(thetas), " * m) * (a - a_star)"),
              paste0(" ~ x2 * (a - a_star)"))

  s <- stringr::str_replace_all(s, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a),
                                               "\\bm\\b" = as.character(m)))

  return(as.formula(s))
}

CDE_boot_function <- function(thetas, treatment, mediator, m, a_star, a, interaction,
                              mreg = "linear", yreg = "linear") {

  if (mreg != "linear" & yreg != "linear")
    cde <- CDE_bin(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  else if (mreg != "linear" & yreg == "linear")
    cde <- CDE_cont(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  else if (mreg == "linear" & yreg != "linear")
    cde <- CDE_bin(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  else if (mreg == "linear" & yreg == "linear")
    cde <- CDE_cont(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  return(cde)
}

CDE_delta_function <- function(thetas, treatment, mediator, m, a_star, a, interaction,
                               mreg = "linear", yreg = "linear") {

  if (mreg != "linear" & yreg != "linear")
    cded <- CDE_bin_delta(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  else if (mreg != "linear" & yreg == "linear")
    cded <- CDE_cont_delta(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  else if (mreg == "linear" & yreg != "linear")
    cded <- CDE_bin_delta(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  else if (mreg == "linear" & yreg == "linear")
    cded <- CDE_cont_delta(thetas = thetas, treatment = treatment, mediator = mediator, interaction = interaction, m = m, a_star = a_star, a = a)

  return(cded)
}

CDE_boot <- function(coef=list(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment

  mediator <- coef$mediator

  interaction <- coef$interaction

  return(CDE_boot_function(coef$thetas, treatment, mediator, m, a_star, a, interaction))
}

#CDE_boot(coef, m = 1, a_star = 0, a = 1)

CDE_delta <- function(coef=list(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment

  mediator <- coef$mediator

  interaction <- coef$interaction

  cde_delta <- CDE_delta_function(coef$thetas, treatment, mediator, m, a_star, a, interaction)

  return(msm::deltamethod(cde_delta, coef$thetas, coef$vcov_thetas))
}

#CDE_delta(coef, m = 1, a_star = 0, a = 1)
