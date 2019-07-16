total_NDE_binbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                             a_star, a, interaction) {
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])

  ORnde <- (exp(thetas[treatment] * a) *
              (1 + exp(thetas[mediator] +
                         interactionTerm * a + betas[1] +
                         betas[treatment] * a + covariatesTerm))) /
    (exp(thetas[treatment] * a_star) *
       (1 + exp(thetas[mediator] + interactionTerm * a_star +
                  betas[1] + betas[treatment] * a + covariatesTerm)))
  unname(ORnde)
}

pure_NDE_binbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                            a_star, a, interaction) {
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])

  ORnde <- (exp(thetas[treatment] * a) *
              (1 + exp(thetas[mediator] +
                         interactionTerm * a + betas[1] + betas[treatment] * a_star +
                         covariatesTerm))) /
    (exp(thetas[treatment] * a_star) *
       (1 + exp(thetas[mediator] + interactionTerm * a_star +
                  betas[1] + betas[treatment] * a_star + covariatesTerm)))
  unname(ORnde)
}

pure_NDE_bincont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                             a_star, a, interaction) {
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])

  nde <- thetas[treatment] * (a - a_star) +
    interactionTerm*(a - a_star) *
    (exp(betas[1] + betas[treatment] * a_star + covariatesTerm) /
       (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)))

  unname(nde)
}

total_NDE_bincont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              a_star, a, interaction) {
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]), 0 , thetas[paste(treatment, mediator, sep = ":")])

  nde <- thetas[treatment] * (a - a_star) +
    interactionTerm*(a - a_star) *
    (exp(betas[1] + betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)))

  unname(nde)
}

pure_NDE_contbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                             variance, a_star, a, interaction) {
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  ORnde <- exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star + covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                 0.5 * interactionTerm ^ 2 * variance*(a^2 - a_star ^ 2))

  unname(ORnde)
}


total_NDE_contbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              variance, a_star, a, interaction){
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])

  ORnde <- exp((thetas[treatment] +
                  interactionTerm * (betas[1] + betas[treatment] * a +
                                       covariatesTerm + thetas[mediator] * variance)) *
                 (a - a_star) + 0.5 * interactionTerm ^ 2 * variance*(a^2 - a_star ^ 2))

  unname(ORnde)
}


pure_NDE_contcont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              a_star, a, interaction){
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])

  nde <- (thetas[treatment] + interactionTerm * betas[1] + interactionTerm * betas[treatment] * a_star + interactionTerm*covariatesTerm)*(a - a_star)
  unname(nde)
}


total_NDE_contcont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                               a_star, a, interaction){
  covariatesTerm <- 0
  for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])

  nde <- (thetas[treatment] + interactionTerm * betas[1] + interactionTerm * betas[treatment] * a + interactionTerm*covariatesTerm)*(a - a_star)
  unname(nde)
}

pure_NDE_contcont_delta <- function(thetas, vecc, interaction = TRUE,
                                    a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j) {
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  f <- "x2"
  if (interaction) {
    f <- paste0("(", f, " + x", k, " * x", k + 1, " + x", k, " * x", k + 2, "*a_star + ")

    fc <- ifelse(j > 0, paste0("x", k + 2 + 1:j, "  * ", "vecc_", 1:j, collapse = " + "), "0")
    fc <- paste0("x", k, " * (", fc, ")")

    f <- paste0(f, fc, ")")
  }
  f <- paste0(" ~ ", f, "*(a-a_star)")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j) {
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}

total_NDE_contcont_delta <- function(thetas, vecc, interaction = TRUE,
                                     a_star, a) {

  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j) {
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  f <- "x2"
  if (interaction) {
    f <- paste0("(", f, " + x", k, " * x", k + 1, " + x", k, "*x", k + 2, "*a + ")

    fc <- ifelse(j > 0, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), "0")
    fc <- paste0("x",k," *(", fc, ")")

    f <- paste0(f, fc, ")")
  }
  f <- paste0(" ~ ", f, "*(a-a_star)")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j) {
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}


pure_NDE_contbin_delta <- function(thetas, vecc, variance, interaction,
                                   a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j) {
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  f <- "exp((x2"
  if (interaction) {
    f <- paste0(f," + x",k," * (x", k + 1," + x", k + 2,"* a_star", " + ",
                paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), " + ",
                "x3 * ", "variance",
                ")",")")
    f <- paste0(f, "*(a-a_star)")
    f <- paste0(f, " + 0.5 * x", k," * x", k," * ", "variance", "* (a * a - a_star * a_star))")
  } else {
    f <- paste0(f,")")
    f <- paste0(f, "*(a-a_star))")
  }

  f <- paste0(" ~ ", f)

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a), "\\bvariance\\b" =
                                                 as.character(variance)))

  if (j > 0) {
    for (i in 1:j) {
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}


total_NDE_contbin_delta <- function(thetas, vecc, variance, interaction,
                                    a_star, a) {

  j <- length(vecc)
  k <- length(thetas)

  for (i in 1:j) {
    assign(paste("vecc", i, sep = "_"), vecc[i])
  }

  f <- "exp((x2"
  if (interaction) {
    f <- paste0(f," + x",k,"*(x", k + 1," + x", k + 2,"*a", " + ",
                paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), " + ",
                "x3*", "variance",
                ")",")")
    f <- paste0(f, "*(a-a_star)")
    f <- paste0(f, " + 0.5*x",k,"*x",k,"*","variance","*(a*a-a_star*a_star))")
  }else{
    f <- paste0(f,")")
    f <- paste0(f, "*(a-a_star))")
  }

  f <- paste0(" ~ ", f)

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a),
                                               "\\bvariance\\b" = as.character(variance)))

  for (i in 1:j) {
    ss <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
  }

  return(as.formula(ss))
}


pure_NDE_bincont_delta <- function(thetas, vecc, interaction = TRUE,
                                   a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  for (i in 1:j) {
    assign(paste("vecc", i, sep = "_"), vecc[i])
  }

  f <- "~ x2 * (a-a_star)"
  if (interaction) {
    f <- paste0(f, " + (x",k," * (a-a_star)")

    # Numerator
    N2 <- paste0("x", k + 1:2, collapse = " + " )
    N2 <- paste0("exp(", N2, "*a_star + ")
    N2 <- paste(N2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")")
    # Denominator
    D2 <- paste0("x", k + 1:2, collapse = " + " )
    D2 <- paste0("1 + exp(", D2, "*a_star + ")
    D2 <- paste(D2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")")
    D2 <- paste0("(", D2, ")")

    F2 <- paste0("(",N2,"/",D2,")")

    f <- paste0(f,"*",F2, ")")

  }

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  for (i in 1:j) {
    ss <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
  }

  return(as.formula(ss))
}


total_NDE_bincont_delta <- function(thetas, vecc, interaction,
                                    a_star, a) {

  j <- length(vecc)
  k <- length(thetas)

  for (i in 1:j) {
    assign(paste("vecc", i, sep = "_"), vecc[i])
  }

  f <- "~ x2 * (a-a_star)"
  if (interaction) {
    f <- paste0(f, " + (x", k, " * (a - a_star)")

    # Numerator
    N2 <- paste0("x", k + 1:2, collapse = " + " )
    N2 <- paste0("exp(", N2, "*a + ")
    N2 <- paste(N2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")")
    # Denominator
    D2 <- paste0("x", k + 1:2, collapse = " + " )
    D2 <- paste0("1 + exp(", D2, "* a + ")
    D2 <- paste(D2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")")
    D2 <- paste0("(", D2, ")")

    F2 <- paste0("(", N2, "/", D2, ")")

    f <- paste0(f, "*", F2, ")")
  }

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  for (i in 1:j) {
    ss <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
  }

  return(as.formula(ss))

}


pure_NDE_binbin_delta <- function(thetas, vecc, interaction,
                                  a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  for (i in 1:j) {
    assign(paste("vecc", i, sep = "_"), vecc[i])
  }

  # First term in nominator
  N1 <- "exp(x2 * a)"

  # Second term in nominator
  if (interaction) {
    N2 <- paste0("(1+exp(x3 + x", k, "* a + x", k + 1," + x", k + 2, "* a_star + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
  } else {
    N2 <- paste0("(1+exp(x3 + x", k + 1," + x", k + 2, "* a_star + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
  }

  # Construct nominator
  N <- paste(N1, N2, sep = "*")

  # First term in denominator
  D1 <- "exp(x2 * a_star)"

  # Second term in denominator
  if (interaction) {
    D2 <- paste0("(1 + exp(x3 + x", k, "* a_star + x", k + 1," + x", k + 2, "* a_star + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
  } else {
    D2 <- paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "* a_star + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
  }

  # Construct denominator
  D <- paste(D1, D2, sep = "*")

  # Construct formula
  f <- paste0(" ~ (", N, "/", D, ")")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j) {
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}

total_NDE_binbin_delta <- function(thetas, vecc, interaction, a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j) {
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  # First term in nominator
  N1 <- "exp(x2 * a)"

  # Second term in nominator
  if (j > 0) {
    if (interaction) {
      N2 <- paste0("(1 + exp(x3 + x", k,"* a + x",k + 1," + x", k + 2, "* a + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
    } else {
      N2 <- paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "* a + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
    }
  } else {
    if (interaction) {
      N2 <- paste0("(1 + exp(x3 + x",k,"* a + x", k + 1," + x", k + 2, "* a", "))")
    } else {
      N2 <- paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "* a", "))")
    }
  }

  # Construct nominator
  N <- paste(N1, N2, sep = "*")

  # First term in denominator
  D1 <- "exp(x2 * a_star)"

  # Second term in denominator
  if (j > 0) {
    if (interaction) {
      D2 <- paste0("(1 + exp(x3 + x", k," * a_star + x", k + 1," + x", k + 2, "* a + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
    }else{
      D2 <- paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "*a + ", paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))")
    }
  } else {
    if (interaction) {
      D2 <- paste0("(1 + exp(x3 + x", k, " * a_star + x", k + 1," + x", k + 2, "*a", "))")
    } else {
      D2 <- paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "*a", "))")
    }
  }

  # Construct denominator
  D <- paste(D1, D2, sep = "*")

  # Construct formula
  f <- paste0(" ~ (",N,"/",D,")")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j) {
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}

NDE_boot_function <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              m, interaction, a_star, a, variance,
                              mreg, yreg) {
  if (mreg != "linear" & yreg != "linear") {
    pnde <- pure_NDE_binbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                            covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnde <- total_NDE_binbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                             covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg != "linear" & yreg == "linear") {
    pnde <- pure_NDE_bincont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                             covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnde <- total_NDE_bincont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                              covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg == "linear" & yreg != "linear") {
    pnde <- pure_NDE_contbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                             covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a, variance = variance)
    tnde <- total_NDE_contbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                              covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a, variance = variance)
  } else if (mreg == "linear" & yreg == "linear") {
    pnde <- pure_NDE_contcont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                              covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnde <- total_NDE_contcont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                               covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  }
  return(list(pnde = pnde, tnde = tnde))
}


NDE_delta_function <- function(thetas, treatment, mediator, m, interaction,
                               vecc, a_star, a, variance,
                               mreg, yreg) {
  if (mreg != "linear" & yreg != "linear") {
    pnded <- pure_NDE_binbin_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnded <- total_NDE_binbin_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg != "linear" & yreg == "linear") {
    pnded <- pure_NDE_bincont_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnded <- total_NDE_bincont_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg == "linear" & yreg != "linear") {
    pnded <- pure_NDE_contbin_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a, variance = variance)
    tnded <- total_NDE_contbin_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a, variance = variance)
  } else if (mreg == "linear" & yreg == "linear") {
    pnded <- pure_NDE_contcont_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnded <- total_NDE_contcont_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  }
  return(list(pnded = pnded, tnded = tnded))
}

NDE_boot <- function(coef=list(), vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment
  mediator <- coef$mediator
  covariates <- coef$covariates
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  return(NDE_boot_function(coef$betas, coef$thetas, treatment, mediator, covariates, vecc,
                    m, interaction, a_star, a, coef$variance, mreg, yreg))
}

#NDE_boot(coef,m=1,a_star=0,a=1)

NDE_delta <- function(coef=list(), vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment
  mediator <- coef$mediator
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  nde_delta <- NDE_delta_function(coef$thetas, treatment, mediator, m, interaction, vecc,
                                  a_star, a, coef$variance, mreg, yreg)

  se_pnde_delta <- msm::deltamethod(nde_delta$pnded, c(coef$thetas, coef$betas), coef$vcov_block)

  se_tnde_delta<- msm::deltamethod(nde_delta$tnded, c(coef$thetas, coef$betas), coef$vcov_block)

  return(list(se_pnde_delta=se_pnde_delta,se_tnde_delta=se_tnde_delta))
}

#NDE_delta(coef1,m=1,a_star=0,a=1)
