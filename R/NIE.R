total_NIE_binbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                             a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates) {
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0, thetas[paste(treatment, mediator, sep = ":")])
  ORnie <- ((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
              (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a + covariatesTerm))) /
    ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) * (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a_star + covariatesTerm)))
  unname(ORnie)
}


pure_NIE_binbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                            a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates) {
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }

  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])
  ORnie <- ((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
              (1 + exp(thetas[mediator] + interactionTerm * a_star + betas[1] + betas[treatment] * a + covariatesTerm))) /
    ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) * (1+ exp(thetas[mediator] + interactionTerm * a_star + betas[1] + betas[treatment] * a_star + covariatesTerm)))
  unname(ORnie)
}


total_NIE_bincont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates){
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }
  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  nie <- (thetas[mediator]+interactionTerm*a) *
    (exp(betas[1] + betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) -
       exp(betas[1] + betas[treatment] * a_star + covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)))
  unname(nie)
}


pure_NIE_bincont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                             a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates){
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }
  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  nie <- (thetas[mediator]+interactionTerm*a_star) *
    (exp(betas[1] + betas[treatment] * a + covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) -
       exp(betas[1] + betas[treatment] * a_star + covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)))
  unname(nie)
}


total_NIE_contbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates) {
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }
  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  ORnie <- exp((thetas[mediator] * betas[treatment] + interactionTerm * betas[treatment] * a) * (a - a_star))
  unname(ORnie)
}


pure_NIE_contbin <- function(betas, thetas, treatment, mediator, covariates, vecc,
                             a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates) {
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }
  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  ORnie <- exp((thetas[mediator] * betas[treatment] + interactionTerm * betas[treatment] * a_star) * (a - a_star))
  unname(ORnie)
}


total_NIE_contcont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                               a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates) {
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }
  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  nie <- (thetas[mediator] * betas[treatment] + interactionTerm * betas[treatment] * a) * (a - a_star)
  unname(nie)
}


pure_NIE_contcont <- function(betas, thetas, treatment, mediator, covariates, vecc,
                              a_star, a, interaction) {
  covariatesTerm <- 0
  if (is.null(vecc)) {
    for (c in covariates) {
      covariatesTerm <- covariatesTerm + betas[c] * apply(df[c], 2, mean, na.rm = TRUE)
    }
  } else {
    for (i in 1:length(covariates)) {
      covariatesTerm <- covariatesTerm + betas[covariates[i]] * vecc[i]
    }
  }
  interactionTerm <- ifelse(is.na(thetas[paste(treatment, mediator, sep = ":")]),
                            0,
                            thetas[paste(treatment, mediator, sep = ":")])

  nie <- (thetas[mediator] * betas[treatment] + interactionTerm * betas[treatment] * a_star) * (a - a_star)
  unname(nie)
}


total_NIE_contcont_delta <- function(thetas, interaction,
                                     a_star, a){
  k <- length(thetas)

  F1 <- paste0("(x3 * x", k+2)
  F2 <- ")"
  F3 <- " * (a-a_star)"

  if(interaction){
    F2 <- paste0(" + x", k, " * x", k+2, " * a)")
  }

  f = paste0(" ~ ", F1, F2, F3)

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))
  return(as.formula(s))
}


pure_NIE_contcont_delta <- function(thetas, interaction,
                                    a_star, a){
  k <- length(thetas)

  F1 <- paste0("(x3 * x", k+2)
  F2 <- ")"
  F3 <- " * (a - a_star)"

  if(interaction){
    F2 <- paste0(" + x", k, " * x", k+2, " * a_star)")
  }

  f = paste0(" ~ ", F1, F2, F3)

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))
  return(as.formula(s))
}


total_NIE_contbin_delta <- function(thetas, interaction,
                                    a_star, a){
  k <- length(thetas)

  F1 <- paste0("exp((x3 * x", k+2)
  F2 <- ")"
  F3 <- " * (a-a_star))"

  if(interaction){
    F2 <- paste0(" + x", k, " * x", k+2, " * a)")
  }

  f = paste0(" ~ ", F1, F2, F3)

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))
  return(as.formula(s))
}


pure_NIE_contbin_delta <- function(thetas, interaction,
                                   a_star, a){
  k <- length(thetas)

  F1 <- paste0("exp((x3 * x", k+2)
  F2 <- ")"
  F3 <- " * (a-a_star))"

  if(interaction){
    F2 <- paste0(" + x", k, " * x", k+2, " * a_star)")
  }

  f = paste0(" ~ ", F1, F2, F3)

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))
  return(as.formula(s))
}


total_NIE_bincont_delta <- function(thetas, vecc, interaction,
                                    a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j){
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  F1 <- paste0("(x3 + x", k, "*a)")

  # First numerator
  N1 <- paste0("x", k+1:2, collapse = " + " )
  N1 <- paste0("exp(", N1, "*a + ")
  N1 <- ifelse(j > 0, paste(N1, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), N1)

  # Second numerator
  N2 <- paste0("x", k+1:2, collapse = " + " )
  N2 <- paste0("exp(", N2, "*a_star + ")
  N2 <- ifelse(j > 0, paste(N2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), N2)

  # First denominator
  D1 <- paste0("x", k+1:2, collapse = " + " )
  D1 <- paste0("1 + exp(", D1, "*a + ")
  D1 <- ifelse(j > 0, paste(D1, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), D1)
  D1 <- paste0("(", D1, ")")

  # Second denominator
  D2 <- paste0("x", k+1:2, collapse = " + " )
  D2 <- paste0("1 + exp(", D2, "*a_star + ")
  D2 <- ifelse(j > 0, paste(D2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), D2)
  D2 <- paste0("(", D2, ")")

  # Construct formula
  f <- paste0(" ~ ", F1, " * ( (", N1, "/", D1, ") - (", N2, "/", D2, ") )")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j){
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}


pure_NIE_bincont_delta <- function(thetas, vecc, interaction,
                                   a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j){
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  F1 <- paste0("(x3 + x", k, "*a_star)")

  # First numerator
  N1 <- paste0("x", k+1:2, collapse = " + " )
  N1 <- paste0("exp(", N1, "*a + ")
  N1 <- ifelse(j > 0,paste(N1, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), N1)

  # Second numerator
  N2 <- paste0("x", k+1:2, collapse = " + " )
  N2 <- paste0("exp(", N2, "*a_star + ")
  N2 <- ifelse(j > 0,paste(N2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), N2)

  # First denominator
  D1 <- paste0("x", k+1:2, collapse = " + " )
  D1 <- paste0("1 + exp(", D1, "*a + ")
  D1 <- ifelse(j > 0,paste(D1, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), D1)
  D1 <- paste0("(", D1, ")")

  # Second denominator
  D2 <- paste0("x", k+1:2, collapse = " + " )
  D2 <- paste0("1 + exp(", D2, "*a_star + ")
  D2 <- ifelse(j > 0, paste(D2, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")"), D2)
  D2 <- paste0("(", D2, ")")

  # Construct formula
  f <- paste0(" ~ ", F1, " * ( (", N1, "/", D1, ") - (", N2, "/", D2, ") )")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j){
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}


total_NIE_binbin_delta <- function(thetas, vecc, interaction,
                                   a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j){
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  X1 <- paste0("x", k + 1:2, collapse = " + ")
  XC <- ifelse(j  > 0, paste0("x", k + 2 + 1:j, "  * ", "vecc_", 1:j, collapse = " + "), "0")
  s1 <- paste0("(1 + exp(", X1, " * a_star + ", XC, "))")
  s3 <- paste0("(1 + exp(", X1, " * a + ", XC, "))")

  if (interaction) {
    X2 <- paste0("x3 +x",k)
    s2 <- paste0("(1 + exp(", X2, " * a + ", X1, " * a +", XC, "))")
    s4 <- paste0("(1 + exp(", X2, " * a + ", X1, " * a_star +", XC, "))")
  } else {
    X2 <- paste0("x3")
    s2 <- paste0("(1 + exp(", X2, " + ", X1, " * a +", XC, "))")
    s4 <- paste0("(1 + exp(", X2, " + ", X1, " * a_star +", XC, "))")
  }

  f <- paste0(" ~ ", "(", s1, "*", s2, ")/(", s3, "*", s4, ")")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j){
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}


pure_NIE_binbin_delta <- function(thetas, vecc, interaction,
                                  a_star, a) {
  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j){
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  X1 <- paste0("x", k + 1:2, collapse = " + ")
  XC <- ifelse(j > 0, paste0("x", k + 2 + 1:j, "  * ", "vecc_", 1:j, collapse = " + "), "0")
  s1 <- paste0("(1 + exp(", X1, " * a_star + ", XC, "))")
  s3 <- paste0("(1 + exp(", X1, " * a + ", XC, "))")

  if (interaction) {
    X2 <- paste0("x3 + x",k)
    s2 <- paste0("(1 + exp(", X2, " * a_star + ", X1, " * a +", XC, "))")
    s4 <- paste0("(1 + exp(", X2, " * a_star + ", X1, " * a_star +", XC, "))")
  } else {
    X2 <- paste0("x3")
    s2 <- paste0("(1 + exp(", X2, " + ", X1, " * a +", XC, "))")
    s4 <- paste0("(1 + exp(", X2, " + ", X1, " * a_star +", XC, "))")
  }

  f <- paste0(" ~ ", "(", s1, "*", s2, ")/(", s3, "*", s4, ")")

  s <- stringr::str_replace_all(f, pattern = c("\\ba_star\\b" = as.character(a_star),
                                               "\\ba\\b" = as.character(a)))

  if (j > 0) {
    for (i in 1:j) {
      s <- stringr::str_replace_all(s, paste("vecc", i, sep = "_"), as.character(vecc[i]))
    }
  }

  return(as.formula(s))
}


NIE_boot_function <- function(betas, thetas, treatment, mediator, covariates, vecc = vecc,
                              m, interaction, a_star, a,
                              mreg = "linear", yreg = "linear") {
  if (mreg != "linear" & yreg != "linear") {
    pnie <- pure_NIE_binbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                            covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnie <- total_NIE_binbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                             covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg != "linear" & yreg == "linear") {
    pnie <- pure_NIE_bincont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                             covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnie <- total_NIE_bincont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                              covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg == "linear" & yreg != "linear") {
    pnie <- pure_NIE_contbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                             covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnie <- total_NIE_contbin(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                              covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg == "linear" & yreg == "linear") {
    pnie <- pure_NIE_contcont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                              covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnie <- total_NIE_contcont(betas = betas, thetas = thetas, treatment = treatment, mediator = mediator,
                               covariates = covariates, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  }
  return(list(pnie = pnie, tnie = tnie))
}


NIE_delta_function <- function(thetas, treatment, mediator, m, vecc, interaction = TRUE,
                               a_star, a,
                               mreg = "linear", yreg = "linear") {
  if (mreg != "linear" & yreg != "linear") {
    pnied <- pure_NIE_binbin_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnied <- total_NIE_binbin_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg != "linear" & yreg == "linear") {
    pnied <- pure_NIE_bincont_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
    tnied <- total_NIE_bincont_delta(thetas = thetas, vecc = vecc, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg == "linear" & yreg != "linear") {
    pnied <- pure_NIE_contbin_delta(thetas = thetas, interaction = interaction, a_star = a_star, a = a)
    tnied <- total_NIE_contbin_delta(thetas = thetas, interaction = interaction, a_star = a_star, a = a)
  } else if (mreg == "linear" & yreg == "linear") {
    pnied <- pure_NIE_contcont_delta(thetas = thetas, interaction = interaction, a_star = a_star, a = a)
    tnied <- total_NIE_contcont_delta(thetas = thetas, interaction = interaction, a_star = a_star, a = a)
  }
  return(list(pnied = pnied, tnied = tnied))
}


NIE_boot <- function(coef, vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment
  mediator <- coef$mediator
  covariates <- coef$covariates
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  NIE_boot_function(coef$betas, coef$thetas, treatment, mediator, covariates, vecc,
                    m, interaction, a_star, a, mreg, yreg)

}

#NIE_boot(coef, m = 1, a_star = 2, a = 3)

NIE_delta <- function(coef = list(), vecc = c(), m = NULL, a_star = NULL, a = NULL) {

  treatment <- coef$treatment
  mediator <- coef$mediator
  interaction <- coef$interaction
  yreg <- coef$outcome_reg
  mreg <- coef$mediator_reg

  nie_delta <- NIE_delta_function(coef$thetas, treatment, mediator, m, vecc, interaction,
                                  a_star, a, mreg, yreg)

  se_pnie_delta <- msm::deltamethod(nie_delta$pnied, c(coef$thetas, coef$betas),coef$vcov_block)

  se_tnie_delta<- msm::deltamethod(nie_delta$tnied, c(coef$thetas, coef$betas), coef$vcov_block)

  list(se_pnie_delta = se_pnie_delta, se_tnie_delta = se_tnie_delta)
}

#NIE_delta(coef, m = 1, a_star = 2, a = 3)
