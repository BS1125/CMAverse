nie_function <- function(betas, thetas, treatment, mediator, covariates, vecc = vecc,
                              m, interaction, a_star, a, mreg, yreg) {

  covariatesTerm <- ifelse(!is.null(vecc), sum(betas[covariates]*t(vecc)), 0)

  interactionTerm <- ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")], 0)

  if (mreg != "linear" & yreg != "linear") {

    ORpnie <- unname(((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
               (1 + exp(thetas[mediator] + interactionTerm * a_star + betas[1] +
                          betas[treatment] * a + covariatesTerm))) /
      ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) *
         (1+ exp(thetas[mediator] + interactionTerm * a_star + betas[1] + betas[treatment] * a_star +
                   covariatesTerm))))

    ORtnie <- unname(((1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm)) *
               (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a +
                covariatesTerm))) / ((1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) *
               (1 + exp(thetas[mediator] + interactionTerm * a + betas[1] + betas[treatment] * a_star +
                covariatesTerm))))

    return(c(ORpnie = ORpnie, ORtnie = ORtnie))

    } else if (mreg != "linear" & yreg == "linear") {

      pnie <- unname((thetas[mediator]+interactionTerm*a_start) * (exp(betas[1] + betas[treatment] * a +
             covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) -
             exp(betas[1] + betas[treatment] * a_star + covariatesTerm) / (1 + exp(betas[1] +
             betas[treatment] * a_star + covariatesTerm))))

      tnie <- unname((thetas[mediator]+interactionTerm*a) * (exp(betas[1] + betas[treatment] * a +
             covariatesTerm) / (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm)) -
             exp(betas[1] + betas[treatment] * a_star + covariatesTerm) / (1 + exp(betas[1] +
             betas[treatment] * a_star + covariatesTerm))))

      return(c(pnie = pnie, tnie = tnie))

    } else if (mreg == "linear" & yreg != "linear") {

      ORpnie <- unname(exp((thetas[mediator] * betas[treatment] +
                     interactionTerm * betas[treatment] * a_star) * (a - a_star)))

      ORtnie <- unname(exp((thetas[mediator] * betas[treatment] +
                     interactionTerm * betas[treatment] * a) * (a - a_star)))

      return(c(ORpnie = ORpnie, ORtnie = ORtnie))

    } else if (mreg == "linear" & yreg == "linear") {

      pnie <- unname((thetas[mediator] * betas[treatment] +
                 interactionTerm * betas[treatment] * a_star) * (a - a_star))

      tnie <- unname((thetas[mediator] * betas[treatment] +
                 interactionTerm * betas[treatment] * a) * (a - a_star))

      return(c(pnie = pnie, tnie = tnie))

      }
  }


nie_se_delta <- function(thetas, betas, vcov_block, treatment, mediator, m, vecc, interaction,
                         a_star, a, mreg, yreg) {

  k <- length(thetas)

  j <- length(vecc)

  if (j > 0) {
    for (i in 1:j){
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  XC <- ifelse(j  > 0, paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), "0")

  if (mreg != "linear" & yreg != "linear") {

    if (interaction) {

      pnie_formula <- paste0(" ~ (1 + exp(", paste0("x", k + 1:2, collapse = " + "), " * a_star + ",
       XC, "))/(1 + exp(", paste0("x", k + 1:2, collapse = " + "), " * a + ", XC, "))*",
       "(1 + exp(x3 +x",k, " * a_star + ", paste0("x", k + 1:2, collapse = " + "), " * a +", XC, "))/",
       "(1 + exp(x3 +x",k, " * a_star + ", paste0("x", k + 1:2, collapse = " + "), " * a_star +", XC, "))")

      tnie_formula <- paste0(" ~ (1 + exp(", paste0("x", k + 1:2, collapse = " + "), " * a_star + ",
       XC, "))/(1 + exp(", paste0("x", k + 1:2, collapse = " + "), " * a + ", XC, "))*",
       "(1 + exp(x3 +x",k, " * a + ", paste0("x", k + 1:2, collapse = " + "), " * a +", XC, "))/",
       "(1 + exp(x3 +x",k, " * a + ", paste0("x", k + 1:2, collapse = " + "), " * a_star +", XC, "))")

    } else { pnie_formula <- tnie_formula <- paste0(" ~ (1 + exp(", paste0("x", k + 1:2, collapse = " + "),
                " * a_star + ", XC, "))/(1 + exp(", paste0("x", k + 1:2, collapse = " + "),
                " * a + ", XC, "))*", "(1 + exp(x3 + ", paste0("x", k + 1:2, collapse = " + "),
                " * a +", XC, "))/", "(1 + exp(x3 + ", paste0("x", k + 1:2, collapse = " + "),
                " * a_star +", XC, "))")
      }

    } else if (mreg != "linear" & yreg == "linear") {

      if (interaction) {

        pnie_formula <- paste0(" ~ (x3 + x", k, " * a_star ) * ( ",
         paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a + "), XC,")"),
         "/( 1 + ", paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a + "), XC,")) - ",
         paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a_star + "), XC,")"),
         "/( 1 + ", paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a_star + "),
                          XC,")))")))

        tnie_formula <- paste0(" ~ (x3 + x", k, " * a ) * ( ",
         paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a + "), XC,")"),
         "/( 1 + ", paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a + "), XC,")) - ",
         paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a_star + "), XC,")"),
         "/( 1 + ", paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a_star + "),
                          XC,")))")))

      } else {

        pnie_formula <- tnie_formula <-  paste0(" ~ x3 * ( ",
         paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a + "), XC,")"),
         "/( 1 + ", paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a + "), XC,")) - ",
         paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a_star + "), XC,")"),
         "/( 1 + ", paste(paste0("exp( ", paste0("x", k+1:2, collapse = " + " ), " * a_star + "),
         XC,")))")))

        }

    } else if (mreg == "linear" & yreg != "linear") {

      if (interaction) {

        pnie_formula <- paste0(" ~ exp((x3 * x", k+2, " + x", k, " * x", k+2, " * a_star) * (a-a_star))")

        tnie_formula <- paste0(" ~ exp((x3 * x", k+2, " + x", k, " * x", k+2, " * a) * (a-a_star))")

      } else { pnie_formula <- tnie_formula <-
        paste0(" ~ exp((x3 * x", k+2, " + x", k, " * x", k+2, " * a) * (a-a_star))") }

    } else if (mreg == "linear" & yreg == "linear") {


      if (interaction) {

        pnie_formula <- paste0(" ~ (x3 * x", k+2, " + x", k, " * x", k+2, " * a_star) * (a-a_star)")

        tnie_formula <- paste0(" ~ (x3 * x", k+2, " + x", k, " * x", k+2, " * a) * (a-a_star)")

      } else { pnie_formula <- tnie_formula <- paste0(" ~ x3 * x", k+2, "* (a-a_star)") }

   }


  pnie_formula <- stringr::str_replace_all(pnie_formula,
                                           pattern = c("\\ba_star\\b" = as.character(a_star),
                                                       "\\ba\\b" = as.character(a)))

  tnie_formula <- stringr::str_replace_all(tnie_formula,
                                           pattern = c("\\ba_star\\b" = as.character(a_star),
                                                       "\\ba\\b" = as.character(a)))

  if (j>0) {
    for (i in 1:j) {

      pnie_formula <- stringr::str_replace_all(pnie_formula,
                                               paste("vecc", i, sep = "_"), as.character(vecc[i]))

      tnie_formula <- stringr::str_replace_all(tnie_formula,
                                               paste("vecc", i, sep = "_"), as.character(vecc[i]))

    }
  }

  pnie_formula <- as.formula(pnie_formula)

  tnie_formula <- as.formula(tnie_formula)

  pnie_se_delta <- msm::deltamethod(pnie_formula, c(thetas, betas), vcov_block)

  tnie_se_delta <- msm::deltamethod(tnie_formula, c(thetas, betas), vcov_block)

  return(c(pnie_se_delta = pnie_se_delta, tnie_se_delta = tnie_se_delta))

}

nie_se_delta(thetas=coef1$thetas, betas=coef1$betas,vcov_block=coef1$vcov_block,treatment, mediator, interaction,
           m, vecc, a_star, a, mreg, yreg)

