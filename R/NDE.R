nde_function <- function(betas, thetas, variance, vcov_block, treatment, mediator,
                         covariates, vecc, interaction, m,  a_star, a, mreg, yreg) {

  covariatesTerm <- sum(betas[covariates]*t(vecc))

  interactionTerm <- ifelse(interaction, thetas[paste(treatment, mediator, sep = ":")], 0)

  if (mreg != "linear" & yreg != "linear") {

    ORpnde <- unname((exp(thetas[treatment] * (a - a_star)) * (1 + exp(thetas[mediator] +
                      interactionTerm * a + betas[1] + betas[treatment] * a_star +
                      covariatesTerm))) /(1 + exp(thetas[mediator] + interactionTerm * a_star +
                      betas[1] + betas[treatment] * a_star + covariatesTerm)))

    ORtnde <- unname((exp(thetas[treatment] * (a - a_star)) * (1 + exp(thetas[mediator] +
                      interactionTerm * a + betas[1] + betas[treatment] * a + covariatesTerm))) /
                     (1 + exp(thetas[mediator] + interactionTerm * a_star +  betas[1] +
                      betas[treatment] * a + covariatesTerm)))

    return(c(ORpnde = ORpnde, ORtnde = ORtnde))

  } else if (mreg != "linear" & yreg == "linear") {

    pnde <- unname(thetas[treatment] * (a - a_star) + interactionTerm*(a - a_star) *
                  (exp(betas[1] + betas[treatment] * a_star + covariatesTerm) /
                  (1 + exp(betas[1] + betas[treatment] * a_star + covariatesTerm))))


    tnde <- unname(thetas[treatment] * (a - a_star) + interactionTerm*(a - a_star) *
                  (exp(betas[1] + betas[treatment] * a + covariatesTerm) /
                  (1 + exp(betas[1] + betas[treatment] * a + covariatesTerm))))

    return(c(pnde = pnde, tnde = tnde))

    } else if (mreg == "linear" & yreg != "linear") {

    ORpnde <- unname(exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star +
                        covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                        0.5 * interactionTerm ^ 2 * variance * (a^2 - a_star ^ 2)))


    ORtnde <- unname(exp((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a +
                        covariatesTerm + thetas[mediator]  * variance)) * (a - a_star) +
                        0.5 * interactionTerm ^ 2 * variance * (a^2 - a_star ^ 2)))

    return(c(ORpnde = ORpnde, ORtnde = ORtnde))

    } else if (mreg == "linear" & yreg == "linear") {

    pnde <- unname((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a_star +
                    covariatesTerm))*(a - a_star))


    tnde <- unname((thetas[treatment] + interactionTerm * (betas[1] + betas[treatment] * a +
                                                             covariatesTerm))*(a - a_star))

    return(c(pnde = pnde, tnde = tnde))

  }
}


nde_se_delta <- function(thetas, betas, vcov_block, treatment, mediator, interaction,
                         m, vecc, a_star, a, variance, mreg, yreg) {

  j <- length(vecc)
  k <- length(thetas)

  if (j > 0) {
    for (i in 1:j) {
      assign(paste("vecc", i, sep = "_"), vecc[i])
    }
  }

  if (mreg != "linear" & yreg != "linear") {

    if (interaction) {

      pnde_formula <- paste0(" ~ exp(x2 * (a - a_star)) * (1+exp(x3 + x", k, "* a + x", k + 1," + x", k + 2,
                             "* a_star + ",
                   paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))/",
                   paste0("(1 + exp(x3 + x", k, "* a_star + x", k + 1," + x", k + 2, "* a_star + ",
                          paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))"))

      tnde_formula <- paste0(" ~ exp(x2 * (a - a_star)) * (1+exp(x3 + x", k, "* a + x", k + 1," + x", k + 2,
                             "* a + ",
                             paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))/",
                             paste0("(1 + exp(x3 + x", k, "* a_star + x", k + 1," + x", k + 2, "* a + ",
                                    paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))"))

    } else {

      pnde_formula <- paste0(" ~ exp(x2 * a) * (1+exp(x3 + x", k + 1," + x", k + 2, "* a_star + ",
                   paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))/",
                   paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "* a_star + ",
                          paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))"))

      tnde_formula <- paste0(" ~ exp(x2 * a) * (1+exp(x3 + x", k + 1," + x", k + 2, "* a + ",
                             paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))/",
                             paste0("(1 + exp(x3 + x", k + 1," + x", k + 2, "* a + ",
                                    paste0("x", k + 2 + 1:j, "*", "vecc_", 1:j, collapse = " + "), "))"))

    }

    } else if (mreg != "linear" & yreg == "linear") {

      if (interaction) {

        pnde_formula <- paste(paste0("~ x2 * (a-a_star) + x",k," * (a-a_star) * "),
                              paste0("exp(", paste0("x", k + 1:2, collapse = " + " ), " * a_star + "),
                              paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")/(1+",
                              paste0("exp(", paste0("x", k + 1:2, collapse = " + " ), " * a_star + "),
                              paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), "))")

        tnde_formula <- paste(paste0("~ x2 * (a-a_star) + x",k," * (a-a_star) * "),paste0("exp(", paste0("x", k + 1:2, collapse = " + " ), " * a + "),
                              paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), ")/(1+",
                              paste0("exp(", paste0("x", k + 1:2, collapse = " + " ), " * a + "),
                              paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), "))")

        } else { pnde_formula <- tnde_formula <- "~ x2 * (a-a_star)" }

      } else if (mreg == "linear" & yreg != "linear") {

        if (interaction) {

          pnde_formula <- paste0(" ~ exp((x2"," + x",k," * (x", k + 1," + x", k + 2,"* a_star", " + ",
                  paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), " + ",
                  "x3 * ", "variance",
                  ")",")) * (a-a_star)"," + 0.5 * x", k,"^2 * ", "variance ", "* (a^2 - a_star^2))")

          tnde_formula <- paste0(" ~ exp((x2"," + x",k," * (x", k + 1," + x", k + 2,"* a", " + ",
                  paste0("x", k + 2 + 1:j, " * ", "vecc_", 1:j, collapse = " + "), " + ",
                  "x3 * ", "variance",
                  ")",")) * (a-a_star)"," + 0.5 * x", k,"^2 * ", "variance ", "* (a^2 - a_star^2))")

          } else {

            pnde_formula <- tnde_formula <-  paste0(" ~ exp(x2 * (a-a_star)"," + 0.5 * x", k,"^2 * ",
                                                "variance ", "* (a^2 - a_star^2))")

      }

    } else if (mreg == "linear" & yreg == "linear") {

      if (interaction) {

        pnde_formula <- paste0(" ~ ", "(x2 + x", k, " * (x", k + 1, " + x", k + 2,"*a_star + ",
                        ifelse(j > 0,paste0("x", k + 2 + 1:j, "  * ", "vecc_", 1:j, collapse = " + "),
                               "0"),")) * (a - a_star)")

        tnde_formula <- paste0(" ~ ", "(x2 + x", k, " * (x", k + 1, " + x", k + 2,"*a + ",
                        ifelse(j > 0,paste0("x", k + 2 + 1:j, "  * ", "vecc_", 1:j, collapse = " + "),
                               "0"),")) * (a - a_star)")

      } else { pnde_formula <- tnde_formula <- paste0(" ~ x2 * (a - a_star)") }

      pnde_formula <- stringr::str_replace_all(pnde_formula,
                                             pattern = c("\\ba_star\\b" = as.character(a_star),
                                                         "\\ba\\b" = as.character(a)))

      tnde_formula <- stringr::str_replace_all(tnde_formula,
                                               pattern = c("\\ba_star\\b" = as.character(a_star),
                                                           "\\ba\\b" = as.character(a)))
    }


  for (i in 1:j) {

        pnde_formula <- stringr::str_replace_all(pnde_formula,
                                               paste("vecc", i, sep = "_"), as.character(vecc[i]))

        tnde_formula <- stringr::str_replace_all(tnde_formula,
                                                 paste("vecc", i, sep = "_"), as.character(vecc[i]))

      }

      pnde_formula <- as.formula(pnde_formula)

      tnde_formula <- as.formula(tnde_formula)

      pnde_se_delta <- msm::deltamethod(pnde_formula, c(thetas, betas), vcov_block)

      tnde_se_delta<- msm::deltamethod(tnde_formula, c(thetas, betas), vcov_block)

      return(c(pnde_se_delta = pnde_se_delta, tnde_se_delta = tnde_se_delta))

}

#nde_se_delta(thetas=coef1$thetas, betas=coef1$betas,vcov_block=coef1$vcov_block,treatment, mediator, interaction,
#            m, vecc, a_star, a, variance, mreg, yreg)
