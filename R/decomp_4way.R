decomp_4way <- function(data, outcome, treatment, mediator, covariates, vecc,
                        interaction, event, mreg, yreg, m_star, a_star, a,
                        type = c("delta", "bootstrap"), nboot = 100, conf = 0.95) {

  if (type == "delta") {
    warning("Delta method not finished")
  } else if (type == "bootstrap") {
    bootstrap_4way(data = data, outcome = outcome, treatment = treatment, mediator = mediator,
                         covariates = covariates, vecc = vecc, interaction = interaction, event = event,
                         mreg = mreg, yreg = yreg, nboot = nboot, conf = conf,
                         m = m, a_star = a_star, a = a)
  }
}

