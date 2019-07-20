cde_function <- function(thetas, treatment, mediator, interaction, m_star, a_star, a, yreg) {

  if (yreg != "linear")
    cde <- unname(exp(thetas[treatment] + ifelse(interaction,
           thetas[paste(treatment, mediator, sep = ":")] * m_star, 0) * (a - a_star)))

  else if (yreg == "linear")
    cde <- unname((thetas[treatment] + ifelse(interaction,
           thetas[paste(treatment, mediator, sep = ":")] * m_star, 0)) * (a - a_star))

  return(cde)

}


cde_se_delta <- function(thetas, vcov_thetas, treatment, mediator, m_star, a_star, a, interaction,
                         yreg ) {

  if (yreg != "linear")
    cde_formula <- as.formula(stringr::str_replace_all(
      ifelse(interaction,
             paste0("~ exp((x2 + x", length(thetas), " * m_star) * (a - a_star))"),
             paste0(" ~exp(x2 * (a - a_star))")),
      pattern = c("\\ba_star\\b" = as.character(a_star),
                  "\\ba\\b" = as.character(a),
                  "\\bm_star\\b" = as.character(m_star))))

  else if (yreg == "linear")
    cde_formula <- as.formula(stringr::str_replace_all(
      ifelse(interaction,
             paste0("~ (x2 + x", length(thetas), " * m_star) * (a - a_star)"),
             paste0(" ~ x2 * (a - a_star)")),
      pattern = c("\\ba_star\\b" = as.character(a_star),
                  "\\ba\\b" = as.character(a),
                  "\\bm_star\\b" = as.character(m_star))))

  cde_se_delta <- msm::deltamethod(cde_formula, thetas, vcov_thetas)

  return(cde_se_delta)

}

cde_se_delta(thetas=coef1$thetas, vcov_thetas=coef1$vcov_thetas,treatment, mediator, interaction,
             m_star, a_star, a, yreg)

cde_function(thetas=coef1$thetas, treatment, mediator, interaction,
             m_star, a_star, a, yreg)
