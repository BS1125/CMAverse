medflex_weight_categorical <- function(formulas = NULL,data = NULL) {

  treatment <- formulas$treatment
  mediator <- formulas$mediator
  outcome_formula <- formulas$outcome_formula
  mediator_formula <- formulas$mediator_formula

  s <- gsub(pattern = treatment,
            replacement = paste0(treatment, "0"),
            x=outcome_formula)
  assign("tmp_data", data, envir = globalenv())
  medflex_formula <- gsub(pattern = mediator, replacement =  paste0(treatment, "1"), s)
  assign("medflex_formula", medflex_formula, envir = globalenv())
  medflex_data <- medflex::neWeight(as.formula(mediator_formula), data = tmp_data)
  assign("medflex_data", medflex_data, envir = globalenv())

  result <- medflex::neModel(as.formula(medflex_formula), expData = medflex_data, se = "robust")

  rm(medflex_formula, tmp_data, medflex_data, envir = globalenv())

  return(summary(result))
}

#medflex_weight_categorical(formulas = formulas, data = data)

medflex_weight_continuous = function(formulas = NULL,data = NULL) {
  treatment <- formulas$treatment
  mediator <- formulas$mediator
  outcome_formula <- formulas$outcome_formula
  mediator_formula <- formulas$mediator_formula

  s <- gsub(pattern = treatment,
            replacement = paste0(treatment, "0"),
            outcome_formula)
  assign("tmp_data", data, envir = globalenv())
  medflex_formula <- gsub(pattern = mediator, replacement =  paste0(treatment, "1"), s)
  assign("medflex_formula", medflex_formula, envir = globalenv())
  medflex_data <- medflex::neWeight(as.formula(mediator_formula), data = tmp_data)
  assign("medflex_data", medflex_data, envir = globalenv())

  result <- medflex::neModel(as.formula(medflex_formula), expData = medflex_data, se = "robust")

  rm(medflex_formula, tmp_data, medflex_data, envir = globalenv())

  return(summary(result))
}

#medflex_weight_continuous(formulas = formulas, data = data)

medflex_imputation_categorical <- function(formulas = NULL,data = NULL) {

  treatment <- formulas$treatment
  mediator <- formulas$mediator
  outcome_formula <- formulas$outcome_formula
  mediator_formula <- formulas$mediator_formula

  s <- gsub(pattern = treatment,
            replacement = paste0(treatment, "0"),
            outcome_formula)
  medflex_formula <- gsub(pattern = mediator, replacement =  paste0(treatment, "1"), s)
  assign("medflex_formula", medflex_formula, envir = globalenv())

  s <- gsub(pattern = treatment,
            replacement = paste0("factor(", treatment, ")"),
            outcome_formula)
  assign("tmp_data", data, envir = globalenv())
  medflex_data <- medflex::neImpute(as.formula(s), data = tmp_data)
  assign("medflex_data", medflex_data, envir = globalenv())

  result <- medflex::neModel(as.formula(medflex_formula), expData = medflex_data, se = "robust")

  rm(medflex_formula, tmp_data, medflex_data, envir = globalenv())

  return(summary(result))
}

#medflex_imputation_categorical(formulas = formulas, data = data)

medflex_imputation_continuous <- function(formulas = NULL,data = NULL) {

  treatment <- formulas$treatment
  mediator <- formulas$mediator
  outcome_formula <- formulas$outcome_formula
  mediator_formula <- formulas$mediator_formula

  s <- gsub(pattern = treatment,
            replacement = paste0(treatment, "0"),
            outcome_formula)
  medflex_formula <- gsub(pattern = mediator, replacement =  paste0(treatment, "1"), s)
  assign("medflex_formula", medflex_formula, envir = globalenv())

  s <- gsub(pattern = treatment,
            replacement = paste0("factor(", treatment, ")"),
            outcome_formula)
  assign("tmp_data", data, envir = globalenv())
  medflex_data <- medflex::neImpute(as.formula(s), data = tmp_data)
  assign("medflex_data", medflex_data, envir = globalenv())

  result <- medflex::neModel(as.formula(medflex_formula), expData = medflex_data, se = "robust")

  rm(medflex_formula, tmp_data, medflex_data, envir = globalenv())

  return(summary(result))
}

#medflex_imputation_continuous(formulas = formulas, data = data)
