create_formulas <- function(data, model,
                            outcome, event, exposure, mediator, EMint, prec, postc,
                            yreg, mreg, ereg, postcreg, wmreg) {

  ########################################Outcome Regression Formula################################

  if (is.character(yreg)) {

    if (EMint) {

      int.terms <- paste(exposure, mediator, sep = "*")

    } else {int.terms <- NULL}


    if (model == "iorw") {

      outcome_formula <- paste0(outcome, "~", paste(c(exposure, prec),
                                                    collapse = "+"))

    } else if (model %in% c("rb", "wb", "ne", "msm")) {

      outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, prec),
                                                    collapse = "+"))

    } else if (model == "g-formula") {

      outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, prec,
                                                      postc),
                                                    collapse = "+"))

    }


    if (yreg %in% c("coxph","aft_exp","aft_weibull")) {

      outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                               strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")

    }

  } else outcome_formula <- NULL

  ########################################Mediator Regression Formula################################

  if (model %in% c("rb", "msm")) {

    mediator_formula <- lapply(1:length(mediator),
                               FUN = function(x) {
                                 if (is.character(mreg[[x]])) {
                                   paste0(mediator[x], "~",
                                          paste(c(exposure, prec),
                                                collapse = "+"))
                                 } else NULL})

  } else if (model == "g-formula") {

    mediator_formula <- lapply(1:length(mediator),
                               FUN = function(x) {
                                 if (is.character(mreg[[x]])) {
                                   paste0(mediator[x], "~",
                                          paste(c(exposure, prec, postc),
                                                collapse = "+"))
                                 } else NULL})

  } else if (model %in% c("wb", "iorw", "ne")) {

    mediator_formula <- NULL

  }

  ###############################Exposure Regression Formula For Weights#############################

  if (!is.null(ereg) && is.character(ereg)) {

    if (model %in% c("wb", "msm")) {

      exposure_formula <- paste0(exposure, "~", paste0(prec,
                                                       collapse = "+"))

    } else if (model == "iorw") {

      exposure_formula <- paste0(exposure, "~", paste0(c(mediator, prec),
                                                       collapse = "+"))

    } else exposure_formula <- NULL

  } else exposure_formula <- NULL

  ###################################Mediator Regression Formula For Weights#########################

  if (model == "msm") {

    wm_nom_formula <- lapply(1:length(mediator),
                             FUN = function(x) paste0(mediator[x], "~",exposure))

    wm_denom_formula <- lapply(1:length(mediator),
                               FUN = function(x) {
                                 if (is.character(wmreg[[x]])) {
                                   paste0(mediator[x], "~",
                                          paste(c(exposure, prec, postc),
                                                collapse = "+"))
                                 } else NULL})

  } else  wm_nom_formula <- wm_denom_formula <- NULL

  ###################################Post-exposure Confounder Regression Formula#####################

  if (!is.null(postc)) {

    postc_formula <- lapply(1:length(postc),
                            FUN = function(x) {
                              if (is.character(postcreg[[x]])) {
                                paste0(postc[x], "~",
                                       paste(c(exposure, prec),
                                             collapse = "+"))
                              } else NULL})

  } else {postc_formula <- NULL}

  #######################################Return Formula List#########################################

  if (model == "rb") {

    formulas <- list(outcome_formula = outcome_formula,
                     mediator_formula = mediator_formula)

  } else if (model %in% c("wb", "iorw")) {

    formulas <- list(outcome_formula = outcome_formula,
                     exposure_formula = exposure_formula)

  } else if (model == "ne") {

    formulas <- list(outcome_formula = outcome_formula)

  } else if (model == "msm") {

    formulas <- list(outcome_formula = outcome_formula,
                     mediator_formula = mediator_formula,
                     exposure_formula = exposure_formula,
                     wm_nom_formula = wm_nom_formula,
                     wm_denom_formula = wm_denom_formula)

  } else if (model == "g-formula") {

    formulas <- list(outcome_formula = outcome_formula,
                     mediator_formula = mediator_formula,
                     postc_formula = postc_formula)
  }

  return(formulas)

}
