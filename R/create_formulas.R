create_formulas <- function(outcome, exposure, mediator, covariates.pre, covariates.post,
                            EMint, MMint, EMMint, EMint.terms, MMint.terms, EMMint.terms,
                            event, mreg, yreg, model) {

  if (length(mediator) == 1) {

    if (EMint) {
      int.terms <- paste(exposure, mediator, sep = "*")
      } else int.terms <- NULL

  } else if (length(mediator) > 1) {

    if (EMint == FALSE) EMint.terms <- NULL
    if (MMint == FALSE) MMint.terms <- NULL
    if (EMMint == FALSE) EMMint.terms <- NULL

    if (!is.null(c(EMint.terms, MMint.terms, EMMint.terms))) {
      if(!all(unlist(strsplit(c(EMint.terms, MMint.terms, EMMint.terms),
                              split = "[*]")) %in% c(exposure, mediator))){
        stop("Unsupported non exposure-mediator interaction")
      }
    }

    if(EMMint == FALSE) {

      if (EMint == TRUE && length(EMint.terms) == 0) {
      stop("Exposure-mediator interaction terms need to be specified")
    } else if (MMint == TRUE && length(MMint.terms) == 0) {
      stop("Mediator-mediator interaction terms need to be specified")
    } else { int.terms <- unique(c(EMint.terms, MMint.terms)) }

    } else if (EMMint == TRUE) {

      if (length(EMMint.terms) == 0) {
      stop("Exposure-mediator-mediator interaction terms need to be specified")
      } else if (length(EMMint.terms) > 0) {

      int.main.effect <- unlist(lapply(1:length(EMMint.terms),
                                       FUN = function(i)
                                         sapply(2:(length(unlist(strsplit(EMMint.terms[i], split = "[*]")))-1),
                                                FUN = function(x)
                                                  sapply(1:length(combn(unlist(strsplit(EMMint.terms[i], split = "[*]")),
                                                                        m = x,simplify = FALSE)),
                                                         FUN = function(y) paste(combn(unlist(strsplit(EMMint.terms[i], split = "[*]")),
                                                                                       m = x,simplify = FALSE)[[y]],collapse = "*")))))

      int.terms <- unique(c(EMint.terms, MMint.terms, EMMint.terms, int.main.effect))

    }
    }
    } else {int.terms <- NULL}

# forumulas for mediator regression

  if (!is.null(mreg)) {

    mediator.nomodel.index <- which(sapply(1:length(mreg), function(x) is.character(mreg[[x]])))

    mediator.nomodel <- mediator[mediator.nomodel.index]

    mreg.nomodel <- mreg[mediator.nomodel.index]

    if (!is.null(mediator.nomodel)) {
      mediator_formula <- lapply(1:length(mediator.nomodel),
           FUN = function(x) paste0(mediator.nomodel[x], "~", paste(c(exposure, covariates.pre),
                                                         collapse = "+")))

     # if (length(mediator_formula) == 1) mediator_formula <- mediator_formula[[1]]

      } else mediator_formula <- NULL

   } else mediator_formula <- NULL

# forumulas for outcome regression
if (is.character(yreg)) {

   if (model == "iorw") {

     outcome_formula <- paste0(outcome, "~", paste(c(exposure, covariates.pre),
                                                   collapse = "+"))

   } else {

     outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, covariates.pre),
                                                 collapse = "+"))

   }

  } else outcome_formula <- NULL


  if (yreg %in% c("coxph","aft_exp","aft_weibull")) {

    outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                             strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
  }

#######################################Weighting-based Approach##############################

  if (model == "wb") {

    exposure_formula <- paste0(exposure, "~", paste0(covariates.pre,
                                                       collapse = "+"))

    formulas <- list(exposure_formula = exposure_formula,
                     outcome_formula = outcome_formula)

  } else if (model == "iorw") {

#######################################Inverse OR Weighting Approach##############################

    exposure_formula <- paste0(exposure, "~", paste0(c(mediator, covariates.pre),
                                                collapse = "+"))

    formulas <- list(exposure_formula = exposure_formula,
                     outcome_formula = outcome_formula)

  } else if (model == "msm") {

###################################Marginal Structural Model#####################################

    wa_denom_formula <- paste0(exposure, "~", paste(covariates.pre, collapse = "+"))

    wz_nom_formula <- lapply(1:length(mediator),
           FUN = function(x) paste0(mediator[x], "~", exposure))

    wz_denom_formula <- lapply(1:length(mediator),
           FUN = function(x) paste0(mediator[x], "~",
                                    paste(c(exposure, mediator[0:(x-1)],
                                            covariates.pre, covariates.post),
                                                            collapse = "+")))

    cde_outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms),
                                                       collapse = "+"))

    formulas <- list(wa_denom_formula = wa_denom_formula,
                     wz_nom_formula = wz_nom_formula,
                     wz_denom_formula = wz_denom_formula,
                     cde_outcome_formula = cde_outcome_formula,
                     mediator_formula = mediator_formula,
                     outcome_formula = outcome_formula)

  } else if (model == "g-formula") {

########################################G-formula Approach#######################################

    if (!is.null(covariates.post)) {
        postcovar_formula <- lapply(1:length(covariates.post),
                                     FUN = function(x) paste0(covariates.post[x], "~",
                                                              paste0(c(exposure,
                                                                       covariates.pre),
                                                                     collapse = "+")))
      } else postcovar_formula <- NULL

      if (!is.null(mediator_formula)) {

        mediator_formula <- lapply(1:length(mediator_formula),
                               FUN = function(x) paste(c(mediator_formula[x], covariates.post), collapse = "+"))

      }

      if (length(mediator_formula) == 1) mediator_formula <- mediator_formula[[1]]

      if (!is.null(outcome_formula)) {

        outcome_formula <- paste0(c(outcome_formula, covariates.post), collapse = "+")

      }

    formulas <- list(postcovar_formula = postcovar_formula,
                     mediator_formula = mediator_formula,
                     outcome_formula = outcome_formula)

  } else {

    formulas <- list(outcome_formula = outcome_formula,
                     mediator_formula = mediator_formula)

  }

  return(formulas)

}
