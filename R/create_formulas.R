create_formulas <- function(outcome, exposure, mediator, covariates.pre, covariates.post,
                            EMint, MMint, EMMint, EMint.terms, MMint.terms, EMMint.terms,
                            event, mreg, yreg, model) {

  if (length(mediator) == 1) {

    if (EMint) {
      int.terms <- paste(exposure, mediator, sep = "*")
      } else int.terms <- NULL

  } else if (length(mediator) > 1) {

    if (!is.null(c(EMint.terms, MMint.terms, EMMint.terms))) {
      if(!all(unlist(strsplit(c(EMint.terms, MMint.terms, EMMint.terms),
                              split = "[*]")) %in% c(exposure, mediator))){
        stop("Variables other than exposure and mediators in interaction terms")
      }
    }

    if (EMint == TRUE && length(EMint.terms) == 0 && EMMint == FALSE) {
      stop("Exposure-mediator interaction terms need to be specified")
    } else if (MMint == TRUE && length(MMint.terms) == 0 && EMMint == FALSE) {
      stop("Mediator-mediator interaction terms need to be specified")
    } else if (EMMint == TRUE && length(EMMint.terms) == 0) {
      stop("Exposure-mediator-mediator interaction terms need to be specified")
    }

    if (EMMint == TRUE && length(EMMint.terms) > 0) {

      int.main.effect <- unlist(lapply(1:length(EMMint.terms),
                                       FUN = function(i)
                                         sapply(2:(length( unlist(strsplit(EMMint.terms[i], split = "[*]")))-1),
                                                FUN = function(x)
                                                  sapply(1:length(combn(unlist(strsplit(EMMint.terms[i], split = "[*]")),
                                                                        m = x,simplify = FALSE)),
                                                         FUN = function(y) paste(combn(unlist(strsplit(EMMint.terms[i], split = "[*]")),
                                                                                       m = x,simplify = FALSE)[[y]],collapse = "*")))))

      int.terms <- unique(c(EMint.terms, MMint.terms, EMMint.terms, int.main.effect))

      } else if (EMMint == FALSE) int.terms <- unique(c(EMint.terms, MMint.terms))
   }

# forumulas for mediator regression

  if (!is.null(mreg)) {

    mediator_formula <- lapply(1:length(mediator),
           FUN = function(x) paste0(mediator[x], "~", paste(c(exposure, covariates.pre),
                                                         collapse = "+")))

   } else mediator_formula <- NULL

# forumulas for outcome regression

   if (model == "iorw") {

     outcome_formula <- paste0(outcome, "~", paste(c(exposure, covariates.pre),
                                                   collapse = "+"))

   } else {

     outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, covariates.pre),
                                                 collapse = "+"))

   }



  if (yreg %in% c("coxph","aft_exp","aft_weibull")) {

    l <- strsplit(outcome_formula, split = " ~ ")

    l[[1]][1] <- paste0("Surv(", outcome, ", ", event, ")")

    outcome_formula <- paste(l[[1]][1], l[[1]][2], sep = " ~ ")
  }

  formulas <- list(outcome_formula = outcome_formula,
                   mediator_formula = mediator_formula)

###################################Marginal Structural Model#####################################

  if (model == "msm") {

    wa_denom_formula <- paste0(exposure, "~", covariates.pre)

    wz_nom_formula <- paste0(mediator, "~", exposure)

    wz_denom_formula <- paste0(mediator, "~", paste(c(exposure, covariates.pre, covariates.post),
                                                     collapse = "+"))

    cde_outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms),
                                                       collapse = "+"))

    formulas <- list(wa_denom_formula = wa_denom_formula,
                     wz_nom_formula = wz_nom_formula,
                     wz_denom_formula = wz_denom_formula,
                     cde_outcome_formula = cde_outcome_formula,
                     mediator_formula = mediator_formula,
                     outcome_formula = outcome_formula)

  }

#######################################Inverse OR Weighting Approach##############################

  if (model == "iorw") {

    exposure_formula <- paste0(exposure, "~", paste0(c(mediator, covariates.pre),
                                                collapse = "+"))

    formulas <- list(exposure_formula = exposure_formula,
                     outcome_formula = outcome_formula)



  }

########################################G-formula Approach#######################################

  if (model == "g-formula") {

    post_covar_formula <- lapply(1:length(covariates.post),
                                 FUN = function(x) paste0(covariates.post[x], "~",
                                    ifelse(x==1,paste0(c(exposure, covariates.pre),
                                                       collapse = "+"),
                                           paste0(c(covariates.post[1:(x-1)],
                                                    exposure, covariates.pre),
                                                  collapse = "+"))))

    mediator_formula <- lapply(1:length(mediator),
                                 FUN = function(x) paste0(mediator[x], "~",
                                                          ifelse(x==1,paste0(c(exposure, covariates.pre),
                                                                             collapse = "+"),
                                                                 paste0(c(mediator[1:(x-1)],
                                                                          exposure, covariates.pre),
                                                                        collapse = "+"))))

    outcome_formula <- paste0(c(outcome_formula, covariates.post), collapse = "+")

    formulas <- list(post_covar_formula = post_covar_formula,
                     mediator_formula = mediator_formula,
                     outcome_formula = outcome_formula)

  }

  return(formulas)

}
