#' Causal Mediation Analysis with the Regression-Based Approach
#'
#' \code{cmest_rb} is used to implement the \emph{the regression-based approach} 
#' by Valeri & VanderWeele (2013) and VanderWeele & Vansteelandt (2014) for causal mediation analysis 
#' with a single exposure, a single outcome, and a single or multiple mediators.
#'
#' @param data a data frame (or object coercible by \link{as.data.frame} to a data frame) 
#' containing the variables in the model. 
#' @param outcome variable name of the outcome.
#' @param event (required when \code{yreg} is \code{coxph}, \code{aft_exp},
#' or \code{aft_weibull}) variable name of the event.
#' @param exposure variable name of the exposure.
#' @param mediator a vector of variable name(s) of mediator(s).
#' @param EMint a logical value. \code{TRUE} indicates there is 
#' exposure-mediator interaction in \code{yreg}. 
#' @param basec a vector of variable names of confounders. See \code{Details}.
#' @param yreg outcome regression model. See \code{Details}.
#' @param mreg a list of mediator regression models following the order in \code{mediator}. See \code{Details}.
#' @param estimation estimation method. \code{paramfunc} and 
#' \code{imputation} are implemented (the first 4 letters are sufficient). Default is \code{imputation}. 
#' See \code{Details}. 
#' @param inference inference method. \code{delta} and 
#' \code{bootstrap} are implemented (the first 4 letters are sufficient). Default is \code{bootstrap}. 
#' See \code{Details}.
#' @param astar the control value of the exposure. 
#' @param a the treatment value of the exposure. 
#' @param mval a list of values at which each mediator is controlled to calculate the \code{cde}, following the order in \code{mediator}.
#' @param basecval (required when \code{estimation} is \code{paramfunc} and \code{EMint} is \code{TRUE}) 
#' a list of values at which each confounder is conditioned on, following the order in \code{basec}. 
#' If \code{NULL}, the mean of each confounder is used.
#' @param yval (required when the outcome is categorical) the level of the outcome at which causal effects on the 
#' risk ratio scale are estimated. If \code{NULL}, the last level is used.
#' @param nboot (used when \code{inference} is \code{bootstrap}) the number of bootstraps applied. 
#' Default is 200.
#' @param boot.ci.type (used when \code{inference} is \code{bootstrap}) the type of bootstrap confidence interval. If \code{per}, percentile bootstrap
#' confidence intervals are estimated; if \code{bca}, bias-corrected and accelerated (BCa) bootstrap 
#' confidence intervals are estimated. Default is \code{per}.
#' @param casecontrol a logical value. \code{TRUE} indicates a case control study in which the
#' first level of the outcome is treated as the control and the second level of the outcome is 
#' treated as the case. Default is \code{FALSE}.
#' @param yrare (used when \code{casecontrol} is \code{TRUE}) a logical value. \code{TRUE} 
#' indicates the case is rare.
#' @param yprevalence (used when \code{casecontrol} is \code{TRUE}) the prevalence of the case.
#' @param multimp a logical value. If 
#' \code{TRUE}, conduct multiple imputations using the \link[mice]{mice} function. Default is 
#' \code{FALSE}.
#' @param args_mice a list of additional arguments passed to the \link[mice]{mice} function. See \link[mice]{mice}
#' for details.
#' @param x an object of class \code{cmest}.
#' @param object an object of class \code{cmest}.
#' @param digits minimal number of significant digits. See \link{print.default}.
#' 
#' 
#' 
#' @details
#' 
#' \strong{Assumptions of the regression-based approach}
#' \itemize{
#' \item{\emph{There is no unmeasured exposure-outcome confounding:} }{given \code{basec} and 
#' \code{postc}, \code{exposure} is independent of \code{outcome}.}
#' \item{\emph{There is no unmeasured mediator-outcome confounding:} }{given \code{exposure} and 
#' \code{basec}, \code{mediator} is independent of \code{outcome}.}
#' \item{\emph{There is no unmeasured exposure-mediator confounding:} }{given \code{basec}, 
#' \code{exposure} is independent of \code{mediator}.}
#' \item{\emph{There is no mediator-outcome confounder affected by the exposure:} }{there is no 
#' variable in \code{basec} affected by \code{exposure}.}
#' }
#' 
#' 
#' \strong{Regression models}
#' 
#' Each regression model in \code{yreg} and \code{mreg} can be specified by a fitted regression 
#' object or the character name of a regression model. 
#' 
#' \emph{The Character Name of a Regression Model:}
#' \itemize{
#' \item{\code{linear}: }{linear regression fitted by \link{glm} with \code{family = gaussian()}}
#' \item{\code{logistic}: }{logistic regression fitted by \link{glm} with \code{family = logit()}}
#' \item{\code{loglinear}: }{loglinear regression fitted by \link{glm} with 
#' \code{family = poisson()}}
#' \item{\code{poisson}: }{poisson regression fitted by \link{glm} with 
#' \code{family = poisson()}}
#' \item{\code{quasipoisson}: }{quasipoisson regression fitted by \link{glm} with 
#' \code{family = quasipoisson()}}
#' \item{\code{negbin}: }{negative binomial regression fitted by \link[MASS]{glm.nb}}
#' \item{\code{multinomial}: }{multinomial regression fitted by \link[nnet]{multinom}}
#' \item{\code{ordinal}: }{ordered logistic regression fitted by \link[MASS]{polr}}
#' \item{\code{coxph}: }{cox proportional hazard model fitted by \link[survival]{coxph}}
#' \item{\code{aft_exp}: }{accelerated failure time model fitted by \link[survival]{survreg}
#' with \code{dist = "exponential"}}
#' \item{\code{aft_weibull}: }{accelerated failure time model fitted by \link[survival]{survreg}
#' with \code{dist = "weibull"}}
#' }
#' \code{coxph}, \code{aft_exp} and \code{aft_weibull} are currently not implemented for \code{mreg}.
#' 
#' If \code{EMint} is \code{TRUE} and \code{yreg} is specified by the character name of a regression 
#' model, \code{yreg} is fitted with the interaction between the exposure and each mediator. 
#' 
#' 
#' \emph{A Fitted Regression Object:} 
#' \itemize{
#' \item{}{Regression objects can be fitted by \link{lm}, \link{glm}, \link{glm.nb}, 
#' \link[mgcv]{gam}, \link[nnet]{multinom}, \link[MASS]{polr}, \link[survival]{coxph} and
#' \link[survival]{survreg}. }
#' \item{}{Regression objects fitted by \link[survival]{coxph} and \link[survival]{survreg} 
#' are currently not supported for \code{mreg}. }
#' \item{}{\code{yreg} should regress \code{outcome} on \code{exposure},
#' \code{mediator} and \code{basec}.}
#' \item{}{For \code{p=1,...,k}, \code{mreg[p]} should regress \code{mediator[p]} on 
#' \code{exposure} and \code{basec}, where \code{k} is the number of mediators.}
#' \item{}{\code{yreg} can't include mediator-mediator interactions when there 
#' are multiple mediators (VanderWeele TJ & Vansteelandt, 2014).}
#' }
#' 
#' 
#' \strong{Estimation Methods} 
#' \itemize{
#' \item{\code{paramfunc}: }{(only available for a single 
#' mediator) \emph{closed-form parameter function estimation} by Valeri & VanderWeele (2013).
#' Each causal effect is estimated by a closed-form formula of regression coefficients. }
#' \item{\code{imputation}: }{\emph{direct counterfactual imputation estimation} by Imai, et al (2010). 
#' Each causal effect is estimated by imputing counterfactuals directly.}
#' }
#' 
#' To use \code{paramfunc}, \code{yreg} and \code{mreg} must be specified by the character name 
#' of a regression model. \code{yreg} can be chosen from \code{linear}, \code{logistic}, \code{loglinear}, 
#' \code{poisson}, \code{quasipoisson}, \code{negbin}, \code{coxph}, \code{aft_exp} and 
#' \code{aft_weibull}. \code{mreg} can be chosen from \code{linear}, \code{logistic} and 
#' \code{multinomial}.
#' 
#' To use \code{paramfunc} with \code{yreg = "logistic"} or \code{yreg = "coxph"}, the outcome must 
#' be rare.
#' 
#' 
#' \strong{Inference Methods}
#' \itemize{
#' \item{\code{delta}: }{(only available when \code{estimation = "paramfunc"}) inferences 
#' about causal effects are obtained by the delta method. }
#' \item{\code{bootstrap}: }{inferences about causal effects are obtained by bootstrapping. }
#' }  
#' 
#' 
#' \strong{Estimated Causal Effects}
#' 
#' For a continuous outcome, causal effects on the difference scale are estimated. For a 
#' categorical, count or survival outcome, causal effects on the ratio scale are estimated. Depending 
#' on the outcome type, the ratio can be risk ratio for a categorical outcome, rate ratio for a count outcome, 
#' hazard ratio for a survival outcome fitted by \link[survival]{coxph}, mean survival ratio for a survival outcome fitted 
#' by \link[survival]{survreg}, etc.
#' 
#' When \code{EMint} is \code{FALSE}, \emph{two-way decomposition} (Valeri & VanderWeele, 2013) is conducted, i.e.,
#' \itemize{
#' \item{for a continuous outcome: }{\code{cde} (controlled direct effect), \code{pnde} (pure natural 
#' direct effect), \code{tnde} (total natural direct effect), \code{pnie} (pure natural indirect 
#' effect), \code{tnie} (total natural indirect effect), \code{te} (total effect), and 
#' \code{pm} (proportion mediated) are estimated. }
#' \item{for a categorical, count or survival outcome: }{\code{Rcde} (\code{cde} ratio), \code{Rpnde} (\code{pnde} ratio), 
#' \code{Rtnde} (\code{tnde} ratio), \code{Rpnie} (\code{pnie} ratio), \code{Rtnie} (\code{tnie} ratio), 
#' \code{Rte} (\code{te} ratio), and \code{pm} are estimated.}
#' }
#' 
#' When \code{EMint} is \code{TRUE}: additional \emph{four-way decomposition} (VanderWeele, 2014) is conducted, i.e., 
#' \itemize{
#' \item{for a continuous outcome: }{ \code{intref} 
#' (reference interaction), \code{intmed} (mediated interaction), 
#' \code{cde(prop)} (proportion \code{cde}), \code{intref(prop)} (proportion 
#' \code{intref}), \code{intmed(prop)} (proportion \code{intmed}), \code{pnie(prop)} 
#' (proportion \code{pnie}), \code{int} (proportion 
#' attributable to interaction), and \code{pe} (proportion eliminated) are estimated.}
#' \item{for a categorical, count or survival outcome: }{\code{ERcde} (excess ratio due to \code{cde}), \code{ERintref} (excess 
#' ratio due to \code{intref}), \code{ERintmed} (excess ratio due to \code{intmed}), \code{ERpnie} 
#' (excess ratio due to \code{pnie}), \code{ERcde(prop)} (proportion \code{ERcde}), 
#' \code{ERintref(prop)} (proportion \code{ERintref}), \code{ERintmed(prop)} (proportion \code{ERintmed}), 
#' \code{ERpnie(prop)} (proportion \code{ERpnie}), \code{int}, and \code{pe} are estimated. }
#' }
#' 
#' When \code{EMint} is \code{TRUE} and \code{estimation} is \code{paramfunc}, 
#' causal effects conditional on \code{basecval} are estimated. 
#' Otherwise, marginal causal effects are estimated.
#' 
#' 
#' @return
#' An object of classes \code{cmest} and \code{cmest_rb} is returned:
#' \item{call}{the function call,}
#' \item{data}{the data frame,}
#' \item{methods}{a list of methods used which may include \code{estimation}, \code{inference}, 
#' \code{nboot}, \code{boot.ci.type}, \code{casecontrol}, \code{yrare}, and \code{yprevalence},}
#' \item{variables}{a list of variables used which may include \code{outcome}, \code{event}, 
#' \code{exposure}, \code{mediator}, \code{EMint}, and \code{basec},}
#' \item{ref}{reference values used which may include \code{astar}, \code{a}, \code{mval}, 
#' \code{basecval} and \code{yval},}
#' \item{reg.input}{a list of regressions input,}
#' \item{reg.output}{a list of regressions output. If \code{multimp} is \code{TRUE}, 
#' reg.output contains regression models fitted by each imputed dataset,}
#' \item{multimp}{a list of arguments used for multiple imputation,}
#' \item{effect.pe}{point estimates of causal effects,}
#' \item{effect.se}{standard errors of causal effects,}
#' \item{effect.ci.low}{lower limits of the 95\% confidence intervals of causal effects,}
#' \item{effect.ci.high}{higher limits of the 95\% confidence intervals of causal effects,}
#' \item{effect.pval}{p-values of causal effects,}
#' ...
#'
#' @seealso \link{cmest_gformula}, \link{cmest_wb}, \link{cmest_iorw}, \link{cmest_msm}, 
#' \link{cmest_multistate}, \link{ggcmest}, \link{cmdag}, \link{cmsens}.
#'
#' @references
#' Valeri L, VanderWeele TJ (2013). Mediation analysis allowing for
#' exposure-mediator interactions and causal interpretation: theoretical assumptions and
#' implementation with SAS and SPSS macros. Psychological Methods. 18(2): 137 - 150.
#' 
#' VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple mediators.
#' Epidemiologic Methods. 2(1): 95 - 115.
#'  
#' VanderWeele TJ (2014). A unification of mediation and interaction: a 4-way decomposition. 
#' Epidemiology. 25(5): 749 - 61.
#' 
#' Imai K, Keele L, Tingley D (2010). A general approach to causal mediation analysis.
#' Psychological Methods. 15(4): 309 - 334.
#' 
#' Schomaker M, Heumann C (2018). Bootstrap inference when using multiple 
#' imputation. Statistics in Medicine. 37(14): 2252 - 2266. 
#' 
#' Efron B (1987). Better Bootstrap Confidence Intervals. Journal of the American Statistical 
#' Association. 82(397): 171-185.
#'
#' @examples
#' 
#' \dontrun{
#' library(CMAverse)
#' 
#' # single-mediator case without exposure-mediator interaction
#' exp1 <- cmest_rb(data = cma2020, outcome = "contY", 
#' exposure = "A", mediator = "M1", basec = c("C1", "C2"), 
#' EMint = FALSE, yreg = "linear", mreg = list("logistic"), 
#' estimation = "paramfunc", inference = "delta", astar = 0, a = 1, mval = list(1))
#' summary(exp1)
#' 
#' # single-mediator case with exposure-mediator interaction
#' exp2 <- cmest_rb(data = cma2020, outcome = "contY", 
#' exposure = "A", mediator = "M2", basec = c("C1", "C2"), 
#' EMint = TRUE, yreg = "linear", mreg = list("multinomial"), 
#' estimation = "paramfunc", inference = "delta", astar = 0, a = 1, mval = list("M2_0"))
#' summary(exp2)
#' 
#' # multiple-mediators case
#' exp3 <- cmest_rb(data = cma2020, outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), EMint = TRUE, basec = c("C1", "C2"), 
#' yreg = "linear", mreg = list("logistic", "multinomial"), 
#' estimation = "imputation", inference = "bootstrap", 
#' astar = 0, a = 1, mval = list(0, "M2_0"), 
#' nboot = 100, boot.ci.type = "bca")
#' summary(exp3)
#' 
#' # specify regression models by fitted regression objects
#' exp4 <- cmest_rb(data = cma2020, outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), EMint = TRUE, basec = c("C1", "C2"), 
#' yreg = lm(contY ~ A + M1 + M2 + C1 + C2, data = cma2020), 
#' mreg = list(glm(M1 ~ A + C1 + C2, data = cma2020, family = binomial()),
#' nnet::multinom(M2 ~ A + C1 + C2, data = cma2020)),
#' estimation = "imputation", inference = "bootstrap", 
#' astar = 0, a = 1, mval = list(0, "M2_0"), 
#' nboot = 100, boot.ci.type = "bca")
#' summary(exp4)
#' 
#' }
#' 
#' @importFrom stats glm binomial poisson as.formula gaussian quasipoisson model.frame printCoefmat 
#' family sd coef vcov sigma predict rbinom rmultinom rnorm rgamma rpois weighted.mean 
#' model.matrix getCall quantile qnorm pnorm lm cov formula update na.pass na.omit
#' @importFrom nnet multinom
#' @importFrom MASS polr glm.nb gamma.shape rnegbin
#' @importFrom survival survreg coxph Surv
#' @importFrom survey svyglm svydesign as.svrepdesign withReplicates
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom EValue evalues.RR
#' @importFrom Matrix bdiag
#' @importFrom msm deltamethod
#' @importFrom SuppDists rinvGauss
#' @importFrom dplyr select left_join count
#' @importFrom medflex neImpute
#' @importFrom boot boot
#' @importFrom mice complete mice
#' @importFrom simex check.mc.matrix
#' @importFrom ggplot2 ggproto ggplot geom_errorbar aes geom_point ylab geom_hline position_dodge2 
#' scale_colour_hue theme element_blank facet_grid 
#' @importFrom predint rqpois
#' 
#' @export

cmest_rb <- function(data = NULL, outcome = NULL, event = NULL, 
                     exposure = NULL, mediator = NULL, EMint = NULL, basec = NULL, 
                     yreg = NULL, mreg = NULL, estimation = "imputation", inference = "bootstrap",
                     astar = NULL, a = NULL, mval = NULL, basecval = NULL, yval = NULL, 
                     nboot = 200, boot.ci.type = "per", casecontrol = FALSE, yrare = NULL, yprevalence = NULL,
                     multimp = FALSE, args_mice = NULL) {
  # function call
  cl <- match.call()
  # output list
  out <- list(call = cl)
  
  ###################################################################################################
  #################################Argument Restrictions And Warnings################################
  ###################################################################################################
  # data
  if (is.null(data)) stop("Unspecified data")
  data <- as.data.frame(data)
  if (!all(c(outcome, event, exposure, mediator, basec) %in% colnames(data))) stop("Variables specified in outcome, event, exposure, mediator and basec not found in data")
  if (sum(is.na(data[, c(outcome, event, exposure, mediator, basec)])) > 0 && !multimp) stop("NAs in outcome, event, exposure, mediator or basec data; delete rows with NAs in these variables from the data or set multimp = TRUE") 
  out$data <- data
  n <- nrow(data)
  
  # estimation, inference, nboot
  if (estimation == "para") estimation <- "paramfunc"
  if (estimation == "impu") estimation <- "imputation"
  if (inference == "delt") inference <- "delta"
  if (inference == "boot") inference <- "bootstrap"
  if (!estimation %in% c("paramfunc", "imputation")) stop("Select estimation from 'paramfunc', 'imputation'")
  if (estimation == "paramfunc") {
    if (!is.character(yreg) | length(yreg) != 1 | 
        !yreg %in% c("linear", "logistic", "loglinear", "poisson",
                     "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) stop(
                       "When estimation = 'paramfunc', select yreg from 'linear', 'logistic', 
                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'")
    if (length(mediator) > 1) stop("The parameter function estimation only supports a single mediator")
    if (!is.character(mreg[[1]]) | length(mreg[[1]]) != 1 | 
        !mreg[[1]] %in% c("linear", "logistic", "multinomial")) stop(
          "When estimation = 'paramfunc', select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
  }
  out$methods$estimation <- estimation
  if (!inference %in% c("delta", "bootstrap")) stop("Select inference from 'delta', 'bootstrap'")
  if (estimation == "imputation" && inference == "delta") stop("Use inference = 'bootstrap' when estimation = 'imputation'")
  out$methods$inference <- inference
  if (inference == "bootstrap") {
    if (!is.numeric(nboot)) stop("nboot should be numeric")
    if (nboot <= 0) stop("nboot must be greater than 0") 
    if (!boot.ci.type %in% c("per", "bca")) stop("Select boot.ci.type from 'per', 'bca'")
    out$methods$nboot <- nboot
    out$methods$boot.ci.type <- boot.ci.type
  }
  
  # casecontrol, yrare, yprevalence
  if (!is.logical(casecontrol)) stop("casecontrol must be TRUE or FALSE")
  out$methods$casecontrol <- casecontrol
  if (casecontrol) {
    if (length(unique(data[, outcome])) != 2) stop("When casecontrol is TRUE, the outcome must be binary")
    if (is.null(yprevalence) && yrare != TRUE) stop("When casecontrol is TRUE and the outcome is not rare, specify yprevalence")
    # imputation-based estimation is biased for case-control studies without specifying yprevalence 
    if (is.null(yprevalence) && estimation != "paramfunc") stop("When casecontrol is TRUE, specify yprevalence or use estimation = 'paramfunc'")
    if (!is.null(yprevalence)) {
      if (!is.numeric(yprevalence)) stop("yprevalence must be numeric")
      if (!(yprevalence > 0 && yprevalence < 1)) stop("yprevalence must be between 0 and 1")
      out$methods$yprevalence <- yprevalence
    } else out$methods$yrare <- yrare
  }
  
  # outcome
  if (length(outcome) == 0) stop("Unspecified outcome")
  if (length(outcome) > 1) stop("length(outcome) > 1")
  out$variables$outcome <- outcome
  
  # event
  if (is.character(yreg)) {
    if (yreg %in% c("coxph", "aft_exp", "aft_weibull") && length(event) == 0) stop("Unspecified event")
    out$variables$event <- event
  } else {
    if (!is.null(event)) warning("event is ignored when yreg is not character")
  }
  
  # exposure
  if (length(exposure) == 0) stop("Unspecified exposure")
  if (length(exposure) > 1) stop("length(exposure) > 1")
  out$variables$exposure <- exposure
  
  # mediator
  if (length(mediator) == 0) stop("Unspecified mediator")
  out$variables$mediator <- mediator
  
  # EMint
  if (!is.logical(EMint)) stop("EMint must be TRUE or FALSE")
  out$variables$EMint <- EMint
  
  # basec
  if (!is.null(basec)) out$variables$basec <- basec
  
  #regs
  out$reg.input <- list(yreg = yreg, mreg = mreg)
  
  # mval
  if (!is.list(mval)) stop("mval should be a list")
  if (length(mval) != length(mediator)) stop("length(mval) != length(mediator)")
  
  # basecval
  if (!(estimation == "paramfunc" && length(basec) != 0) && !is.null(basecval)) warning("basecval is ignored")
  
  # multimp
  out$multimp <- list(multimp = multimp)
  if (!is.logical(multimp)) stop("multimp must be TRUE or FALSE")
  # args_mice
  if (multimp) {
    if (!is.null(args_mice) && !is.list(args_mice)) stop("args_mice must be a list")
    if (is.null(args_mice)) args_mice <- list()
    args_mice$print <- FALSE
    if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
    args_mice$data <- data
    out$multimp$args_mice <- args_mice
  }
  
  # run regressions
  environment(regrun_rb) <- environment()
  regs <- regrun_rb()
  yreg <- regs$yreg
  mreg <- regs$mreg
  
  ###################################################################################################
  ############################################Estimation and Inference###############################
  ###################################################################################################
  # add a progress bar for bootstrap inference
  if (inference == "bootstrap") {
    env <- environment()
    counter <- 0
    progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
  }
  # estimation and inference of causal effects
  environment(estinf_rb) <- environment()
  out <- c(out, estinf_rb())
  class(out) <- c("cmest", "cmest_rb")
  return(out)
}


# This function runs regressions for the regression-based approach
regrun_rb <- function() {
  # Mediator Regression: a mediator regression is required for each mediator
  if (is.null(mreg)) stop("mreg is required")
  if (!is.list(mreg)) stop("mreg must be a list")
  if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
  for (p in 1:length(mreg)) {
    if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
    if (is.character(mreg[[p]])) {
      if (mreg[[p]] == "loglinear" && length(unique(na.omit(data[, mediator[p]]))) != 2) stop(paste0("When mreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
      if (!mreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                            "negbin", "multinomial", "ordinal")) stop(
                              paste0("Select character mreg[[", p, "]] from 'linear', 'logistic',
                                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
      # regress each mediator on the exposure and confounders
      mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, basec), collapse = "+"))
      switch(mreg[[p]],
             linear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = gaussian(), data = .(data)))),
             logistic = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = binomial(), data = .(data)))),
             loglinear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
             poisson = mreg[[p]]  <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
             quasipoisson = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = quasipoisson(), data = .(data)))),
             negbin = mreg[[p]] <- eval(bquote(MASS::glm.nb(.(as.formula(mediator_formula)), data = .(data)))),
             multinomial = mreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(mediator_formula)), data = .(data), trace = FALSE))),
             ordinal = mreg[[p]] <- eval(bquote(MASS::polr(.(as.formula(mediator_formula)), data = .(data)))))
    } else if (!(identical(class(mreg[[p]]), "lm") | identical(class(mreg[[p]]), c("glm", "lm")) |
                 identical(class(mreg[[p]]), c("negbin", "glm", "lm")) | identical(class(mreg[[p]]), c("multinom", "nnet")) |
                 identical(class(mreg[[p]]), c("gam", "glm", "lm")) | identical(class(mreg[[p]]), "polr"))) {
      stop(paste0("Fit mreg[[", p, "]] by lm, glm, glm.nb, gam, multinom, polr"))
    } else if (inherits(mreg[[p]], "lm") | inherits(mreg[[p]], "glm")) {
      family_mreg <- family(mreg[[p]])
      if (!(family_mreg$family %in% 
            c("gaussian", "inverse.gaussian", "poisson", "quasipoisson", 
              "Gamma", "binomial", "multinom") |
            startsWith(family_mreg$family, "Negative Binomial") |
            startsWith(family_mreg$family, "Ordered Categorical"))) stop(paste0("Unsupported mreg[[", p, "]]"))
    } else {
      mreg_formula <- formula(mreg[[p]])
      d_var <- unique(all.vars(mreg_formula[[2]]))
      ind_var <- unique(all.vars(mreg_formula[[3]]))
      if ((mediator[p] != d_var) | !all(ind_var %in% c(exposure, basec))) stop(
        paste0("For mreg[[", p, "]], regress mediator[", p, "] on variables in c(exposure, basec)"))
      rm(mreg_formula, d_var, ind_var)
    }
  }
  
  # Outcome Regression
  if (is.null(yreg)) stop("yreg is required")
  if (is.character(yreg)) {
    if (yreg == "loglinear" && length(unique(na.omit(data[, outcome]))) != 2) stop("When yreg is 'loglinear', outcome should be binary")
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                     "negbin", "multinomial", "ordinal", "coxph", "aft_exp",
                     "aft_weibull")) stop(
                       paste0("Select character yreg from 'linear', 'logistic',
                              'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal',
                              'coxph', 'aft_exp', 'aft_weibull'"))
    if (estimation == "paramfunc" && yreg %in% c("logistic", "coxph")) warning("When estimation is 'paramfunc' and yreg is 'logistic' or 'coxph', the outcome must be rare; ignore this warning if the outcome is rare")
    int.terms <- switch(EMint + 1, "1" = NULL, "2" = paste(exposure, mediator, sep = "*"))
    # regress the outcome on the exposure, mediators and confounders
    outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+"))
    if (yreg %in% c("coxph","aft_exp","aft_weibull")) {
      if (!is.null(event)) {
        outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                                 strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
      } else outcome_formula <- paste(paste0("Surv(", outcome, ")"),
                                      strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
    }
    switch(yreg,
           linear = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                            family = gaussian(), data = .(data)))),
           logistic = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                              family = binomial(), data = .(data)))),
           loglinear = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                               family = poisson(), data = .(data)))),
           poisson = yreg  <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                              family = poisson(), data = .(data)))),
           quasipoisson = yreg <- eval(bquote(glm(formula = .(as.formula(outcome_formula)),
                                                  family = quasipoisson(), data = .(data)))),
           negbin = yreg <- eval(bquote(MASS::glm.nb(formula = .(as.formula(outcome_formula)),
                                                     data = .(data)))),
           multinomial = yreg <- eval(bquote(nnet::multinom(formula = .(as.formula(outcome_formula)),
                                                            data = .(data), trace = FALSE))),
           ordinal = yreg <- eval(bquote(MASS::polr(formula = .(as.formula(outcome_formula)),
                                                    data = .(data)))),
           coxph = yreg <- eval(bquote(survival::coxph(formula = .(as.formula(outcome_formula)),
                                                       data = .(data)))),
           aft_exp = yreg <- eval(bquote(survival::survreg(formula = .(as.formula(outcome_formula)),
                                                           dist = "exponential", data = .(data)))),
           aft_weibull = yreg <- eval(bquote(survival::survreg(formula = .(as.formula(outcome_formula)),
                                                               dist = "weibull", data = .(data)))))
  } else if (!(identical(class(yreg), "lm") | identical(class(yreg), c("glm", "lm")) |
               identical(class(yreg), c("negbin", "glm", "lm")) | identical(class(yreg), c("multinom", "nnet")) |
               identical(class(yreg), c("gam", "glm", "lm")) | identical(class(yreg), "polr") |
               identical(class(yreg), "coxph") | identical(class(yreg), "survreg"))) {
    stop("Fit yreg by lm, glm, glm.nb, gam, multinom, polr, coxph, survreg")
  } else if (inherits(yreg, "lm") | inherits(yreg, "glm")) {
    family_yreg <- family(yreg)
    if (!(family_yreg$family %in% 
          c("gaussian", "inverse.gaussian", "quasi", "poisson", "quasipoisson", 
            "Gamma", "binomial", "quasibinomial", "multinom", "ziplss") |
          startsWith(family_yreg$family, "Negative Binomial") |
          startsWith(family_yreg$family, "Zero inflated Poisson") |
          startsWith(family_yreg$family, "Ordered Categorical"))) stop("Unsupported yreg")
  } else {
    warning("Make sure there is no mediator-mediator interaction in yreg; ignore this warning if there isn't")
    yreg_formula <- formula(yreg)
    d_var <- unique(all.vars(yreg_formula[[2]]))
    ind_var <- unique(all.vars(yreg_formula[[3]]))
    if (!(outcome %in% d_var) | !all(ind_var %in% c(exposure, mediator, basec))) stop(
      "For yreg, regress outcome on variables in c(exposure, mediator, basec)")
    rm(yreg_formula, d_var, ind_var)
  }
  
  
  # for delta method inference, use survey regressions for yreg and mreg when weights are applied
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    yreg <- suppressWarnings(eval(bquote(survey::svyglm(formula = .(formula(yreg)), family = .(family(yreg)),
                                                        design = survey::svydesign(~ 1, data = .(data))))))
  }
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    if (inherits(mreg[[1]], "glm")) mreg[[1]] <- suppressWarnings(eval(bquote(survey::svyglm(formula = .(formula(mreg[[1]])), 
                                                                                             family = .(family(mreg[[1]])),
                                                                                             design = survey::svydesign(~ 1, data = .(data))))))
    if (inherits(mreg[[1]], "multinom")) mreg[[1]] <- suppressWarnings(eval(bquote(svymultinom(formula = .(formula(mreg[[1]])), 
                                                                                               data = .(data)))))
  }
  out <- list(yreg = yreg, mreg = mreg)
  return(out)
}



# This function estimates causal effects and make inferences using the regression-based approach
estinf_rb <- function() {
  # restrict data types of variables
  allvar <- c(outcome, event, exposure, mediator, basec)
  for (i in 1:length(allvar))
    if (!(is.numeric(data[, allvar[i]]) | is.factor(data[, allvar[i]]) |
          is.character(data[, allvar[i]]))) stop(paste0("The variable ", allvar[i], " must be numeric, factor or character"))
  
  # output list
  out <- list()
  
  # obtain regression calls
  call_yreg <- getCall(yreg)
  call_mreg <- lapply(1:length(mreg), function(x) getCall(mreg[[x]]))
  
  # obtain outcome regression classes
  if (inherits(yreg, "rcreg") | inherits(yreg, "simexreg")) {
    yreg_mid <- yreg$NAIVEreg
  } else yreg_mid <- yreg
  is_lm_yreg <- inherits(yreg_mid, "lm")
  is_glm_yreg <- inherits(yreg_mid, "glm")
  if (is_lm_yreg | is_glm_yreg) family_yreg <- family(yreg_mid)
  is_svyglm_yreg <- inherits(yreg_mid, "svyglm")
  is_gam_yreg <- inherits(yreg_mid, "gam")
  is_multinom_yreg <- inherits(yreg_mid, "multinom")
  is_svymultinom_yreg <- inherits(yreg_mid, "svymultinom")
  is_polr_yreg <- inherits(yreg_mid, "polr")
  is_survreg_yreg <- inherits(yreg_mid, "survreg")
  is_coxph_yreg <- inherits(yreg_mid, "coxph")
  rm(yreg_mid)
  
  # obtain outcome regression weights
  if (is_svyglm_yreg) {
    weights_yreg <- yreg$data$.survey.prob.weights
  } else {
    weights_yreg <- model.frame(yreg)$'(weights)'
  }
  
  # obtain mediator regression classes 
  mreg_mid <- lapply(1:length(mreg), function(x) {
    if (inherits(mreg[[x]], "rcreg") | inherits(mreg[[x]], "simexreg")) {
      mreg[[x]]$NAIVEreg
    } else mreg[[x]] 
  })
  is_lm_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "lm"))
  is_glm_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "glm"))
  family_mreg <- lapply(1:length(mreg_mid), function(x) if (is_lm_mreg[x] | is_glm_mreg[x]) family(mreg_mid[[x]]))
  is_svyglm_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "svyglm"))
  is_gam_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "gam"))
  is_multinom_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "multinom"))
  is_svymultinom_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "svymultinom"))
  is_polr_mreg <- sapply(1:length(mreg_mid), function(x) inherits(mreg_mid[[x]], "polr"))
  rm(mreg_mid)
  
  # obtain mediator regression weights
  weights_mreg <- lapply(1:length(mreg), function(x) {
    if (is_svyglm_mreg[x]) { 
      mreg[[x]]$data$.survey.prob.weights
    } else model.frame(mreg[[x]])$'(weights)'
  })
  
  # reference values for the exposure
  if (is.factor(data[, exposure]) | is.character(data[, exposure])) {
    a_lev <- levels(droplevels(as.factor(data[, exposure])))
    if (is.null(a) | !a %in% a_lev) {
      a <- a_lev[length(a_lev)]
      warning(paste0("a is not a value of the exposure; ", a, " is used"))
    }
    if (is.null(astar) | !astar %in% a_lev) {
      astar <- a_lev[1]
      warning(paste0("astar is not a value of the exposure; ", astar, " is used"))
    }
  } else{
    if (is.null(a)) {
      a <- quantile(data[, exposure], probs = 0.75, names = F)
      warning(paste0("a is missing; the 3rd quartile is used"))
    }
    if (is.null(astar)) {
      astar <- quantile(data[, exposure], probs = 0.25, names = F)
      warning(paste0("astar is missing; the 1st quartile is used"))
    }
  }
  out$ref$a <- a
  out$ref$astar <- astar
  
  # yval: the reference level for a categorical outcome
  if ((is_glm_yreg && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                       startsWith(family_yreg$family, "Ordered Categorical"))) |
      is_multinom_yreg | is_polr_yreg) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    # if yval is not provided or yval provided is not a value of the outcome, use the last level of the outcome
    if (is.null(yval)) {
      yval <- y_lev[length(y_lev)]
      warning(paste0("yval is not specified; ", yval, " is used"))
    }
    if (!yval %in% y_lev) {
      yval <- y_lev[length(y_lev)]
      warning(paste0("yval is not a value of the outcome; ", yval, " is used"))
    }
    out$ref$yval <- yval
  }
  
  # reference values for the mediators
  for (i in 1:length(mediator)) {
    if ((is_glm_mreg[i] && ((family_mreg[[i]]$family %in% c("binomial", "multinom")) |
                            startsWith(family_mreg[[i]]$family, "Ordered Categorical")))|
        is_multinom_mreg[i] | is_polr_mreg[i]) {
      m_lev <- levels(droplevels(as.factor(data[, mediator[i]])))
      if (is.numeric(data[, mediator[i]])) m_lev <- as.numeric(m_lev)
      if (is.na(mval[[i]]) | !mval[[i]] %in% m_lev) {
        mval[[i]] <- m_lev[length(m_lev)]
        warning(paste0("mval[[", i, "]] is not a value of mediator[", i, "]; ", mval[[i]], " is used"))
      }
    } else {
      if (is.na(mval[[i]])) {
        mval[[i]] <- mean(data[, mediator[i]])
        warning(paste0("mval[[", i, "]] is missing, the mean value is used"))
      }
    }
  }
  out$ref$mval <- mval
  
  # get the level of the case and the level of the control
  if (casecontrol) {
    y_lev <- levels(droplevels(as.factor(data[, outcome])))
    y_control <- y_lev[1]
    y_case <- y_lev[2]
  }
  
  # closed-form parameter function estimation
  if (estimation == "paramfunc") {
    # create a list of covariate values to calculate conditional causal effects
    if (length(basec) != 0) {
      if (!is.null(basecval)) {
        if (!is.list(basecval)) stop("basecval should be a list")
        if (length(basecval) != length(basec)) stop("length(basecval) != length(basec)")
      }
      if (is.null(basecval)) basecval <- rep(list(NULL), length(basec))
      # if NULL, set basecval[[i]] to be the mean value of basec[i]
      for (i in 1:length(basec)) {
        if (is.factor(data[, basec[i]]) | is.character(data[, basec[i]])) {
          c_lev <- levels(droplevels(as.factor(data[, basec[i]])))
          if (is.null(basecval[[i]])) {
            c_data <- data[, basec[i], drop = FALSE]
            c_data[, basec[i]] <- factor(c_data[, basec[i]], levels = c_lev)
            # set basecval[[i]] to be the mean values of dummy variables
            basecval[[i]] <- unname(colMeans(as.matrix(model.matrix(as.formula(paste0("~", basec[i])),
                                                                    data = model.frame(~., data = c_data, 
                                                                                       na.action = na.pass))[, -1]), 
                                             na.rm = TRUE))
            rm(c_data)
            # dummy values of basecval[[i]]
          } else basecval[[i]] <- as.numeric(c_lev == basecval[[i]])[-1]
          rm(c_lev)
        } else if (is.numeric(data[, basec[i]])) {
          if (is.null(basecval[[i]])) {
            # set basecval[[i]] to be the mean value of basec[i]
            basecval[[i]] <- mean(data[, basec[i]], na.rm = TRUE)
          } else basecval[[i]] <- basecval[[i]]
        } 
      }
      out$ref$basecval <- basecval
    }
  }
  
  # estimation and inference
  environment(est_rb) <- environment()
  if (!multimp) {
    # point estimates of causal effects
    est <- est_rb(data = data, indices = NULL, outReg = TRUE)
    effect.pe <- est$est
    n_effect <- length(effect.pe)
    out$reg.output <- est$reg.output
    if (inference == "bootstrap") {
      # bootstrap results
      boots <- boot(data = data, statistic = est_rb, R = nboot, outReg = FALSE)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else if (inference == "delta") {
      yreg <- est$reg.output$yreg
      mreg <- est$reg.output$mreg[[1]]
      # standard errors by the delta method
      environment(inf.delta) <- environment()
      effect.se <- inf.delta(data = data, yreg = yreg, mreg = mreg)
      # critical value
      z0 <- qnorm(0.975)
      z <- effect.pe/effect.se
      # delta method CIs
      effect.ci.low <- effect.pe - z0 * effect.se
      effect.ci.high <- effect.pe + z0 * effect.se
      # delta method p-values
      effect.pval <- 2 * (1 - pnorm(abs(z)))
    }
  } else {
    # imputed data sets
    data_imp <- complete(do.call(mice, args_mice), action = "all")
    m <- length(data_imp)
    # estimate causal effects for each imputed data set
    est_imp <- lapply(1:m, function(x)
      est_rb(data = data_imp[[x]], indices = NULL, outReg = TRUE))
    est_imp_df <- do.call(rbind, lapply(1:m, function(x) est_imp[[x]]$est))
    effect.pe <- colMeans(est_imp_df)
    n_effect <- length(effect.pe)
    out$reg.output <- lapply(1:m, function(x) est_imp[[x]]$reg.output)
    
    if (inference == "bootstrap") {
      boot.step <- function(data = NULL, indices = NULL) {
        data_boot <- data[indices, ]
        args_mice$data <- data_boot
        data_imp <- complete(do.call(mice, args_mice), action = "all")
        curVal <- get("counter", envir = env)
        assign("counter", curVal + 1, envir = env)
        setTxtProgressBar(get("progbar", envir = env), curVal + 1)
        return(colMeans(do.call(rbind, lapply(1:m, function(x)
          est_rb(data = data_imp[[x]], outReg = FALSE)))))
      }
      environment(boot.step) <- environment()
      # bootstrap results
      boots <- boot(data = data, statistic = boot.step, R = nboot)
      # bootstrap CIs
      environment(boot.ci) <- environment()
      effect.ci <- boot.ci(boots = boots)
      effect.ci.low <- effect.ci[, 1]
      effect.ci.high <- effect.ci[, 2]
      # bootstrap p-values
      effect.pval <- sapply(1:n_effect, function(x) boot.pval(boots = boots$t[, x], pe = effect.pe[x]))
    } else if (inference == "delta") {
      environment(inf.delta) <- environment()
      # standard errors by the delta method
      se_imp <- do.call(rbind, lapply(1:m, function(x)
        inf.delta(data = data_imp[[x]], yreg = est_imp[[x]]$reg.output$yreg, mreg = est_imp[[x]]$reg.output$mreg[[1]])))
      # pool the results by Rubin's rule
      var_within <- colMeans(se_imp ^ 2)
      var_between <- colSums((est_imp_df - t(replicate(m, effect.pe)))^2)/(m - 1)
      effect.se <- sqrt(var_within + var_between * (m + 1) / m)
      z0 <- qnorm(0.975)
      z <- effect.pe/effect.se
      effect.ci.low <- effect.pe - z0 * effect.se
      effect.ci.high <- effect.pe + z0 * effect.se
      effect.pval <- 2 * (1 - pnorm(abs(z)))
    }
  }
  
  if ((is_lm_yreg | is_glm_yreg) &&
      (family_yreg$family %in% c("gaussian", "inverse.gaussian", "Gamma", "quasi"))) {
    # standard errors by bootstrapping
    if (inference == "bootstrap") effect.se <- sapply(1:n_effect, function(x) sd(boots$t[, x], na.rm = TRUE))
    # effect names
    if (EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", 
                                "intref", "intmed", "cde(prop)", "intref(prop)", "intmed(prop)", "pnie(prop)",
                                "pm", "int", "pe")
    if (!EMint) effect_name <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
  } else {
    # transform standard errors of effects in log scale
    if (inference == "bootstrap") effect.se <- sapply(1:n_effect, function(x)
      ifelse(x <= 6, sd(exp(boots$t[, x])), sd(boots$t[, x]))) #
    if (inference == "delta") effect.se[1:6] <- sapply(1:6, function(x)
      deltamethod(as.formula("~exp(x1)"), effect.pe[x], effect.se[x]^2))
    # transform effects in log ratio scale into effects in ratio scale
    effect.pe[1:6] <- exp(effect.pe[1:6])
    effect.ci.low[1:6] <- exp(effect.ci.low[1:6])
    effect.ci.high[1:6] <- exp(effect.ci.high[1:6])
    # effect names
    if (EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", 
                                "ERcde", "ERintref", "ERintmed", "ERpnie",
                                "ERcde(prop)", "ERintref(prop)", "ERintmed(prop)", "ERpnie(prop)",
                                "pm", "int", "pe")
    if (!EMint) effect_name <- c("Rcde", "Rpnde", "Rtnde", "Rpnie", "Rtnie", "Rte", "pm")
  }
  
  names(effect.pe) <- names(effect.se) <- names(effect.ci.low) <- names(effect.ci.high) <-
    names(effect.pval) <- effect_name
  out$effect.pe <- effect.pe
  out$effect.se <- effect.se
  out$effect.ci.low <- effect.ci.low
  out$effect.ci.high <- effect.ci.high
  out$effect.pval <- effect.pval
  
  return(out)
}




