#' Causal Mediation Analysis
#'
#' \code{cmest} is used to conduct causal mediation analysis via six statistical approaches 
#' including \emph{the regression-based approach} by Valeri et al. (2013) and VanderWeele 
#' et al. (2014), \emph{the weighting-based approach} by VanderWeele et al. (2014), 
#' \emph{the inverse odd-ratio weighting approach} by Tchetgen Tchetgen (2013), 
#' \emph{the natural effect model} by Vansteelandt et al. (2012), \emph{the marginal structural 
#' model} by VanderWeele et al. (2017), and \emph{the g-formula approach} by Lin et al. (2017).
#'
#' @param data dataset for causal mediation analysis.
#' @param model causal mediation analysis approach. "\code{rb}", "\code{wb}", "\code{iorw}", 
#' "\code{ne}", "\code{msm}", "\code{gformula}" are implemented. See \code{Details}.
#' @param full a logical value. If \code{TRUE}, output a full list of causal effects; if 
#' \code{FALSE}, output a reduced list of causal effects. Default is \code{TRUE}. See 
#' \code{Details}.
#' @param casecontrol a logical value. \code{TRUE} indicates a case control study in which the
#' first level of the outcome is treated as the control and the second level of the outcome is 
#' treated as the case. Default is \code{FALSE}.
#' @param yrare a logical value (used when \code{casecontrol} is \code{TRUE}). \code{TRUE} 
#' indicates the case is rare.
#' @param yprevalence the prevalence of the case (used when \code{casecontrol} is \code{TRUE}).
#' @param estimation estimation method for causal effects. \code{paramfunc} and 
#' \code{imputation} are implemented (the first 4 letters are enough). Default is \code{imputation}. 
#' See \code{Details}. 
#' @param inference inference method for causal effects. \code{delta} and \code{bootstrap} are 
#' implemented (the first 4 letters are enough). Default is \code{bootstrap}. See 
#' \code{Details}.
#' @param outcome variable name of the outcome.
#' @param event variable name of the event (used when \code{yreg} is \code{coxph}, \code{aft_exp},
#' or \code{aft_weibull}).
#' @param exposure variable name of the exposure.
#' @param mediator a vector of variable name(s) of the mediator(s).
#' @param EMint a logical value (used when \code{yreg} is character). If \code{TRUE}, the 
#' outcome regression formula includes the interaction(s) between the exposure and each of the 
#' mediator(s). Default is \code{FALSE}.
#' @param basec a vector of variable name(s) of the baseline confounder(s)
#' @param postc a vector of variable name(s) of mediator-outcome confounder(s) affected 
#' by the exposure
#' @param yreg outcome regression. See \code{Details}.
#' @param mreg a list of regression(s) specifying a regression model for each of the 
#' mediator(s) (used when \code{model} is \code{rb}, \code{msm} or \code{gformula}). See 
#' \code{Details}.
#' @param wmreg a list of regression(s) specifying a regression model for calculating 
#' weights with respect to each of the mediator(s) (used when \code{model} is \code{msm}). 
#' See \code{Details}.
#' @param ereg exposure regression for calculating weights with respect to the exposure (used 
#' when \code{model} is \code{wb} or \code{msm} with \code{length(basec)} != 0 or when 
#' \code{model} is \code{iorw}). See \code{Details}.
#' @param postcreg a list of regression(s) specifying a regression model for each variable
#' in \code{postc} (used when \code{model} is \code{gformula}). See \code{Details}.
#' @param astar the first reference value for the exposure. Default is \code{0}.
#' @param a the second reference value for the exposure. Default is \code{1}.
#' @param mval a list specifying a reference value for each of the mediator(s) (used when 
#' \code{model} is \code{rb}, \code{wb}, \code{ne}, \code{msm} or \code{gformula}).
#' @param yref reference value for the outcome (used when the outcome is categorical).
#' @param basecval a list specifying a conditional value for each of the baseline 
#' confounder(s) for estimating conditional causal effects (used when \code{estimation} is 
#' \code{paramfunc}). If \code{NULL}, mean value(s) of the baseline confounder(s) are used.
#' @param nboot the number of boots applied (used when \code{inference} is \code{bootstrap}). 
#' Default is 200. 
#' @param nRep number of replications or hypothetical values of the exposure to sample for 
#' each observation unit (used when \code{model} is \code{ne}). See \link[medflex]{neImpute.default}
#' for details. Default is \code{5}.
#' @param multimp a logical value (used when \code{data} contains missing values). If 
#' \code{TRUE}, conduct multiple imputation using the \link[mice]{mice} function. Default is 
#' \code{FALSE}.
#' @param x an object of class 'cmest'
#' @param object an object of class 'cmest'
#' @param digits minimal number of significant digits. See \link{print.default}.
#' @param ... Additional arguments passed to the \link[mice]{mice} function. See \link[mice]{mice}
#' for details.
#'
#' @details
#' 
#' \strong{Regressions}
#'
#' Each regression in \code{yreg}, \code{mreg}, \code{wmreg},
#' \code{ereg} and \code{postcreg} can be specified by a user-defined regression 
#' object or the character name of a regression. 
#' 
#' \emph{The Character Name of A Regression}
#' 
#' \itemize{
#' \item{"\code{linear}": }{linear regression fitted by \link{glm} with \code{family = gaussian()}}
#' \item{"\code{logistic}": }{logistic regression fitted by \link{glm} with \code{family = logit()}}
#' \item{"\code{loglinear}": }{log linear regression fitted by \link{glm} with 
#' \code{family = poisson()} for a binary response}
#' \item{"\code{poisson}": }{poisson regression fitted by \link{glm} with 
#' \code{family = poisson()} for a count response}
#' \item{"\code{quasipoisson}": }{quasipoisson regression fitted by \link{glm} with 
#' \code{family = quasipoisson()}}
#' \item{"\code{negbin}": }{negative binomial regression fitted by \link[MASS]{glm.nb}}
#' \item{"\code{multinomial}": }{multinomial regression fitted by \link[nnet]{multinom}}
#' \item{"\code{ordinal}": }{ordered logistic regression fitted by \link[MASS]{polr}}
#' \item{"\code{coxph}": }{cox proportional hazard model fitted by \link[survival]{coxph}}
#' \item{"\code{aft_exp}": }{accelerated failure time model fitted by \link[survival]{survreg}
#' with \code{dist = "exponential"}}
#' \item{"\code{aft_weibull}": }{accelerated failure time model fitted by \link[survival]{survreg}
#' with \code{dist = "weibull"}}
#' }
#' \code{coxph}, \code{aft_exp} and \code{aft_weibull} are currently not implemented for 
#' \code{mreg}, \code{wmreg}, \code{ereg} and \code{postcreg}.
#' 
#' \emph{The User-defined Regression Object} 
#' 
#' A user-defined regression object can be fitted by \link{lm}, \link{glm}, \link{glm.nb}, 
#' \link[mgcv]{gam}, \link[nnet]{multinom}, \link[MASS]{polr}, \link[survival]{coxph} and
#' \link[survival]{survreg}. Objects fitted by \link[survival]{coxph} and \link[survival]{survreg} 
#' are currently not supported for \code{mreg}, \code{wmreg}, 
#' \code{ereg} and \code{postcreg}.
#' 
#' The \code{cmest} function calculates weights for regressions when weighting is required. If a 
#' user-defined regression object is fitted with prior weights, the final weights for this 
#' regression are constructed by multiplying the prior weights with the weights calculated by 
#' the \code{cmest} function. 
#' 
#' 
#' \strong{Causal Mediation Analysis Approaches}
#' 
#' Let \code{Y} denote the outcome, \code{A} denote the exposure, \code{M=(M_1,...,M_k)^T} 
#' denote the set of mediator(s) following the temporal order, \code{C} denote the set of 
#' variables in \code{basec}, \code{L=(L_1,...,L_s)^T} denote the set of variables in 
#' \code{postc} following the temporal order.
#' 
#' \itemize{
#'     \item{\code{rb}: }{\emph{the regression-based approach} by Valeri et al. (2013) and
#'       VanderWeele et al. (2014). \code{yreg} and \code{mreg} are required. If specified as
#'       a user defined regression object, \code{yreg} should regress \code{Y} on \code{A},
#' \code{M} and \code{C} and \code{mreg[p]} should regress \code{M_p} on \code{A},
#' \code{M_1}, ..., \code{M_{p-1}} and \code{C}} for \code{p=1,...,k}.
#' 
#'     \item{\code{wb}: }{\emph{the weighting-based approach} by VanderWeele et al. (2014).
#'     \code{yreg} is required. When \code{basec} is not empty, \code{ereg} is also required 
#'     and \code{A} must be categorical. If specified as a user defined regression object, 
#'     \code{yreg} should regress \code{Y} on \code{A}, \code{M} and \code{C} and \code{ereg} 
#'     should regress \code{A} on \code{C}.}
#'     
#'     \item{\code{iorw}: }{\emph{the inverse odd-ratio weighting approach} by
#'       Tchetgen Tchetgen (2013). \code{yreg} and \code{ereg} are required and 
#'       \code{A} must be categorical. If specified as a user defined regression object, 
#'       \code{yreg} should regress \code{Y} on \code{A} and \code{C} and \code{ereg} should 
#'       regress \code{A} on \code{M} and \code{C}.}
#'       
#'     \item{\code{ne}: }{\emph{the natural effect model} by Vansteelandt et al. (2012). 
#'     \code{yreg} is required. If specified as a user defined regression object, \code{yreg}
#'     should regress \code{Y} on \code{A}, \code{M} and \code{C}. The variables in the 
#'     formula of \code{yreg} should follow the order of \code{A}, \code{M} and \code{C}, i.e., 
#'     the first variable needs to point to the exposure, the variable(s) right after the 
#'     exposure need to point to the mediator(s), e.g., \code{Y ~ A + M_1 + M_2 + A*M_1 + C}.}
#'       
#'     \item{\code{msm}: }{\emph{the marginal structural model} by VanderWeele et al. (2017).
#'     \code{yreg}, \code{mreg} and \code{wmreg} are required and all 
#'     mediators must be categorical. When \code{basec} is not empty, \code{ereg} is also 
#'     required and \code{A} must be categorical. If specified as a user defined regression 
#'     object, \code{yreg} should regress \code{Y} on \code{A} and \code{M}, 
#'     \code{mreg[p]} should regress \code{M_p} on \code{A}, \code{M_1}, ..., and 
#'     \code{M_{p-1}} for \code{p=1,...,k}, 
#'     \code{wmreg[p]} should regress \code{M_p} on \code{A}, \code{M_1}, ..., \code{M_{p-1}},
#'     \code{C} and \code{L} for \code{p=1,...,k}, and \code{ereg} should regress \code{A} on 
#'     \code{C}.}
#'     
#'     \item{\code{gformula}: }{\emph{the g-formula approach} by Lin et al. (2017). 
#'     \code{yreg}, \code{mreg} are required. \code{postcreg} is also required when \code{postc}
#'     is not empty. If specified as a user defined regression object, \code{yreg} should 
#'     regress \code{Y} on \code{A}, \code{M}, \code{C} and \code{L}, \code{mreg[p]} should 
#'     regress \code{M_p} on \code{A}, \code{M_1}, ..., \code{M_{p-1}}, \code{C} and \code{L},
#'     for \code{p=1,...,k}, \code{postcreg[p]} should regress \code{L_q} on \code{A}, 
#'     \code{L_1}, ..., \code{L_{q-1}} and \code{C} for \code{q=1,...,s}.}
#'   }
#'   
#' When \code{postc} is not empty, \code{msm} and \code{gformula} can be used, while \code{rb}, 
#' \code{wb}, \code{iorw} and \code{ne} can't be used.
#'   
#' 
#' \strong{Estimation Methods} 
#'   
#' \itemize{
#' \item{\code{paramfunc}: }{closed-form parameter function estimation (only available when 
#' \code{model = "rb"} with \code{length(mediator) = 1}). The point estimate of each causal 
#' effect is obtained by its closed-form formula of regression coefficients. Effects 
#' conditional on pre-specified values of \code{C}, i.e., \code{basecval}, are estimated.}
#' \item{\code{imputation}: }{direct counterfactual imputation estimation. The point estimate 
#' of each causal effect is obtained by imputing counterfactuals directly.}
#' }
#' 
#' To use closed-form parameter function estimation, \code{yreg} and \code{mreg} must be 
#' specified by the character name of the regression. \code{yreg} can be chosen from 
#' \code{linear}, \code{logistic}, \code{loglinear}, \code{poisson}, \code{quasipoisson}, 
#' \code{negbin}, \code{coxph}, \code{aft_exp} and \code{aft_weibull}. \code{mreg} can be 
#' chosen from \code{linear}, \code{logistic} and \code{multinomial}.
#' 
#' To use \code{yreg = "logistic"} in closed-form parameter function estimation, the binary 
#' outcome must be rare. For common binary outcomes, use \code{yreg = "loglinear"}.
#' 
#' 
#' \strong{Inference Methods}
#' 
#' \itemize{
#' \item{\code{delta}: }{delta method (only available when \code{estimation = "paramfunc"}). 
#' The standard errors of causal effects are obtained by the delta method. The confidence 
#' intervals of causal effects are obtained by normal distribution approximation.}
#' \item{\code{bootstrap}: }{bootstrapping. The standard errors of causal effects are 
#' obtained by the standard deviations of bootstrapped results. The confidence intervals of 
#' causal effects are obtained by percentiles of bootstrapped results.}
#' }  
#' 
#' 
#' \strong{Estimated Causal Effects}
#' 
#' When there are multiple mediators, the joint effect of them is studied.
#' 
#' For a continuous outcome, the causal effects on the difference scale are estimated. For a 
#' categorical, count or survival outcome, the causal effects on the ratio scale are estimated.
#' 
#' \emph{Continuous Outcome}
#' 
#' When \code{model = "rb", "wb", "ne", "msm" or "gformula"} with \code{length(postc) = 0} and 
#' \code{full = TRUE}, \code{cde} (controlled direct effect), \code{pnde} (pure natural direct 
#' effect), \code{tnde} (total natural direct effect), \code{pnie} (pure natural indirect effect), 
#' \code{tnie} (total natural indirect effect), \code{te} (total effect),  \code{intref} 
#' (reference interaction), \code{intmed} (mediated interaction), 
#' \code{cde(prop)} (proportion \code{cde}), \code{intref(prop)} (proportion 
#' \code{intref}), \code{intmed(prop)} (proportion \code{intmed}), \code{pnie(prop)} 
#' (proportion \code{pnie}), \code{pm} (proportion mediated), \code{int} (proportion 
#' attributable to interaction) and \code{pe} (proportion eliminated) are estimated. 
#' 
#' When \code{model = "rb", "wb", "ne", "msm" or "gformula"} with \code{length(postc) = 0} and
#' \code{full = FALSE}, \code{cde}, \code{pnde}, \code{tnde}, \code{pnie}, \code{tnie} and 
#' \code{te} are estimated.
#' 
#' When \code{model = "msm" or "gformula"} with \code{length(postc) != 0} and \code{full = TRUE}, 
#' \code{pnde}, \code{tnde}, \code{pnie}, \code{tnie},  \code{intref}, \code{intmed}, 
#' \code{intref(prop)}, \code{intmed(prop)}, \code{pie(prop)}, \code{pm}, \code{int} and 
#' \code{pe} are replaced by their randomized analogues \code{rpnde}, \code{rtnde}, 
#' \code{rpnie}, \code{rtnie}, \code{rintref}, \code{rintmed}, \code{rintref(prop)}, 
#' \code{rintmed(prop)}, \code{rpie(prop)}, \code{rpm}, \code{rint} and \code{rpe}.
#' 
#' When \code{model = "msm" or "gformula"} with \code{length(postc) != 0} and
#' \code{full = FALSE}, \code{cde}, \code{rpnde}, \code{rtnde}, \code{rpnie}, \code{rtnie} and 
#' \code{te} are estimated.
#' 
#' When \code{model = "iorw"} and \code{full = TRUE}, \code{te}, \code{pnde}, \code{tnie} and
#' \code{pm} are estimated. When \code{model = "iorw"} and \code{full = FALSE}, \code{te}, 
#' \code{pnde} and \code{tnie} are estimated.  
#' 
#' \emph{Categorical, Count or Survival Outcome}
#' 
#' When \code{model = "rb", "wb", "ne", "msm" or "gformula"} with \code{length(postc) = 0} and 
#' \code{full = TRUE}, \code{RRcde} (\code{cde} risk rario or rate ratio), \code{RRpnde} 
#' (\code{pnde} risk rario or rate ratio), \code{RRtnde} (\code{tnde} risk rario or rate ratio), 
#' \code{RRpnie} (\code{pnie} risk rario or rate ratio), \code{RRtnie} (\code{tnie} risk rario 
#' or rate ratio), \code{RRte} (\code{te} risk rario or rate ratio), \code{ERRcde} 
#' (excess risk rario or rate ratio due to \code{cde}), \code{ERRintref} (excess risk rario or 
#' rate ratio due to \code{intref}), \code{ERRintmed} (excess risk rario or rate ratio due to 
#' \code{intmed}), \code{ERRpie} (excess risk rario or rate ratio due to \code{pie}), 
#' \code{ERRcde(prop)} (proportion \code{ERRcde}), \code{ERRintref(prop)} (proportion 
#' \code{ERRintref}), \code{ERRintmed(prop)} (proportion \code{ERRintmed}), \code{ERRpie(prop)} 
#' (proportion \code{ERRpie}), \code{pm}, \code{int} (proportion attributable to 
#' interaction) and \code{pe} (proportion eliminated) are estimated. 
#' 
#' When \code{model = "rb", "wb", "ne", "msm" or "gformula"} with \code{length(postc) = 0} and
#' \code{full = FALSE}, \code{RRcde}, \code{RRpnde}, \code{RRtnde}, \code{RRpnie}, \code{RRtnie} 
#' and \code{RRte} are estimated.
#' 
#' When \code{model = "msm" or "gformula"} with \code{length(postc) != 0} and \code{full = TRUE}, 
#' \code{RRpnde}, \code{RRtnde}, \code{RRpnie}, \code{RRtnie},  \code{ERRintref}, 
#' \code{ERRintmed}, \code{ERRpie}, \code{ERRintref(prop)}, \code{ERRintmed(prop)}, 
#' \code{ERRpie(prop)}, \code{pm}, \code{int} and \code{pe} are 
#' replaced by their randomized analogues \code{rRRpnde}, \code{rRRtnde}, \code{rRRpnie}, 
#' \code{rRRtnie}, \code{rERRintref}, \code{rERRintmed}, \code{rERRpie}, 
#' \code{rERRintref(prop)}, \code{rERRintmed(prop)}, \code{rERRpie(prop)}, 
#' \code{rpm}, \code{rint} and \code{rpe}.
#' 
#' When \code{model = "msm" or "gformula"} with \code{length(postc) != 0} and
#' \code{full = FALSE}, \code{RRcde}, \code{rRRpnde}, \code{rRRtnde}, \code{rRRpnie}, 
#' \code{rRRtnie} and \code{RRte} are estimated.
#' 
#' When \code{model = "iorw"} and \code{full = TRUE}, \code{RRte}, \code{RRpnde}, 
#' \code{RRtnie} and \code{pm} are estimated. When \code{model = "iorw"} and 
#' \code{full = FALSE}, \code{RRte}, \code{RRpnde} and \code{RRtnie} are estimated.  
#' 
#' @return
#' An object of class 'cmest' is returned:
#' \item{call}{the function call,}
#' \item{data}{the dataset,}
#' \item{methods}{a list of methods used which may include \code{model}, \code{full}, 
#' \code{casecontrol}, \code{yprevalence}, \code{yrare}, \code{estimation}, \code{inference} 
#' and \code{nboot},}
#' \item{variables}{a list of variables used which may include \code{outcome}, \code{event}, 
#' \code{exposure}, \code{mediator}, \code{EMint}, \code{basec} and \code{postc},}
#' \item{reg.input}{a list of regressions input,}
#' \item{multimp}{a list of arguments used for multiple imputation,}
#' \item{ref}{reference values used,}
#' \item{reg.output}{a list of regressions output,}
#' \item{effect.pe}{point estimates of causal effects,}
#' \item{effect.se}{standard errors of causal effects,}
#' \item{effect.ci.low}{the lower limits of confidence intervals of causal effects,}
#' \item{effect.ci.high}{the higher limits of confidence intervals of causal effects,}
#' \item{effect.pe}{point estimates of causal effects,}
#' \item{effect.pval}{p-values of causal effects,}
#' ...
#'
#' @seealso \code{\link{ggcmest}}, \code{\link{cmdag}}, 
#' \code{\link{cmsens}}.
#'
#' @references
#' Valeri L, VanderWeele TJ (2013). Mediation analysis allowing for
#' exposure-mediator interactions and causal interpretation: theoretical assumptions and
#' implementation with SAS and SPSS macros. Psychological Methods. 18(2): 137 - 150.
#' 
#' VanderWeele TJ, Vansteelandt S (2014). Mediation analysis with multiple mediators.
#' Epidemiologic Methods. 2(1): 95 - 115.
#' 
#' Tchetgen Tchetgen EJ (2013). Inverse odds ratio-weighted estimation for causal
#' mediation analysis. Statistics in medicine. 32: 4567 - 4580.
#' 
#' Nguyen QC, Osypuk TL, Schmidt NM, Glymour MM, Tchetgen Tchetgen EJ. Practical guidance
#' for conducting mediation analysis with multiple mediators using inverse odds ratio
#' weighting (2015). American Journal of Epidemiology. 181(5): 349 - 356.
#' 
#' VanderWeele TJ, Tchetgen Tchetgen EJ (2017). Mediation analysis with time varying
#' exposures and mediators. Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology). 79(3): 917 - 938.
#' 
#' Lin SH, Young J, Logan R, Tchetgen Tchetgen EJ, VanderWeele TJ (2017). Parametric
#' mediational g-formula approach to mediation analysis with time-varying exposures,
#' mediators, and confounders. Epidemiology. 28: 266 - 274.
#' 
#' Vansteelandt S, Bekaert M, Lange T. (2012). Imputation Strategies for the Estimation 
#' of Natural Direct and Indirect Effects. Epidemiologic Methods. 1(1): 131 - 158.
#' 
#' Steen J, Loeys T, Moerkerke B, Vansteelandt S (2017). Medflex: an R package for
#' flexible mediation analysis using natural effect models. Journal of Statistical
#' Software. 76(11).
#' 
#' VanderWeele TJ. A unification of mediation and interaction: a 4-way decomposition (2014). 
#' Epidemiology. 25(5): 749 - 61.
#' 
#' Imai K, Keele L, Tingley D. A general approach to causal mediation analysis (2010).
#' Psychological Methods. 15(4): 309 - 334.
#' 
#' Schomaker M, Heumann C. Bootstrap inference when using multiple 
#' imputation (2018). Statistics in Medicine. 37(14): 2252 - 2266. 
#'
#' @examples
#' 
#' \dontrun{
#' library(CMAverse)
#' 
#' # single-mediator case with rb
#' exp1 <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = "M2", basec = c("C1", "C2"), 
#' EMint = TRUE, mreg = list("multinomial"), yreg = "linear", 
#' astar = 0, a = 1, mval = list("M2_0"), estimation = "paramfunc", 
#' inference = "delta")
#' summary(exp1)
#' 
#' # multiple-mediator case with rb
#' # 10 boots are used for illustration
#' exp2 <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), basec = c("C1", "C2"), 
#' EMint = TRUE, mreg = list("logistic", "multinomial"), 
#' yreg = "linear", astar = 0, a = 1, mval = list(0, "M2_0"), 
#' estimation = "imputation", inference = "bootstrap", nboot = 10)
#' 
#' # multiple-mediator case with ne
#' # 10 boots are used for illustration
#' exp3 <- cmest(data = cma2020, model = "ne", outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), basec = c("C1", "C2"), 
#' yreg = glm(contY ~ A + M1 + M2 + A*M1 + A*M2 + C1 + C2, family = gaussian, data = cma2020), 
#' astar = 0, a = 1, mval = list(0, "M2_0"), estimation = "imputation", 
#' inference = "bootstrap", nboot = 10)
#' 
#' # case control study with msm
#' # 10 boots are used for illustration
#' exp4 <- cmest(data = cma2020, model = "msm", casecontrol = TRUE, 
#' yrare = TRUE, outcome = "binY", exposure = "A", 
#' mediator = c("M1", "M2"), EMint = TRUE, basec = c("C1", "C2"), yreg = "logistic", 
#' ereg = "logistic", mreg = list(glm(M1 ~ A, family = binomial, 
#' data = cma2020), nnet::multinom(M2 ~ A + M1, data = cma2020, trace = FALSE)), 
#' wmreg = list(glm(M1 ~ A + C1 + C2, family = binomial, data = cma2020), 
#' nnet::multinom(M2 ~ A + M1 + C1 + C2, data = cma2020, trace = FALSE)), astar = 0, a = 1, 
#' mval = list(0, "M2_0"), estimation = "imputation", 
#' inference = "bootstrap", nboot = 10)
#' }
#' 
#' @importFrom stats glm binomial poisson as.formula gaussian quasipoisson model.frame printCoefmat 
#' family sd coef vcov sigma predict rbinom rmultinom rnorm rgamma rpois weighted.mean 
#' model.matrix getCall quantile qnorm pnorm lm cov formula update
#' @importFrom nnet multinom
#' @importFrom MASS polr glm.nb gamma.shape rnegbin
#' @importFrom survival survreg coxph
#' @importFrom survey svyglm svydesign svysurvreg svycoxph as.svrepdesign withReplicates
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
#'
#' @export

cmest <- function(data = NULL, model = NULL,
                  full = TRUE, casecontrol = FALSE, yrare = NULL, yprevalence = NULL,
                  estimation = "imputation", inference = "bootstrap",
                  outcome = NULL, event = NULL,
                  exposure = NULL, mediator = NULL, EMint = NULL, basec = NULL, postc = NULL,
                  yreg = NULL, mreg = NULL, wmreg = NULL, ereg = NULL, postcreg = NULL,
                  astar = 0, a = 1, mval = NULL, yref = NULL, basecval = NULL,
                  nboot = 200, nRep = 5, multimp = FALSE, ...) {
  
  cl <- match.call()
  
  n <- nrow(data)
  # output list
  out <- list(call = cl)
  
  ###################################################################################################
  #################################Argument Restrictions And Warnings################################
  ###################################################################################################
  
  # data
  if (is.null(data)) stop("Unspecified data")
  data <- as.data.frame(data)
  out$data <- data
  
  # model
  if (!model %in% c("rb", "wb", "iorw", "ne", "gformula", "msm")) stop("Select model from 'rb', 'wb', 'iorw', 'ne', 'gformula', 'msm'")
  out$methods$model <- model
  
  # full
  if (!is.logical(full)) stop("full should be TRUE or FALSE")
  out$methods$full <- full
  
  # casecontrol, yrare, yprevalence
  if (!is.logical(casecontrol)) stop("casecontrol should be TRUE or FALSE")
  out$methods$casecontrol <- casecontrol
  if (!casecontrol) {
    if (!is.null(yrare)) warning("When casecontrol is FALSE, yrare is ignored")
    if (!is.null(yprevalence)) warning("When casecontrol is FALSE, yprevalence is ignored ")
  } else {
    if (is.null(yprevalence) && !yrare == TRUE) stop("When casecontrol is TRUE, specify yprevalence or set yrare to be TRUE")
    if (!is.null(yprevalence)) {
      if (!is.numeric(yprevalence)) stop("yprevalence should be numeric")
      if (!is.null(yrare)) out$methods$yrare <- yrare
      out$methods$yprevalence <- yprevalence
    } else out$methods$yrare <- yrare
  }
  
  # outcome
  if (length(outcome) == 0) stop("Unspecified outcome")
  if (length(outcome) > 1) stop("length(outcome) > 1")
  out$variables$outcome <- outcome
  # event
  if (is.character(yreg) && !is.null(event)) out$variables$event <- event
  if (!is.character(yreg) && !is.null(event)) warning("event is ignored when yreg is not character")
  # exposure
  if (length(exposure) == 0) stop("Unspecified exposure")
  if (length(exposure) > 1) stop("length(exposure) > 1")
  out$variables$exposure <- exposure
  # mediator
  if (length(mediator) == 0) stop("Unspecified mediator")
  out$variables$mediator <- mediator
  # EMint
  if (is.character(yreg) && model != "iorw") {
    if (!is.logical(EMint)) stop("EMint should be TRUE or FALSE")
    out$variables$EMint <- EMint
  } else if (!is.null(EMint)) warning("EMint is ignored")
  if (!is.null(basec)) out$variables$basec <- basec
  # postc
  if (length(postc) != 0) {
    if (!model %in% c("msm", "gformula")) stop("When postc is not empty, select model from 'msm' and 'gformula'")
    out$variables$postc <- postc
  }
  
  #regs
  if (model == "rb") out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "wb") out$reg.input <- list(yreg = yreg, ereg = ereg)
  if (model == "gformula") out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "msm") out$reg.input <- list(yreg = yreg, mreg = mreg, ereg = ereg, wmreg = wmreg)
  if (length(postc) != 0 && model %in% c("msm", "gformula")) out$reg.input$postcreg <- postcreg
  if (model == "iorw") out$reg.input <- list(yreg = yreg, ereg = ereg)
  if (model == "ne") out$reg.input <- list(yreg = yreg)
  
  # a, astar
  if (is.null(a) | is.null(astar)) stop("Unspecified a or astar")
  # mval
  if (model %in% c("rb", "wb", "msm", "gformula", "ne")) {
    if (!is.list(mval)) stop("mval should be a list")
    if (length(mval) != length(mediator)) stop("length(mval) != length(mediator)")
    # for rb, wb, ne, msm, gformula, a reference value for each mediator is required for estimating cd
    for (i in 1:length(mval)) if (is.null(mval[[i]])) stop(paste0("Unspecified mval[[", i, "]]"))
  } else if (!is.null(mval)) warning("mval is ignored when model = 'iorw'")
  # basecval
  if (!(model == "rb" && estimation == "paramfunc" && length(basec) != 0) && !is.null(basecval)) warning("basecval is ignored")
  
  # estimation and inference
  if (model == "rb" && !estimation %in% c("paramfunc", "imputation")) stop("When model = 'rb', select estimation from 'paramfunc', 'imputation'")
  if (model != "rb" && !estimation == "imputation") stop("Use estimation = 'imputation'")
  if (estimation == "paramfunc") {
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson",
                     "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) stop(
                       paste0("When estimation = 'paramfunc', select yreg from 'linear', 'logistic',
                 'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'"))
    if (model == "rb" && !mreg[[1]] %in% c("linear", "logistic", "multinomial")) stop(
      "When model = rb and estimation = paramfunc, select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
  }
  out$methods$estimation <- estimation
  if (estimation == "imputation" && inference == "delta") stop("Use inference = 'bootstrap' when estimation = 'imputation'")
  out$methods$inference <- inference
  if (inference == "bootstrap") {
    if (!is.numeric(nboot)) stop("nboot should be numeric")
    out$methods$nboot <- nboot
  }
  
  if (model == "ne") out$methods$nRep <- nRep
  
  # multimp
  if (!is.logical(multimp)) stop("multimp should be TRUE or FALSE")
  # args_mice
  args_mice <- list(...)
  args_mice$print <- FALSE
  
  ###################################################################################################
  ##########################################Run Regressions##########################################
  ###################################################################################################
  
  # Y: outcome; M: mediator A: exposure; C: baseline confounder; L: post-exposure confounder
  # the variable used to calculate weights is required to be categorical
  
  ####################################Exposure Regression For Weighting##############################
  
  # for wb and msm, the exposure regression is required for calculating weights if basec is not empty, w_{a,i}=P(A=A_i)/P(A=A_i|C_i)
  # for iorw, the exposure regression is required for calculating weights, w_{a,i}=P(A=0|M_i,C_i)/P(A=A_i|M_i,C_i)
  if ((model %in% c("wb", "msm") && length(basec) > 0) | model == "iorw") {
    if (is.null(ereg)) stop("ereg is required when model is 'wb' or 'msm' with length(basec) > 0 or model is 'iorw'")
    if (is.character(ereg)) {
      # fit glm with family = poisson() rather than family = binomial("log") for "loglinear"
      if (ereg == "loglinear" && length(unique(data[, exposure])) != 2) stop("When ereg is 'loglinear', exposure should be binary")
      if (!ereg %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop("Select character ereg from 'logistic', 'loglinear', 'multinomial', 'ordinal'")
      exposure_formula <- switch((model == "iorw") + 1, "1" = paste0(exposure, "~", paste0(basec, collapse = "+")),
                                 "2" = paste0(exposure, "~", paste0(c(mediator, basec), collapse = "+")))
      switch(ereg,
             logistic = ereg <- eval(bquote(glm(.(as.formula(exposure_formula)), family = binomial(), data = .(data)))),
             loglinear = ereg <- eval(bquote(glm(.(as.formula(exposure_formula)), family = poisson(), data = .(data)))),
             multinomial = ereg <- eval(bquote(nnet::multinom(.(as.formula(exposure_formula)), data = .(data), trace = FALSE))),
             ordinal = ereg <- eval(bquote(polr(.(as.formula(exposure_formula)), data = .(data)))))
    }
  } else {
    if (!is.null(ereg)) warning("ereg is ignored when model is 'wb' or 'msm' with length(basec) = 0 or model is 'rb', 'ne' or 'gformula'")
    ereg <- NULL
  }
  
  ####################################Mediator Regression For Weighting##############################
  
  # for msm, a mediator regression for weighting is required for each mediator
  if (model == "msm") {
    if (!is.list(wmreg)) stop("wmreg should be a list")
    if (length(wmreg) != length(mediator)) stop("length(wmreg) != length(mediator)")
    for (p in 1:length(wmreg)) {
      if (is.null(wmreg[[p]])) stop(paste0("Unspecified wmreg[[", p, "]]"))
      if (is.character(wmreg[[p]])) {
        if (wmreg[[p]] == "loglinear" && length(unique(data[, mediator[p]])) != 2) stop(paste0("When wmreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!wmreg[[p]] %in% c("logistic", "loglinear", "multinomial", "ordinal")) stop(paste0("Select character wmreg[[", p, "]] from 'logistic', 'loglinear', 'multinomial', 'ordinal'"))
        # w_{m_p,i}=P(M_p=M_{p,i}|A=A_i)/P(M_p=M_{p,i}|A=A_i,C=C_i,L=L_i,M_1=M_{1,i},...,M_{p-1}=M_{p-1,i})
        wmreg_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec, postc), collapse = "+"))
        # regression for the denominator of w_{m_p,i}
        switch(wmreg[[p]],
               logistic = wmreg[[p]] <- eval(bquote(glm(.(as.formula(wmreg_formula)), family = binomial(), data = .(data)))),
               loglinear = wmreg[[p]] <- eval(bquote(glm(.(as.formula(wmreg_formula)), family = poisson(), data = .(data)))),
               multinomial = wmreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(wmreg_formula)), data = .(data), trace = FALSE))),
               ordinal = wmreg[[p]] <- eval(bquote(polr(.(as.formula(wmreg_formula)), data = .(data)))))
      }
    }
  } else {
    if (!is.null(wmreg)) warning("wmreg is ignored when model is not 'msm'")
    wmreg <- NULL
  }
  
  ###########################################Mediator Regression#####################################
  
  # for rb, msm and gformula, a mediator regression is required for each mediator
  if (model %in% c("rb", "msm", "gformula")) {
    
    if (!is.list(mreg)) stop("mreg should be a list")
    if (length(mreg) != length(mediator)) stop("length(mreg) != length(mediator)")
    
    for (p in 1:length(mreg)) {
      if (is.null(mreg[[p]])) stop(paste0("Unspecified mreg[[", p, "]]"))
      if (is.character(mreg[[p]])) {
        if (mreg[[p]] == "loglinear" && length(unique(data[, mediator[p]])) != 2) stop(paste0("When mreg[[", p, "]] is 'loglinear', mediator[[", p, "]] should be binary"))
        if (!mreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                              "negbin", "multinomial", "ordinal")) stop(
                                paste0("Select character mreg[[", p, "]] from 'linear', 'logistic',
'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
        
        # for rb, regress each mediator on A, C and previous mediators
        # for msm, regress each mediator on A and previous mediators
        # for gformula, regress each mediator on A, C, L and previous mediators
        switch(model,
               rb = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec), collapse = "+")),
               msm = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)]), collapse = "+")),
               gformula = mediator_formula <- paste0(mediator[p], "~", paste(c(exposure, mediator[0:(p-1)], basec, postc), collapse = "+")))
        switch(mreg[[p]],
               linear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = gaussian(), data = .(data)))),
               logistic = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = binomial(), data = .(data)))),
               loglinear = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
               poisson = mreg[[p]]  <- eval(bquote(glm(.(as.formula(mediator_formula)), family = poisson(), data = .(data)))),
               quasipoisson = mreg[[p]] <- eval(bquote(glm(.(as.formula(mediator_formula)), family = quasipoisson(), data = .(data)))),
               negbin = mreg[[p]] <- eval(bquote(glm.nb(.(as.formula(mediator_formula)), data = .(data)))),
               multinomial = mreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(mediator_formula)), data = .(data), trace = FALSE))),
               ordinal = mreg[[p]] <- eval(bquote(polr(.(as.formula(mediator_formula)), data = .(data)))))
      }
    }
  } else {
    if (!is.null(mreg)) warning("mreg is ignored when model is 'wb', 'iorw' or 'ne'")
    mreg <- NULL
  }
  
  ####################################Post-exposure Confounder Regression############################
  
  # for gformula, a regression is required for each post-exposure confounder
  if (model == "gformula" && length(postc) > 0) {
    if (!is.list(postcreg)) stop("postcreg should be a list")
    if (length(postcreg) != length(postc)) stop("length(postcreg) != length(postc)")
    for (p in 1:length(postcreg)) {
      if (is.null(postcreg[[p]])) stop(paste0("Unspecified postcreg[[", p, "]]"))
      if (is.character(postcreg[[p]])) {
        if (postcreg[[p]] == "loglinear" && length(unique(data[, postc[p]])) != 2) stop(paste0("When postcreg[[", p, "]] is 'loglinear', postc[[", p, "]] should be binary"))
        if (!postcreg[[p]] %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                                  "negbin", "multinomial", "ordinal"))  stop(
                                    paste0("Select character postcreg[[", p, "]] from 'linear', 'logistic',
                                           'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal'"))
        # regress each post-exposure confounder on A, C and previous post-exposure confounders
        postc_formula <- paste0(postc[p], "~", paste(c(exposure, basec, postc[0:(p-1)]), collapse = "+"))
        switch(postcreg[[p]],
               linear = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = gaussian(), data = .(data)))),
               logistic = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = binomial(), data = .(data)))),
               loglinear = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = poisson(), data = .(data)))),
               poisson = postcreg[[p]]  <- eval(bquote(glm(.(as.formula(postc_formula)), family = poisson(), data = .(data)))),
               quasipoisson = postcreg[[p]] <- eval(bquote(glm(.(as.formula(postc_formula)), family = quasipoisson(), data = .(data)))),
               negbin = postcreg[[p]] <- eval(bquote(glm.nb(.(as.formula(postc_formula)), data = .(data)))),
               multinomial = postcreg[[p]] <- eval(bquote(nnet::multinom(.(as.formula(postc_formula)), data = .(data), trace = FALSE))),
               ordinal = postcreg[[p]] <- eval(bquote(polr(.(as.formula(postc_formula)), data = .(data)))))
      }
    }
  } else {
    if (!is.null(postcreg)) warning("postcreg is ignored when model is not 'gformula' or length(postc) = 0")
    postcreg <- NULL
  }
  
  ###########################################Outcome Regression######################################
  
  if (is.null(yreg)) stop("yreg is required")
  if (is.character(yreg)) {
    if (yreg == "loglinear" && length(unique(data[, outcome])) != 2) stop("When yreg is 'loglinear', outcome should be binary")
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson", "quasipoisson",
                     "negbin", "multinomial", "ordinal", "coxph", "aft_exp",
                     "aft_weibull")) stop(
                       paste0("Select character yreg from 'linear', 'logistic',
                              'loglinear', 'poisson', 'quasipoisson', 'negbin', 'multinomial', 'ordinal',
                              'coxph', 'aft_exp', 'aft_weibull'"))
    if (model != "iorw") {
      out$variables$EMint <- EMint
      int.terms <- switch((model != "iorw" && EMint) + 1, "1" = NULL, "2" = paste(exposure, mediator, sep = "*"))
    }
    
    # for rb, wb and ne, regress Y on A, M and C
    # for iorw, regress Y on A and C
    # for msm, regress Y on A and M
    # for gformula, regress Y on A, M, C and L
    switch(model,
           rb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           wb = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           ne = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec), collapse = "+")),
           iorw = outcome_formula <- paste0(outcome, "~", paste(c(exposure, basec), collapse = "+")),
           msm = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms), collapse = "+")),
           gformula = outcome_formula <- paste0(outcome, "~", paste(c(exposure, mediator, int.terms, basec, postc), collapse = "+")))
    
    if (yreg %in% c("coxph","aft_exp","aft_weibull")) {
      if (!is.null(event)) {
        out$variables$event <- event
        outcome_formula <- paste(paste0("Surv(", outcome, ", ", event, ")"),
                                 strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
      } else outcome_formula <- paste(paste0("Surv(", outcome, ")"),
                                      strsplit(outcome_formula, split = "~")[[1]][2], sep = " ~ ")
    }
    
    # run outcome regression
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
           negbin = yreg <- eval(bquote(glm.nb(formula = .(as.formula(outcome_formula)),
                                               data = .(data)))),
           multinomial = yreg <- eval(bquote(nnet::multinom(formula = .(as.formula(outcome_formula)),
                                                            data = .(data), trace = FALSE))),
           ordinal = yreg <- eval(bquote(polr(formula = .(as.formula(outcome_formula)),
                                              data = .(data)))),
           coxph = yreg <- eval(bquote(coxph(formula = .(as.formula(outcome_formula)),
                                             data = .(data)))),
           aft_exp = yreg <- eval(bquote(survreg(formula = .(as.formula(outcome_formula)),
                                                 dist = "exponential", data = .(data)))),
           aft_weibull = yreg <- eval(bquote(survreg(formula = .(as.formula(outcome_formula)),
                                                     dist = "weibull", data = .(data)))))
  }
  
  # for delta method inference, use survey regressions for yreg and mreg when weights are applied
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    if (inherits(yreg, "glm")) yreg <- eval(bquote(svyglm(formula = .(formula(yreg)), family = .(family(yreg)),
                                                          design = svydesign(~ 1, data = .(data)))))
    if (inherits(yreg, "survreg")) yreg <- eval(bquote(svysurvreg(formula = .(formula(yreg)),
                                                                  design = svydesign(~ 1, weights = .(model.frame(yreg)$'(weights)'), data = .(data)))))
    if (inherits(yreg, "coxph")) yreg <- eval(bquote(svycoxph(formula = .(formula(yreg)),
                                                              design = svydesign(~ 1, data = .(data)))))
  }
  if (inference == "delta" && casecontrol && !is.null(yprevalence)) {
    if (inherits(mreg[[1]], "glm")) mreg[[1]] <- eval(bquote(svyglm(formula = .(formula(mreg[[1]])), 
                                                                    family = .(family(mreg[[1]])),
                                                                    design = svydesign(~ 1, data = .(data)))))
    if (inherits(mreg[[1]], "multinom")) mreg[[1]] <- eval(bquote(svymultinom(formula = .(formula(mreg[[1]])), 
                                                                              data = .(data))))
  }
  
  ###################################################################################################
  ############################################Estimation and Inference###############################
  ###################################################################################################
  
  # add a progress bar for bootstrap inference
  if (inference == "bootstrap") {
    env <- environment()
    counter <- 0
    progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
  }
  
  environment(estinf) <- environment()
  out <- c(out, estinf())
  class(out) <- "cmest"
  return(out)
}

#' @describeIn cmest Print the results of cmest nicely
#' @export
print.cmest <- function(x, ...) {
  # print regression models used
  if (!x$multimp$multimp) {
    regnames <- names(x$reg.output)
    for (name in regnames) {
      if (name == "yreg") {
        if (!is.null(x$reg.output$yreg)) {
          cat("# Outcome Regression: \n")
          x$reg.output$yreg$call <- update(x$reg.output$yreg,data=getCall(x$reg.output$yreg)$data,
                                           weights=getCall(x$reg.output$yreg)$weights, evaluate = FALSE)
          print(x$reg.output$yreg)
        }
      }
      if (name == "ereg") {
        if (!is.null(x$reg.output$ereg)) {
          cat("# Exposure Regression for Weighting: \n")
          x$reg.output$ereg$call <- update(x$reg.output$ereg,data=getCall(x$reg.output$ereg)$data,
                                           weights=getCall(x$reg.output$ereg)$weights, evaluate = FALSE)
          print(x$reg.output$ereg)
        }
      }
      if (name == "mreg") {
        if (!is.null(x$reg.output$mreg)) {
          cat("# Mediator Regression: \n")
          for (i in 1:length(x$reg.output$mreg)) {
            x$reg.output$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[i]], 
                                                              data=getCall(x$reg.output$mreg[[.(i)]])$data, 
                                                              weights=getCall(x$reg.output$mreg[[.(i)]])$weights,
                                                              evaluate = FALSE)))
            print(x$reg.output$mreg[[i]])
            if (i < length(x$reg.output$mreg)) cat("\n")
          }
        }
      }
      if (name == "wmreg") {
        if (!is.null(x$reg.output$wmreg)) {
          cat("# Mediator Regression for Weighting: \n")
          for (i in 1:length(x$reg.output$wmreg)) {
            x$reg.output$wmreg[[i]]$call <- eval(bquote(update(x$reg.output$wmreg[[i]], 
                                                               data=getCall(x$reg.output$wmreg[[.(i)]])$data, 
                                                               weights=getCall(x$reg.output$wmreg[[.(i)]])$weights,
                                                               evaluate = FALSE)))
            print(x$reg.output$wmreg[[i]])
            if (i < length(x$reg.output$wmreg)) cat("\n")
          }
        }
      }
      if (name == "postcreg") {
        if (!is.null(x$reg.output$postcreg)) {
          cat("# Post-exposure Confounder Regression: \n")
          for (i in 1:length(x$reg.output$postcreg)) {
            x$reg.output$postcreg[[i]]$call <- eval(bquote(update(x$reg.output$postcreg[[i]], 
                                                                  data=getCall(x$reg.output$postcreg[[.(i)]])$data, 
                                                                  weights=getCall(x$reg.output$postcreg[[.(i)]])$weights,
                                                                  evaluate = FALSE)))
            print(x$reg.output$postcreg[[i]])
            if (i < length(x$reg.output$postcreg)) cat("\n")
          }
        }
      }
      cat("\n")
    }
  } else {
    for (m in 1:length(x$reg.output)){ 
      cat(paste("# Regressions with Imputed Dataset", m, "\n\n"))
      regnames <- names(x$reg.output[[m]])
      for (name in regnames) {
        if (name == "yreg") {
          if (!is.null(x$reg.output[[m]]$yreg)) {
            cat("## Outcome Regression: \n")
            x$reg.output[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,
                                                              data=getCall(x$reg.output[[.(m)]]$yreg)$data,
                                                              weights=getCall(x$reg.output[[.(m)]]$yreg)$weights, 
                                                              evaluate = FALSE)))
            print(x$reg.output[[m]]$yreg)
          }
        }
        if (name == "ereg") {
          if (!is.null(x$reg.output[[m]]$ereg)) {
            cat("## Exposure Regression for Weighting: \n")
            x$reg.output[[m]]$ereg$call <- eval(bquote(update(x$reg.output[[.(m)]]$ereg,
                                                              data=getCall(x$reg.output[[.(m)]]$ereg)$data,
                                                              weights=getCall(x$reg.output[[.(m)]]$ereg)$weights, 
                                                              evaluate = FALSE)))
            print(x$reg.output[[m]]$ereg)
          }
        }
        if (name == "mreg") {
          if (!is.null(x$reg.output[[m]]$mreg)) {
            cat("## Mediator Regression: \n")
            for (i in 1:length(x$reg.output[[m]]$mreg)) {
              x$reg.output[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[i]], 
                                                                     data=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$data, 
                                                                     weights=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$weights, 
                                                                     evaluate = FALSE)))
              print(x$reg.output[[m]]$mreg[[i]])
              if (i < length(x$reg.output[[m]]$mreg)) cat("\n")
            }
          }
        }
        if (name == "wmreg") {
          if (!is.null(x$reg.output[[m]]$wmreg)) {
            cat("## Mediator Regression for Weighting: \n")
            for (i in 1:length(x$reg.output[[m]]$wmreg)) {
              x$reg.output[[m]]$wmreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmreg[[i]], 
                                                                      data=getCall(x$reg.output[[.(m)]]$wmreg[[.(i)]])$data, 
                                                                      weights=getCall(x$reg.output[[.(m)]]$wmreg[[.(i)]])$weights, 
                                                                      evaluate = FALSE)))
              print(x$reg.output[[m]]$wmreg[[i]])
              if (i < length(x$reg.output[[m]]$wmreg)) cat("\n")
            }
          }
        }
        if (name == "postcreg") {
          if (!is.null(x$reg.output[[m]]$postcreg)) {
            cat("## Post-exposure Confounder Regression: \n")
            for (i in 1:length(x$reg.output[[m]]$postcreg)) {
              x$reg.output[[m]]$postcreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$postcreg[[i]], 
                                                                         data=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$data, 
                                                                         weights=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$weights, 
                                                                         evaluate = FALSE)))
              print(x$reg.output[[m]]$postcreg[[i]])
              if (i < length(x$reg.output[[m]]$postcreg)) cat("\n")
            }
          }
        }
        cat("\n")
      }
    }
  }
  # print causal mediation analysis results
  if (x$methods$model == "rb") model_str <- "Regression-based Approach"
  if (x$methods$model == "wb") model_str <- "Weighting-based Approach"
  if (x$methods$model == "ne") model_str <- "Natural Effect Model"
  if (x$methods$model == "iorw") model_str <- "Inverse Odds Ratio Weighting Approach"
  if (x$methods$model == "msm") model_str <- "Marginal Structural Model"
  if (x$methods$model == "gformula") model_str <- "G-formula Approach"
  if (x$multimp$multimp) model_str <- paste(model_str, "with Multiple Imputation")
  if (x$methods$estimation == "paramfunc") est_str <- "Closed-form parameter function estimation"
  if (x$methods$estimation == "imputation") est_str <- "Direct counterfactual imputation estimation"
  if (x$methods$inference == "delta") inf_str <- "delta method standard errors, confidence intervals and p-values"
  if (x$methods$inference == "bootstrap") inf_str <- "bootstrap standard errors, percentile confidence intervals and p-values"
  if (x$methods$model != "ne" && (x$methods$casecontrol)) cat("# Causal Mediation Analysis for a Case Control Study via the ")
  if (!(x$methods$model != "ne" && (x$methods$casecontrol))) cat("# Causal Mediation Analysis via the ")
  cat(model_str)
  cat("\n \n")
  cat(est_str)
  cat(paste(" with \n", inf_str, "\n \n"))
  print(x$effect.pe)
  cat("\n")
  cat("Reference values: \n")
  print(x$ref)
}


#' @describeIn cmest Summarize the results of cmest nicely
#' @export
summary.cmest <- function(object, ...) {
  # summarize regression models used
  regsumm <- list()  
  if (!object$multimp$multimp) {
    regnames <- names(object$reg.output)
    for (name in regnames) {
      if (name == "yreg") {
        if (!is.null(object$reg.output$yreg)) {
          regsumm$yreg <- summary(object$reg.output$yreg)
        }
      }
      if (name == "ereg") {
        if (!is.null(object$reg.output$ereg)) {
          regsumm$ereg <- summary(object$reg.output$ereg)
        }
      }
      if (name == "mreg") {
        if (!is.null(object$reg.output$mreg)) {
          regsumm$mreg <- lapply(1:length(object$reg.output$mreg), function(i)
            summary(object$reg.output$mreg[[i]]))
        }
      }
      if (name == "wmreg") {
        if (!is.null(object$reg.output$wmreg)) {
          regsumm$wmreg <- lapply(1:length(object$reg.output$wmreg), function(i)
            summary(object$reg.output$wmreg[[i]]))
        }
      }
      if (name == "postcreg") {
        if (!is.null(object$reg.output$postcreg)) {
          regsumm$postcreg <- lapply(1:length(object$reg.output$postcreg), function(i)
            summary(object$reg.output$postcreg[[i]]))
        }
      }
    }
  } else {
    for (m in 1:length(object$reg.output)){ 
      regsumm[[m]] <- list()
      regnames <- names(object$reg.output[[m]])
      for (name in regnames) {
        if (name == "yreg") {
          if (!is.null(object$reg.output[[m]]$yreg)) {
            regsumm[[m]]$yreg <- summary(object$reg.output[[m]]$yreg)
          }
        }
        if (name == "ereg") {
          if (!is.null(object$reg.output[[m]]$ereg)) {
            regsumm[[m]]$ereg <- summary(object$reg.output[[m]]$ereg)
          }
        }
        if (name == "mreg") {
          if (!is.null(object$reg.output[[m]]$mreg)) {
            regsumm[[m]]$mreg <- lapply(1:length(object$reg.output[[m]]$mreg), function(i)
              summary(object$reg.output[[m]]$mreg[[i]]))
          }
        }
        if (name == "wmreg") {
          if (!is.null(object$reg.output[[m]]$wmreg)) {
            regsumm[[m]]$wmreg <- lapply(1:length(object$reg.output[[m]]$wmreg), function(i)
              summary(object$reg.output[[m]]$wmreg[[i]]))
          }
        }
        if (name == "postcreg") {
          if (!is.null(object$reg.output[[m]]$postcreg)) {
            regsumm[[m]]$postcreg <- lapply(1:length(object$reg.output[[m]]$postcreg), function(i)
              summary(object$reg.output[[m]]$postcreg[[i]]))
          }
        }
      }
    }
  }
  # summarize causal mediation analysis results
  summarydf <- data.frame(object$effect.pe, object$effect.se, object$effect.ci.low, 
                          object$effect.ci.high, object$effect.pval)
  colnames(summarydf) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "P.val")
  out <- list(call = object$call, data = object$data, methods = object$methods,
              variables = object$variables, multimp = object$multimp, ref = object$ref, 
              regsumm = regsumm, summarydf = summarydf)
  class(out) <- c("summary.cmest")
  return(out)
}

#' @describeIn cmest Print the summary of cmest nicely
#' @export
print.summary.cmest <- function(x, digits = 4, ...) {
  # print summary of causal mediation analysis results
  if (x$methods$model == "rb") model_str <- "Regression-based Approach"
  if (x$methods$model == "wb") model_str <- "Weighting-based Approach"
  if (x$methods$model == "ne") model_str <- "Natural Effect Model"
  if (x$methods$model == "iorw") model_str <- "Inverse Odds Ratio Weighting Approach"
  if (x$methods$model == "msm") model_str <- "Marginal Structural Model"
  if (x$methods$model == "gformula") model_str <- "G-formula Approach"
  if (x$multimp$multimp) model_str <- paste(model_str, "with Multiple Imputation")
  if (x$methods$estimation == "paramfunc") est_str <- "Closed-form parameter function estimation"
  if (x$methods$estimation == "imputation") est_str <- "Direct counterfactual imputation estimation"
  if (x$methods$inference == "delta") inf_str <- "delta method standard errors, confidence intervals and p-values"
  if (x$methods$inference == "bootstrap") inf_str <- "bootstrap standard errors, percentile confidence intervals and p-values"
  if (x$methods$model != "ne" && (x$methods$casecontrol)) cat("# Causal Mediation Analysis for a Case Control Study via the ")
  if (!(x$methods$model != "ne" && (x$methods$casecontrol))) cat("# Causal Mediation Analysis via the ")
  cat(model_str)
  cat("\n \n")
  cat(est_str)
  cat(paste(" with \n", inf_str, "\n \n"))
  printCoefmat(x$summarydf, digits = digits, has.Pvalue = TRUE)
  cat("\n")
  cat("Reference values: \n")
  print(x$ref)
}

#' Plotting Point Estimates and Confidence IntervalsCausal of Causal Effect Estimates
#' 
#' This function is used to plot the results of \code{cmest} nicely via plotting functions
#' in the \code{ggplot2} package. Additional layers can be added to this plot using other 
#' plotting functions in the \code{ggplot2} package.
#' 
#' @param x an object of class \code{cmest}
#' @param errorbar.width width of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.size size of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.colour colour of errorbars for confidence intervals. Default is \code{black}.
#' @param point.size size of points for point estimates. Default is \code{1}.
#' @param point.colour colour of points for point estimates. Default is \code{black}.
#' @param refline a logical value. If \code{true}, include a reference line at 
#' \code{y = 0} when effects are on the difference scale and include a reference line at 
#' \code{y = 1} when effects are on the ratio scale. Default is \code{TRUE}.
#' @param refline.colour colour of the reference line. Default is \code{red}.
#' @param refline.size size of the reference line. Default is \code{0.3}.
#' 
#' @seealso \code{\link{cmest}}, \code{\link{ggplot2}}.
#' 
#' @examples
#' 
#' library(CMAverse)
#' library(ggplot2)
#' 
#' x <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = "M2", basec = c("C1", "C2"), 
#' EMint = TRUE, mreg = list("multinomial"), yreg = "linear", 
#' astar = 0, a = 1, mval = list("M2_0"), estimation = "paramfunc", 
#' inference = "delta")
#' 
#' ggcmest(x) +
#' theme(axis.text.x = element_text(angle = 45))
#' 
#' ggcmest(x) +
#' coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#' 
#' @export
#' 
ggcmest <- function(x, errorbar.width = 0.3, errorbar.size = 0.3, errorbar.colour = "black",
                    point.size = 1, point.colour = "blue", 
                    refline = TRUE, refline.colour = "red", refline.size = 0.3) {
  # create a data frame for results of cmest
  effect_df <- data.frame(Effect = factor(names(x$effect.pe), levels = names(x$effect.pe)),
                          PE = x$effect.pe, CIlower = x$effect.ci.low,
                          CIupper = x$effect.ci.high)
  # reference line
  if (refline) {
    if (!x$multimp$multimp) {
      if ((inherits(x$reg.output$yreg, "lm") | inherits(x$reg.output$yreg, "glm")) &&
          (family(x$reg.output$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else ref <- 1
    } else {
      if ((inherits(x$reg.output[[1]]$yreg, "lm") | inherits(x$reg.output[[1]]$yreg, "glm")) &&
          (family(x$reg.output[[1]]$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else ref <- 1
    }
  } else ref <- NULL
  # plot
  ggplot() +
    geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper),
                  width = errorbar.width, size = errorbar.size, colour = errorbar.colour,
                  data = effect_df) +
    geom_point(aes(x = Effect, y = PE),
               size = point.size, colour = point.colour, data = effect_df) +
    ylab("Point Estimate and 95% CI") +
    geom_hline(yintercept = ref, color = refline.colour, size = refline.size)
}
