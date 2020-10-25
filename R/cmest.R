#' Causal Mediation Analysis
#'
#' \code{cmest} is used to implement six causal mediation analysis approaches including 
#' \emph{the regression-based approach} by Valeri et al. (2013) and VanderWeele 
#' et al. (2014), \emph{the weighting-based approach} by VanderWeele et al. (2014), 
#' \emph{the inverse odd-ratio weighting approach} by Tchetgen Tchetgen (2013), 
#' \emph{the natural effect model} by Vansteelandt et al. (2012), \emph{the marginal structural 
#' model} by VanderWeele et al. (2017), and \emph{the g-formula approach} by Robins (1986).
#'
#' @param data dataset
#' @param model causal mediation analysis approach. \code{rb}, \code{wb}, \code{iorw}, 
#' \code{ne}, \code{msm} and \code{gformula} are implemented. Default is \code{rb}. 
#' See \code{Details}.
#' @param full a logical value. If \code{TRUE}, output a full list of causal effects; if 
#' \code{FALSE}, output a reduced list of causal effects. Default is \code{TRUE}. 
#' @param casecontrol a logical value. \code{TRUE} indicates a case control study in which the
#' first level of the outcome is treated as the control and the second level of the outcome is 
#' treated as the case. Default is \code{FALSE}.
#' @param yrare a logical value (used when \code{casecontrol} is \code{TRUE}). \code{TRUE} 
#' indicates the case is rare.
#' @param yprevalence the prevalence of the case (used when \code{casecontrol} is \code{TRUE}).
#' @param estimation method for estimating causal effects. \code{paramfunc} and 
#' \code{imputation} are implemented (the first 4 letters are enough). Default is \code{imputation}. 
#' See \code{Details}. 
#' @param inference method for estimating standard errors of causal effects. \code{delta} and 
#' \code{bootstrap} are implemented (the first 4 letters are enough). Default is \code{bootstrap}. 
#' See \code{Details}.
#' @param outcome variable name of the outcome.
#' @param event variable name of the event (used when \code{yreg} is \code{coxph}, \code{aft_exp},
#' or \code{aft_weibull}).
#' @param exposure variable name of the exposure.
#' @param mediator a vector of variable name(s) of the mediator(s) following the temporal order.
#' @param EMint a logical value (used when \code{yreg} is character and \code{model} is not 
#' \code{iorw}). If \code{TRUE}, the outcome regression formula includes the interaction(s) 
#' between the exposure and each of the mediator(s).
#' @param basec a vector of variable name(s) of the exposure-outcome confounder(s), exposure-mediator 
#' confounder(s) and mediator-outcome confounder(s) not affected by the exposure
#' @param postc a vector of variable name(s) of the mediator-outcome confounder(s) affected 
#' by the exposure following the temporal order
#' @param yreg outcome regression model. See \code{Details}.
#' @param mreg a list specifying a regression model for each variable in \code{mediator} (used 
#' when \code{model} is \code{rb}, \code{msm} or \code{gformula}). The order of regression models 
#' must follow the order of variables in \code{mediator}. See \code{Details}.
#' @param wmnomreg a list specifying a regression model for calculating the nominators of 
#' weights with respect to each variable in \code{mediator} (used when \code{model} is \code{msm}). 
#' The order of regression models must follow the order of variables in \code{mediator}. See 
#' \code{Details}.
#' @param wmdenomreg a list specifying a regression model for calculating the denominators of 
#' weights with respect to each variable in \code{mediator} (used when \code{model} is \code{msm}). 
#' The order of regression models must follow the order of variables in \code{mediator}. See 
#' \code{Details}.
#' @param ereg exposure regression model for calculating weights with respect to the exposure (used 
#' when \code{model} is \code{wb} or \code{msm} with a non-empty \code{basec} or when 
#' \code{model} is \code{iorw}). See \code{Details}.
#' @param postcreg a list specifying a regression model for each variable in \code{postc} (used 
#' when \code{model} is \code{gformula}). The order of regression models must follow the order of 
#' variables in \code{postc}. See \code{Details}.
#' @param astar the control value for the exposure. Default is \code{0}.
#' @param a the active value for the exposure. Default is \code{1}.
#' @param mval a list specifying a value for each variable in \code{mediator} at which the variable is 
#' controlled (used when \code{model} is \code{rb}, \code{wb}, \code{ne}, \code{msm} or \code{gformula}).
#' @param yval the value of the outcome at which causal effects on the risk/odds ratio 
#' scale are estimated (used when the outcome is categorical).
#' @param basecval a list specifying a conditional value for each variable in \code{basec} conditional 
#' on which causal effects are estimated (used when \code{estimation} is \code{paramfunc}). The order 
#' of conditional values must follow the order of variables in \code{basec}. If \code{NULL}, 
#' mean values of variable(s) in \code{basec} are used.
#' @param nboot the number of boots applied (used when \code{inference} is \code{bootstrap}). 
#' Default is 200.
#' @param boot.ci.type the type of bootstrap confidence interval. If \code{per}, percentile bootstrap
#' confidence intervals are estimated; if \code{bca}, bias-corrected and accelerated (BCa) bootstrap 
#' confidence intervals are estimated. Default is \code{per}.
#' @param nRep number of replications or hypothetical values of the exposure to sample for 
#' each observation unit (used when \code{model} is \code{ne}). Default is \code{5}.
#' @param multimp a logical value (used when \code{data} contains missing values). If 
#' \code{TRUE}, conduct multiple imputations using the \link[mice]{mice} function. Default is 
#' \code{FALSE}.
#' @param ... Additional arguments passed to the \link[mice]{mice} function. See \link[mice]{mice}
#' for details.
#' @param x an object of class \code{cmest}
#' @param object an object of class \code{cmest}
#' @param digits minimal number of significant digits. See \link{print.default}.
#' 
#' @details
#' 
#' \strong{Regressions}
#'
#' Each regression in \code{yreg}, \code{mreg}, \code{wmnomreg}, \code{wmdenomreg},
#' \code{ereg} and \code{postcreg} can be specified by a user-defined regression 
#' object or the character name of the regression. 
#' 
#' \emph{The Character Name of A Regression}
#' 
#' \itemize{
#' \item{\code{linear}: }{linear regression fitted by \link{glm} with \code{family = gaussian()}}
#' \item{\code{logistic}: }{logistic regression fitted by \link{glm} with \code{family = logit()}}
#' \item{\code{loglinear}: }{log linear regression fitted by \link{glm} with 
#' \code{family = poisson()} for a binary response}
#' \item{\code{poisson}: }{poisson regression fitted by \link{glm} with 
#' \code{family = poisson()} for a count response}
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
#' \code{coxph}, \code{aft_exp} and \code{aft_weibull} are currently not implemented for 
#' \code{mreg}, \code{wmnomreg}, \code{wmdenomreg}, \code{ereg} and \code{postcreg}.
#' 
#' \emph{The User-defined Regression Object} 
#' 
#' A user-defined regression object can be fitted by \link{lm}, \link{glm}, \link{glm.nb}, 
#' \link[mgcv]{gam}, \link[nnet]{multinom}, \link[MASS]{polr}, \link[survival]{coxph} and
#' \link[survival]{survreg}. Objects fitted by \link[survival]{coxph} and \link[survival]{survreg} 
#' are currently not supported for \code{mreg}, \code{wmnomreg}, \code{wmdenomreg}, 
#' \code{ereg} and \code{postcreg}.
#' 
#' The \code{cmest} function calculates weights for regressions when weighting is required. If a 
#' user-defined regression object is fitted with prior weights, the final weights for this 
#' regression object are constructed by multiplying the prior weights and the weights calculated 
#' inside the \code{cmest} function. 
#' 
#' 
#' \strong{Causal Mediation Analysis Approaches}
#' 
#' Let \code{Y} denote the outcome, \code{A} denote the exposure, \code{M=(M_1,...,M_k)^T} 
#' denote \code{mediator}, \code{C} denote \code{basec}, \code{L=(L_1,...,L_s)^T} denote \code{postc}.
#' 
#' \itemize{
#'     \item{\code{rb}: }{\emph{the regression-based approach} by Valeri et al. (2013) and
#'       VanderWeele et al. (2014). \code{yreg} and \code{mreg} are required. If specified as
#'       a user-defined regression object, \code{yreg} should regress \code{Y} on \code{A},
#' \code{M} and \code{C} and \code{mreg[p]} should regress \code{M_p} on \code{A} and \code{C}} 
#' for \code{p=1,...,k}.
#' 
#'     \item{\code{wb}: }{\emph{the weighting-based approach} by VanderWeele et al. (2014).
#'     \code{yreg} is required. When \code{basec} is not empty, \code{ereg} is also required 
#'     and \code{A} must be categorical. If specified as a user-defined regression object, 
#'     \code{yreg} should regress \code{Y} on \code{A}, \code{M} and \code{C} and \code{ereg} 
#'     should regress \code{A} on \code{C}.}
#'     
#'     \item{\code{iorw}: }{\emph{the inverse odd-ratio weighting approach} by
#'       Tchetgen Tchetgen (2013). \code{yreg} and \code{ereg} are required and 
#'       \code{A} must be categorical. If specified as a user-defined regression object, 
#'       \code{yreg} should regress \code{Y} on \code{A} and \code{C} and \code{ereg} should 
#'       regress \code{A} on \code{M} and \code{C}.}
#'       
#'     \item{\code{ne}: }{\emph{the natural effect model} by Vansteelandt et al. (2012). 
#'     \code{yreg} is required. If specified as a user-defined regression object, \code{yreg}
#'     should regress \code{Y} on \code{A}, \code{M} and \code{C}. The variables in the 
#'     formula of \code{yreg} must follow the order of \code{A}, \code{M} and \code{C}, i.e., 
#'     the first variable must point to the exposure, the variable(s) right after the 
#'     exposure must point to the mediator(s), e.g., \code{Y ~ A + M_1 + M_2 + A*M_1 + C}.}
#'       
#'     \item{\code{msm}: }{\emph{the marginal structural model} by VanderWeele et al. (2017).
#'     \code{yreg}, \code{mreg}, \code{wmnomreg} and \code{wmdenomreg} are required and all 
#'     mediators must be categorical. When \code{basec} is not empty, \code{ereg} is also 
#'     required and \code{A} must be categorical. If specified as a user-defined regression 
#'     object, \code{yreg} should regress \code{Y} on \code{A} and \code{M}; \code{mreg[p]} 
#'     should regress \code{M_p} on \code{A} for \code{p=1,...,k}; \code{wmnomreg[p]} should 
#'     regress \code{M_p} on \code{A}, \code{M_1}, ..., \code{M_{p-1}} for \code{p=1,...,k};
#'     \code{wmdenomreg[p]} should regress \code{M_p} on \code{A}, \code{M_1}, ..., \code{M_{p-1}},
#'     \code{C} and \code{L} for \code{p=1,...,k}; and \code{ereg} should regress \code{A} on 
#'     \code{C}.}
#'     
#'     \item{\code{gformula}: }{\emph{the g-formula approach} by Robins (1986). 
#'     \code{yreg}, \code{mreg} are required. \code{postcreg} is also required when \code{postc}
#'     is not empty. If specified as a user-defined regression object, \code{yreg} should 
#'     regress \code{Y} on \code{A}, \code{M}, \code{C} and \code{L}, \code{mreg[p]} should 
#'     regress \code{M_p} on \code{A}, \code{C} and \code{L}
#'     for \code{p=1,...,k}, \code{postcreg[q]} should regress \code{L_q} on \code{A} and \code{C} 
#'     for \code{q=1,...,s}.}
#'   }
#'   
#' When \code{postc} is not empty, only \code{msm} and \code{gformula} can be used.
#'   
#' 
#' \strong{Estimation Methods} 
#'   
#' \itemize{
#' \item{\code{paramfunc}: }{closed-form parameter function estimation (only available when 
#' \code{model = "rb"} and \code{length(mediator) = 1}). The point estimate of each causal 
#' effect is obtained by a closed-form formula of regression coefficients. Effects conditional 
#' on \code{basecval} are estimated.}
#' \item{\code{imputation}: }{direct counterfactual imputation estimation. The point estimate 
#' of each causal effect is obtained by imputing counterfactuals directly.}
#' }
#' 
#' To use \code{paramfunc}, \code{yreg} and \code{mreg} must be specified by the character name 
#' of the regression. \code{yreg} can be chosen from \code{linear}, \code{logistic}, \code{loglinear}, 
#' \code{poisson}, \code{quasipoisson}, \code{negbin}, \code{coxph}, \code{aft_exp} and 
#' \code{aft_weibull}. \code{mreg} can be chosen from \code{linear}, \code{logistic} and 
#' \code{multinomial}.
#' 
#' To use \code{paramfunc} with \code{yreg = "logistic"}, the outcome must be rare. To use 
#' \code{paramfunc} with \code{yreg = "coxph"}, the outcome at the end of follow-up must be rare.
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
#' For a continuous outcome, the causal effects on the difference scale are estimated. For a 
#' categorical, count or survival outcome, the causal effects on the ratio scale are estimated. The 
#' interpretation of the ratio depends on the type of the outcome and it can be risk 
#' ratio for a categorical outcome, rate ratio for a count outcome, hazard ratio for a survival 
#' outcome fitted by \link[survival]{coxph}, mean survival ratio for a survival outcome fitted 
#' by \link[survival]{survreg}, etc.
#' 
#' \emph{Continuous Outcome}
#' 
#' When \code{model = "rb", "wb", "ne", "msm" or "gformula"} with an empty \code{postc}, 
#' \code{cde} (controlled direct effect), \code{pnde} (pure natural direct effect), 
#' \code{tnde} (total natural direct effect), \code{pnie} (pure natural indirect effect), 
#' \code{tnie} (total natural indirect effect), \code{te} (total effect),  \code{intref} 
#' (reference interaction), \code{intmed} (mediated interaction), 
#' \code{cde(prop)} (proportion \code{cde}), \code{intref(prop)} (proportion 
#' \code{intref}), \code{intmed(prop)} (proportion \code{intmed}), \code{pnie(prop)} 
#' (proportion \code{pnie}), \code{pm} (proportion mediated), \code{int} (proportion 
#' attributable to interaction) and \code{pe} (proportion eliminated) are estimated. 
#' 
#' When \code{postc} is not empty, \code{pnde}, \code{tnde}, \code{pnie}, \code{tnie},  
#' \code{intref}, \code{intmed}, \code{intref(prop)}, \code{intmed(prop)}, \code{pnie(prop)}, 
#' \code{pm}, \code{int} and \code{pe} are replaced by their randomized analogues \code{rpnde}, 
#' \code{rtnde}, \code{rpnie}, \code{rtnie}, \code{rintref}, \code{rintmed}, \code{rintref(prop)}, 
#' \code{rintmed(prop)}, \code{rpnie(prop)}, \code{rpm}, \code{rint} and \code{rpe}.
#' 
#' When \code{model = "iorw"}, \code{te}, \code{pnde}, \code{tnie} and \code{pm} are estimated. 
#' 
#' \emph{Categorical, Count or Survival Outcome}
#' 
#' When \code{model = "rb", "wb", "ne", "msm" or "gformula"} with an empty \code{postc}, 
#' \code{Rcde} (\code{cde} ratio), \code{Rpnde} (\code{pnde} ratio), \code{Rtnde} (\code{tnde} 
#' ratio), \code{Rtnie} (\code{tnie} ratio), \code{Rte} (\code{te} ratio), \code{ERcde} (excess 
#' ratio due to \code{cde}), \code{ERintref} (excess ratio due to \code{intref}), 
#' \code{ERintmed} (excess ratio due to \code{intmed}), \code{ERpnie} (excess ratio due to 
#' \code{pnie}), \code{ERcde(prop)} (proportion \code{ERcde}), \code{ERintref(prop)} (proportion 
#' \code{ERintref}), \code{ERintmed(prop)} (proportion \code{ERintmed}), 
#' \code{ERpnie(prop)} (proportion \code{ERpnie}), \code{pm}, \code{int} and \code{pe} are estimated. 
#' 
#' When \code{model = "msm" or "gformula"} with a non-empty \code{postc}, \code{Rpnde}, \code{Rtnde}, 
#' \code{Rpnie}, \code{Rtnie},  \code{ERintref}, \code{ERintmed}, \code{ERpnie}, 
#' \code{ERintref(prop)}, \code{ERintmed(prop)}, \code{ERpnie(prop)}, \code{pm}, \code{int} and 
#' \code{pe} are replaced by their randomized analogues \code{rRpnde}, \code{rRtnde}, \code{rRpnie}, 
#' \code{rRtnie}, \code{rERintref}, \code{rERintmed}, \code{rERpnie}, \code{rERintref(prop)}, 
#' \code{rERintmed(prop)}, \code{rERpnie(prop)}, \code{rpm}, \code{rint} and \code{rpe}.
#' 
#' When \code{model = "iorw"}, \code{Rte}, \code{Rpnde}, \code{Rtnie} and \code{pm} are estimated. 
#' 
#' @return
#' An object of class \code{cmest} is returned:
#' \item{call}{the function call,}
#' \item{data}{the dataset,}
#' \item{methods}{a list of methods used which may include \code{model}, \code{full}, 
#' \code{casecontrol}, \code{yprevalence}, \code{yrare}, \code{estimation}, \code{inference}, 
#' \code{nboot}, \code{boot.ci.type} and \code{nRep},}
#' \item{variables}{a list of variables used which may include \code{outcome}, \code{event}, 
#' \code{exposure}, \code{mediator}, \code{EMint}, \code{basec} and \code{postc},}
#' \item{reg.input}{a list of regressions input,}
#' \item{multimp}{a list of arguments used for multiple imputation,}
#' \item{ref}{reference values used which may include \code{a}, \code{astar}, \code{mval}, 
#' \code{yval} and \code{basecval},}
#' \item{reg.output}{a list of regressions output. If \code{multimp} is \code{TRUE}, 
#' reg.output contains regressions fitted by each of the imputed dataset,}
#' \item{effect.pe}{point estimates of causal effects,}
#' \item{effect.se}{standard errors of causal effects,}
#' \item{effect.ci.low}{the lower limits of 95\% confidence intervals of causal effects,}
#' \item{effect.ci.high}{the higher limits of 95\% confidence intervals of causal effects,}
#' \item{effect.pval}{p-values of causal effects,}
#' ...
#'
#' @seealso \link{ggcmest}, \link{cmdag}, \link{cmsens}, \link{svymultinom}.
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
#' Nguyen QC, Osypuk TL, Schmidt NM, Glymour MM, Tchetgen Tchetgen EJ (2015). Practical guidance
#' for conducting mediation analysis with multiple mediators using inverse odds ratio
#' weighting. American Journal of Epidemiology. 181(5): 349 - 356.
#' 
#' VanderWeele TJ, Tchetgen Tchetgen EJ (2017). Mediation analysis with time varying
#' exposures and mediators. Journal of the Royal Statistical Society: Series B (Statistical
#' Methodology). 79(3): 917 - 938.
#' 
#' Robins JM (1986). A new approach to causal inference in mortality studies with a sustained 
#' exposure period-Application to control of the healthy worker survivor effect. Mathematical 
#' Modelling. 7: 1393 - 1512.
#' 
#' Vansteelandt S, Bekaert M, Lange T (2012). Imputation Strategies for the Estimation 
#' of Natural Direct and Indirect Effects. Epidemiologic Methods. 1(1): 131 - 158.
#' 
#' Steen J, Loeys T, Moerkerke B, Vansteelandt S (2017). Medflex: an R package for
#' flexible mediation analysis using natural effect models. Journal of Statistical
#' Software. 76(11).
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
#' estimation = "imputation", inference = "bootstrap", nboot = 10,
#' boot.ci.type = "bca")
#' 
#' # multiple-mediator case with ne
#' exp3 <- cmest(data = cma2020, model = "ne", outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), basec = c("C1", "C2"), 
#' yreg = glm(contY ~ A + M1 + M2 + A*M1 + A*M2 + C1 + C2, family = gaussian, data = cma2020), 
#' astar = 0, a = 1, mval = list(0, "M2_0"), estimation = "imputation", 
#' inference = "bootstrap", nboot = 10)
#' 
#' # case control study with msm
#' exp4 <- cmest(data = cma2020, model = "msm", casecontrol = TRUE, 
#' yrare = TRUE, outcome = "binY", exposure = "A", 
#' mediator = c("M1", "M2"), EMint = TRUE, basec = c("C1", "C2"), yreg = "logistic", 
#' ereg = "logistic", mreg = list(glm(M1 ~ A, family = binomial, 
#' data = cma2020), nnet::multinom(M2 ~ A, data = cma2020, trace = FALSE)), 
#' wmnomreg = list(glm(M1 ~ A, family = binomial, data = cma2020), 
#' nnet::multinom(M2 ~ A + M1, data = cma2020, trace = FALSE)),
#' wmdenomreg = list(glm(M1 ~ A + C1 + C2, family = binomial, data = cma2020), 
#' nnet::multinom(M2 ~ A + M1 + C1 + C2, data = cma2020, trace = FALSE)), astar = 0, a = 1, 
#' mval = list(0, "M2_0"), estimation = "imputation", 
#' inference = "bootstrap", nboot = 10)
#' }
#' 
#' @importFrom stats glm binomial poisson as.formula gaussian quasipoisson model.frame printCoefmat 
#' family sd coef vcov sigma predict rbinom rmultinom rnorm rgamma rpois weighted.mean 
#' model.matrix getCall quantile qnorm pnorm lm cov formula update na.pass
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
#'
#' @export

cmest <- function(data = NULL, model = "rb",
                  full = TRUE, casecontrol = FALSE, yrare = NULL, yprevalence = NULL,
                  estimation = "imputation", inference = "bootstrap",
                  outcome = NULL, event = NULL,
                  exposure = NULL, mediator = NULL, EMint = NULL, basec = NULL, postc = NULL,
                  yreg = NULL, mreg = NULL, wmnomreg = NULL, wmdenomreg = NULL, ereg = NULL, 
                  postcreg = NULL,
                  astar = 0, a = 1, mval = NULL, yval = NULL, basecval = NULL,
                  nboot = 200, boot.ci.type = "per", nRep = 5, multimp = FALSE, ...) {
  # function call
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
  if (casecontrol) {
    if (length(unique(data[, outcome])) != 2) stop("When casecontrol is TRUE, the outcome must be binary")
    if (is.null(yprevalence) && yrare != TRUE) stop("When casecontrol is TRUE, specify yprevalence or set yrare to be TRUE")
    if (!is.null(yprevalence)) {
      if (!is.numeric(yprevalence)) stop("yprevalence should be numeric")
      out$methods$yprevalence <- yprevalence
    } else out$methods$yrare <- yrare
  }
  # estimation, inference, nboot
  if (estimation == "para") estimation <- "paramfunc"
  if (estimation == "impu") estimation <- "imputation"
  if (inference == "delt") inference <- "delta"
  if (inference == "boot") inference <- "bootstrap"
  if (model == "rb" && !estimation %in% c("paramfunc", "imputation")) stop("When model = 'rb', select estimation from 'paramfunc', 'imputation'")
  if (model != "rb" && !estimation == "imputation") stop("Use estimation = 'imputation'")
  if (estimation == "paramfunc") {
    if (!yreg %in% c("linear", "logistic", "loglinear", "poisson",
                     "quasipoisson", "negbin", "coxph", "aft_exp", "aft_weibull")) stop(
                       "When estimation = 'paramfunc', select yreg from 'linear', 'logistic', 
                       'loglinear', 'poisson', 'quasipoisson', 'negbin', 'coxph', 'aft_exp', 'aft_weibull'")
    if (length(mediator) > 1) stop("'paramfunc' only supports a single mediator")
    if (!mreg[[1]] %in% c("linear", "logistic", "multinomial")) stop(
      "When estimation = 'paramfunc', select mreg[[1]] from 'linear', 'logistic', 'multinomial'")
  }
  out$methods$estimation <- estimation
  if (!inference %in% c("delta", "bootstrap")) stop("Select inference from 'delta', 'bootstrap'")
  if (estimation == "imputation" && inference == "delta") stop("Use inference = 'bootstrap' when estimation = 'imputation'")
  out$methods$inference <- inference
  if (inference == "bootstrap") {
    if (!is.numeric(nboot)) stop("nboot should be numeric")
    if (!boot.ci.type %in% c("per", "bca")) stop("Select boot.ci.type from 'per', 'bca'")
    out$methods$nboot <- nboot
    out$methods$boot.ci.type <- boot.ci.type
  }
  # outcome
  if (length(outcome) == 0) stop("Unspecified outcome")
  if (length(outcome) > 1) stop("length(outcome) > 1")
  out$variables$outcome <- outcome
  # event
  if (yreg %in% c("coxph", "aft_exp", "aft_weibull") && !is.null(event)) out$variables$event <- event
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
  # basec
  if (!is.null(basec)) out$variables$basec <- basec
  # postc
  if (length(postc) != 0) {
    if (!model %in% c("msm", "gformula")) stop("When postc is not empty, select model from 'msm' and 'gformula'")
    out$variables$postc <- postc
  }
  #regs
  if (model == "rb") out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "wb") out$reg.input <- list(yreg = yreg)
  if (model == "gformula") out$reg.input <- list(yreg = yreg, mreg = mreg)
  if (model == "msm") out$reg.input <- list(yreg = yreg, mreg = mreg, wmnomreg = wmnomreg, wmdenomreg = wmdenomreg)
  if (model == "iorw") out$reg.input <- list(yreg = yreg, ereg = ereg)
  if (model == "ne") out$reg.input <- list(yreg = yreg) 
  if (length(basec) != 0 && model %in% c("wb", "msm")) out$reg.input$ereg <- ereg
  if (length(postc) != 0 && model == "gformula") out$reg.input$postcreg <- postcreg
  # a, astar
  if (is.null(a) | is.null(astar)) stop("Unspecified a or astar")
  # mval
  if (model != "iorw") {
    if (!is.list(mval)) stop("mval should be a list")
    if (length(mval) != length(mediator)) stop("length(mval) != length(mediator)")
    for (p in 1:length(mval)) if (is.null(mval[[p]])) stop(paste0("Unspecified mval[[", p, "]]"))
  } 
  # basecval
  if (!(model == "rb" && estimation == "paramfunc" && length(basec) != 0) && !is.null(basecval)) warning("basecval is ignored")
  # nRep
  if (model == "ne") {
    if (!is.numeric(nRep)) stop("nRep should be numeric")
    out$methods$nRep <- nRep
  }
  # multimp
  out$multimp <- list(multimp = multimp)
  if (!is.logical(multimp)) stop("multimp should be TRUE or FALSE")
  # args_mice
  if (multimp) {
  args_mice <- list(...)
  args_mice$print <- FALSE
  if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
  args_mice$data <- data
  out$multimp$args_mice <- args_mice
  }
  if (!multimp && model %in% c("wb", "iorw", "ne", "msm")) {
    if (sum(is.na(data)) > 0) stop("Selected model doesn't support missing values; use multimp = TRUE")
  }
  
  # run regressions
  environment(regrun) <- environment()
  regs <- regrun()
  yreg <- regs$yreg
  ereg <- regs$ereg
  mreg <- regs$mreg
  wmnomreg <- regs$wmnomreg
  wmdenomreg <- regs$wmdenomreg
  postcreg <- regs$postcreg
  
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
  environment(estinf) <- environment()
  out <- c(out, estinf())
  class(out) <- "cmest"
  return(out)
}


#' @describeIn cmest Print the results of \code{cmest} nicely
#' @export
print.cmest <- function(x, ...) {
  # print regression models used
  if (!x$multimp$multimp) {
    regnames <- names(x$reg.output)
    for (name in regnames) {
      if (name == "yreg") {
        cat("# Outcome Regression:\n")
        if (inherits(x$reg.output$yreg, "svyglm")) {
          x$reg.output$yreg$call <- update(x$reg.output$yreg,design = getCall(x$reg.output$yreg)$design,
                                           family = getCall(x$reg.output$yreg)$family, evaluate = FALSE)
          x$reg.output$yreg$survey.design$call <- as.call(update(summary(x$reg.output$yreg)$survey.design,
                                                                 data = getCall(summary(x$reg.output$yreg)$survey.design)$data,
                                                                 weights = getCall(summary(x$reg.output$yreg)$survey.design)$weights, 
                                                                 evaluate = FALSE))
          print(x$reg.output$yreg)
        } else {
          x$reg.output$yreg$call <- update(x$reg.output$yreg,data=getCall(x$reg.output$yreg)$data,
                                           weights=getCall(x$reg.output$yreg)$weights, evaluate = FALSE)
          print(x$reg.output$yreg)
        }
      }
      if (name == "yregTot") {
        cat("# Outcome Regression for the Total Effect: \n")
        x$reg.output$yregTot$call <- update(x$reg.output$yregTot,data=getCall(x$reg.output$yregTot)$data,
                                            weights=getCall(x$reg.output$yregTot)$weights, evaluate = FALSE)
        print(x$reg.output$yregTot)
      }
      if (name == "yregDir") {
        cat("# Outcome Regression for the Direct Effect: \n")
        x$reg.output$yregDir$call <- update(x$reg.output$yregDir,data=getCall(x$reg.output$yregDir)$data,
                                            weights=getCall(x$reg.output$yregDir)$weights, evaluate = FALSE)
        print(x$reg.output$yregDir)
      }
      if (name == "ereg") {
        cat("# Exposure Regression for Weighting: \n")
        x$reg.output$ereg$call <- update(x$reg.output$ereg,data=getCall(x$reg.output$ereg)$data,
                                         weights=getCall(x$reg.output$ereg)$weights, evaluate = FALSE)
        print(x$reg.output$ereg)
      }
      if (name == "mreg") {
        cat("# Mediator Regressions: \n")
        for (i in 1:length(x$reg.output$mreg)) {
          if (inherits(x$reg.output$mreg[[i]], "svyglm")) {
            x$reg.output$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]],
                                                         design = getCall(x$reg.output$mreg[[.(i)]])$design,
                                                         family = getCall(x$reg.output$mreg[[.(i)]])$family, evaluate = FALSE)))
            x$reg.output$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(summary(x$reg.output$mreg[[.(i)]])$survey.design,
                                                                   data = getCall(summary(x$reg.output$mreg[[.(i)]])$survey.design)$data,
                                                                   weights = getCall(summary(x$reg.output$mreg[[.(i)]])$survey.design)$weights, 
                                                                   evaluate = FALSE))))
            print(x$reg.output$mreg[[i]])
          } else {
            x$reg.output$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]], 
                                                              data=getCall(x$reg.output$mreg[[.(i)]])$data, 
                                                              weights=getCall(x$reg.output$mreg[[.(i)]])$weights,
                                                              evaluate = FALSE)))
            print(x$reg.output$mreg[[i]])
          }
          if (i < length(x$reg.output$mreg)) cat("\n")
        }
      }
      if (name == "wmdenomreg") {
        cat("# Mediator Regressions for Weighting (Denominator): \n")
        for (i in 1:length(x$reg.output$wmdenomreg)) {
          x$reg.output$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmdenomreg[[i]], 
                                                                  data=getCall(x$reg.output$wmdenomreg[[.(i)]])$data, 
                                                                  weights=getCall(x$reg.output$wmdenomreg[[.(i)]])$weights,
                                                                  evaluate = FALSE)))
          print(x$reg.output$wmdenomreg[[i]])
          if (i < length(x$reg.output$wmdenomreg)) cat("\n")
        }
      }
      if (name == "wmnomreg") {
        cat("# Mediator Regressions for Weighting (Nominator): \n")
        for (i in 1:length(x$reg.output$wmnomreg)) {
          x$reg.output$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmnomreg[[i]], 
                                                                data=getCall(x$reg.output$wmnomreg[[.(i)]])$data, 
                                                                weights=getCall(x$reg.output$wmnomreg[[.(i)]])$weights,
                                                                evaluate = FALSE)))
          print(x$reg.output$wmnomreg[[i]])
          if (i < length(x$reg.output$wmnomreg)) cat("\n")
        }
      }
      if (name == "postcreg") {
        cat("# Regressions for Mediator-outcome Confounders Affected by the Exposure: \n")
        for (i in 1:length(x$reg.output$postcreg)) {
          x$reg.output$postcreg[[i]]$call <- eval(bquote(update(x$reg.output$postcreg[[i]], 
                                                                data=getCall(x$reg.output$postcreg[[.(i)]])$data, 
                                                                weights=getCall(x$reg.output$postcreg[[.(i)]])$weights,
                                                                evaluate = FALSE)))
          print(x$reg.output$postcreg[[i]])
          if (i < length(x$reg.output$postcreg)) cat("\n")
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
          cat("## Outcome Regression: \n")
          if (inherits(x$reg.output[[m]]$yreg, "svyglm")) {
            x$reg.output[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,design = getCall(x$reg.output[[.(m)]]$yreg)$design,
                                             family = getCall(x$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
            x$reg.output[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(summary(x$reg.output[[.(m)]]$yreg)$survey.design,
                                                                   data = getCall(summary(x$reg.output[[.(m)]]$yreg)$survey.design)$data,
                                                                   weights = getCall(summary(x$reg.output[[.(m)]]$yreg)$survey.design)$weights, 
                                                                   evaluate = FALSE))))
            print(x$reg.output[[m]]$yreg)
          } else {
            x$reg.output[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,
                                                              data = getCall(x$reg.output[[.(m)]]$yreg)$data,
                                                              weights = getCall(x$reg.output[[.(m)]]$yreg)$weights, 
                                                              evaluate = FALSE)))
            print(x$reg.output[[m]]$yreg)
          }
        }
        if (name == "yregTot") {
          cat("## Outcome Regression for the Total Effect: \n")
          x$reg.output[[m]]$yregTot$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregTot,
                                                               data=getCall(x$reg.output[[.(m)]]$yregTot)$data,
                                                               weights=getCall(x$reg.output[[.(m)]]$yregTot)$weights, 
                                                               evaluate = FALSE)))
          print(x$reg.output[[m]]$yregTot)
        }
        if (name == "yregDir") {
          cat("## Outcome Regression for the Direct Effect: \n")
          x$reg.output[[m]]$yregDir$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregDir,
                                                               data=getCall(x$reg.output[[.(m)]]$yregDir)$data,
                                                               weights=getCall(x$reg.output[[.(m)]]$yregDir)$weights, 
                                                               evaluate = FALSE)))
          print(x$reg.output[[m]]$yregDir)
        }
        if (name == "ereg") {
          cat("## Exposure Regression for Weighting: \n")
          x$reg.output[[m]]$ereg$call <- eval(bquote(update(x$reg.output[[.(m)]]$ereg,
                                                            data=getCall(x$reg.output[[.(m)]]$ereg)$data,
                                                            weights=getCall(x$reg.output[[.(m)]]$ereg)$weights, 
                                                            evaluate = FALSE)))
          print(x$reg.output[[m]]$ereg)
        }
        if (name == "mreg") {
          cat("## Mediator Regressions: \n")
          for (i in 1:length(x$reg.output[[m]]$mreg)) {
            if (inherits(x$reg.output[[m]]$mreg[[i]], "svyglm")) {
              x$reg.output[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]],
                                                                design = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$design,
                                                                family = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$family, evaluate = FALSE)))
              x$reg.output[[m]]$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(summary(x$reg.output[[.(m)]]$mreg[[.(i)]])$survey.design,
                                                                                      data = getCall(summary(x$reg.output[[.(m)]]$mreg[[.(i)]])$survey.design)$data,
                                                                                      weights = getCall(summary(x$reg.output[[.(m)]]$mreg[[.(i)]])$survey.design)$weights, 
                                                                                      evaluate = FALSE))))
              print(x$reg.output[[m]]$mreg[[i]])
            } else {
              x$reg.output[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]], 
                                                                     data=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$data, 
                                                                     weights=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$weights, 
                                                                     evaluate = FALSE)))
              print(x$reg.output[[m]]$mreg[[i]])
            }
            if (i < length(x$reg.output[[m]]$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          if (!is.null(x$reg.output[[m]]$wmdenomreg)) {
            cat("## Mediator Regressions for Weighting (Denominator): \n")
            for (i in 1:length(x$reg.output[[m]]$wmdenomreg)) {
              x$reg.output[[m]]$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmdenomreg[[i]], 
                                                                           data=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$data, 
                                                                           weights=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$weights, 
                                                                           evaluate = FALSE)))
              print(x$reg.output[[m]]$wmdenomreg[[i]])
              if (i < length(x$reg.output[[m]]$wmdenomreg)) cat("\n")
            }
          }
        }
        if (name == "wmnomreg") {
          if (!is.null(x$reg.output[[m]]$wmnomreg)) {
            cat("## Mediator Regressions for Weighting (Nominator): \n")
            for (i in 1:length(x$reg.output[[m]]$wmnomreg)) {
              x$reg.output[[m]]$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmnomreg[[i]], 
                                                                         data=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$data, 
                                                                         weights=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$weights, 
                                                                         evaluate = FALSE)))
              print(x$reg.output[[m]]$wmnomreg[[i]])
              if (i < length(x$reg.output[[m]]$wmnomreg)) cat("\n")
            }
          }
        }
        if (name == "postcreg") {
          if (!is.null(x$reg.output[[m]]$postcreg)) {
            cat("## Regressions for Mediator-outcome Confounders Affected by the Exposure: \n")
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
  if (x$methods$inference == "bootstrap") {
    if (x$methods$boot.ci.type == "per") inf_str <- "bootstrap standard errors, percentile confidence intervals and p-values"
    if (x$methods$boot.ci.type == "bca") inf_str <- "bootstrap standard errors, bias-corrected and accelerated confidence intervals and p-values"
  }
  if (x$methods$model != "ne" && (x$methods$casecontrol)) cat("Causal Mediation Analysis for a Case Control Study via the ")
  if (!(x$methods$model != "ne" && (x$methods$casecontrol))) cat("Causal Mediation Analysis via the ")
  cat(model_str)
  cat("\n \n")
  cat(est_str)
  cat(paste(" with \n", inf_str, "\n \n"))
  print(x$effect.pe)
  cat("\n")
  cat("Relevant parameter values: \n")
  print(x$ref)
}


#' @describeIn cmest Summarize the results of \code{cmest} nicely
#' @export
summary.cmest <- function(object, ...) {
  # summarize regressions
  # print regression models used
  out <- object
  out$reg.output.summary <- out$reg.output
  if (!object$multimp$multimp) {
    regnames <- names(object$reg.output)
    for (name in regnames) {
      if (name == "yreg") out$reg.output.summary$yreg <- summary(object$reg.output$yreg)
      if (name == "yregTot") out$reg.output.summary$yregTot <- summary(object$reg.output$yregTot)
      if (name == "yregDir") out$reg.output.summary$yregDir <- summary(object$reg.output$yregDir)
      if (name == "ereg") out$reg.output.summary$ereg <- summary(object$reg.output$ereg)
      if (name == "mreg") out$reg.output.summary$mreg <- lapply(1:length(object$reg.output$mreg), function(i) 
          summary(object$reg.output$mreg[[i]]))
      if (name == "wmnomreg") out$reg.output.summary$wmnomreg <- lapply(1:length(object$reg.output$wmnomreg), function(i) 
        summary(object$reg.output$wmnomreg[[i]]))
      if (name == "wmdenomreg") out$reg.output.summary$wmdenomreg <- lapply(1:length(object$reg.output$wmdenomreg), function(i) 
        summary(object$reg.output$wmdenomreg[[i]]))
      if (name == "postcreg") out$reg.output.summary$postcreg <- lapply(1:length(object$reg.output$postcreg), function(i) 
        summary(object$reg.output$postcreg[[i]]))
    }
  } else {
    for (m in 1:length(object$reg.output)){ 
      regnames <- names(object$reg.output[[m]])
      for (name in regnames) {
        if (name == "yreg") out$reg.output.summary[[m]]$yreg <- summary(object$reg.output[[m]]$yreg)
        if (name == "yregTot") out$reg.output.summary[[m]]$yregTot <- summary(object$reg.output[[m]]$yregTot)
        if (name == "yregDir") out$reg.output.summary[[m]]$yregDir <- summary(object$reg.output[[m]]$yregDir)
        if (name == "ereg") out$reg.output.summary[[m]]$ereg <- summary(object$reg.output[[m]]$ereg)
        if (name == "mreg") out$reg.output.summary[[m]]$mreg <- lapply(1:length(object$reg.output[[m]]$mreg), function(i) 
          summary(object$reg.output[[m]]$mreg[[i]]))
        if (name == "wmnomreg") out$reg.output.summary[[m]]$wmnomreg <- lapply(1:length(object$reg.output[[m]]$wmnomreg), function(i) 
          summary(object$reg.output[[m]]$wmnomreg[[i]]))
        if (name == "wmdenomreg") out$reg.output.summary[[m]]$wmdenomreg <- lapply(1:length(object$reg.output[[m]]$wmdenomreg), function(i) 
          summary(object$reg.output[[m]]$wmdenomreg[[i]]))
        if (name == "postcreg") out$reg.output.summary[[m]]$postcreg <- lapply(1:length(object$reg.output[[m]]$postcreg), function(i) 
          summary(object$reg.output[[m]]$postcreg[[i]]))
      }
    }
  }
  
  # summarize causal mediation analysis results
  summarydf <- data.frame(object$effect.pe, object$effect.se, object$effect.ci.low, 
                          object$effect.ci.high, object$effect.pval)
  colnames(summarydf) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "P.val")
  out$summarydf <- summarydf
  class(out) <- c("summary.cmest")
  return(out)
}


#' @describeIn cmest Print the summary of \code{cmest} nicely
#' @export
print.summary.cmest <- function(x, digits = 4, ...) {
  # print summary of regression models used
  if (!x$multimp$multimp) {
    regnames <- names(x$reg.output.summary)
    for (name in regnames) {
      if (name == "yreg") {
        cat("# Outcome Regression:\n")
        if (inherits(x$reg.output$yreg, "svyglm")) {
          x$reg.output.summary$yreg$call <- update(x$reg.output$yreg,design = getCall(x$reg.output$yreg)$design,
                                                   family = getCall(x$reg.output$yreg)$family, evaluate = FALSE)
          x$reg.output.summary$yreg$survey.design$call <- as.call(update(x$reg.output.summary$yreg$survey.design,
                                                                         data = getCall(x$reg.output.summary$yreg$survey.design)$data,
                                                                         weights = getCall(x$reg.output.summary$yreg$survey.design)$weights, 
                                                                         evaluate = FALSE)) 
          print(x$reg.output.summary$yreg)
        } else {
          x$reg.output.summary$yreg$call <- update(x$reg.output$yreg,data=getCall(x$reg.output$yreg)$data,
                                           weights=getCall(x$reg.output$yreg)$weights, evaluate = FALSE)
          print(x$reg.output.summary$yreg)
        }
      }
      if (name == "yregTot") {
        cat("# Outcome Regression for the Total Effect: \n")
        x$reg.output.summary$yregTot$call <- update(x$reg.output$yregTot,data=getCall(x$reg.output$yregTot)$data,
                                                 weights=getCall(x$reg.output$yregTot)$weights, evaluate = FALSE)
        print(x$reg.output.summary$yregTot)
      }
      if (name == "yregDir") {
        cat("# Outcome Regression for the Direct Effect: \n")
        x$reg.output.summary$yregDir$call <- update(x$reg.output$yregDir,data=getCall(x$reg.output$yregDir)$data,
                                                    weights=getCall(x$reg.output$yregDir)$weights, evaluate = FALSE)
        print(x$reg.output.summary$yregDir)
      }
      if (name == "ereg") {
        cat("# Exposure Regression for Weighting: \n")
        x$reg.output.summary$ereg$call <- update(x$reg.output$ereg,data=getCall(x$reg.output$ereg)$data,
                                                    weights=getCall(x$reg.output$ereg)$weights, evaluate = FALSE)
        print(x$reg.output.summary$ereg)
      }
      if (name == "mreg") {
        cat("# Mediator Regressions: \n")
        for (i in 1:length(x$reg.output.summary$mreg)) {
          if (inherits(x$reg.output$mreg[[i]], "svyglm")) {
            x$reg.output.summary$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]],design = getCall(x$reg.output$mreg[[.(i)]])$design,
                                                     family = getCall(x$reg.output$mreg[[.(i)]])$family, evaluate = FALSE)))
            x$reg.output.summary$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(x$reg.output.summary$mreg[[.(i)]]$survey.design,
                                                                           data = getCall(x$reg.output.summary$mreg[[.(i)]]$survey.design)$data,
                                                                           weights = getCall(x$reg.output.summary$mreg[[.(i)]]$survey.design)$weights, 
                                                                           evaluate = FALSE))))
            print(x$reg.output.summary$mreg[[i]])
          } else {
            x$reg.output.summary$mreg[[i]]$call <- eval(bquote(update(x$reg.output$mreg[[.(i)]], 
                                                              data=getCall(x$reg.output$mreg[[.(i)]])$data, 
                                                              weights=getCall(x$reg.output$mreg[[.(i)]])$weights,
                                                              evaluate = FALSE)))
            print(x$reg.output.summary$mreg[[i]])
          }
          if (i < length(x$reg.output$mreg)) cat("\n")
        }
      }
      if (name == "wmdenomreg") {
        cat("# Mediator Regressions for Weighting (Denominator): \n")
        for (i in 1:length(x$reg.output.summary$wmdenomreg)) {
          x$reg.output.summary$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmdenomreg[[i]], 
                                                                  data=getCall(x$reg.output$wmdenomreg[[.(i)]])$data, 
                                                                  weights=getCall(x$reg.output$wmdenomreg[[.(i)]])$weights,
                                                                  evaluate = FALSE)))
          print(x$reg.output.summary$wmdenomreg[[i]])
          if (i < length(x$reg.output.summary$wmdenomreg)) cat("\n")
        }
      }
      if (name == "wmnomreg") {
        cat("# Mediator Regressions for Weighting (Nominator): \n")
        for (i in 1:length(x$reg.output.summary$wmnomreg)) {
          x$reg.output.summary$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output$wmnomreg[[i]], 
                                                                data=getCall(x$reg.output$wmnomreg[[.(i)]])$data, 
                                                                weights=getCall(x$reg.output$wmnomreg[[.(i)]])$weights,
                                                                evaluate = FALSE)))
          print(x$reg.output.summary$wmnomreg[[i]])
          if (i < length(x$reg.output.summary$wmnomreg)) cat("\n")
        }
      }
      if (name == "postcreg") {
        cat("# Regressions for Mediator-outcome Confounders Affected by the Exposure: \n")
        for (i in 1:length(x$reg.output.summary$postcreg)) {
          x$reg.output.summary$postcreg[[i]]$call <- eval(bquote(update(x$reg.output$postcreg[[i]], 
                                                                data=getCall(x$reg.output$postcreg[[.(i)]])$data, 
                                                                weights=getCall(x$reg.output$postcreg[[.(i)]])$weights,
                                                                evaluate = FALSE)))
          print(x$reg.output.summary$postcreg[[i]])
          if (i < length(x$reg.output.summary$postcreg)) cat("\n")
        }
      }
      cat("\n")
    }
  } else {
    for (m in 1:length(x$reg.output.summary)){ 
      cat(paste("# Regressions with Imputed Dataset", m, "\n\n"))
      regnames <- names(x$reg.output.summary[[m]])
      for (name in regnames) {
        if (name == "yreg") {
          cat("## Outcome Regression: \n")
          if (inherits(x$reg.output[[m]]$yreg, "svyglm")) {
            x$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,design = getCall(x$reg.output[[.(m)]]$yreg)$design,
                                                     family = getCall(x$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
            x$reg.output.summary[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(x$reg.output.summary[[.(m)]]$yreg$survey.design,
                                                                           data = getCall(x$reg.output.summary[[.(m)]]$yreg$survey.design)$data,
                                                                           weights = getCall(x$reg.output.summary[[.(m)]]$yreg$survey.design)$weights, 
                                                                           evaluate = FALSE))))
            print(x$reg.output.summary[[m]]$yreg)
          } else {
            x$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$reg.output[[.(m)]]$yreg,
                                                              data = getCall(x$reg.output[[.(m)]]$yreg)$data,
                                                              weights = getCall(x$reg.output[[.(m)]]$yreg)$weights, 
                                                              evaluate = FALSE)))
            print(x$reg.output.summary[[m]]$yreg)
          }
        }
        if (name == "yregTot") {
          cat("## Outcome Regression for the Total Effect: \n")
          x$reg.output.summary[[m]]$yregTot$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregTot,
                                                               data=getCall(x$reg.output[[.(m)]]$yregTot)$data,
                                                               weights=getCall(x$reg.output[[.(m)]]$yregTot)$weights, 
                                                               evaluate = FALSE)))
          print(x$reg.output.summary[[m]]$yregTot)
        }
        if (name == "yregDir") {
          cat("## Outcome Regression for the Direct Effect: \n")
          x$reg.output.summary[[m]]$yregDir$call <- eval(bquote(update(x$reg.output[[.(m)]]$yregDir,
                                                               data=getCall(x$reg.output[[.(m)]]$yregDir)$data,
                                                               weights=getCall(x$reg.output[[.(m)]]$yregDir)$weights, 
                                                               evaluate = FALSE)))
          print(x$reg.output.summary[[m]]$yregDir)
        }
        if (name == "ereg") {
          cat("## Exposure Regression for Weighting: \n")
          x$reg.output.summary[[m]]$ereg$call <- eval(bquote(update(x$reg.output[[.(m)]]$ereg,
                                                            data=getCall(x$reg.output[[.(m)]]$ereg)$data,
                                                            weights=getCall(x$reg.output[[.(m)]]$ereg)$weights, 
                                                            evaluate = FALSE)))
          print(x$reg.output.summary[[m]]$ereg)
        }
        if (name == "mreg") {
          cat("## Mediator Regressions: \n")
          for (i in 1:length(x$reg.output.summary[[m]]$mreg)) {
            if (inherits(x$reg.output[[m]]$mreg[[i]], "svyglm")) {
              x$reg.output.summary[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]],
                                                                     design = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$design,
                                                                     family = getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$family, evaluate = FALSE)))
              x$reg.output.summary[[m]]$mreg[[i]]$survey.design$call <- eval(bquote(as.call(update(x$reg.output.summary[[.(m)]]$mreg[[.(i)]]$survey.design,
                                                                                           data = getCall(x$reg.output.summary[[.(m)]]$mreg[[.(i)]]$survey.design)$data,
                                                                                           weights = getCall(x$reg.output.summary[[.(m)]]$mreg[[.(i)]]$survey.design)$weights, 
                                                                                           evaluate = FALSE))))
              print(x$reg.output.summary[[m]]$mreg[[i]])
            } else {
              x$reg.output.summary[[m]]$mreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$mreg[[.(i)]], 
                                                                     data=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$data, 
                                                                     weights=getCall(x$reg.output[[.(m)]]$mreg[[.(i)]])$weights, 
                                                                     evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$mreg[[i]])
            }
            if (i < length(x$reg.output.summary[[m]]$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          if (!is.null(x$reg.output.summary[[m]]$wmdenomreg)) {
            cat("## Mediator Regressions for Weighting (Denominator): \n")
            for (i in 1:length(x$reg.output.summary[[m]]$wmdenomreg)) {
              x$reg.output.summary[[m]]$wmdenomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmdenomreg[[i]], 
                                                                           data=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$data, 
                                                                           weights=getCall(x$reg.output[[.(m)]]$wmdenomreg[[.(i)]])$weights, 
                                                                           evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$wmdenomreg[[i]])
              if (i < length(x$reg.output.summary[[m]]$wmdenomreg)) cat("\n")
            }
          }
        }
        if (name == "wmnomreg") {
          if (!is.null(x$reg.output.summary[[m]]$wmnomreg)) {
            cat("## Mediator Regressions for Weighting (Nominator): \n")
            for (i in 1:length(x$reg.output.summary[[m]]$wmnomreg)) {
              x$reg.output.summary[[m]]$wmnomreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$wmnomreg[[i]], 
                                                                         data=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$data, 
                                                                         weights=getCall(x$reg.output[[.(m)]]$wmnomreg[[.(i)]])$weights, 
                                                                         evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$wmnomreg[[i]])
              if (i < length(x$reg.output.summary[[m]]$wmnomreg)) cat("\n")
            }
          }
        }
        if (name == "postcreg") {
          if (!is.null(x$reg.output.summary[[m]]$postcreg)) {
            cat("## Regressions for Mediator-outcome Confounders Affected by the Exposure: \n")
            for (i in 1:length(x$reg.output.summary[[m]]$postcreg)) {
              x$reg.output.summary[[m]]$postcreg[[i]]$call <- eval(bquote(update(x$reg.output[[.(m)]]$postcreg[[i]], 
                                                                         data=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$data, 
                                                                         weights=getCall(x$reg.output[[.(m)]]$postcreg[[.(i)]])$weights, 
                                                                         evaluate = FALSE)))
              print(x$reg.output.summary[[m]]$postcreg[[i]])
              if (i < length(x$reg.output.summary[[m]]$postcreg)) cat("\n")
            }
          }
        }
        cat("\n")
      }
    }
  }
  
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
  if (x$methods$inference == "bootstrap") {
    if (x$methods$boot.ci.type == "per") inf_str <- "bootstrap standard errors, percentile confidence intervals and p-values"
    if (x$methods$boot.ci.type == "bca") inf_str <- "bootstrap standard errors, bias-corrected and accelerated confidence intervals and p-values"
  }
  if (x$methods$model != "ne" && (x$methods$casecontrol)) cat("Causal Mediation Analysis for a Case Control Study via the ")
  if (!(x$methods$model != "ne" && (x$methods$casecontrol))) cat("Causal Mediation Analysis via the ")
  cat(model_str)
  cat("\n \n")
  cat(est_str)
  cat(paste(" with \n", inf_str, "\n \n"))
  printCoefmat(x$summarydf, digits = digits, has.Pvalue = TRUE)
  cat("\n")
  cat("Relevant parameter values: \n")
  print(x$ref)
}


#' Plotting Point Estimates and Confidence Intervals of Causal Effects
#' 
#' \code{ggcmest} is used to plot results of \code{cmest} nicely with plotting functions
#' in the \link{ggplot2} package. Additional layers can be added to this plot using other 
#' plotting functions in the \link{ggplot2} package.
#' 
#' @param x an object of class \code{cmest}
#' @param errorbar.width width of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.size size of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.colour colour of errorbars for confidence intervals. Default is \code{black}.
#' @param point.size size of points for point estimates. Default is \code{1}.
#' @param point.colour colour of points for point estimates. Default is \code{blue}.
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
      } else {
        if (x$methods$full) ref <- c(0, 1)
        if (!x$methods$full) ref <- 1
      }
    } else {
      if ((inherits(x$reg.output[[1]]$yreg, "lm") | inherits(x$reg.output[[1]]$yreg, "glm")) &&
          (family(x$reg.output[[1]]$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else {
        if (x$methods$full) ref <- c(0, 1)
        if (!x$methods$full) ref <- 1
      }
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
