#' Sensitivity Analysis For Unmeasured Confounding and Measurement Error
#'
#' \code{cmsens} is used to conduct sensitivity analysis for unmeasured confounding via 
#' the \emph{E-value} approach by Vanderweele et al. (2017) and Smith et al. (2019), and 
#' sensitivity analysis for measurement error via \emph{regression calibration} by Carroll 
#' et al. (1995) and \emph{SIMEX} by Cook et al. (1994) and Küchenhoff et al. (2006).
#'
#' @param object an object of class \code{cmsens}.
#' @param sens sensitivity analysis for unmeasured confounding or measurement error. 
#' \code{uc} represents unmeasured confounding and \code{me} represents measurement error.
#' See \code{Details}.
#' @param MEmethod method for measurement error correction. \code{rc} represents regression
#' calibration and \code{simex} represents SIMEX. See \code{Details}.
#' @param MEvariable variable measured with error.
#' @param MEvartype type of the variable measured with error. Can be \code{continuous} or 
#' \code{categorical} (first 3 letters are enough).
#' @param MEerror a vector of standard deviations of the measurement error (when \code{MEvartype}
#' is \code{continuous}) or a list of misclassification matrices (when \code{MEvartype}
#' is \code{categorical}). 
#' @param lambda a vector of lambdas for SIMEX. Default is \code{c(0.5, 1, 1.5, 2)}. 
#' @param B number of simulations for SIMEX. Default is \code{200}.
#' @param nboot.rc number of boots for estimating the var-cov matrix of coefficients 
#' with regression calibration. Default is \code{400}.
#' @param x an object of class \code{cmsens}
#' @param digits minimal number of significant digits. See \link{print.default}.
#' @param ... other arguments.
#' 
#' @details
#' 
#' For unmeasured confounding, all E-values are on the (risk or rate) ratio scale. If the causal 
#' effects are estimated on the difference scale, the point estimates are transformed into 
#' ratios using the transformation described by Vanderweele et al. (2017).
#' 
#' Sensitivity analysis for measurement error only supports a single variable measured with error.
#' Regression calibration only supports an independent continuous variable in a regression 
#' formula measured with error. SIMEX supports a continuous or categorical variable measured 
#' with error.
#' 
#' @return
#' If \code{sens} is \code{uc}, an object of class 'cmsens.uc' is returned:
#' \item{call}{the function call,}
#' \item{evalues}{the data frame in which the first three columns are the point estimates, 
#' lower limits of confidence intervals, upper limits of confidence intervals of causal effects on 
#' the ratio scale and the last three columns are the E-values on the ratio scale for them,}
#' If \code{sens} is \code{me}, an object of class 'cmsens.me' is returned:
#' \item{call}{the function call,}
#' \item{ME}{a list which might contain MEvariable, MEvartype, MEerror, reliability ratio, lambda and B,}
#' \item{naive}{naive causal mediation analysis results}
#' \item{sens}{a list of causal mediation analysis results after correcting errors in 
#' \code{MEerror},}
#' 
#' ...
#'
#' @seealso \code{\link{ggcmsens}}, \code{\link{cmdag}}, \code{\link{cmest}}
#'
#' @references
#' VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the 
#' E-Value (2017). Annals of Internal Medicine. 167(4): 268 - 274.
#' 
#' Smith LH, VanderWeele TJ. Mediational E-values: Approximate sensitivity analysis for 
#' unmeasured mediator-outcome confounding (2019). Epidemiology. 30(6): 835 - 837.
#' 
#' Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C. Measurement Error in Nonlinear Models: 
#' A Modern Perspective, Second Edition (2006). London: Chapman & Hall.
#' 
#' Cook JR, Stefanski LA. Simulation-extrapolation estimation in parametric measurement error 
#' models (1994). Journal of the American Statistical Association, 89(428): 1314 - 1328.
#' 
#' Küchenhoff H, Mwalili SM, Lesaffre E. A general method for dealing with misclassification 
#' in regression: the misclassification SIMEX (2006). Biometrics. 62(1): 85 - 96.
#' 
#' Stefanski LA, Cook JR. Simulation-extrapolation: the measurement error jackknife (1995). 
#' Journal of the American Statistical Association. 90(432): 1247 - 56.
#' 
#' Valeri L, Lin X, VanderWeele TJ. Mediation analysis when a continuous mediator is measured 
#' with error and the outcome follows a generalized linear model (2014). Statistics in 
#' medicine, 33(28): 4875–4890. 
#'
#' @examples
#' 
#' \dontrun{
#' library(CMAverse)
#' 
#' # 10 boots are used for illustration
#' naive <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), 
#' basec = c("C1", "C2"), EMint = TRUE,
#' mreg = list("logistic", "multinomial"), yreg = "linear",
#' astar = 0, a = 1, mval = list(0, "M2_0"),
#' estimation = "imputation", inference = "bootstrap", nboot = 10)
#' 
#' exp1 <- cmsens(object = naive, sens = "uc")
#' exp2 <- cmsens(object = naive, sens = "me", MEmethod = "rc", 
#' MEvariable = "C1", MEvartype = "con",
#' MEerror = c(0.1, 0.2))
#' summary(exp2)
#' 
#' # B = 10 is used for illustration
#' exp3 <- cmsens(object = naive, sens = "me", MEmethod = "simex", 
#' MEvariable = "M1", MEvartype = "cat",
#' MEerror = list(matrix(c(0.95,0.05,0.05,0.95), nrow = 2), 
#' matrix(c(0.9,0.1,0.1,0.9), nrow = 2)), B = 10)
#' summary(exp3)
#' }
#' 
#' @importFrom stats glm binomial poisson as.formula gaussian quasipoisson model.frame printCoefmat 
#' family sd coef vcov sigma predict rbinom rmultinom rnorm rgamma rpois weighted.mean 
#' model.matrix getCall quantile qnorm pnorm lm cov formula update
#' @importFrom nnet multinom 
#' @importFrom MASS polr glm.nb gamma.shape rnegbin
#' @importFrom survival survreg coxph
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

cmsens <- function(object = NULL, sens = "uc", MEmethod = "simex",
                   MEvariable = NULL, MEvartype = NULL, MEerror = NULL,
                   lambda = c(0.5, 1, 1.5, 2), B = 200, nboot.rc = 400) {
  
  if (is.null(object)) stop("Unspecified object")
  
  cl <- match.call()
  out <- list(call = cl)
  
  data <- object$data
  model <- object$methods$model
  reg.output <- object$reg.output
  effect.pe <- object$effect.pe
  effect.se <- object$effect.se
  outcome <- object$variables$outcome
  
  #################################################################################################
  ############################Sensitivity Analysis for Unmeasured Confounding#######################
  #################################################################################################
  
  if (sens == "uc") {
    if (model %in% c("rb", "wb", "msm", "g-formula", "ne")) {
      out.index <- 1:6
    } else if (model == "iorw") {
      out.index <- 1:3
    } 
    
    effect.pe <- effect.pe[out.index]
    effect.se <- effect.se[out.index]
    
    yreg <- switch((model == "iorw") + 1, "1" = reg.output$yreg,
                   "2" = reg.output$yregTot)
    yreg4class <- switch(("rcreg" %in% class(yreg) | "simexreg" %in% class(yreg)) + 1,
                         "1" = yreg, "2" = yreg$NAIVEreg)
    if (inherits(yreg4class, "lm")) yreg_family <- family(yreg)$family
    
    evalues <- c()
    if (inherits(yreg4class, "lm") && 
        (yreg_family %in% c("gaussian","Gamma","inverse.gaussian","quasi", "gaulss", "gevlss") |
         startsWith(yreg_family, "Tweedie") | startsWith(yreg_family, "Beta regression") |
         startsWith(yreg_family, "Scaled t"))) {
      
      # transform effect on the difference scale to RR scale
      outcome.se <- sd(data[, outcome], na.rm = TRUE)
      d <- unname(effect.pe / outcome.se)
      sd <- unname(effect.se / outcome.se)
      for (i in out.index) {
        evalues_mid <- evalues.RR(est = exp(0.91 * d[i]), lo = exp(0.91 * d[i] - 1.78 * sd[i]),
                                          hi = exp(0.91 * d[i] + 1.78 * sd[i]))
        evalues_mid <- c(evalues_mid[1, ], evalues_mid[2, ])
        names(evalues_mid) <- c("estRR", "lowerRR", "upperRR", "Evalue.estRR", "Evalue.lowerRR", "Evalue.upperRR")
        evalues <- rbind(evalues, evalues_mid)
      }
    } else {
      ci_lo <- object$effect.ci.low
      ci_up <- object$effect.ci.high
      for (i in out.index) {
        evalues_mid <- evalues.RR(est = effect.pe[i], lo = ci_lo[i], hi = ci_up[i])
        evalues_mid <- c(evalues_mid[1, ], evalues_mid[2, ])
        names(evalues_mid) <- c("estRR", "lowerRR", "upperRR", "Evalue.estRR", "Evalue.lowerRR", "Evalue.upperRR")
        evalues <- rbind(evalues, evalues_mid)
      }
    }
    rownames(evalues) <- names(effect.pe)
    out$evalues <- evalues
    class(out) <- "cmsens.uc"
    
    #################################################################################################
    ############################Sensitiviti Analysis for Measurement Error###########################
    #################################################################################################
    
  } else if (sens == "me") {
    
    if(!model %in% c("rb", "gformula")) stop("Measurement error correction currently only supports model = 'rb' or 'gformula'")
    if (length(MEvariable) == 0) stop("Unspecified MEvariable")
    if (length(MEvariable) != 1) stop("Currently only supports one variable measured with error")
    if (length(MEvartype) != length(MEvariable)) stop("length(MEvartype) != length(MEvariable)")
    
    if (MEvartype == "con") MEvartype <- "continuous"
    if (MEvartype == "cat") MEvartype <- "categorical"
    if (MEvartype == "continuous") {
      if (!is.vector(MEerror)) stop("MEerror should be a vector for a continuous variable measured with error")
    } else if (MEvartype == "categorical") {
      if (!is.list(MEerror)) stop("MEerror should be a list for a categorical variable measured with error")
    }
    
    if (MEmethod == "rc") out$ME <- list(MEvariable = MEvariable, MEvartype = MEvartype, MEerror = MEerror)
    if (MEmethod == "simex") out$ME <- list(MEvariable = MEvariable, MEvartype = MEvartype, MEerror = MEerror,
                                            lambda = lambda, B = B)
    if (MEvartype == "continuous") out$ME$reliabilityRatio <- 
        1 - MEerror/sd(data[, MEvariable], na.rm = TRUE)
    out$naive <- object
    
    n <- nrow(data)
    estimation <- object$methods$estimation
    inference <- object$methods$inference
    casecontrol <- object$methods$casecontrol
    yrare <- object$methods$yrare
    yprevalence <- object$methods$yprevalence
    full <- object$methods$full
    outcome <- object$variables$outcome
    event <- object$variables$event
    exposure <- object$variables$exposure
    mediator <- object$variables$mediator
    EMint <- object$variables$EMint
    basec <- object$variables$basec
    postc <- object$variables$postc
    multimp <- object$multimp$multimp
    a <- object$ref$a
    astar <- object$ref$asta
    mval <- object$ref$mval
    yref <- object$ref$yref
    basecval <- object$ref$basecval
    nboot <- object$methods$nboot
    nRep <- object$methods$nRep
    args_mice <- object$multimp$args_mice
    
    if (inference == "delta") variance <- TRUE
    if (inference == "bootstrap") variance <- FALSE
    sens <- list()
    
    for (i in 1:length(MEerror)) {
      
      reg.output.mid <- reg.output
      
      if (MEmethod == "rc") {
        if (MEvartype != "continuous") stop("Regression calibration only supports a continuous variable measured with error")
        for (r in 1:length(reg.output.mid)) {
          if (inherits(reg.output.mid[[r]], "list")) {
            reg.output.mid[[r]] <- lapply(1:length(reg.output.mid[[r]]), function(x)
              eval(bquote(rcreg(reg = .(reg.output.mid[[r]][[x]]), data = .(data), 
                                MEvariable = .(MEvariable), MEerror = .(MEerror[[i]]), 
                                variance = .(variance), nboot = .(nboot.rc)))))
          } else {
            reg.output.mid[[r]] <- eval(bquote(rcreg(reg = .(reg.output.mid[[r]]), 
                                                          data = .(data),
                                                          MEvariable = .(MEvariable), 
                                                          MEerror = .(MEerror[[i]]), 
                                                          variance = .(variance), 
                                                          nboot = .(nboot.rc))))
          }
        }
      } else if (MEmethod == "simex") {
        for (r in 1:length(reg.output.mid)) {
          if (inherits(reg.output.mid[[r]], "list")) {
            reg.output.mid[[r]] <- lapply(1:length(reg.output.mid[[r]]), function(x)
              eval(bquote(simexreg(reg = .(reg.output.mid[[r]][[x]]), data = .(data), 
                                   MEvariable = .(MEvariable), MEvartype = .(MEvartype), 
                                   MEerror = .(MEerror[[i]]), variance = .(variance), 
                                   lambda = .(lambda), B = .(B)))))
          } else reg.output.mid[[r]] <- eval(bquote(simexreg(reg = .(reg.output.mid[[r]]), 
                                                             data = .(data), 
                                                         MEvariable = .(MEvariable), 
                                                         MEvartype = .(MEvartype), 
                                                         MEerror = .(MEerror[[i]]), 
                                                         variance = .(variance), 
                                                         lambda = .(lambda), B = .(B))))
        }
      } else stop("Unsupported MEmethod; use 'rc' or 'simex'")
      
      yreg <- reg.output.mid$yreg
      mreg <- reg.output.mid$mreg
      ereg <- reg.output.mid$ereg
      wmreg <- reg.output.mid$wmreg
      postcreg <- reg.output.mid$postcreg
      
      # add a progress bar for bootstrap inference
      if (inference == "bootstrap") {
        env <- environment()
        counter <- 0
        progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
      }
      environment(estinf) <- environment()
      sens[[i]] <- estinf()[c("effect.pe", "effect.se", "effect.ci.low", 
                               "effect.ci.high", "effect.pval")]
      sens[[i]]$reg.output <- reg.output.mid
      
    }
    
    out$sens <- sens
    class(out) <- "cmsens.me"
    
  } else stop("Unsupported sens; use 'uc' or 'me'")
  
  return(out)
  
}

#' @describeIn cmsens Print the results of cmsens.uc nicely
#' @export
print.cmsens.uc <- function(x, ...) {
  cat("Sensitivity Analysis For Unmeasured Confounding \n")
  cat("\nEvalues on the ratio scale: \n")
  print(x$evalues)
}

#' @describeIn cmsens Print the results of cmsens.me nicely
#' @export
print.cmsens.me <- function(x, ...) {
  cat("Sensitivity Analysis For Measurement Error \n \n")
  cat("The variable measured with error: ")
  cat(x$ME$MEvariable)
  cat("\nType of the variable measured with error: ")
  cat(x$ME$MEvartype)
  cat("\n")
  for (i in 1:length(x$sens)) {
    cat(paste0("\nMeasurement error ", i, ": \n"))
    if (x$ME$MEvartype == "continuous") print(x$ME$MEerror[i])
    if (x$ME$MEvartype == "categorical") print(x$ME$MEerror[[i]])
    cat(paste0("\nMeasurement error correction for measurement error ", i, ": \n"))
    out <- data.frame(x$sens[[i]]$effect.pe, x$sens[[i]]$effect.se, 
                      x$sens[[i]]$effect.ci.low, x$sens[[i]]$effect.ci.high, 
                      x$sens[[i]]$effect.pval)
    colnames(out) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "P.val")
    print(out)
    cat("----------------------------------------------------------------\n")
  }
}

#' @describeIn cmsens Summarize the results of cmsens.me nicely
#' @export
summary.cmsens.me <- function(object, ...) {
  summarydf <- list()
  for (i in 1:length(object$sens)) {
    summarydf[[i]] <- data.frame(object$sens[[i]]$effect.pe, object$sens[[i]]$effect.se, 
                                  object$sens[[i]]$effect.ci.low, object$sens[[i]]$effect.ci.high, 
                                  object$sens[[i]]$effect.pval)
    colnames(summarydf[[i]]) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "P.val")
  }
  out <- c(object, list(summarydf = summarydf))
  class(out) <- c("summary.cmsens.me")
  return(out)
}

#' @describeIn cmsens Print the summary of cmsens.me nicely
#' @export
print.summary.cmsens.me <- function(x, digits = 4, ...) {
  cat("Sensitivity Analysis For Measurement Error \n \n")
  cat("The variable measured with error: ")
  cat(x$ME$MEvariable)
  cat("\nType of the variable measured with error: ")
  cat(x$ME$MEvartype)
  cat("\n")
  for (i in 1:length(x$sens)) {
    cat(paste0("\nMeasurement error ", i, ": \n"))
    if (x$ME$MEvartype == "continuous") print(x$ME$MEerror[i])
    if (x$ME$MEvartype == "categorical") print(x$ME$MEerror[[i]])
    cat(paste0("\nMeasurement error correction for measurement error ", i, ": \n"))
    printCoefmat(x$summarydf[[i]], digits = digits, has.Pvalue = TRUE)
    cat("----------------------------------------------------------------\n")
  }
}

#' Plotting Results of Sensitivity Analysis for Measurement Error
#' 
#' This function is used to plot the results of \code{cmsens} nicely via plotting functions
#' in the \code{ggplot2} package. Additional layers can be added to this plot using other 
#' plotting functions in the \code{ggplot2} package.
#' 
#' @param x an object of class \code{cmsens.me}
#' @param errorbar.width width of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.size size of errorbars for confidence intervals. Default is \code{0.3}.
#' @param errorbar.position position adjustment for confidence intervals, either as a string, 
#' or the result of a call to a position adjustment function. Default is 
#' \code{position_dodge2(width = 0.5)}. See \link[ggplot2]{geom_errorbar} for details.
#' @param point.size size of points for point estimates. Default is \code{1}.
#' @param point.position position adjustment for point estimates, either as a string, or 
#' the result of a call to a position adjustment function. Default is 
#' \code{position_dodge2(width = 0.3)}. See \link[ggplot2]{geom_errorbar} for details.
#' @param refline a logical value. If \code{true}, include a reference line at 
#' \code{y = 0} when effects are on the difference scale and include a reference line at 
#' \code{y = 1} when effects are on the ratio scale. Default is \code{TRUE}.
#' @param refline.colour colour of the reference line. Default is \code{red}.
#' @param refline.size size of the reference line. Default is \code{0.3}.
#' 
#' @seealso \code{\link{cmsens}}, \code{\link{ggplot2}}.
#' 
#' @examples
#' 
#' library(CMAverse)
#' library(ggplot2)
#' 
#' naive <- cmest(data = cma2020, model = "rb", outcome = "contY", 
#' exposure = "A", mediator = c("M1", "M2"), 
#' basec = c("C1", "C2"), EMint = TRUE,
#' mreg = list("logistic", "multinomial"), yreg = "linear",
#' astar = 0, a = 1, mval = list(0, "M2_0"),
#' estimation = "imputation", inference = "bootstrap", nboot = 10)
#' 
#' x <- cmsens(object = naive, sens = "me", MEmethod = "rc", 
#' MEvariable = "C1", MEvartype = "con",
#' MEerror = c(0.1, 0.2))
#' 
#' ggcmsens(x) +
#' ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45))
#' 
#' ggcmsens(x) +
#' ggplot2::coord_flip(xlim = NULL, ylim = NULL, expand = TRUE, clip = "on")
#' 
#' @export
#' 
ggcmsens <- function(x, errorbar.width = 0.3, errorbar.size = 0.3, 
                     errorbar.position = position_dodge2(width = 0.5),
                     point.size = 1, point.position = position_dodge2(width = 0.3),
                     refline = TRUE, refline.colour = "red", refline.size = 0.3) {
  # reference line
  if (refline) {
    if (!x$naive$multimp$multimp) {
      if ((inherits(x$naive$reg.output$yreg, "lm") | inherits(x$naive$reg.output$yreg, "glm")) &&
          (family(x$naive$reg.output$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else ref <- 1
    } else {
      if ((inherits(x$naive$reg.output[[1]]$yreg, "lm") | inherits(x$naive$reg.output[[1]]$yreg, "glm")) &&
          (family(x$naive$reg.output[[1]]$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else ref <- 1
    }
  } else ref <- NULL
  # naive results
  naive.pe <- x$naive$effect.pe
  naive.ci.low <- x$naive$effect.ci.low
  naive.ci.high <- x$naive$effect.ci.high
  if (x$ME$MEvartype == "continuous") {
    effect_df <- data.frame(Effect = factor(names(naive.pe), levels = names(naive.pe)),
                            Point = naive.pe,
                            CIlower = naive.ci.low,
                            CIupper = naive.ci.high,
                            ReliabilityRatio = rep(1, length(naive.pe)))
    for (i in 1:length(x$sens)) {
      pe.mid <- x$sens[[i]]$effect.pe
      ci.low.mid <- x$sens[[i]]$effect.ci.low
      ci.high.mid <- x$sens[[i]]$effect.ci.high
      effect_df <- rbind(effect_df,
                         data.frame(Effect = factor(names(pe.mid), levels = names(pe.mid)),
                                    Point = pe.mid,
                                    CIlower = ci.low.mid,
                                    CIupper = ci.high.mid,
                                    ReliabilityRatio = rep(factor(round(x$ME$reliabilityRatio[i], 2)), 
                                                           length(pe.mid))))
      
    }
    ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper,
                        colour = ReliabilityRatio), data = effect_df, 
                    width = errorbar.width, size = errorbar.size, 
                    position = errorbar.position) +
      geom_point(aes(x = Effect, y = Point, colour = ReliabilityRatio),
                 data = effect_df, size = point.size,
                 position = point.position) +
      ylab("Point Estimate and 95% CI") +
      scale_colour_hue() +
      geom_hline(yintercept = ref, color = refline.colour, size = refline.size) +
      theme(legend.position = "bottom")
  } else if (x$ME$MEvartype == "categorical") {
    pe.mid <- x$sens[[1]]$effect.pe
    ci.low.mid <- x$sens[[1]]$effect.ci.low
    ci.high.mid <- x$sens[[1]]$effect.ci.high
    effect_df <- data.frame(Effect = rep(factor(names(pe.mid), levels = names(pe.mid)), 2),
                            Point = c(pe.mid, naive.pe),
                            CIlower = c(ci.low.mid, naive.ci.low),
                            CIupper = c(ci.high.mid, naive.ci.high),
                            MC = factor(c(rep("SIMEX", length(pe.mid)), 
                                          rep("Naive", length(naive.pe)))),
                            MisclassificationMAtrix = factor(rep("MEerror[1]", 
                                                                 length(pe.mid)+length(naive.pe))))
    if (length(x$sens) > 1) {
      for (i in 2:length(x$sens)) {
        pe.mid <- x$sens[[i]]$effect.pe
        ci.low.mid <- x$sens[[i]]$effect.ci.low
        ci.high.mid <- x$sens[[i]]$effect.ci.high
        effect_df <- rbind(effect_df,
                           data.frame(Effect = rep(factor(names(pe.mid), levels = names(pe.mid)), 2),
                                      Point = c(pe.mid, naive.pe),
                                      CIlower = c(ci.low.mid, naive.ci.low),
                                      CIupper = c(ci.high.mid, naive.ci.high),
                                      MC = factor(c(rep("SIMEX", length(pe.mid)), 
                                                    rep("Naive", length(naive.pe)))),
                                      MisclassificationMAtrix = factor(rep(paste0("MEerror[", i, "]"), 
                                                                           length(pe.mid)+length(naive.pe)))))
        
      }
    }
    ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper,
                        colour = MC), data = effect_df, 
                    width = errorbar.width, size = errorbar.size, 
                    position = errorbar.position) +
      geom_point(aes(x = Effect, y = Point, colour = MC),
                 data = effect_df, size = point.size,
                 position = point.position) +
      facet_grid(MisclassificationMAtrix~.) +
      ylab("Point Estimate and 95% CI") +
      scale_colour_hue() +
      geom_hline(yintercept = ref, color = refline.colour, size = refline.size) +
      theme(legend.position = "bottom", legend.title = element_blank())
  }
}

