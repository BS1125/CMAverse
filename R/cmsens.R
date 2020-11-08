#' Sensitivity Analysis For Unmeasured Confounding and Measurement Error
#'
#' \code{cmsens} is used to conduct sensitivity analysis for unmeasured confounding via 
#' the \emph{E-value} approach by Vanderweele et al. (2017) and Smith et al. (2019), and 
#' sensitivity analysis for measurement error via \emph{regression calibration} by Carroll 
#' et al. (1995) and \emph{SIMEX} by Cook et al. (1994) and Küchenhoff et al. (2006).
#'
#' @param object an object of class \code{cmest}.
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
#' @param nboot.rc number of boots for correcting the var-cov matrix of coefficients 
#' with regression calibration. Default is \code{400}.
#' @param x an object of class \code{cmsens}
#' @param digits minimal number of significant digits. See \link{print.default}.
#' @param ... other arguments.
#' 
#' @details
#' 
#' \strong{Sensitivity Analysis for Unmeasured Confounding}
#' 
#' Currently, sensitivity analysis for unmeasured confounding are available when the outcome 
#' regression model is fitted by \link{lm}, \link{glm}, \link{glm.nb}, \link[mgcv]{gam}, 
#' \link[nnet]{multinom}, \link[MASS]{polr}.
#' 
#' All E-values are reported on the risk or rate ratio scale. If the causal effects are estimated 
#' on the difference scale (i.e., the outcome is continuous), they are transformed into risk
#' ratios using the transformation described by Vanderweele et al. (2017).
#' 
#' \strong{Sensitivity Analysis for Measurement Error}
#' 
#' Currently, sensitivity analysis for measurement error are available:
#' 
#' 1) when the regression 
#' model involving the variable measured with error is fitted by \link{lm}, 
#' \link{glm} (with family \code{gaussian}, \code{binomial} or \code{poisson}), \link[nnet]{multinom}, 
#' \link[MASS]{polr}, \link[survival]{coxph} or \link[survival]{survreg} and \code{model} is 
#' \code{rb} or \code{gformula};
#' 2) when \code{estimation} is \code{paramfunc}.
#' 
#' Sensitivity analysis for measurement error only supports a single variable measured with error.
#' Regression calibration requires that the variable measured with error be an independent 
#' continuous variable in the regression it's involved in. SIMEX supports a continuous or 
#' categorical variable measured with error. Quadratic extrapolation method is implemented for
#' SIMEX.
#' 
#' @return
#' If \code{sens} is \code{uc}, an object of class \code{cmsens.uc} is returned:
#' \item{call}{the function call,}
#' \item{evalues}{a data frame in which the first three columns are point estimates, lower limits 
#' of 95\% confidence intervals and upper limits of 95\% confidence intervals of causal effects on 
#' the risk or rate ratio scale and the last three columns are E-values on the risk or rate ratio scale,}
#' 
#' If \code{sens} is \code{me}, an object of class \code{cmsens.me} is returned:
#' \item{call}{the function call,}
#' \item{ME}{a list which might contain \code{MEmethod}, \code{MEvariable}, \code{MEvartype}, 
#' \code{MEerror}, \code{lambda}, \code{B}, \code{nboot.rc} and reliability ratio (which is calculated 
#' by \code{1 - MEerror[i]/sd(data[, MEvariable])} for \code{i=1,...,length(MEerror)} when 
#' \code{MEvartype} is \code{continuous}),}
#' \item{naive}{naive causal mediation analysis results,}
#' \item{sens}{a list of causal mediation analysis results after correcting errors in \code{MEerror},}
#' ...
#'
#' @seealso \link{ggcmsens}, \link{cmdag}, \link{cmest}, \link{rcreg}, \link{simexreg}.
#'
#' @references
#' VanderWeele TJ, Ding P (2017). Sensitivity analysis in observational research: introducing the 
#' E-Value. Annals of Internal Medicine. 167(4): 268 - 274.
#' 
#' Smith LH, VanderWeele TJ (2019). Mediational E-values: Approximate sensitivity analysis for 
#' unmeasured mediator-outcome confounding. Epidemiology. 30(6): 835 - 837.
#' 
#' Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C (2006). Measurement Error in Nonlinear Models: 
#' A Modern Perspective, Second Edition. London: Chapman & Hall.
#' 
#' Cook JR, Stefanski LA (1994). Simulation-extrapolation estimation in parametric measurement error 
#' models. Journal of the American Statistical Association, 89(428): 1314 - 1328.
#' 
#' Küchenhoff H, Mwalili SM, Lesaffre E (2006). A general method for dealing with misclassification 
#' in regression: the misclassification SIMEX. Biometrics. 62(1): 85 - 96.
#' 
#' Stefanski LA, Cook JR (1995). Simulation-extrapolation: the measurement error jackknife. 
#' Journal of the American Statistical Association. 90(432): 1247 - 56.
#' 
#' Valeri L, Lin X, VanderWeele TJ (2014). Mediation analysis when a continuous mediator is measured 
#' with error and the outcome follows a generalized linear model. Statistics in 
#' medicine, 33(28): 4875 – 4890. 
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
#' @importFrom stats model.frame printCoefmat family sd getCall
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom EValue evalues.RR
#' @importFrom simex check.mc.matrix
#' @importFrom ggplot2 ggplot geom_errorbar aes geom_point ylab geom_hline position_dodge2 
#' scale_colour_hue theme element_blank facet_grid 
#' 
#' @export

cmsens <- function(object = NULL, sens = "uc", MEmethod = "simex",
                   MEvariable = NULL, MEvartype = NULL, MEerror = NULL,
                   lambda = c(0.5, 1, 1.5, 2), B = 200, nboot.rc = 400) {
  if (is.null(object)) stop("Unspecified object")
  cl <- match.call()
  out <- list(call = cl)
  # extract objects from cmest object
  data <- object$data
  model <- object$methods$model
  reg.input <- object$reg.input
  reg.output <- object$reg.output
  effect.pe <- object$effect.pe
  effect.se <- object$effect.se
  outcome <- object$variables$outcome
  
  #################################################################################################
  ############################Sensitivity Analysis for Unmeasured Confounding#######################
  #################################################################################################
  if (sens == "uc") {
    if (inherits(reg.input$yreg, "coxph") | inherits(reg.input$yreg, "survreg") | 
        (is.character(reg.input$yreg) && reg.input$yreg %in% c("coxph", "aft_exp", "aft_weibull"))
        ) stop("sensitivity analysis for unmeasured confounding currently doesn't support survival outcomes")
    # report evalues for cde, pnde, tnde, pnie, tnie, te when model is not iorw
    if (model %in% c("rb", "wb", "msm", "g-formula", "ne")) out.index <- 1:6
    # report evalues for te, pnde, tnie when model is iorw
    if (model == "iorw") out.index <- 1:3
    
    effect.pe <- effect.pe[out.index]
    effect.se <- effect.se[out.index]
    yreg <- switch((model == "iorw") + 1, "1" = reg.output$yreg, "2" = reg.output$yregTot)
    yreg4class <- switch(("rcreg" %in% class(yreg) | "simexreg" %in% class(yreg)) + 1,
                         "1" = yreg, "2" = yreg$NAIVEreg)
    if (inherits(yreg4class, "lm")) yreg_family <- family(yreg)$family
    
    evalues <- c()
    if (inherits(yreg4class, "lm") && 
        (yreg_family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
      # transform effect on the difference scale to RR scale
      outcome.se <- sd(data[, outcome], na.rm = TRUE)
      d <- unname(effect.pe / outcome.se)
      sd <- unname(effect.se / outcome.se)
      for (i in out.index) {
        evalues_mid <- evalues.RR(est = exp(0.91 * d[i]), lo = exp(0.91 * d[i] - 1.78 * sd[i]),
                                  hi = exp(0.91 * d[i] + 1.78 * sd[i]))
        evalues_mid <- c(evalues_mid[1, ], evalues_mid[2, ])
        evalues <- rbind(evalues, evalues_mid)
      }
    } else {
      ci_lo <- object$effect.ci.low
      ci_hi <- object$effect.ci.high
      for (i in out.index) {
        evalues_mid <- evalues.RR(est = effect.pe[i], lo = ci_lo[i], hi = ci_hi[i])
        evalues_mid <- c(evalues_mid[1, ], evalues_mid[2, ])
        evalues <- rbind(evalues, evalues_mid)
      }
    }
    colnames(evalues) <- c("estRR", "lowerRR", "upperRR", "Evalue.estRR", "Evalue.lowerRR", "Evalue.upperRR")
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
      if (!is.vector(MEerror)) stop("MEerror should be a vector of error standard deviations for a continuous variable measured with error")
    } else if (MEvartype == "categorical") {
      if (!is.list(MEerror)) stop("MEerror should be a list of misclassification matrices for a categorical variable measured with error")
    }
    
    if (MEmethod == "rc") out$ME <- list(MEmethod = MEmethod, MEvariable = MEvariable, MEvartype = MEvartype, MEerror = MEerror)
    if (MEmethod == "simex") out$ME <- list(MEmethod = MEmethod, MEvariable = MEvariable, MEvartype = MEvartype, MEerror = MEerror,
                                            lambda = lambda, B = B)
    if (MEvartype == "continuous") out$ME$reliabilityRatio <- 
        1 - MEerror/sd(data[, MEvariable], na.rm = TRUE)
    out$naive <- object
    
    n <- nrow(data)
    full <- object$methods$full   
    casecontrol <- object$methods$casecontrol
    yrare <- object$methods$yrare
    yprevalence <- object$methods$yprevalence    
    estimation <- object$methods$estimation
    inference <- object$methods$inference
    outcome <- object$variables$outcome
    event <- object$variables$event
    exposure <- object$variables$exposure
    mediator <- object$variables$mediator
    EMint <- object$variables$EMint
    basec <- object$variables$basec
    postc <- object$variables$postc
    yreg <- reg.input$yreg
    ereg <- reg.input$ereg
    mreg <- reg.input$mreg
    wmnomreg <- reg.input$wmnomreg
    wmdenomreg <- reg.input$wmdenomreg
    postcreg <- reg.input$postcreg
    a <- object$ref$a
    astar <- object$ref$astar
    mval <- object$ref$mval
    yval <- object$ref$yval
    basecval <- object$ref$basecval
    nboot <- object$methods$nboot
    boot.ci.type <- object$methods$boot.ci.type
    nRep <- object$methods$nRep    
    multimp <- object$multimp$multimp
    args_mice <- object$multimp$args_mice
    
    # run regressions
    environment(regrun) <- environment()
    regs <- regrun()
    
    if (inference == "delta") {
      variance <- TRUE
      if (MEmethod == "rc") out$ME$nboot.rc <- nboot.rc
    } else variance <- FALSE
    
    sens <- list()
    for (i in 1:length(MEerror)) {
      regs_mid <- regs
      if (MEmethod == "rc") {
        if (MEvartype != "continuous") stop("Regression calibration only supports a continuous variable measured with error")
        for (r in 1:length(regs_mid)) {
          if (!is.null(regs_mid[[r]])) {
            if (inherits(regs_mid[[r]], "list")) {
              regs_mid[[r]] <- lapply(1:length(regs_mid[[r]]), function(x)
                eval(bquote(rcreg(reg = .(regs_mid[[r]][[x]]), formula = .(formula(regs_mid[[r]][[x]])),
                                  data = .(data), 
                                  weights = .(model.frame(regs_mid[[r]][[x]])$'(weights)'),
                                  MEvariable = .(MEvariable), MEerror = .(MEerror[[i]]), 
                                  variance = .(variance), nboot = .(nboot.rc)))))
            } else {
              regs_mid[[r]] <- eval(bquote(rcreg(reg = .(regs_mid[[r]]), 
                                                 formula = .(formula(regs_mid[[r]])), data = .(data),
                                                 weights = .(model.frame(regs_mid[[r]])$'(weights)'),
                                                 MEvariable = .(MEvariable), MEerror = .(MEerror[[i]]), 
                                                 variance = .(variance), nboot = .(nboot.rc))))
            }
          }
        }
      } else if (MEmethod == "simex") {
        for (r in 1:length(regs_mid)) {
          if (!is.null(regs_mid[[r]])) {
            if (inherits(regs_mid[[r]], "list")) {
              regs_mid[[r]] <- lapply(1:length(regs_mid[[r]]), function(x)
                eval(bquote(simexreg(reg = .(regs_mid[[r]][[x]]), 
                                     formula = .(formula(regs_mid[[r]][[x]])), data = .(data), 
                                     weights = .(model.frame(regs_mid[[r]][[x]])$'(weights)'),
                                     MEvariable = .(MEvariable), MEvartype = .(MEvartype), 
                                     MEerror = .(MEerror[[i]]), variance = .(variance), 
                                     lambda = .(lambda), B = .(B)))))
            } else regs_mid[[r]] <- eval(bquote(simexreg(reg = .(regs_mid[[r]]), 
                                                         formula = .(formula(regs_mid[[r]])), data = .(data), 
                                                         weights = .(model.frame(regs_mid[[r]])$'(weights)'),
                                                         MEvariable = .(MEvariable), MEvartype = .(MEvartype), 
                                                         MEerror = .(MEerror[[i]]), variance = .(variance), 
                                                         lambda = .(lambda), B = .(B))))
          }
        }
      } else stop("Unsupported MEmethod; use 'rc' or 'simex'")
      
      yreg <- regs_mid$yreg
      mreg <- regs_mid$mreg
      ereg <- regs_mid$ereg
      wmnomreg <- regs_mid$wmnomreg
      wmdenomreg <- regs_mid$wmdenomreg
      postcreg <- regs_mid$postcreg
      # add a progress bar for bootstrap inference
      if (inference == "bootstrap") {
        env <- environment()
        counter <- 0
        progbar <- txtProgressBar(min = 0, max = nboot, style = 3)
      }
      environment(estinf) <- environment()
      estinf_res <- estinf()
      sens[[i]] <- estinf_res[c("effect.pe", "effect.se", "effect.ci.low", "effect.ci.high", "effect.pval")]
      sens[[i]]$reg.output <- estinf_res$reg.output
    }
    
    out$sens <- sens
    class(out) <- "cmsens.me"
  } else stop("Unsupported sens; use 'uc' or 'me'")
  
  return(out)
  
}


#' @describeIn cmsens Print results of \code{cmsens.uc} nicely
#' @export
print.cmsens.uc <- function(x, ...) {
  cat("Sensitivity Analysis For Unmeasured Confounding \n")
  cat("\nEvalues on the risk or rate ratio scale: \n")
  print(x$evalues)
}


#' @describeIn cmsens Print results of \code{cmsens.me} nicely
#' @export
print.cmsens.me <- function(x, ...) {
  cat("Sensitivity Analysis For Measurement Error \n \n")
  cat("The variable measured with error: ")
  cat(x$ME$MEvariable)
  cat("\nType of the variable measured with error: ")
  cat(x$ME$MEvartype)
  cat("\n\n")
  
  for (i in 1:length(x$sens)) {
    cat(paste0("# Measurement error ", i, ": \n"))
    if (x$ME$MEvartype == "continuous") print(x$ME$MEerror[i])
    if (x$ME$MEvartype == "categorical") print(x$ME$MEerror[[i]])
    cat(paste0("\n## Error-corrected regressions for measurement error ", i, ": \n\n"))
    # print error corrected regression models
    if (!x$naive$multimp$multimp) {
      regnames <- names(x$sens[[i]]$reg.output)
      for (name in regnames) {
        if (name == "yreg") {
          cat("### Outcome regression:\n")
          if (inherits(x$sens[[i]]$reg.output$yreg, "svyglm")) {
            x$sens[[i]]$reg.output$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yreg,design = getCall(x$sens[[.(i)]]$reg.output$yreg)$design,
                                                                   family = getCall(x$sens[[.(i)]]$reg.output$yreg)$family, evaluate = FALSE)))
            x$sens[[i]]$reg.output$yreg$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output$yreg)$survey.design,
                                                                                         data = getCall(summary(x$sens[[.(i)]]$reg.output$yreg)$survey.design)$data,
                                                                                         weights = getCall(summary(x$sens[[.(i)]]$reg.output$yreg)$survey.design)$weights, 
                                                                                         evaluate = FALSE))))
            print(x$sens[[i]]$reg.output$yreg)
          } else if (inherits(x$sens[[i]]$reg.output$yreg, "rcreg")|inherits(x$sens[[i]]$reg.output$yreg, "simexreg")) {
            x$sens[[i]]$reg.output$yreg$call <- eval(bquote(update(x$sens[[i]]$reg.output$yreg, 
                                                                   reg = getCall(x$sens[[.(i)]]$reg.output$yreg)$reg,
                                                                   data=getCall(x$sens[[.(i)]]$reg.output$yreg)$data,
                                                                   weights=getCall(x$sens[[.(i)]]$reg.output$yreg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yreg)
          } else {
            x$sens[[i]]$reg.output$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yreg,
                                                                   data=getCall(x$sens[[.(i)]]$reg.output$yreg)$data,
                                                                   weights=getCall(x$sens[[.(i)]]$reg.output$yreg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yreg)
          }
        }
        if (name == "yregTot") {
          cat("### Outcome regression for the total effect: \n")
          if (inherits(x$sens[[i]]$reg.output$yregTot, "rcreg")|inherits(x$sens[[i]]$reg.output$yregTot, "simexreg")) {
            x$sens[[i]]$reg.output$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregTot, reg = getCall(x$sens[[.(i)]]$reg.output$yregTot)$reg,
                                                                      data=getCall(x$sens[[.(i)]]$reg.output$yregTot)$data,
                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$yregTot)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregTot)
          } else {
            x$sens[[i]]$reg.output$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregTot,data=getCall(x$sens[[.(i)]]$reg.output$yregTot)$data,
                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$yregTot)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregTot)
          }
        }
        if (name == "yregDir") {
          cat("### Outcome regression for the direct effect: \n")
          if (inherits(x$sens[[i]]$reg.output$yregDir, "rcreg")|inherits(x$sens[[i]]$reg.output$yregDir, "simexreg")) {
            x$sens[[i]]$reg.output$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregDir, reg = getCall(x$sens[[.(i)]]$reg.output$yregDir)$reg,
                                                                      data=getCall(x$sens[[.(i)]]$reg.output$yregDir)$data,
                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$yregDir)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregDir)
          } else {
            x$sens[[i]]$reg.output$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregDir,data=getCall(x$sens[[.(i)]]$reg.output$yregDir)$data,
                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$yregDir)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$yregDir)
          }
        }
        if (name == "ereg") {
          cat("### Exposure regression for weighting: \n")
          if (inherits(x$sens[[i]]$reg.output$ereg, "rcreg")|inherits(x$sens[[i]]$reg.output$ereg, "simexreg")) {
            x$sens[[i]]$reg.output$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$ereg, reg = getCall(x$sens[[.(i)]]$reg.output$ereg)$reg,
                                                                   data=getCall(x$sens[[.(i)]]$reg.output$ereg)$data,
                                                                   weights=getCall(x$sens[[.(i)]]$reg.output$ereg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$ereg)
          } else {
            x$sens[[i]]$reg.output$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$ereg,data=getCall(x$sens[[.(i)]]$reg.output$ereg)$data,
                                                                   weights=getCall(x$sens[[.(i)]]$reg.output$ereg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output$ereg)
          }
        }
        if (name == "mreg") {
          cat("### Mediator regressions: \n")
          for (j in 1:length(x$sens[[i]]$reg.output$mreg)) {
            if (inherits(x$sens[[i]]$reg.output$mreg[[j]], "svyglm")) {
              x$sens[[i]]$reg.output$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]],
                                                                          design = getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$design,
                                                                          family = getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$family, evaluate = FALSE)))
              x$sens[[i]]$reg.output$mreg[[j]]$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$survey.design,
                                                                                                data = getCall(summary(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$survey.design)$data,
                                                                                                weights = getCall(summary(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$survey.design)$weights, 
                                                                                                evaluate = FALSE))))
              print(x$sens[[i]]$reg.output$mreg[[j]])
            } else if (inherits(x$sens[[i]]$reg.output$mreg[[j]], "rcreg") |
                       inherits(x$sens[[i]]$reg.output$mreg[[j]], "simexreg")){
              x$sens[[i]]$reg.output$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], 
                                                                          reg=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$reg,
                                                                          data=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$data, 
                                                                          weights=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$weights,
                                                                          evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$mreg[[j]])
            } else {
              x$sens[[i]]$reg.output$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], 
                                                                          data=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$data, 
                                                                          weights=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$weights,
                                                                          evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$mreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          cat("### Mediator regressions for weighting (denominator): \n")
          for (j in 1:length(x$sens[[i]]$reg.output$wmdenomreg)) {
            if (inherits(x$sens[[i]]$reg.output$wmdenomreg[[j]], "rcreg") |
                inherits(x$sens[[i]]$reg.output$wmdenomreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]], 
                                                                                reg=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$reg,
                                                                                data=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$data, 
                                                                                weights=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$weights,
                                                                                evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmdenomreg[[j]])
            } else {
              x$sens[[i]]$reg.output$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]], 
                                                                                data=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$data, 
                                                                                weights=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$weights,
                                                                                evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmdenomreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$wmdenomreg)) cat("\n")
          }
        }
        if (name == "wmnomreg") {
          cat("### Mediator regressions for weighting (nominator): \n")
          for (j in 1:length(x$sens[[i]]$reg.output$wmnomreg)) {
            if (inherits(x$sens[[i]]$reg.output$wmnomreg[[j]], "rcreg") |
                inherits(x$sens[[i]]$reg.output$wmnomreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]], 
                                                                              reg=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$reg,
                                                                              data=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$data, 
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$weights,
                                                                              evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmnomreg[[j]])
            } else {
              x$sens[[i]]$reg.output$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]], 
                                                                              data=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$data, 
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$weights,
                                                                              evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$wmnomreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$wmnomreg)) cat("\n")
          }
        }
        if (name == "postcreg") {
          cat("### Regressions for mediator-outcome confounders affected by the exposure: \n")
          for (j in 1:length(x$sens[[i]]$reg.output$postcreg)) {
            if (inherits(x$sens[[i]]$reg.output$postcreg[[j]], "rcreg") |
                inherits(x$sens[[i]]$reg.output$postcreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]], 
                                                                              reg=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$reg,
                                                                              data=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$data, 
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$weights,
                                                                              evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$postcreg[[j]])
            } else {
              x$sens[[i]]$reg.output$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]], 
                                                                              data=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$data, 
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$weights,
                                                                              evaluate = FALSE)))
              print(x$sens[[i]]$reg.output$postcreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output$postcreg)) cat("\n")
          }
        }
        cat("\n")
      }
    } else {
      for (m in 1:length(x$sens[[i]]$reg.output)) { 
        cat(paste("### Regressions with imputed dataset", m, "\n\n"))
        regnames <- names(x$sens[[i]]$reg.output[[m]])
        for (name in regnames) {
          if (name == "yreg") {
            cat("#### Outcome regression: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yreg, "svyglm")) {
              x$sens[[i]]$reg.output[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg,
                                                                          design = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$design,
                                                                          family = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
              x$sens[[i]]$reg.output[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$survey.design,
                                                                                                data = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$survey.design)$data,
                                                                                                weights = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$survey.design)$weights, 
                                                                                                evaluate = FALSE))))
              print(x$sens[[i]]$reg.output[[m]]$yreg)
            } else if (inherits(x$sens[[i]]$reg.output[[m]]$yreg, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yreg, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg,
                                                                          reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$reg,
                                                                          data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$data,
                                                                          weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$weights, 
                                                                          evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yreg)
            } else {
              x$sens[[i]]$reg.output[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg,
                                                                          data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$data,
                                                                          weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$weights, 
                                                                          evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yreg)
            }
          }
          if (name == "yregTot") {
            cat("#### Outcome regression for the total effect: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yregTot, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yregTot, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot,
                                                                             reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$reg,
                                                                             data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$data,
                                                                             weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$weights, 
                                                                             evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregTot)
            } else {
              x$sens[[i]]$reg.output[[m]]$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot,
                                                                             data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$data,
                                                                             weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$weights, 
                                                                             evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregTot)
            }
          }
          if (name == "yregDir") {
            cat("#### Outcome regression for the direct effect: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yregDir, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yregDir, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir,
                                                                             reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$reg,
                                                                             data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$data,
                                                                             weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$weights, 
                                                                             evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregDir)
            } else {
              x$sens[[i]]$reg.output[[m]]$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir,
                                                                             data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$data,
                                                                             weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$weights, 
                                                                             evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$yregDir)
            }
          }
          if (name == "ereg") {
            cat("#### Exposure regression for weighting: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$ereg, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$ereg, "simexreg")){
              x$sens[[i]]$reg.output[[m]]$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg,
                                                                          reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$reg,
                                                                          data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$data,
                                                                          weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$weights, 
                                                                          evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$ereg)
            } else {
              x$sens[[i]]$reg.output[[m]]$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg,
                                                                          data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$data,
                                                                          weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$weights, 
                                                                          evaluate = FALSE)))
              print(x$sens[[i]]$reg.output[[m]]$ereg)
            }
          }
          if (name == "mreg") {
            cat("#### Mediator regressions: \n")
            for (j in 1:length(x$sens[[i]]$reg.output[[m]]$mreg)) {
              if (inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "svyglm")) {
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]],
                                                                                 design = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$design,
                                                                                 family = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$family, evaluate = FALSE)))
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$survey.design$call <- eval(bquote(as.call(update(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$survey.design,
                                                                                                       data = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$survey.design)$data,
                                                                                                       weights = getCall(summary(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$survey.design)$weights, 
                                                                                                       evaluate = FALSE))))
                print(x$sens[[i]]$reg.output[[m]]$mreg[[j]])
              } else if (inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "simexreg")){
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], 
                                                                                 reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$reg,
                                                                                 data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$data, 
                                                                                 weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$weights, 
                                                                                 evaluate = FALSE)))
                print(x$sens[[i]]$reg.output[[m]]$mreg[[j]])
              } else {
                x$sens[[i]]$reg.output[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], 
                                                                                 data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$data, 
                                                                                 weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$weights, 
                                                                                 evaluate = FALSE)))
                print(x$sens[[i]]$reg.output[[m]]$mreg[[j]])
              }
              if (j < length(x$sens[[i]]$reg.output[[m]]$mreg)) cat("\n")
            }
          }
          if (name == "wmdenomreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) {
              cat("#### Mediator regressions for weighting (denominator): \n")
              for (j in 1:length(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]], 
                                                                                         reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$reg,
                                                                                         data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$data, 
                                                                                         weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$weights, 
                                                                                         evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]])
                } else {
                  x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]], 
                                                                                         data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$data, 
                                                                                         weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$weights, 
                                                                                         evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) cat("\n")
              }
            }
          }
          if (name == "wmnomreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$wmnomreg)) {
              cat("#### Mediator regressions for weighting (nominator): \n")
              for (j in 1:length(x$sens[[i]]$reg.output[[m]]$wmnomreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]], 
                                                                                       reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$reg,
                                                                                       data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$data, 
                                                                                       weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$weights, 
                                                                                       evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]])
                } else {
                  x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]], 
                                                                                       data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$data, 
                                                                                       weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$weights, 
                                                                                       evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output[[m]]$wmnomreg)) cat("\n")
              }
            }
          }
          if (name == "postcreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$postcreg)) {
              cat("#### Regressions for mediator-outcome confounders affected by the exposure: \n")
              for (j in 1:length(x$sens[[i]]$reg.output[[m]]$postcreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$postcreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$postcreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output[[m]]$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]], 
                                                                                       reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$reg,
                                                                                       data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$data, 
                                                                                       weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$weights, 
                                                                                       evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$postcreg[[j]])
                } else {
                  x$sens[[i]]$reg.output[[m]]$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]], 
                                                                                       data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$data, 
                                                                                       weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$weights, 
                                                                                       evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output[[m]]$postcreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output[[m]]$postcreg)) cat("\n")
              }
            }
          }
          cat("\n")
        }
      }
    }
    
    # scale and legend
    full <- x$naive$methods$full
    model <- x$naive$methods$model
    EMint <- x$naive$variables$EMint
    if (model == "iorw") {
      if (x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output[[1]]$yregTot
      if (!x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output$yregTot
    } else {
      if (x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output[[1]]$yreg
      if (!x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output$yreg
    }
    if (inherits(yreg_mid, "rcreg") | inherits(yreg_mid, "simexreg")) yreg_mid <- yreg_mid$NAIVEreg
    is_lm <- inherits(yreg_mid, "lm")
    is_glm <- inherits(yreg_mid, "glm")
    is_svyglm <- inherits(yreg_mid, "svyglm")
    is_gam <- inherits(yreg_mid, "gam")
    if (is_lm | is_glm) family_yreg <- family(yreg_mid)
    is_multinom <- inherits(yreg_mid, "multinom")
    is_svymultinom <- inherits(yreg_mid, "svymultinom")
    is_polr <- inherits(yreg_mid, "polr")
    is_survreg <- inherits(yreg_mid, "survreg")
    is_coxph <- inherits(yreg_mid, "coxph")
    if ((is_lm | is_glm) && (family_yreg$family %in% c("gaussian", "inverse.gaussian", "quasi", "Gamma"))) {
      scale <- "mean difference scale"
      if (model == "iorw") {
        if (full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect; pm: proportion mediated)"   
        if (!full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect)"
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rintref: randomized analogue of reference interaction; rintmed: randomized analogue of mediated interaction; cde(prop): proportion cde; rintref(prop): proportion rintref; rintmed(prop): proportion rintmed; rpnie(prop): proportion rpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect)"
      } else {
        if (full) {
          if (EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; intref: reference interaction; intmed: mediated interaction; cde(prop): proportion cde; intref(prop): proportion intref; intmed(prop): proportion intmed; pnie(prop): proportion pnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; pm: overall proportion mediated)"
        } else legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect)"
      }
    } else if ((is_lm | is_glm) && (family_yreg$family %in% c("poisson", "quasipoisson", "ziplss") |
                                    startsWith(family_yreg$family, "Negative Binomial") |
                                    startsWith(family_yreg$family, "Zero inflated Poisson"))) {
      scale <- "rate ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio)" 
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; rERintref: randomized analogue of excess relative rate due to reference interaction; rERintmed: randomized analogue of excess relative rate due to mediated interaction; rERpnie: randomized analogue of excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; ERintref: excess relative rate due to reference interaction; ERintmed: excess relative rate due to mediated interaction; ERpnie: excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio)"
      }
    } else if (((is_lm | is_glm) && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                                     startsWith(family_yreg$family, "Ordered Categorical"))) |
               is_multinom | is_polr) {
      scale <- "risk ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio)" 
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; rERintref: randomized analogue of excess relative risk due to reference interaction; rERintmed: randomized analogue of excess relative risk due to mediated interaction; rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; ERintref: excess relative risk due to reference interaction; ERintmed: excess relative risk due to mediated interaction; ERpnie: excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      }
    } else if (is_coxph) {
      scale <- "hazard ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio)" 
      }  else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; rERintref: randomized analogue of excess relative hazard due to reference interaction; rERintmed: randomized analogue of excess relative hazard due to mediated interaction; rERpnie: randomized analogue of excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; ERintref: excess relative hazard due to reference interaction; ERintmed: excess relative hazard due to mediated interaction; ERpnie: excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
      }
    } else if (is_survreg) {
      scale <- "mean survival scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio)" 
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; rERintref: randomized analogue of excess mean survival ratio due to reference interaction; rERintmed: randomized analogue of excess mean survival ratio due to mediated interaction; rERpnie: randomized analogue of excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; ERintref: excess mean survival ratio due to reference interaction; ERintmed: excess mean survival ratio due to mediated interaction; ERpnie: excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
      }
    }
    
    # print error-corrected results
    cat(paste0(paste("\n## Error-corrected causal effects on the", scale, "for measurement error "), i, ": \n"))
    print(x$sens[[i]]$effect.pe)
    cat("----------------------------------------------------------------\n")
  }
  cat("\n")
  cat(legend)
  cat("\n\nRelevant variable values: \n")
  print(x$naive$ref)
}


#' @describeIn cmsens Summarize results of \code{cmsens.me} nicely
#' @export
summary.cmsens.me <- function(object, ...) {
  out <- object
  for (i in 1:length(object$sens)) { 
    out$sens[[i]]$reg.output.summary <- out$sens[[i]]$reg.output
    if (!object$naive$multimp$multimp) {
      regnames <- names(object$sens[[i]]$reg.output)
      for (name in regnames) {
        if (name == "yreg") out$sens[[i]]$reg.output.summary$yreg <- summary(object$sens[[i]]$reg.output$yreg)
        if (name == "yregTot") out$sens[[i]]$reg.output.summary$yregTot <- summary(object$sens[[i]]$reg.output$yregTot)
        if (name == "yregDir") out$sens[[i]]$reg.output.summary$yregDir <- summary(object$sens[[i]]$reg.output$yregDir)
        if (name == "ereg") out$sens[[i]]$reg.output.summary$ereg <- summary(object$sens[[i]]$reg.output$ereg)
        if (name == "mreg") out$sens[[i]]$reg.output.summary$mreg <- lapply(1:length(object$sens[[i]]$reg.output$mreg), function(i) 
          summary(object$sens[[i]]$reg.output$mreg[[i]]))
        if (name == "wmnomreg") out$sens[[i]]$reg.output.summary$wmnomreg <- lapply(1:length(object$sens[[i]]$reg.output$wmnomreg), function(i) 
          summary(object$sens[[i]]$reg.output$wmnomreg[[i]]))
        if (name == "wmdenomreg") out$sens[[i]]$reg.output.summary$wmdenomreg <- lapply(1:length(object$sens[[i]]$reg.output$wmdenomreg), function(i) 
          summary(object$sens[[i]]$reg.output$wmdenomreg[[i]]))
        if (name == "postcreg") out$sens[[i]]$reg.output.summary$postcreg <- lapply(1:length(object$sens[[i]]$reg.output$postcreg), function(i) 
          summary(object$sens[[i]]$reg.output$postcreg[[i]]))
      }
    } else {
      for (m in 1:length(object$sens[[i]]$reg.output)){ 
        regnames <- names(object$sens[[i]]$reg.output[[m]])
        for (name in regnames) {
          if (name == "yreg") out$sens[[i]]$reg.output.summary[[m]]$yreg <- summary(object$sens[[i]]$reg.output[[m]]$yreg)
          if (name == "yregTot") out$sens[[i]]$reg.output.summary[[m]]$yregTot <- summary(object$sens[[i]]$reg.output[[m]]$yregTot)
          if (name == "yregDir") out$sens[[i]]$reg.output.summary[[m]]$yregDir <- summary(object$sens[[i]]$reg.output[[m]]$yregDir)
          if (name == "ereg") out$sens[[i]]$reg.output.summary[[m]]$ereg <- summary(object$sens[[i]]$reg.output[[m]]$ereg)
          if (name == "mreg") out$sens[[i]]$reg.output.summary[[m]]$mreg <- lapply(1:length(object$sens[[i]]$reg.output[[m]]$mreg), function(i) 
            summary(object$sens[[i]]$reg.output[[m]]$mreg[[i]]))
          if (name == "wmnomreg") out$sens[[i]]$reg.output.summary[[m]]$wmnomreg <- lapply(1:length(object$sens[[i]]$reg.output[[m]]$wmnomreg), function(i) 
            summary(object$sens[[i]]$reg.output[[m]]$wmnomreg[[i]]))
          if (name == "wmdenomreg") out$sens[[i]]$reg.output.summary[[m]]$wmdenomreg <- lapply(1:length(object$sens[[i]]$reg.output[[m]]$wmdenomreg), function(i) 
            summary(object$sens[[i]]$reg.output[[m]]$wmdenomreg[[i]]))
          if (name == "postcreg") out$sens[[i]]$reg.output.summary[[m]]$postcreg <- lapply(1:length(object$sens[[i]]$reg.output[[m]]$postcreg), function(i) 
            summary(object$sens[[i]]$reg.output[[m]]$postcreg[[i]]))
        }
      }
    }
  }
  
  summarydf <- list()
  for (i in 1:length(object$sens)) {
    summarydf[[i]] <- data.frame(object$sens[[i]]$effect.pe, object$sens[[i]]$effect.se, 
                                 object$sens[[i]]$effect.ci.low, object$sens[[i]]$effect.ci.high, 
                                 object$sens[[i]]$effect.pval)
    colnames(summarydf[[i]]) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "P.val")
  }
  out$summarydf <- summarydf
  class(out) <- c("summary.cmsens.me")
  return(out)
}


#' @describeIn cmsens Print the summary of \code{cmsens.me} nicely
#' @export
print.summary.cmsens.me <- function(x, digits = 4, ...) {
  cat("Sensitivity Analysis For Measurement Error \n \n")
  cat("The variable measured with error: ")
  cat(x$ME$MEvariable)
  cat("\nType of the variable measured with error: ")
  cat(x$ME$MEvartype)
  cat("\n")
  for (i in 1:length(x$sens)) {
    cat(paste0("\n# Measurement error ", i, ": \n"))
    if (x$ME$MEvartype == "continuous") print(x$ME$MEerror[i])
    if (x$ME$MEvartype == "categorical") print(x$ME$MEerror[[i]])
    cat(paste0("\n## Error-corrected regressions for measurement error ", i, ": \n\n"))
    # print error corrected regression models
    if (!x$naive$multimp$multimp) {
      regnames <- names(x$sens[[i]]$reg.output)
      for (name in regnames) {
        if (name == "yreg") {
          cat("### Outcome regression:\n")
          if (inherits(x$sens[[i]]$reg.output$yreg, "svyglm")) {
            x$sens[[i]]$reg.output.summary$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yreg,design = getCall(x$sens[[.(i)]]$reg.output$yreg)$design,
                                                                           family = getCall(x$sens[[.(i)]]$reg.output$yreg)$family, evaluate = FALSE)))
            x$sens[[i]]$reg.output.summary$yreg$survey.design$call <- eval(bquote(as.call(update(x$sens[[.(i)]]$reg.output.summary$yreg$survey.design,
                                                                                                 data = getCall(x$sens[[.(i)]]$reg.output.summary$yreg$survey.design)$data,
                                                                                                 weights = getCall(x$sens[[.(i)]]$reg.output.summary$yreg$survey.design)$weights, 
                                                                                                 evaluate = FALSE))))
            print(x$sens[[i]]$reg.output.summary$yreg)
          } else if (inherits(x$sens[[i]]$reg.output$yreg, "rcreg")|inherits(x$sens[[i]]$reg.output$yreg, "simexreg")) {
            x$sens[[i]]$reg.output.summary$yreg$call <- eval(bquote(update(x$sens[[i]]$reg.output$yreg, 
                                                                           reg = getCall(x$sens[[.(i)]]$reg.output$yreg)$reg,
                                                                           data=getCall(x$sens[[.(i)]]$reg.output$yreg)$data,
                                                                           weights=getCall(x$sens[[.(i)]]$reg.output$yreg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$yreg)
          } else {
            x$sens[[i]]$reg.output.summary$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yreg,
                                                                           data=getCall(x$sens[[.(i)]]$reg.output$yreg)$data,
                                                                           weights=getCall(x$sens[[.(i)]]$reg.output$yreg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$yreg)
          }
        }
        if (name == "yregTot") {
          cat("### Outcome regression for the total effect: \n")
          if (inherits(x$sens[[i]]$reg.output$yregTot, "rcreg")|inherits(x$sens[[i]]$reg.output$yregTot, "simexreg")) {
            x$sens[[i]]$reg.output.summary$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregTot, reg = getCall(x$sens[[.(i)]]$reg.output$yregTot)$reg,
                                                                              data=getCall(x$sens[[.(i)]]$reg.output$yregTot)$data,
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$yregTot)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$yregTot)
          } else {
            x$sens[[i]]$reg.output.summary$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregTot,data=getCall(x$sens[[.(i)]]$reg.output$yregTot)$data,
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$yregTot)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$yregTot)
          }
        }
        if (name == "yregDir") {
          cat("### Outcome regression for the direct effect: \n")
          if (inherits(x$sens[[i]]$reg.output$yregDir, "rcreg")|inherits(x$sens[[i]]$reg.output$yregDir, "simexreg")) {
            x$sens[[i]]$reg.output.summary$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregDir, reg = getCall(x$sens[[.(i)]]$reg.output$yregDir)$reg,
                                                                              data=getCall(x$sens[[.(i)]]$reg.output$yregDir)$data,
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$yregDir)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$yregDir)
          } else {
            x$sens[[i]]$reg.output.summary$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$yregDir,data=getCall(x$sens[[.(i)]]$reg.output$yregDir)$data,
                                                                              weights=getCall(x$sens[[.(i)]]$reg.output$yregDir)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$yregDir)
          }
        }
        if (name == "ereg") {
          cat("### Exposure regression for weighting: \n")
          if (inherits(x$sens[[i]]$reg.output$ereg, "rcreg")|inherits(x$sens[[i]]$reg.output$ereg, "simexreg")) {
            x$sens[[i]]$reg.output.summary$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$ereg, reg = getCall(x$sens[[.(i)]]$reg.output$ereg)$reg,
                                                                           data=getCall(x$sens[[.(i)]]$reg.output$ereg)$data,
                                                                           weights=getCall(x$sens[[.(i)]]$reg.output$ereg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$ereg)
          } else {
            x$sens[[i]]$reg.output.summary$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$ereg,data=getCall(x$sens[[.(i)]]$reg.output$ereg)$data,
                                                                           weights=getCall(x$sens[[.(i)]]$reg.output$ereg)$weights, evaluate = FALSE)))
            print(x$sens[[i]]$reg.output.summary$ereg)
          }
        }
        if (name == "mreg") {
          cat("### Mediator regressions: \n")
          for (j in 1:length(x$sens[[i]]$reg.output.summary$mreg)) {
            if (inherits(x$sens[[i]]$reg.output$mreg[[j]], "svyglm")) {
              x$sens[[i]]$reg.output.summary$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]],
                                                                                  design = getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$design,
                                                                                  family = getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$family, evaluate = FALSE)))
              x$sens[[i]]$reg.output.summary$mreg[[j]]$survey.design$call <- eval(bquote(as.call(update(x$sens[[.(i)]]$reg.output.summary$mreg[[.(j)]]$survey.design,
                                                                                                        data = getCall(x$sens[[.(i)]]$reg.output.summary$mreg[[.(j)]]$survey.design)$data,
                                                                                                        weights = getCall(x$sens[[.(i)]]$reg.output.summary$mreg[[.(j)]]$survey.design)$weights, 
                                                                                                        evaluate = FALSE))))
              print(x$sens[[i]]$reg.output.summary$mreg[[j]])
            } else if (inherits(x$sens[[i]]$reg.output$mreg[[j]], "rcreg") |
                       inherits(x$sens[[i]]$reg.output$mreg[[j]], "simexreg")){
              x$sens[[i]]$reg.output.summary$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], 
                                                                                  reg=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$reg,
                                                                                  data=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$data, 
                                                                                  weights=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$weights,
                                                                                  evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$mreg[[j]])
            } else {
              x$sens[[i]]$reg.output.summary$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$mreg[[.(j)]], 
                                                                                  data=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$data, 
                                                                                  weights=getCall(x$sens[[.(i)]]$reg.output$mreg[[.(j)]])$weights,
                                                                                  evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$mreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output.summary$mreg)) cat("\n")
          }
        }
        if (name == "wmdenomreg") {
          cat("### Mediator regressions for weighting (denominator): \n")
          for (j in 1:length(x$sens[[i]]$reg.output.summary$wmdenomreg)) {
            if (inherits(x$sens[[i]]$reg.output$wmdenomreg[[j]], "rcreg") |
                inherits(x$sens[[i]]$reg.output$wmdenomreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output.summary$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]], 
                                                                                        reg=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$reg,
                                                                                        data=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$data, 
                                                                                        weights=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$weights,
                                                                                        evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$wmdenomreg[[j]])
            } else {
              x$sens[[i]]$reg.output.summary$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]], 
                                                                                        data=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$data, 
                                                                                        weights=getCall(x$sens[[.(i)]]$reg.output$wmdenomreg[[.(j)]])$weights,
                                                                                        evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$wmdenomreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output.summary$wmdenomreg)) cat("\n")
          }
        }
        if (name == "wmnomreg") {
          cat("### Mediator regressions for weighting (nominator): \n")
          for (j in 1:length(x$sens[[i]]$reg.output.summary$wmnomreg)) {
            if (inherits(x$sens[[i]]$reg.output$wmnomreg[[j]], "rcreg") |
                inherits(x$sens[[i]]$reg.output$wmnomreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output.summary$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]], 
                                                                                      reg=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$reg,
                                                                                      data=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$data, 
                                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$weights,
                                                                                      evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$wmnomreg[[j]])
            } else {
              x$sens[[i]]$reg.output.summary$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]], 
                                                                                      data=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$data, 
                                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$wmnomreg[[.(j)]])$weights,
                                                                                      evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$wmnomreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output.summary$wmnomreg)) cat("\n")
          }
        }
        if (name == "postcreg") {
          cat("### Regressions for mediator-outcome confounders affected by the exposure: \n")
          for (j in 1:length(x$sens[[i]]$reg.output.summary$postcreg)) {
            if (inherits(x$sens[[i]]$reg.output$postcreg[[j]], "rcreg") |
                inherits(x$sens[[i]]$reg.output$postcreg[[j]], "simexreg")) {
              x$sens[[i]]$reg.output.summary$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]], 
                                                                                      reg=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$reg,
                                                                                      data=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$data, 
                                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$weights,
                                                                                      evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$postcreg[[j]])
            } else {
              x$sens[[i]]$reg.output.summary$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]], 
                                                                                      data=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$data, 
                                                                                      weights=getCall(x$sens[[.(i)]]$reg.output$postcreg[[.(j)]])$weights,
                                                                                      evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary$postcreg[[j]])
            }
            if (j < length(x$sens[[i]]$reg.output.summary$postcreg)) cat("\n")
          }
        }
        cat("\n")
      }
    } else {
      for (m in 1:length(x$sens[[i]]$reg.output)) { 
        cat(paste("### Regressions with imputed dataset", m, "\n\n"))
        regnames <- names(x$sens[[i]]$reg.output[[m]])
        for (name in regnames) {
          if (name == "yreg") {
            cat("#### Outcome regression: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yreg, "svyglm")) {
              x$sens[[i]]$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg,
                                                                                  design = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$design,
                                                                                  family = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$family, evaluate = FALSE)))
              x$sens[[i]]$reg.output.summary[[m]]$yreg$survey.design$call <- eval(bquote(as.call(update(x$sens[[.(i)]]$reg.output.summary[[.(m)]]$yreg$survey.design,
                                                                                                        data = getCall(x$sens[[.(i)]]$reg.output.summary[[.(m)]]$yreg$survey.design)$data,
                                                                                                        weights = getCall(x$sens[[.(i)]]$reg.output.summary[[.(m)]]$yreg$survey.design)$weights, 
                                                                                                        evaluate = FALSE))))
              print(x$sens[[i]]$reg.output.summary[[m]]$yreg)
            } else if (inherits(x$sens[[i]]$reg.output[[m]]$yreg, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yreg, "simexreg")){
              x$sens[[i]]$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg,
                                                                                  reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$reg,
                                                                                  data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$data,
                                                                                  weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$weights, 
                                                                                  evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$yreg)
            } else {
              x$sens[[i]]$reg.output.summary[[m]]$yreg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg,
                                                                                  data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$data,
                                                                                  weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yreg)$weights, 
                                                                                  evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$yreg)
            }
          }
          if (name == "yregTot") {
            cat("#### Outcome regression for the total effect: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yregTot, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yregTot, "simexreg")){
              x$sens[[i]]$reg.output.summary[[m]]$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot,
                                                                                     reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$reg,
                                                                                     data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$data,
                                                                                     weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$weights, 
                                                                                     evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$yregTot)
            } else {
              x$sens[[i]]$reg.output.summary[[m]]$yregTot$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot,
                                                                                     data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$data,
                                                                                     weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregTot)$weights, 
                                                                                     evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$yregTot)
            }
          }
          if (name == "yregDir") {
            cat("#### Outcome regression for the direct effect: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$yregDir, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$yregDir, "simexreg")){
              x$sens[[i]]$reg.output.summary[[m]]$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir,
                                                                                     reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$reg,
                                                                                     data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$data,
                                                                                     weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$weights, 
                                                                                     evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$yregDir)
            } else {
              x$sens[[i]]$reg.output.summary[[m]]$yregDir$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir,
                                                                                     data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$data,
                                                                                     weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$yregDir)$weights, 
                                                                                     evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$yregDir)
            }
          }
          if (name == "ereg") {
            cat("#### Exposure regression for weighting: \n")
            if (inherits(x$sens[[i]]$reg.output[[m]]$ereg, "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$ereg, "simexreg")){
              x$sens[[i]]$reg.output.summary[[m]]$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg,
                                                                                  reg = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$reg,
                                                                                  data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$data,
                                                                                  weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$weights, 
                                                                                  evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$ereg)
            } else {
              x$sens[[i]]$reg.output.summary[[m]]$ereg$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg,
                                                                                  data = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$data,
                                                                                  weights = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$ereg)$weights, 
                                                                                  evaluate = FALSE)))
              print(x$sens[[i]]$reg.output.summary[[m]]$ereg)
            }
          }
          if (name == "mreg") {
            cat("#### Mediator regressions: \n")
            for (j in 1:length(x$sens[[i]]$reg.output.summary[[m]]$mreg)) {
              if (inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "svyglm")) {
                x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]],
                                                                                         design = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$design,
                                                                                         family = getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$family, evaluate = FALSE)))
                x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]]$survey.design$call <- eval(bquote(as.call(update(x$sens[[.(i)]]$reg.output.summary[[.(m)]]$mreg[[.(j)]]$survey.design,
                                                                                                               data = getCall(x$sens[[.(i)]]$reg.output.summary[[.(m)]]$mreg[[.(j)]]$survey.design)$data,
                                                                                                               weights = getCall(x$sens[[.(i)]]$reg.output.summary[[.(m)]]$mreg[[.(j)]]$survey.design)$weights, 
                                                                                                               evaluate = FALSE))))
                print(x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]])
              } else if (inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$mreg[[j]], "simexreg")){
                x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], 
                                                                                         reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$reg,
                                                                                         data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$data, 
                                                                                         weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$weights, 
                                                                                         evaluate = FALSE)))
                print(x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]])
              } else {
                x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]], 
                                                                                         data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$data, 
                                                                                         weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$mreg[[.(j)]])$weights, 
                                                                                         evaluate = FALSE)))
                print(x$sens[[i]]$reg.output.summary[[m]]$mreg[[j]])
              }
              if (j < length(x$sens[[i]]$reg.output.summary[[m]]$mreg)) cat("\n")
            }
          }
          if (name == "wmdenomreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$wmdenomreg)) {
              cat("#### Mediator regressions for weighting (denominator): \n")
              for (j in 1:length(x$sens[[i]]$reg.output.summary[[m]]$wmdenomreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$wmdenomreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output.summary[[m]]$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]], 
                                                                                                 reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$reg,
                                                                                                 data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$data, 
                                                                                                 weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$weights, 
                                                                                                 evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output.summary[[m]]$wmdenomreg[[j]])
                } else {
                  x$sens[[i]]$reg.output.summary[[m]]$wmdenomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]], 
                                                                                                 data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$data, 
                                                                                                 weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmdenomreg[[.(j)]])$weights, 
                                                                                                 evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output.summary[[m]]$wmdenomreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output.summary[[m]]$wmdenomreg)) cat("\n")
              }
            }
          }
          if (name == "wmnomreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$wmnomreg)) {
              cat("#### Mediator regressions for weighting (nominator): \n")
              for (j in 1:length(x$sens[[i]]$reg.output.summary[[m]]$wmnomreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$wmnomreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output.summary[[m]]$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]], 
                                                                                               reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$reg,
                                                                                               data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$data, 
                                                                                               weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$weights, 
                                                                                               evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output.summary[[m]]$wmnomreg[[j]])
                } else {
                  x$sens[[i]]$reg.output.summary[[m]]$wmnomreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]], 
                                                                                               data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$data, 
                                                                                               weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$wmnomreg[[.(j)]])$weights, 
                                                                                               evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output.summary[[m]]$wmnomreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output.summary[[m]]$wmnomreg)) cat("\n")
              }
            }
          }
          if (name == "postcreg") {
            if (!is.null(x$sens[[i]]$reg.output[[m]]$postcreg)) {
              cat("#### Regressions for mediator-outcome confounders affected by the exposure: \n")
              for (j in 1:length(x$sens[[i]]$reg.output.summary[[m]]$postcreg)) {
                if (inherits(x$sens[[i]]$reg.output[[m]]$postcreg[[j]], "rcreg")|inherits(x$sens[[i]]$reg.output[[m]]$postcreg[[j]], "simexreg")){
                  x$sens[[i]]$reg.output.summary[[m]]$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]], 
                                                                                               reg=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$reg,
                                                                                               data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$data, 
                                                                                               weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$weights, 
                                                                                               evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output.summary[[m]]$postcreg[[j]])
                } else {
                  x$sens[[i]]$reg.output.summary[[m]]$postcreg[[j]]$call <- eval(bquote(update(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]], 
                                                                                               data=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$data, 
                                                                                               weights=getCall(x$sens[[.(i)]]$reg.output[[.(m)]]$postcreg[[.(j)]])$weights, 
                                                                                               evaluate = FALSE)))
                  print(x$sens[[i]]$reg.output.summary[[m]]$postcreg[[j]])
                }
                if (j < length(x$sens[[i]]$reg.output.summary[[m]]$postcreg)) cat("\n")
              }
            }
          }
          cat("\n")
        }
      }
    }
    
    # scale and legend
    full <- x$naive$methods$full
    model <- x$naive$methods$model
    EMint <- x$naive$variables$EMint
    if (model == "iorw") {
      if (x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output[[1]]$yregTot
      if (!x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output$yregTot
    } else {
      if (x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output[[1]]$yreg
      if (!x$naive$multimp$multimp) yreg_mid <- x$sens[[1]]$reg.output$yreg
    }
    if (inherits(yreg_mid, "rcreg") | inherits(yreg_mid, "simexreg")) yreg_mid <- yreg_mid$NAIVEreg
    is_lm <- inherits(yreg_mid, "lm")
    is_glm <- inherits(yreg_mid, "glm")
    is_svyglm <- inherits(yreg_mid, "svyglm")
    is_gam <- inherits(yreg_mid, "gam")
    if (is_lm | is_glm) family_yreg <- family(yreg_mid)
    is_multinom <- inherits(yreg_mid, "multinom")
    is_svymultinom <- inherits(yreg_mid, "svymultinom")
    is_polr <- inherits(yreg_mid, "polr")
    is_survreg <- inherits(yreg_mid, "survreg")
    is_coxph <- inherits(yreg_mid, "coxph")
    if ((is_lm | is_glm) && (family_yreg$family %in% c("gaussian", "inverse.gaussian", "quasi", "Gamma"))) {
      scale <- "mean difference scale"
      if (model == "iorw") {
        if (full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect; pm: proportion mediated)"   
        if (!full) legend <- "(te: total effect; pnde: pure natural direct effect; tnie: total natural indirect effect)"
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rintref: randomized analogue of reference interaction; rintmed: randomized analogue of mediated interaction; cde(prop): proportion cde; rintref(prop): proportion rintref; rintmed(prop): proportion rintmed; rpnie(prop): proportion rpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(cde: controlled direct effect; rpnde: randomized analogue of pure natural direct effect; rtnde: randomized analogue of total natural direct effect; rpnie: randomized analogue of pure natural indirect effect; rtnie: randomized analogue of total natural indirect effect; te: total effect)"
      } else {
        if (full) {
          if (EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; intref: reference interaction; intmed: mediated interaction; cde(prop): proportion cde; intref(prop): proportion intref; intmed(prop): proportion intmed; pnie(prop): proportion pnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect; pm: overall proportion mediated)"
        } else legend <- "(cde: controlled direct effect; pnde: pure natural direct effect; tnde: total natural direct effect; pnie: pure natural indirect effect; tnie: total natural indirect effect; te: total effect)"
      }
    } else if ((is_lm | is_glm) && (family_yreg$family %in% c("poisson", "quasipoisson", "ziplss") |
                                    startsWith(family_yreg$family, "Negative Binomial") |
                                    startsWith(family_yreg$family, "Zero inflated Poisson"))) {
      scale <- "rate ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnie: total natural indirect effect rate ratio)" 
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; rERintref: randomized analogue of excess relative rate due to reference interaction; rERintmed: randomized analogue of excess relative rate due to mediated interaction; rERpnie: randomized analogue of excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect rate ratio; rRpnde: randomized analogue of pure natural direct effect rate ratio; rRtnde: randomized analogue of total natural direct effect rate ratio; rRpnie: randomized analogue of pure natural indirect effect rate ratio; rRtnie: randomized analogue of total natural indirect effect rate ratio; Rte: total effect rate ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; ERcde: excess relative rate due to controlled direct effect; ERintref: excess relative rate due to reference interaction; ERintmed: excess relative rate due to mediated interaction; ERpnie: excess relative rate due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect rate ratio; Rpnde: pure natural direct effect rate ratio; Rtnde: total natural direct effect rate ratio; Rpnie: pure natural indirect effect rate ratio; Rtnie: total natural indirect effect rate ratio; Rte: total effect rate ratio)"
      }
    } else if (((is_lm | is_glm) && (family_yreg$family %in% c("binomial", "quasibinomial", "multinom") |
                                     startsWith(family_yreg$family, "Ordered Categorical"))) |
               is_multinom | is_polr) {
      scale <- "risk ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnie: total natural indirect effect risk ratio)" 
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; rERintref: randomized analogue of excess relative risk due to reference interaction; rERintmed: randomized analogue of excess relative risk due to mediated interaction; rERpnie: randomized analogue of excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; rRpnde: randomized analogue of pure natural direct effect risk ratio; rRtnde: randomized analogue of total natural direct effect risk ratio; rRpnie: randomized analogue of pure natural indirect effect risk ratio; rRtnie: randomized analogue of total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; ERcde: excess relative risk due to controlled direct effect; ERintref: excess relative risk due to reference interaction; ERintmed: excess relative risk due to mediated interaction; ERpnie: excess relative risk due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect risk ratio; Rpnde: pure natural direct effect risk ratio; Rtnde: total natural direct effect risk ratio; Rpnie: pure natural indirect effect risk ratio; Rtnie: total natural indirect effect risk ratio; Rte: total effect risk ratio)"
      }
    } else if (is_coxph) {
      scale <- "hazard ratio scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnie: total natural indirect effect hazard ratio)" 
      }  else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; rERintref: randomized analogue of excess relative hazard due to reference interaction; rERintmed: randomized analogue of excess relative hazard due to mediated interaction; rERpnie: randomized analogue of excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect hazard ratio; rRpnde: randomized analogue of pure natural direct effect hazard ratio; rRtnde: randomized analogue of total natural direct effect hazard ratio; rRpnie: randomized analogue of pure natural indirect effect hazard ratio; rRtnie: randomized analogue of total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; ERcde: excess relative hazard due to controlled direct effect; ERintref: excess relative hazard due to reference interaction; ERintmed: excess relative hazard due to mediated interaction; ERpnie: excess relative hazard due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect hazard ratio; Rpnde: pure natural direct effect hazard ratio; Rtnde: total natural direct effect hazard ratio; Rpnie: pure natural indirect effect hazard ratio; Rtnie: total natural indirect effect hazard ratio; Rte: total effect hazard ratio)"
      }
    } else if (is_survreg) {
      scale <- "mean survival scale"
      if (model == "iorw") {
        if (full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; pm: proportion mediated)"   
        if (!full) legend <- "(Rte: total effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio)" 
      } else if (length(x$naive$variables$postc) != 0) {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; rERintref: randomized analogue of excess mean survival ratio due to reference interaction; rERintmed: randomized analogue of excess mean survival ratio due to mediated interaction; rERpnie: randomized analogue of excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; rERintref(prop): proportion rERintref; rERintmed(prop): proportion rERintmed; rERpnie(prop): proportion rERpnie; rpm: randomized analogue of overall proportion mediated; rint: randomized analogue of overall proportion attributable to interaction; rpe: randomized analogue of overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; rpm: randomized analogue of overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect mean survival ratio; rRpnde: randomized analogue of pure natural direct effect mean survival ratio; rRtnde: randomized analogue of total natural direct effect mean survival ratio; rRpnie: randomized analogue of pure natural indirect effect mean survival ratio; rRtnie: randomized analogue of total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
      } else {
        if (full) {
          if (EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; ERcde: excess mean survival ratio due to controlled direct effect; ERintref: excess mean survival ratio due to reference interaction; ERintmed: excess mean survival ratio due to mediated interaction; ERpnie: excess mean survival ratio due to pure natural indirect effect; ERcde(prop): proportion ERcde; ERintref(prop): proportion ERintref; ERintmed(prop): proportion ERintmed; ERpnie(prop): proportion ERpnie; pm: overall proportion mediated; int: overall proportion attributable to interaction; pe: overall proportion eliminated)"
          if (!EMint) legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio; pm: overall proportion mediated)"
        } else legend <- "(Rcde: controlled direct effect mean survival ratio; Rpnde: pure natural direct effect mean survival ratio; Rtnde: total natural direct effect mean survival ratio; Rpnie: pure natural indirect effect mean survival ratio; Rtnie: total natural indirect effect mean survival ratio; Rte: total effect mean survival ratio)"
      }
    }
    
    # print error-corrected results
    cat(paste0(paste("\n## Error-corrected causal effects on the", scale, "for measurement error "), i, ": \n"))
    printCoefmat(x$summarydf[[i]], digits = digits, has.Pvalue = TRUE)
    cat("----------------------------------------------------------------\n")
  }
  cat("\n")
  cat(legend)
  cat("\n\nRelevant variable values: \n")
  print(x$naive$ref)
}


#' Plotting Results of Sensitivity Analysis for Measurement Error
#' 
#' \code{ggcmsens} is used to plot results of \code{cmsens} nicely with plotting functions
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
      } else {
        if (x$naive$methods$full) ref <- c(0, 1)
        if (!x$naive$methods$full) ref <- 1
      }
    } else {
      if ((inherits(x$naive$reg.output[[1]]$yreg, "lm") | inherits(x$naive$reg.output[[1]]$yreg, "glm")) &&
          (family(x$naive$reg.output[[1]]$yreg)$family %in% c("gaussian","Gamma","inverse.gaussian","quasi"))) {
        ref <- 0
      } else {
        if (x$naive$methods$full) ref <- c(0, 1)
        if (!x$naive$methods$full) ref <- 1
      }
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
                                    Point = pe.mid, CIlower = ci.low.mid, CIupper = ci.high.mid,
                                    ReliabilityRatio = rep(factor(round(x$ME$reliabilityRatio[i], 2)), 
                                                           length(pe.mid))))
    }
    ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper, colour = ReliabilityRatio), 
                    data = effect_df, width = errorbar.width, size = errorbar.size, position = errorbar.position) +
      geom_point(aes(x = Effect, y = Point, colour = ReliabilityRatio),
                 data = effect_df, size = point.size, position = point.position) +
      ylab("Point Estimate and 95% CI") +
      scale_colour_hue() +
      geom_hline(yintercept = ref, color = refline.colour, size = refline.size) +
      theme(legend.position = "bottom")
  } else if (x$ME$MEvartype == "categorical") {
    pe.mid <- x$sens[[1]]$effect.pe
    ci.low.mid <- x$sens[[1]]$effect.ci.low
    ci.high.mid <- x$sens[[1]]$effect.ci.high
    effect_df <- data.frame(Effect = rep(factor(names(pe.mid), levels = names(pe.mid)), 2),
                            Point = c(pe.mid, naive.pe), CIlower = c(ci.low.mid, naive.ci.low),
                            CIupper = c(ci.high.mid, naive.ci.high),
                            MC = factor(c(rep("SIMEX", length(pe.mid)), rep("Naive", length(naive.pe)))),
                            MisclassificationMAtrix = factor(rep("MEerror[1]", length(pe.mid)+length(naive.pe))))
    if (length(x$sens) > 1) {
      for (i in 2:length(x$sens)) {
        pe.mid <- x$sens[[i]]$effect.pe
        ci.low.mid <- x$sens[[i]]$effect.ci.low
        ci.high.mid <- x$sens[[i]]$effect.ci.high
        effect_df <- rbind(effect_df,
                           data.frame(Effect = rep(factor(names(pe.mid), levels = names(pe.mid)), 2),
                                      Point = c(pe.mid, naive.pe), CIlower = c(ci.low.mid, naive.ci.low),
                                      CIupper = c(ci.high.mid, naive.ci.high),
                                      MC = factor(c(rep("SIMEX", length(pe.mid)), rep("Naive", length(naive.pe)))),
                                      MisclassificationMAtrix = factor(rep(paste0("MEerror[", i, "]"), 
                                                                           length(pe.mid)+length(naive.pe)))))
        
      }
    }
    ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper, colour = MC), data = effect_df, 
                    width = errorbar.width, size = errorbar.size, position = errorbar.position) +
      geom_point(aes(x = Effect, y = Point, colour = MC), data = effect_df, size = point.size, position = point.position) +
      facet_grid(MisclassificationMAtrix~.) +
      ylab("Point Estimate and 95% CI") +
      scale_colour_hue() +
      geom_hline(yintercept = ref, color = refline.colour, size = refline.size) +
      theme(legend.position = "bottom", legend.title = element_blank())
  }
}

