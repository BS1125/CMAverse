#' Simulation and Extrapolation for Measurement Error Correction
#'
#' \code{simexreg} is used to correct a regression object with a variable measured with 
#' error via \emph{SIMEX} by Cook et al. (1994) and Küchenhoff et al. (2006).
#'
#' @param reg naive regression object. See \code{Details}.
#' @param formula regression formula
#' @param data new dataset for \code{reg}
#' @param weights new weights for \code{reg}
#' @param MEvariable variable measured with error
#' @param MEvartype type of the variable measured with error. Can be \code{continuous} or 
#' \code{categorical} (first 3 letters are enough).
#' @param MEerror the standard deviation of the measurement error (when \code{MEvartype}
#' is \code{continuous}) or the misclassification matrix (when \code{MEvartype}
#' is \code{categorical}). 
#' @param variance a logical value. If \code{TRUE}, estimate the var-cov matrix of
#' coefficients through Jackknife. Default is \code{FALSE}.
#' @param lambda a vector of lambdas for \emph{SIMEX}. Default is \code{c(0.5, 1, 1.5, 2)}. 
#' @param B number of simulations for \emph{SIMEX}. Default is \code{200}.
#' @param x an object of class \code{simexreg}
#' @param object an object of class \code{simexreg}
#' @param digits minimal number of significant digits. See \link{print.default}.
#' @param evaluate a logical value. If \code{TRUE}, the updated call is evaluated. Default
#' is \code{TRUE}.
#' @param ... additional arguments
#' 
#' @details
#' 
#' \code{reg} fitted by \link{lm}, \link{glm} (with family \code{gaussian}, \code{binomial} or
#' \code{poisson}), \link[nnet]{multinom}, \link[MASS]{polr}, \link[survival]{coxph} or
#' \link[survival]{survreg} is supported.
#'  
#' @return
#' If \code{MEvariable} is not in the regression formula, \code{reg} is returned. If 
#' \code{MEvariable} is in the regression formula, an object of class \code{simexreg} is returned:
#' \item{call}{the function call,}
#' \item{NAIVEreg}{the naive regression object,}
#' \item{ME}{a list of \code{MEvariable}, \code{MEvartype}, \code{MEerror}, \code{variance},
#' \code{lambda} and \code{B},}
#' \item{RCcoef}{coefficient estimates corrected by SIMEX,}
#' \item{RCsigma}{the residual standard deviation of a linear regression object corrected by 
#' SIMEX,}
#' \item{RCvcov}{the var-cov matrix of coefficients corrected by SIMEX,}
#' ...
#'
#' @seealso \code{\link{rcreg}}, \code{\link{cmsens}}, \code{\link{cmest}}.
#'
#' @references
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
#' @examples
#' 
#' \dontrun{
#' rm(list=ls())
#' library(CMAverse)
#' 
#' # lm
#' n <- 1000
#' x1 <- rnorm(n, mean = 5, sd = 3)
#' x2_true <- rnorm(n, mean = 2, sd = 1)
#' error1 <- rnorm(n, mean = 0, sd = 0.5)
#' x2_error <- x2_true + error1
#' x3 <- rbinom(n, size = 1, prob = 0.4)
#' y <- 1 + 2 * x1 + 4 * x2_true + 2 * x3  + rnorm(n, mean = 0, sd = 2)
#' data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
#'                    x3 = x3, y = y)
#' reg_naive <- lm(y ~ x1 + x2_error + x3, data = data)
#' reg_true <- lm(y ~ x1 + x2_true + x3, data = data)
#' reg_simex <- simexreg(reg = reg_naive, data = data, 
#' MEvariable = "x2_error", MEvartype = "con", MEerror = 0.5, variance = TRUE)
#' coef(reg_simex)
#' vcov(reg_simex)
#' sigma(reg_simex)
#' formula(reg_simex)
#' family(reg_simex)
#' predict(reg_simex, newdata = data[1, ])
#' reg_simex_model <- model.frame(reg_simex)
#' reg_simex_update <- update(reg_simex, data = data, weights = rep(1, n))
#' reg_simex_summ <- summary(reg_simex)
#'                 
#' # glm
#' n <- 1000
#' x1 <- rnorm(n, mean = 5, sd = 3)
#' x2_true <- sample(x = c(1:3), size = n, prob = c(0.2,0.3,0.5), replace = TRUE)
#' MEerror <- matrix(c(0.8,0.1,0.1,0.2,0.7,0.1,0.05,0.25,0.7), nrow = 3)
#' x2_error <- x2_true
#' for (j in 1:3) {
#'   x2_error[which(x2_error == c(1:3)[j])] <-
#'     sample(x = c(1:3), size = length(which(x2_error == c(1:3)[j])),
#'            prob = MEerror[, j], replace = TRUE)
#' }
#' x2_true <- as.factor(x2_true)
#' x2_error <- as.factor(x2_error)
#' x3 <- rnorm(n, mean = 2, sd = 1)
#' linearpred <- 1 + 0.3 * x1 - 1.5*(x2_true == 2) - 2.5*(x2_true == 3) - 0.2 * x3
#' py <- exp(linearpred) / (1 + exp(linearpred))
#' y <- rbinom(n, size = 1, prob = py)
#' data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
#'                    x3 = x3, y = y)
#' reg_naive <- glm(y ~ x1 + x2_error + x3, data = data, family = binomial("logit"))
#' reg_true <- glm(y ~ x1 + x2_true + x3, data = data, family = binomial("logit"))
#' reg_simex <- simexreg(reg = reg_naive, data = data, 
#' MEvariable = "x2_error", MEerror = MEerror, variance = TRUE, MEvartype = "cat")
#' }
#'                   
#' @importFrom stats as.formula model.frame family coef predict model.matrix getCall cov 
#' formula vcov pt
#' @importFrom boot boot
#' 
#' @export
#' 
simexreg <- function (reg = NULL, formula = NULL, data = NULL, weights = NULL, 
                      MEvariable = NULL, MEvartype = NULL, MEerror = NULL,
                      variance = FALSE, lambda = c(0.5, 1, 1.5, 2), B = 200) {
  
  cl <- match.call()
  
  # assign svyglm the global environment
  assign2glob <- function(key, val, pos) assign(key, val, envir = as.environment(pos))
  assign2glob("svyglm", survey::svyglm, 1L)
  
  # the vector of variable names in the regression formula
  if (is.null(formula)) formula <- as.formula(formula)
  if (!is.null(formula)) formula <- formula(reg)
  var_vec <- unique(all.vars(formula))
  
  if (length(MEvariable) == 0 | !MEvariable %in% var_vec) {
    out <- reg
  } else {
    if (is.null(data)) stop("Unspecified data")
    if (!all(lambda>0)) stop("lambda should be positive")
    if (length(MEvariable) > 1) stop("Currently only supports one variable measured with error")
    if (length(MEvariable) != length(MEvartype)) stop("length(MEvariable) != length(MEvartype)")
    if (MEvartype == "con") MEvartype <- "continuous"
    if (MEvartype == "cat") MEvartype <- "categorical"
    if (MEvartype == "continuous") {
      if (length(MEvariable) != length(MEerror)) stop("length(MEvariable) != length(MEerror)")
      if (length(MEerror) != 0 && MEerror < 0) stop("MEerror should be >= 0")
    } else if (MEvartype == "categorical") {
      if (length(MEvariable) != 0 && !is.matrix(MEerror)) stop("MEerror should be a matrix")
      if (dim(MEerror)[1] != dim(MEerror)[2] |
          dim(MEerror)[1] != length(unique(data[, MEvariable]))) stop("Incorrect dimension of MEerror")
      MEvar_lev <- levels(as.factor(data[, MEvariable]))
      dimnames(MEerror) <- list(MEvar_lev, MEvar_lev)
      if (!check.mc.matrix(list(MEerror))) stop("MEerror may contain negative values for exponents smaller than 1")
    } else stop("Unsupported MEvartype; use 'continuous' or 'categorical'")
    
    regCall <- getCall(reg)
    regClass <- class(reg)
    if (!(identical(regClass, c("svyglm", "glm", "lm")) | identical(regClass, c("svymultinom")) |
          identical(regClass, "lm") | 
          (identical(regClass, c("glm", "lm")) && 
           family(reg)$family %in% c("gaussian", "binomial", "poisson")) | 
          identical(regClass, c("multinom", "nnet")) | identical(regClass, "polr") | 
          identical(regClass, "coxph") | identical(regClass, "survreg"))) stop(
            "SIMEX applied to unsupported regression object")
    
    # update reg with data and weights
    if (identical(regClass, c("svyglm", "glm", "lm"))) {
      designCall <- getCall(reg$survey.design)
      designCall$data <- data
      designCall$weights <- weights
      regCall$design <- eval.parent(designCall)
      regCall$formula <- formula
    } else {
      regCall$formula <- formula
      regCall$data <- data
      regCall$weights <- weights
      if (inherits(reg, "multinom")) regCall$trace <- FALSE
      if (inherits(reg, "polr")) regCall$Hess <- TRUE
      if (inherits(reg, "coxph")) regCall$model <- TRUE
    }
    reg <- eval.parent(regCall)
    
    n <- nrow(data)
    # output list
    out <- list(call = cl, NAIVEreg = reg,
                ME = list(MEvariable = MEvariable, MEvartype = MEvartype,
                          MEerror = MEerror, variance = variance, lambda = lambda, B = B))
    
    if (!((((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
            identical(regClass, "lm")) && MEvartype == "categorical") | 
          inherits(reg, "svyglm") | inherits(reg, "svymultinom") | 
          inherits(reg, "multinom") | inherits(reg, "survreg") | inherits(reg, "coxph"))) {
      if (MEvartype == "continuous") {
        SIMEX <- simex::simex(model = reg, SIMEXvariable = MEvariable, measurement.error = MEerror, 
                              lambda = lambda, B = B, fitting.method = "quadratic",
                              jackknife.estimation = "quadratic", asymptotic = FALSE)
        if (!inherits(reg, "polr")) out$SIMEXcoef <- SIMEX$coefficients
        if (inherits(reg, "polr")) out$SIMEXcoef <- c(SIMEX$coefficients, SIMEX$zeta)
        out$SIMEXvcov <- SIMEX$variance.jackknife
        if ((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
            identical(regClass, "lm")) {
          reg_fit <- reg
          reg_fit$coefficients <- out$SIMEXcoef
          if (MEvariable == all.vars(formula(reg))[1]) SIMEXsigma <- sqrt(sigma(reg)^2 - MEerror^2)
          if (MEvariable != all.vars(formula(reg))[1]) SIMEXsigma <- 
            sqrt(sum((model.frame(reg_fit)[, 1] - predict(reg_fit, newdata = data)) ^ 2) / 
                   (n - length(out$SIMEXcoef)) - (out$SIMEXcoef[MEvariable] * MEerror) ^ 2)
          out$SIMEXsigma <- unname(SIMEXsigma)
        }
      } else if (MEvartype == "categorical") {
        if (!is.factor(data[, MEvariable])) {
        data[, MEvariable] <- as.factor(data[, MEvariable])
        regCall$data <- data
        reg <- eval.parent(regCall)
        }
        SIMEX <- simex::mcsimex(model = reg, SIMEXvariable = MEvariable, mc.matrix = MEerror, 
                                lambda = lambda, B = B, fitting.method = "quadratic",
                                jackknife.estimation = "quadratic", asymptotic = FALSE)
        if (!inherits(reg, "polr")) out$SIMEXcoef <- SIMEX$coefficients
        if (inherits(reg, "polr")) out$SIMEXcoef <- c(SIMEX$coefficients, SIMEX$zeta)
        out$SIMEXvcov <- SIMEX$variance.jackknife
      }
    } else {  
      lambda <- c(0, lambda)
      
      # coefficient estimates
      if (identical(class(reg), "polr")) {
        SIMcoef <- c(coef(reg), reg$zeta)
      } else if (identical(class(reg), c("multinom", "nnet"))) {
        SIMcoef <- coef(reg)
        if (length(reg$lev) > 2) {
          coefnames <- dimnames(SIMcoef)
          SIMcoef <- as.vector(t(SIMcoef))
          coefnames <- as.vector(outer(coefnames[[2]], coefnames[[1]], function(name2, name1)
            paste(name1, name2, sep = ":")))
          names(SIMcoef) <- coefnames
        }
      } else if (inherits(reg, "svymultinom")) {
        SIMcoef <- coef(reg)
        if (length(reg$NAIVEreg$lev) > 2) {
          coefnames <- dimnames(SIMcoef)
          SIMcoef <- as.vector(t(SIMcoef))
          coefnames <- as.vector(outer(coefnames[[2]], coefnames[[1]], function(name2, name1)
            paste(name1, name2, sep = ":")))
          names(SIMcoef) <- coefnames
        }
      } else SIMcoef <- coef(reg)
      ncoef <- length(SIMcoef)
      
      # var-cov matrix of coefficient estimates
      if (variance) SIMvcov <- list(vcov(reg)[1:ncoef, 1:ncoef])
      
      # residual standard deviation for a simple linear regression
      if (((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
           identical(regClass, "lm")) && MEvartype == "categorical") SIMsigma <- sigma(reg)
      
      for (i in 2:length(lambda)) {
        SIMcoef_mid <- matrix(nrow = B, ncol = ncoef)
        if (variance) SIMvcov_mid <- matrix(rep(0, ncol(SIMvcov[[1]]) ^ 2), nrow = ncol(SIMvcov[[1]]))
        SIMsigma_mid <- c()
        for (b in 1:B) {
          SIMdata <- data
          if (MEvartype == "continuous") {
            # randomly generate true MEvariable
            Z <- rnorm(n = n)
            SIMdata[, MEvariable] <- SIMdata[, MEvariable] + (sqrt(lambda[i]) * MEerror  *  Z )
          } else if (MEvartype == "categorical") {
            MEvar_lev <- levels(as.factor(SIMdata[, MEvariable]))
            category <- unique(data[, MEvariable])[match(MEvar_lev, unique(data[, MEvariable]))]
            mcmdecomp <- eigen(MEerror)
            SIMmcm <- mcmdecomp$vectors %*% diag(mcmdecomp$values)^lambda[i] %*% solve(mcmdecomp$vectors)
            for (j in 1:length(category)) {
              SIMdata[, MEvariable][which(SIMdata[, MEvariable] == MEvar_lev[j])] <-
                sample(x = category, size = length(which(SIMdata[, MEvariable] == category[j])),
                       prob = SIMmcm[, j], replace = TRUE)
            }
          }
          if (identical(regClass, c("svyglm", "glm", "lm"))) {
            designCall <- getCall(reg$survey.design)
            designCall$data <- SIMdata
            regCall$design <- eval.parent(designCall)
          } else regCall$data <- SIMdata
          SIMreg <- eval.parent(regCall)
          
          if (identical(class(SIMreg), "polr")) {
            SIMcoef_mid[b, ] <- c(coef(SIMreg), SIMreg$zeta)
          } else if (identical(class(SIMreg), c("multinom", "nnet")) | inherits(SIMreg, "svymultinom")) {
            SIMcoef_mid[b, ]  <- as.vector(t(coef(SIMreg)))
          } else {SIMcoef_mid[b, ]  <- c(coef(SIMreg))}
          if (variance) SIMvcov_mid <- SIMvcov_mid + vcov(SIMreg)[1:ncoef, 1:ncoef]/B
          if (((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
               identical(regClass, "lm")) && MEvartype == "categorical") SIMsigma_mid <- 
            c(SIMsigma_mid, sigma(SIMreg))
        }
        SIMcoef <- rbind(SIMcoef, colMeans(SIMcoef_mid))
        if (variance) SIMvcov[[i]]<- SIMvcov_mid[1:ncoef, 1:ncoef] - cov(SIMcoef_mid)
        if (((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
             identical(regClass, "lm")) && MEvartype == "categorical") SIMsigma <- 
          c(SIMsigma, mean(SIMsigma_mid))
      }
      # extrapolation
      extrap_coef <- lm(SIMcoef ~ lambda + I(lambda^2))
      SIMEXcoef <- as.vector(predict(extrap_coef, newdata = data.frame(lambda = -1)))
      names(SIMEXcoef) <- colnames(SIMcoef)
      out$SIMEXcoef <- SIMEXcoef
      if (variance) {
        SIMvcov.element <- matrix(unlist(SIMvcov), nrow = length(lambda) , byrow = TRUE)
        extrap_vcov <- lm(SIMvcov.element ~ lambda + I(lambda^2))
        SIMEXvcov <- matrix(predict(extrap_vcov, newdata = data.frame(lambda = -1)), nrow = length(SIMEXcoef))
        dimnames(SIMEXvcov) <- list(names(SIMEXcoef), names(SIMEXcoef))
        out$SIMEXvcov <- SIMEXvcov
      }
      
      if (((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
           identical(regClass, "lm")) && MEvartype == "categorical") {
        extrap_sigma <- lm(SIMsigma ~ lambda + I(lambda^2))
        SIMEXsigma <- as.vector(predict(extrap_sigma, newdata = data.frame(lambda = -1)))
        out$SIMEXsigma <- SIMEXsigma
      }
    }  
    
    # remove svyglm from the global environment
    rm(svyglm, envir = .GlobalEnv)
    
    class(out) <- "simexreg"
  }
  return(out)
}


#' @describeIn simexreg Extract coefficients corrected by \code{simexreg}
#' @export
coef.simexreg <- function(object, ...) {
  return(object$SIMEXcoef)
}


#' @describeIn simexreg Extract the var-cov matrix of coefficients corrected by 
#' \code{simexreg}
#' @export
vcov.simexreg <- function(object, ...) {
  return(object$SIMEXvcov)
}


#' @describeIn simexreg Extract the residual standard deviation of a linear regression object 
#' corrected by \code{simexreg}
#' @export
sigma.simexreg <- function(object, ...) {
  return(object$SIMEXsigma)
}


#' @describeIn simexreg Extract the regression formula
#' @export
formula.simexreg <- function(x, ...) {
  return(formula(x$NAIVEreg))
}


#' @describeIn simexreg Extract the family of a regression of class \code{lm} or \code{glm}
#' @export
family.simexreg <- function(object, ...) {
  if (inherits(object$NAIVEreg, "lm") | inherits(object$NAIVEreg, "glm")) {
    return(family(object$NAIVEreg, ...))
  } else return(NULL)
}


#' @describeIn simexreg Predict with new data
#' @export
predict.simexreg <- function(object, ...){
  reg <- object$NAIVEreg
  if (identical(class(reg), c("multinom", "nnet"))) {
    if(length(reg$lev) == 2) {
      coef_index <- 1+(1:length(reg$vcoefnames))
    } else {
      coef_index <- as.vector(t(matrix(1:length(reg$wts), nrow = reg$n[3],
                                       byrow=TRUE)[, 1+(1:length(reg$vcoefnames)),
                                                   drop=FALSE][-1, , drop=FALSE]))
    }
    reg$wts[coef_index] <- object$SIMEXcoef
  } else if (inherits(reg, "svymultinom")) {
    if(length(reg$NAIVEreg$lev) == 2) {
      coef_index <- 1+(1:length(reg$NAIVEreg$vcoefnames))
    } else {
      coef_index <- as.vector(t(matrix(1:length(reg$NAIVEreg$wts), nrow = reg$NAIVEreg$n[3],
                                       byrow=TRUE)[, 1+(1:length(reg$NAIVEreg$vcoefnames)),
                                                   drop=FALSE][-1, , drop=FALSE]))
    }
    reg$NAIVEreg$wts[coef_index] <- object$SIMEXcoef
  } else if (identical(class(reg), "polr")) {
    reg$coefficients <- object$SIMEXcoef[1:length(coef(reg))]
    reg$zeta <- object$SIMEXcoef[(length(coef(reg)) + 1):length(object$SIMEXcoef)]
  } else reg$coefficients <- object$SIMEXcoef
  out <- predict(reg, ...)
  return(out)
}


#' @describeIn simexreg Extract the model frame
#' @export
model.frame.simexreg <- function(formula, ...) {
  return(model.frame(formula$NAIVEreg, ...))
}


#' @describeIn simexreg Print results of \code{simexreg} nicely
#' @export
print.simexreg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(paste("\nNaive coefficient estimates: \n"))
  print(coef(x$NAIVEreg))
  cat("\nVariable measured with error:\n")
  cat(x$ME$MEvariable)
  cat("\nMeasurement error:\n")
  print(x$ME$MEerror)
  cat("\nError-corrected coefficient estimates:\n")
  print(x$SIMEXcoef)
}


#' @describeIn simexreg Summarize results of \code{simexreg} nicely
#' @export
summary.simexreg <- function(object, ...) {
  if (!object$ME$variance) stop("Set variance = TRUE when fitting simexreg")
  RCcoef <- coef(object)
  RCse <- sqrt(diag(vcov(object)))
  RCt <- RCcoef/RCse
  df <- nrow(model.frame(object)) - length(RCcoef)
  RCp <- 2 * pt(-abs(RCt), df)
  summarydf <- data.frame(RCcoef, RCse, RCt, RCp)
  colnames(summarydf) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(summarydf) <- names(RCcoef)
  out <- c(object, list(summarydf = summarydf))
  class(out) <- "summary.simexreg"
  return(out)
}


#' @describeIn simexreg Print summary of \code{simexreg} nicely
#' @export
print.summary.simexreg <- function(x, digits = 4, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNaive coefficient estimates: \n")
  print(coef(x$NAIVEreg))
  cat("\nNaive var-cov estimates: \n")
  print(vcov(x$NAIVEreg))
  cat("\nVariable measured with error:\n")
  cat(x$ME$MEvariable)
  cat("\nMeasurement error:\n")
  print(x$ME$MEerror)
  cat("\nError-corrected results:\n")
  printCoefmat(x$summarydf, digits = digits)
}


#' @describeIn simexreg Update \code{simexreg}
#' @export
update.simexreg <- function (object, ..., evaluate = TRUE) {
  simexreg_call <- getCall(object)
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(simexreg_call)))
    for (a in names(extras)[existing]) simexreg_call[[a]] <- extras[[a]]
    if (any(!existing)) {
      simexreg_call <- c(as.list(simexreg_call), extras[!existing])
      simexreg_call <- as.call(simexreg_call)
    }
  }
  if (evaluate) 
    eval(simexreg_call, parent.frame())
  else simexreg_call
}