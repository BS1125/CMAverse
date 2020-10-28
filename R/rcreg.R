#' Regression Calibration for Measurement Error Correction
#'
#' \code{rcreg} is used to correct a regression object with a continuous independent variable
#' measured with error via \emph{regression calibration} by Carroll et al. (1995).
#'
#' @param reg naive regression object. See \code{Details}.
#' @param formula regression formula
#' @param data new dataset for \code{reg}
#' @param weights new weights for \code{reg}
#' @param MEvariable variable measured with error
#' @param MEerror standard deviation of the measurement error
#' @param variance a logical value. If \code{TRUE}, correct the var-cov matrix of
#' coefficients by bootstrapping. Default is \code{FALSE}.
#' @param nboot number of boots for correcting the var-cov matrix of coefficients. Default 
#' is \code{400}.
#' @param x an object of class \code{rcreg}
#' @param object an object of class \code{rcreg}
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
#' \code{MEvariable} is a continuous independent variable in the regression formula, an 
#' object of class \code{rcreg} is returned:
#' \item{call}{the function call,}
#' \item{NAIVEreg}{the naive regression object,}
#' \item{ME}{a list of \code{MEvariable}, \code{MEerror}, \code{variance} and \code{nboot},}
#' \item{RCcoef}{coefficient estimates corrected by regression calibration,}
#' \item{RCsigma}{the residual standard deviation of a linear regression object corrected by 
#' regression calibration,}
#' \item{RCvcov}{the var-cov matrix of coefficients corrected by regression calibration,}
#' ...
#'
#' @seealso \code{\link{simexreg}}, \code{cmsens}, \code{\link{cmest}}
#'
#' @references
#' 
#' Carrol RJ, Ruppert D, Stefanski LA, Crainiceanu C (2006). Measurement Error in Nonlinear Models: 
#' A Modern Perspective, Second Edition. London: Chapman & Hall.
#'
#' @examples
#' 
#' \dontrun{
#' rm(list=ls())
#' library(CMAverse)
#' 
#' # 2 boots are used for illustration
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
#' reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEerror = 0.5, 
#' variance = TRUE, nboot = 2)
#' coef(reg_rc)
#' vcov(reg_rc)
#' sigma(reg_rc)
#' formula(reg_rc)
#' family(reg_rc)
#' predict(reg_rc, newdata = data[1, ])
#' reg_rc_model <- model.frame(reg_rc)
#' reg_rc_update <- update(reg_rc, data = data, weights = rep(1, n))
#' reg_rc_summ <- summary(reg_rc)
#' 
#' #glm
#' n <- 1000
#' x1 <- rnorm(n, mean = 0, sd = 1)
#' x2_true <- rnorm(n, mean = 1, sd = 1)
#' error1 <- rnorm(n, mean = 0, sd = 0.5)
#' x2_error <- x2_true + error1
#' x3 <- rbinom(n, size = 1, prob = 0.4)
#' linearpred <- 1 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
#' py <- exp(linearpred) / (1 + exp(linearpred))
#' y <- rbinom(n, size = 1, prob = py)
#' data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
#'                    x3 = x3, y = y)
#' reg_naive <- glm(y ~ x1 + x2_error + x3, data = data, family = binomial("logit"))
#' reg_true <- glm(y ~ x1 + x2_true + x3, data = data, family = binomial("logit"))
#' reg_rc <- rcreg(reg = reg_naive, data = data, 
#' MEvariable = "x2_error", MEerror = 0.5, variance = TRUE, nboot = 2)
#' 
#' # multinom
#' n <- 1000
#' x1 <- rnorm(n, mean = 0, sd = 1)
#' x2_true <- rnorm(n, mean = 1, sd = 1)
#' error1 <- rnorm(n, mean = 0, sd = 0.5)
#' x2_error <- x2_true + error1
#' x3 <- rbinom(n, size = 1, prob = 0.4)
#' linearpred1 <- 1 + 0.3 * x1 - 0.5 * x2_true - 0.2 * x3
#' linearpred2 <- 2 + 1 * x1 - 2 * x2_true - 1 * x3
#' py2 <- exp(linearpred1) / (1 + exp(linearpred1) + exp(linearpred2))
#' py3 <- exp(linearpred2) / (1 + exp(linearpred1) + exp(linearpred2))
#' py1 <- 1 - py2 - py3
#' y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
#' data <- data.frame(x1 = x1, x2_true = x2_true, x2_error = x2_error,
#'                    x3 = x3, y = y)
#' reg_naive <- nnet::multinom(factor(y) ~ x1 + x2_error + x3, data = data)
#' reg_true <- nnet::multinom(factor(y) ~ x1 + x2_true + x3, data = data)
#' reg_rc <- rcreg(reg = reg_naive, data = data, MEvariable = "x2_error", MEerror = 0.5, 
#' variance = TRUE, nboot = 2)
#' }                
#'
#' @importFrom stats as.formula model.frame family coef predict model.matrix getCall cov 
#' formula vcov pt na.pass
#' @importFrom boot boot
#' 
#' @export
#' 
rcreg <- function(reg = NULL, formula = NULL, data = NULL, weights = NULL,
                  MEvariable = NULL, MEerror = NULL, variance = FALSE, nboot = 400) {
  
  cl <- match.call()
  
  # assign svyglm the global environment
  assign2glob <- function(key, val, pos) assign(key, val, envir = as.environment(pos))
  assign2glob("svyglm", survey::svyglm, 1L)
  
  if (is.null(formula)) formula <- as.formula(formula)
  if (!is.null(formula)) formula <- formula(reg)
  # the vector of names of all variables in the regression formula 
  var_vec <- unique(all.vars(formula))
  # the vector of names of all independent variables in the regression formula
  ind_var <- unique(all.vars(formula[[3]]))
  
  if (length(MEvariable) == 0 | !MEvariable %in% var_vec) {
    out <- reg
  } else if (!MEvariable %in% ind_var) {
    stop("Regression calibration only supports an independent variable in the regression formula measured with error")
  } else if (length(ind_var) == 1) {
    stop("No independent variable measured without error in the regression formula")
  } else {
    if (is.null(data)) stop("Unspecified data")
    if (length(MEvariable) > 1) stop("Currently only supports one variable measured with error")
    if (length(MEvariable) != length(MEerror)) stop("length(MEvariable) != length(MEerror)")
    if (length(MEerror) != 0 && MEerror < 0) stop("MEerror should be >= 0")
    regCall <- getCall(reg)
    regClass <- class(reg)
    if (!(identical(regClass, c("svyglm", "glm", "lm")) | identical(regClass, c("svymultinom")) |
          identical(regClass, "lm") | 
          (identical(regClass, c("glm", "lm")) && 
           family(reg)$family %in% c("gaussian", "binomial", "poisson")) | 
          identical(regClass, c("multinom", "nnet")) | identical(regClass, "polr") | 
          identical(regClass, "coxph") | identical(regClass, "survreg"))) stop(
            "Regression calibration applied to unsupported regression object")
    
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
    }
    reg <- eval.parent(regCall)
    out <- list(call = cl, NAIVEreg = reg,
                ME = list(MEvariable = MEvariable, MEerror = MEerror, variance = variance))
    
    MEvar_index <- which(ind_var %in% MEvariable)
    # code categorical variables into dummy variables
    ind_data <- model.matrix(as.formula(paste0("~", paste(c(MEvariable, ind_var[-MEvar_index]), collapse = "+"))),
                             model.frame(~., data = data[, ind_var], na.action = na.pass))[, -1]
    
    rc_step <- function(data = NULL, ind_data = NULL, indices = NULL, outreg = FALSE) {
      data <- data[indices, ]
      ind_data <- ind_data[indices, ]
      # mean values of independent variables in the regression formula
      mean <- colMeans(ind_data, na.rm = TRUE)
      W <- ind_data - t(mean)[rep(1, n), ]
      # var-cov matrix of the variable measured with error with all independent variables
      covmat <- cov(ind_data[, MEvariable], ind_data, use = "complete.obs")
      # true var-cov matrix of the variable measured with error with all independent variables
      covmat[1, 1] <- covmat[1, 1] - MEerror ^ 2
      EX <- as.vector(rep(mean[1], n) + W %*% t(covmat %*% solve(cov(ind_data, use = "complete.obs"))))
      rcdata <- data
      rcdata[, MEvariable] <- EX
      
      # error-corrected regression object
      RCreg_call <- getCall(reg)
      if (identical(regClass, c("svyglm", "glm", "lm"))) {
        designCall <- getCall(reg$survey.design)
        designCall$data <- rcdata
        RCreg_call$design <- eval.parent(designCall)
      } else {
        RCreg_call$data <- rcdata
        if (inherits(reg, "multinom")) RCreg_call$trace <- FALSE
        if (inherits(reg, "polr")) RCreg_call$Hess <- TRUE
      }
      RCreg <- eval.parent(RCreg_call)
      
      if (identical(class(reg), "polr")) {
        RCcoef <- c(coef(RCreg), RCreg$zeta)
      } else if (identical(class(reg), c("multinom", "nnet"))) {
        RCcoef <- coef(RCreg)
        if (length(reg$lev) > 2) {
          coefnames <- dimnames(RCcoef)
          RCcoef <- as.vector(t(RCcoef))
          coefnames <- as.vector(outer(coefnames[[2]], coefnames[[1]], function(name2, name1)
            paste(name1, name2, sep = ":")))
          names(RCcoef) <- coefnames
        }
      } else if (inherits(reg, "svymultinom")) {
        RCcoef <- coef(RCreg)
        if (length(reg$NAIVEreg$lev) > 2) {
          coefnames <- dimnames(RCcoef)
          RCcoef <- as.vector(t(RCcoef))
          coefnames <- as.vector(outer(coefnames[[2]], coefnames[[1]], function(name2, name1)
            paste(name1, name2, sep = ":")))
          names(RCcoef) <- coefnames
        }
      } else RCcoef <- coef(RCreg)
      if (!outreg) out <- RCcoef
      if (outreg) out <- list(RCreg = RCreg, RCcoef = RCcoef)
      return(out)
    }
    
    n <- nrow(data)
    RC <- rc_step(data = data, ind_data = ind_data, indices = 1:n, outreg = TRUE)
    RCcoef <- RC$RCcoef
    RCreg <- RC$RCreg
    out$RCcoef <- RCcoef
    n_coef <- length(RCcoef)
    if (identical(class(reg), "lm") |
        (identical(class(reg), c("glm", "lm")) && family(reg)$family == "gaussian")) {
      RCsigma <- sqrt(sum((model.frame(RCreg)[, 1] - predict(RCreg, newdata = data)) ^ 2, na.rm = TRUE) /
                        (n-n_coef) - (RCcoef[MEvariable]* MEerror) ^ 2 )
      out$RCsigma <- unname(RCsigma)
    }
    
    if (variance) {
      boots <- boot(data = data, ind_data = ind_data, statistic = rc_step, R = nboot)
      mean_boots <- apply(boots$t, 2, function(x) mean(x, na.rm = TRUE))
      RCvcov <- Reduce("+", lapply(1:nboot, function(x)
        as.matrix(boots$t[x,] - mean_boots)%*%t(as.matrix(boots$t[x,] - mean_boots)))) / (nboot - 1)
      dimnames(RCvcov) <- list(names(RCcoef), names(RCcoef))
      out$RCvcov <- RCvcov
      out$ME$nboot <- nboot
    }
    
    # remove svyglm from the global environment
    rm(svyglm, envir = .GlobalEnv)
    
    class(out) <- c("rcreg")
  }
  return(out)
}


#' @describeIn rcreg Extract coefficients corrected by \code{rcreg}
#' @export
coef.rcreg <- function(object, ...) {
  return(object$RCcoef)
}


#' @describeIn rcreg Extract the var-cov matrix of coefficients corrected by \code{rcreg}
#' @export
vcov.rcreg <- function(object, ...) {
  return(object$RCvcov)
}


#' @describeIn rcreg Extract the residual standard deviation of a linear regression object 
#' corrected by \code{rcreg}
#' @export
sigma.rcreg <- function(object, ...) {
  return(object$RCsigma)
}


#' @describeIn rcreg Extract the regression formula
#' @export
formula.rcreg <- function(x, ...) {
  return(formula(x$NAIVEreg, ...))
}


#' @describeIn rcreg Extract the family of a regression of class \code{lm} or \code{glm}
#' @export
family.rcreg <- function(object, ...) {
  if (inherits(object$NAIVEreg, "lm") | inherits(object$NAIVEreg, "glm")) {
    return(family(object$NAIVEreg, ...))
  } else return(NULL)
}


#' @describeIn rcreg Predict with new data
#' @export
predict.rcreg <- function(object, ...) {
  reg <- object$NAIVEreg
  if (identical(class(reg), c("multinom", "nnet"))) {
    if(length(reg$lev) == 2) {
      coef_index <- 1+(1:length(reg$vcoefnames))
    } else {
      coef_index <- as.vector(t(matrix(1:length(reg$wts), nrow = reg$n[3],
                                       byrow=TRUE)[, 1+(1:length(reg$vcoefnames)),
                                                   drop=FALSE][-1, , drop=FALSE]))
    }
    reg$wts[coef_index] <- object$RCcoef
  } else if (inherits(reg, "svymultinom")) {
    if(length(reg$NAIVEreg$lev) == 2) {
      coef_index <- 1+(1:length(reg$NAIVEreg$vcoefnames))
    } else {
      coef_index <- as.vector(t(matrix(1:length(reg$NAIVEreg$wts), nrow = reg$NAIVEreg$n[3],
                                       byrow=TRUE)[, 1+(1:length(reg$NAIVEreg$vcoefnames)),
                                                   drop=FALSE][-1, , drop=FALSE]))
    }
    reg$NAIVEreg$wts[coef_index] <- object$RCcoef
  } else if (identical(class(reg), "polr")) {
    reg$coefficients <- object$RCcoef[1:length(coef(reg))]
    reg$zeta <- object$RCcoef[(length(coef(reg)) + 1):length(object$RCcoef)]
  } else reg$coefficients <- object$RCcoef
  out <- predict(reg, ...)
  return(out)
}


#' @describeIn rcreg Extract the model frame
#' @export
model.frame.rcreg <- function(formula, ...) {
  return(model.frame(formula$NAIVEreg, ...))
}


#' @describeIn rcreg Print results of \code{rcreg} nicely
#' @export
print.rcreg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNaive coefficient estimates: \n")
  print(coef(x$NAIVEreg))
  cat("\nVariable measured with error:\n")
  cat(x$ME$MEvariable)
  cat("\nMeasurement error:\n")
  cat(x$ME$MEerror)
  cat("\nError-corrected coefficient estimates:\n")
  print(x$RCcoef)
}

#' @describeIn rcreg Summarize results of \code{rcreg} nicely
#' @export
summary.rcreg <- function(object, ...) {
  if (!object$ME$variance) stop("Set variance = TRUE when fitting rcreg")
  RCcoef <- coef(object)
  RCse <- sqrt(diag(vcov(object)))
  RCt <- RCcoef/RCse
  df <- nrow(model.frame(object)) - length(RCcoef)
  RCp <- 2 * pt(-abs(RCt), df)
  summarydf <- data.frame(RCcoef, RCse, RCt, RCp)
  colnames(summarydf) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(summarydf) <- names(RCcoef)
  out <- c(object, list(summarydf = summarydf))
  class(out) <- "summary.rcreg"
  return(out)
}


#' @describeIn rcreg Print summary of \code{rcreg} nicely
#' @export
print.summary.rcreg <- function(x, digits = 4, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\nNaive coefficient estimates: \n")
  print(coef(x$NAIVEreg))
  cat("\nNaive var-cov estimates: \n")
  print(vcov(x$NAIVEreg))
  cat("\nVariable measured with error:\n")
  cat(x$ME$MEvariable)
  cat("\nMeasurement error:\n")
  cat(x$ME$MEerror)
  cat("\nError-corrected results:\n")
  printCoefmat(x$summarydf, digits = digits)
}


#' @describeIn rcreg Update \code{rcreg}
#' @export
update.rcreg <- function (object, ..., evaluate = TRUE) {
  rcreg_call <- getCall(object)
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(rcreg_call)))
    for (a in names(extras)[existing]) rcreg_call[[a]] <- extras[[a]]
    if (any(!existing)) {
      rcreg_call <- c(as.list(rcreg_call), extras[!existing])
      rcreg_call <- as.call(rcreg_call)
    }
  }
  if (evaluate) 
    eval(rcreg_call, parent.frame())
  else rcreg_call
}
