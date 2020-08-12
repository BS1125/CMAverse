#' @export
rcreg <- function(reg = NULL, data = NULL, weights = NULL,
                  MEvariable = NULL, MEerror = NULL, variance = FALSE, nboot = 400) {

  cl <- match.call()
  
  reg_formula <- formula(reg)
  # the vector of names of all variables in the regression formula 
  var_vec <- unique(all.vars(reg_formula))
  # the vector of names of all independent variables in the regression formula
  ind_var <- unique(all.vars(reg_formula[[3]]))

  if (length(MEvariable) == 0 | !MEvariable %in% var_vec) {
    out <- reg
  } else if (!MEvariable %in% ind_var) {
    stop("Regression calibration only supports an independent variable in the regression formula measured with error")
  } else if (length(ind_var) == 1) {
    stop("No independent variable measured without error in the regression formula")
  } else {
    if (is.null(data)) stop("Unspecified data")
    if (length(MEvariable) > 1) stop("length(MEvariable) > 1")
    if (length(MEvariable) != length(MEerror)) stop("length(MEvariable) != length(MEerror)")
    if (length(MEerror) != 0 && MEerror < 0) stop("MEerror should be >= 0")
    regCall <- getCall(reg)
    regClass <- class(reg)
    if (!(identical(regClass, "rcreg") | identical(regClass, "lm") | 
          (identical(regClass, c("glm", "lm")) && 
           family(reg)$family %in% c("gaussian", "binomial", "poisson")) | 
          identical(regClass, c("multinom", "nnet")) | identical(regClass, "polr") | 
          identical(regClass, "coxph") | identical(regClass, "survreg"))) stop(
            "Measurement error correction applied to unsupported regression object")

    # update reg with data and weights
    regCall$data <- data
    regCall$weights <- weights
    if (inherits(reg, "multinom")) regCall$trace <- FALSE
    if (inherits(reg, "polr")) regCall$Hess <- TRUE
    
    reg <- eval.parent(regCall)
    out <- list(call = cl, NAIVEreg = reg,
                ME = list(MEvariable = MEvariable, MEerror = MEerror, variance = variance))

    MEvar_index <- which(ind_var %in% MEvariable)
    # code categorical variables into dummy variables
    ind_data <- model.matrix(as.formula(paste0("~", paste(c(MEvariable, ind_var[-MEvar_index]),
                                                          collapse = "+"))),
                             data = data[, ind_var])[, -1]

    rc_step <- function(data = NULL, ind_data = NULL, indices = NULL, outreg = FALSE) {
      data <- data[indices, ]
      ind_data <- ind_data[indices, ]
      n <- nrow(ind_data)
      # mean values of independent variables in the regression formula
      mean <- colMeans(ind_data)
      W <- ind_data - t(mean)[rep(1, n), ]
      # var-cov matrix of variables measured with error with all independent variables
      covmat <- cov(ind_data[, MEvariable], ind_data)
      # true var-cov matrix of variables measured with error with all independent variables
      covmat[1, 1] <- covmat[1, 1] - MEerror ^ 2
      EX <- as.vector(rep(mean[1], n) + W %*% t(covmat %*% solve(cov(ind_data))))
      rcdata <- data
      rcdata[, MEvariable] <- EX
      # error-corrected regression object
      RCreg_call <- getCall(reg)
      RCreg_call$data <- rcdata
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
      } else RCcoef <- coef(RCreg)
      if (!outreg) out <- RCcoef
      if (outreg) out <- list(RCreg = RCreg, RCcoef = RCcoef)
      return(out)
    }

    RC <- rc_step(data = data, ind_data = ind_data, indices = 1:nrow(data), outreg = TRUE)
    RCcoef <- RC$RCcoef
    RCreg <- RC$RCreg
    out$RCcoef <- RCcoef
    n <- nrow(data)
    n_coef <- length(RCcoef)
    if (identical(class(reg), "lm") |
        (identical(class(reg), c("glm", "lm")) && family(reg)$family == "gaussian")) {
      RCsigma <- sqrt(sum((model.frame(RCreg)[, 1] - predict(RCreg, newdata = data)) ^ 2) /
                        (n-n_coef) - (RCcoef[MEvariable]* MEerror) ^ 2 )
      out$RCsigma <- unname(RCsigma)
    }

    if (variance) {
      boots <- boot(data = data, ind_data = ind_data, statistic = rc_step, R = nboot)
      mean_boots <- apply(boots$t, 2, mean)
      RCvcov <- Reduce("+", lapply(1:nboot, function(x)
        as.matrix(boots$t[x,] - mean_boots)%*%t(as.matrix(boots$t[x,] - mean_boots)))) /
        (nboot - 1)
      dimnames(RCvcov) <- list(names(RCcoef), names(RCcoef))
      out$RCvcov <- RCvcov
      out$ME$nboot <- nboot
    }
  }

  class(out) <- c("rcreg")
  return(out)

}

#' @export
coef.rcreg <- function(object, ...) {
  return(object$RCcoef)
}

#' @export
vcov.rcreg <- function(object, ...) {
  return(object$RCvcov)
}

#' @export
sigma.rcreg <- function(object, ...) {
  return(object$RCsigma)
}

#' @export
formula.rcreg <- function(x, ...) {
  return(formula(x$NAIVEreg, ...))
}

#' @export
family.rcreg <- function(object, ...) {
  if (inherits(object$NAIVEreg, "lm") | inherits(object$NAIVEreg, "glm")) {
    return(family(object$NAIVEreg, ...))
  } else return(NULL)
}

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
  } else if (identical(class(reg), "polr")) {
    reg$coefficients <- object$RCcoef[1:length(coef(reg))]
    reg$zeta <- object$RCcoef[(length(coef(reg)) + 1):length(object$RCcoef)]
  } else reg$coefficients <- object$RCcoef
  out <- predict(reg, ...)
  return(out)
}

#' @export
model.frame.rcreg <- function(formula, ...) {
  return(model.frame(formula$NAIVEreg))
}

#' @export
print.rcreg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(paste("\nNaive regression object: \n"))
  print(x$NAIVEreg, ...)
  cat("\nVariable measured with error:\n")
  cat(x$ME$MEvariable)
  cat("\nMeasurement error:\n")
  cat(x$ME$MEerror)
  cat("\nError-corrected coefficient estimates:\n")
  print(x$RCcoef)
}

