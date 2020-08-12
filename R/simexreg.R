#' @export
simexreg <- function (reg = NULL, data = NULL, weights = NULL, model = TRUE,
                      MEvariable = NULL, MEvartype = NULL, MEerror = NULL,
                      variance = FALSE, lambda = c(0.5, 1, 1.5, 2), B = 200) {

  cl <- match.call()
  
  regCall <- getCall(reg)
  regClass <- class(reg)
  if (!(identical(regClass, "lm") | 
        (identical(regClass, c("glm", "lm")) && 
         family(reg)$family %in% c("gaussian", "binomial", "poisson")) | 
        identical(regClass, c("multinom", "nnet")) | identical(regClass, "polr") | 
        identical(regClass, "coxph") | identical(regClass, "survreg"))) stop(
          "Measurement error correction applied to unsupported regression object")
  
  if (length(MEvariable) > 1) stop("length(MEvariable) > 1")
  if (length(MEvariable) != length(MEvartype)) stop("length(MEvariable) != length(MEvartype)")
  
  if (MEvartype == "continuous") {
    if (length(MEvariable) != length(MEerror)) stop("length(MEvariable) != length(MEerror)")
    if (length(MEerror) != 0 && MEerror < 0) stop("MEerror should be >= 0")
  } else if (MEvartype == "categorical") {
    if (length(MEvariable) != 0 && !is.matrix(MEerror)) stop("MEerror should be a matrix")
    if (dim(MEerror)[1] != dim(MEerror)[2] |
        dim(MEerror)[1] != length(unique(data[, MEvariable]))) stop("Incorrect dimension of MEerror")
    if (!check.mc.matrix(list(MEerror))) stop("MEerror may contain negative values for exponents smaller than 1")
  } else stop("Unsupported MEvartype; use 'continuous' or 'categorical'")
  if (is.null(data)) stop("Unspecified data")
  if (!all(lambda>0)) stop("lambda should be positive")
  
  # update reg with data and weights
  regCall$data <- data
  regCall$weights <- weights
  if (inherits(reg, "multinom")) regCall$trace <- FALSE
  if (inherits(reg, "polr")) regCall$Hess <- TRUE
  reg <- eval.parent(regCall)
  
  # the vector of variable names in the regression formula
  reg_formula <- formula(reg)
  var_vec <- unique(all.vars(reg_formula))

  if (length(MEvariable) == 0 | !MEvariable %in% var_vec) {
    out <- reg
  } else {
    n <- nrow(data)
    # output list
    out <- list(call = cl, NAIVEreg = reg,
                ME = list(MEvariable = MEvariable, MEvartype = MEvartype,
                          MEerror = MEerror, variance = variance, lambda = lambda, B = B))
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
            MEvar_lev <- levels(droplevels(as.factor(SIMdata[, MEvariable])))
            category <- unique(data[, MEvariable])[match(MEvar_lev, unique(data[, MEvariable]))]
            mcmdecomp <- eigen(MEerror)
            SIMmcm <- mcmdecomp$vectors %*% diag(mcmdecomp$values)^lambda[i] %*% solve(mcmdecomp$vectors)
            for (j in 1:length(category)) {
              SIMdata[, MEvariable][which(SIMdata[, MEvariable] == MEvar_lev[j])] <-
                sample(x = category, size = length(which(SIMdata[, MEvariable] == category[j])),
                       prob = SIMmcm[, j], replace = TRUE)
              }
          }
        regCall$data <- SIMdata
        SIMreg <- eval.parent(regCall)
        if (identical(class(SIMreg), "polr")) {
          SIMcoef_mid[b, ] <- c(coef(SIMreg), SIMreg$zeta)
        }  else if (identical(class(SIMreg), c("multinom", "nnet"))) {
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
    
    if (((identical(regClass, c("glm", "lm")) && family(reg)$family == "gaussian") |
         identical(regClass, "lm")) && MEvartype == "continuous") {
      reg_fit <- reg
      reg_fit$coefficients <- SIMEXcoef
      SIMEXsigma <- sqrt(sum((model.frame(reg_fit)[, 1] - predict(reg_fit, newdata = data)) ^ 2) / 
                           (n - ncoef) - (SIMEXcoef[MEvariable] * MEerror) ^ 2)
      out$SIMEXsigma <- SIMEXsigma
    }

    class(out) <- "simexreg"

  }

  return(out)

}

#' @export
coef.simexreg <- function(object, ...) {
  return(object$SIMEXcoef)
}

#' @export
vcov.simexreg <- function(object, ...) {
  return(object$SIMEXvcov)
}

#' @export
sigma.simexreg <- function(object, ...) {
  return(object$SIMEXsigma)
}

#' @export
formula.simexreg <- function(x, ...) {
  return(formula(x$NAIVEreg))
}

#' @export
family.simexreg <- function(object, ...) {
  if (inherits(object$NAIVEreg, "lm") | inherits(object$NAIVEreg, "glm")) {
    return(family(object$NAIVEreg, ...))
  } else return(NULL)
}

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
  } else if (identical(class(reg), "polr")) {
    reg$coefficients <- object$SIMEXcoef[1:length(coef(reg))]
    reg$zeta <- object$SIMEXcoef[(length(coef(reg)) + 1):length(object$SIMEXcoef)]
  } else reg$coefficients <- object$SIMEXcoef
  out <- predict(reg, ...)
  return(out)
}

#' @export
model.frame.simexreg <- function(formula, ...) {
  return(model.frame(formula$NAIVEreg))
}

#' @export
print.simexreg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(paste("\nNaive regression object: \n"))
  print(x$NAIVEreg, ...)
  cat("\nVariable measured with error:\n")
  cat(x$ME$MEvariable)
  cat("\nMeasurement error:\n")
  cat(x$ME$MEerror)
  cat("\nError-corrected coefficient estimates:\n")
  print(x$SIMEXcoef)
}

