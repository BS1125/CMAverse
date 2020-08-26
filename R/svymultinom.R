#' Multinomial Regression with Survey Data
#'
#' \code{svymultinom} uses the \link[survey]{withReplicates} function to compute the replicate-based 
#' estimate of the variance-covariance matrix of coefficients for a multinomial regression fitted by 
#' \link[nnet]{multinom}.
#'
#' @param formula regression formula
#' @param data dataset for regression
#' @param weights weights for regression
#' @param x an object of class 'svymultinom'
#' @param object an object of class 'svymultinom'
#' @param digits minimal number of significant digits. See \link{print.default}.
#' @param evaluate a logical value. If \code{TRUE}, the updated call is evaluated. Default
#' is \code{TRUE}.
#' @param ... additional arguments
#' 
#' @return
#' An object of class 'svymultinom' is returned:
#' \item{call}{the function call,}
#' \item{NAIVEreg}{the naive multinomial regression object,}
#' \item{vcov}{the replicate-based estimate of the variance-covariance matrix of coefficients,}
#' ...
#'
#' @seealso \link[nnet]{multinom}, \link[survey]{withReplicates}
#' 
#' @examples
#' 
#' rm(list=ls())
#' library(CMAverse)
#' 
#' # multinom
#' n <- 1000
#' x1 <- rnorm(n, mean = 0, sd = 1)
#' x2 <- rnorm(n, mean = 1, sd = 1)
#' x3 <- rbinom(n, size = 1, prob = 0.4)
#' linearpred1 <- 1 + 0.3 * x1 - 0.5 * x2 - 0.2 * x3
#' linearpred2 <- 2 + 1 * x1 - 2 * x2 - 1 * x3
#' py2 <- exp(linearpred1) / (1 + exp(linearpred1) + exp(linearpred2))
#' py3 <- exp(linearpred2) / (1 + exp(linearpred1) + exp(linearpred2))
#' py1 <- 1 - py2 - py3
#' y <- sapply(1:n, function(x) sample(size = 1, c(1:3), prob = c(py1[x], py2[x], py3[x])))
#' w <- ifelse(x3 == 0, 0.4, 0.6)
#' data <- data.frame(x1 = x1, x2 = x2, x3 = x3, y = y)
#' 
#' reg <- svymultinom(y ~ x1 + x2 + x3, weights = w, data = data)
#' coef(reg)
#' vcov(reg)
#' formula(reg)
#' predict(reg, newdata = data[1, ], type = "probs")
#' model.frame(reg)
#' summary(reg)
#' update(reg, weights = w[1:500], data = data[1:500, ])
#'  
#' @importFrom stats model.frame family coef predict model.matrix getCall cov 
#' formula vcov pt
#' @importFrom nnet multinom
#' @importFrom survey as.svrepdesign svydesign withReplicates
#' 
#' @export
#' 
svymultinom <- function(formula = NULL, weights = NULL, data = NULL) {
  cl <- match.call()
  out <- list()  
  out$call <- cl
  environment(formula) <- environment()
  reg <- multinom(formula = formula, data = data, weights = weights, trace = FALSE)
  out$NAIVEreg <- reg
  coefnames <- dimnames(coef(reg))
  coefnames <- switch((!is.null(coefnames[[1]])) + 1, "1" =  coefnames[[2]],
                      "2" = as.vector(outer(coefnames[[2]], coefnames[[1]], function(name2, name1)
                        paste(name1, name2, sep = ":"))))
  # survey design
  if (is.null(weights)) svydes <- as.svrepdesign(svydesign(~1, data = data), type = "JK1")
  if (!is.null(weights)) svydes <- as.svrepdesign(svydesign(~1, weights = ~weights, 
                                                            data = data), type = "JK1")
  # survey vcov
  out$vcov <- vcov(withReplicates(svydes, function(w, data) {
    environment(formula) <- environment()
    as.vector(t(coef(eval(bquote(nnet::multinom(
      formula = formula, weights = w, data = data, trace = FALSE))))))
  }))
  attributes(out$vcov)$means <- NULL
  dimnames(out$vcov) <- list(coefnames, coefnames)
  class(out) <- "svymultinom"
  return(out)
}

#' @describeIn svymultinom Extract coefficients
#' @export
coef.svymultinom <- function(object, ...) {
  return(coef(object$NAIVEreg, ...))
}

#' @describeIn svymultinom Extract the var-cov matrix of coefficients
#' @export
vcov.svymultinom <- function(object, ...) {
  return(object$vcov)
}

#' @describeIn svymultinom Extract the regression formula
#' @export
formula.svymultinom <- function(x, ...) {
  return(formula(x$NAIVEreg, ...))
}

#' @describeIn svymultinom Predict with new data
#' @export
predict.svymultinom <- function(object, ...) {
  return(predict(object$NAIVEreg, ...))
}

#' @describeIn svymultinom Extract the model frame
#' @export
model.frame.svymultinom <- function(formula, ...) {
  return(model.frame(formula$NAIVEreg, ...))
}

#' @describeIn svymultinom Print the results of \code{svymultinom} nicely
#' @export
print.svymultinom <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat(paste("\n# Naive regression object: \n"))
  print(x$NAIVEreg)
  cat("\n# Var-cov matrix of coefficients:\n")
  print(x$vcov)
}

#' @describeIn svymultinom Summarize results of \code{svymultinom} nicely
#' @export
summary.svymultinom <- function(object, ...) {
  coef <- as.vector(t(coef(object)))
  se <- sqrt(diag(vcov(object)))
  t <- coef/se
  df <- nrow(model.frame(object)) - length(coef)
  p <- 2 * pt(-abs(t), df)
  summarydf <- data.frame(coef, se, t, p)
  colnames(summarydf) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
  rownames(summarydf) <- names(se)
  out <- c(object, list(summarydf = summarydf))
  class(out) <- "summary.svymultinom"
  return(out)
}

#' @describeIn svymultinom Print summary of \code{svymultinom} nicely
#' @export
print.summary.svymultinom <- function(x, digits = 4, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$summarydf, digits = digits)
}

#' @describeIn svymultinom Update \code{svymultinom}
#' @export
update.svymultinom <- function (object, ..., evaluate = TRUE) {
  svymultinom_call <- getCall(object)
  extras <- match.call(expand.dots = FALSE)$...
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(svymultinom_call)))
    for (a in names(extras)[existing]) svymultinom_call[[a]] <- extras[[a]]
    if (any(!existing)) {
      svymultinom_call <- c(as.list(svymultinom_call), extras[!existing])
      svymultinom_call <- as.call(svymultinom_call)
    }
  }
  if (evaluate) 
    eval(svymultinom_call, parent.frame())
  else svymultinom_call
}


