# multinomial regression with survey data
svymultinom <- function(formula = NULL, weights = NULL, data = NULL) {
  cl <- match.call()
  out <- list()
  environment(formula) <- environment()
  reg <- nnet::multinom(formula = formula, data = data,
                  weights = weights, trace = FALSE)
  out$model <- model.frame(reg)
  out$coef <- coef(reg)
  coefnames <- dimnames(out$coef)
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
  dimnames(out$vcov) <- list(coefnames, coefnames)
  out$call <- cl
  class(out) <- c("svymultinom", "multinom", "nnet")
  return(out)
}

coef.svymultinom <- function(object, ...) {
  object$coef
}

vcov.svymultinom <- function(x, ...) {
  x$vcov
}

model.frame.svymultinom <- function(x, ...) {
  x$model
}

print.svymultinom <- function(x, ...) {
  cat("Call: \n")
  print(x$call)
  cat("\n Coefficients: \n")
  print(x$coef)
}

