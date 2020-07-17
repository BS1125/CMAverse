# multinomial regression with survey data
svymultinom <- function(formula = NULL, weights = NULL, data = NULL) {
  cl <- match.call()
  .formula <- formula
  out <- list()
  reg <- nnet::multinom(formula = as.formula(.formula), data = data,
                        weights = weights, trace = FALSE, model = TRUE)
  out$model <- model.frame(reg)
  out$coef <- coef(reg)
  coefnames <- dimnames(out$coef)
  coefnames <- switch((!is.null(coefnames[[1]])) + 1, "1" =  coefnames[[2]],
                      "2" = as.vector(outer(coefnames[[2]], coefnames[[1]], function(name2, name1)
                        paste(name1, name2, sep = ":"))))
  if (is.null(weights)) {
    svydes <- survey::as.svrepdesign(survey::svydesign(id = ~1, data = data), type = "JK1")
    out$vcov <- vcov(survey::withReplicates(svydes, quote(as.vector(t(coef(nnet::multinom(
      formula = as.formula(.formula), weights = .weights, data = data, trace = FALSE)))))))
    dimnames(out$vcov) <- list(coefnames, coefnames)
  } else {
    svydes <- survey::as.svrepdesign(survey::svydesign(id = ~1, weights = ~weights, data = data), type = "JK1")
    out$vcov <- vcov(survey::withReplicates(svydes, quote(as.vector(t(coef(nnet::multinom(
      formula = as.formula(.formula), weights = .weights, data = data, trace = FALSE)))))))
    dimnames(out$vcov) <- list(coefnames, coefnames)
  }
  out$call <- cl
  class(out) <- c("svymultinom", "multinom", "nnet")
  return(out)
}

coef.svymultinom <- function(svymultinom) {
  svymultinom$coef
}

vcov.svymultinom <- function(svymultinom) {
  svymultinom$vcov
}

model.frame.svymultinom <- function(svymultinom) {
  svymultinom$model
}

print.svymultinom <- function(svymultinom) {
  cat("Call: \n")
  print(svymultinom$call)
  cat("\n Coefficients: \n")
  print(svymultinom$coef)
}

# y <- factor(sample(c("a","b","c"),100,prob = c(0.2,0.3,0.5),replace = T),
#             level = c("a", "b", "c", "d"))
# x <- rnorm(100)
# data <- data.frame(y,x)
# reg=nnet::multinom(y~x,data=data)
# coef(reg)
# vcov(reg)
# lalala="y~x"
# try=eval(bquote((svymultinom(formula = "y~x", weights = NULL, data = data))))
# vcov(try)
# coef(try)
# model.frame(try)

