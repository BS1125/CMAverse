simex_reg <- function (reg, data, MEvariable, MEvariable.type, measurement.error,
                         lambda = c(0.5, 1, 1.5, 2), B = 100) {

  if (!all(lambda>0)) {
    stop("lambda should be positive")
  }

  lambda <- c(0, lambda)

  n <- nrow(data)

  if (identical(class(reg), "polr")) {
    SIMcoef <- c(coef(reg), reg$zeta)
  } else if (identical(class(reg), c("multinom", "nnet"))) {
    SIMcoef <- as.vector(t(coef(reg)))
    names(SIMcoef) <- colnames(vcov(reg))
  } else {SIMcoef <- coef(reg)}

  if (identical(class(reg), "lm")) SIMsigma <- sigma(reg)

  SIMvar <- list(vcov(reg))

  ncoef <- length(SIMcoef)

  if (MEvariable.type == "continuous") {

    if (!length(measurement.error) == 1) {
      stop("measurement.error should be a single value for a continuous variable measured with error")
    }

    error <- rep(measurement.error, n)

    for (i in 2:length(lambda)) {

      SIMcoef_mid <- c()

      if (identical(class(reg), "lm")) SIMsigma_mid <- c()

      SIMvar_mid <- matrix(rep(0, ncoef^2), nrow = ncoef)

      for (b in 1:B) {

        Z <- rnorm(n=n)

        SIMdata <- data

        SIMdata[, MEvariable] <- SIMdata[, MEvariable] + (sqrt(lambda[i]) *  Z * measurement.error)

        SIMreg <- update(reg, data = SIMdata)

        if (identical(class(SIMreg), "polr")) {
          SIMcoef_mid <- rbind(SIMcoef_mid, c(coef(SIMreg), SIMreg$zeta))
        }  else if (identical(class(SIMreg), c("multinom", "nnet"))) {
          SIMcoef_mid <- rbind(SIMcoef_mid, as.vector(t(coef(SIMreg))))
        } else {SIMcoef_mid <- rbind(SIMcoef_mid, c(coef(SIMreg)))}

        if (identical(class(reg), "lm")) SIMsigma_mid <- c(SIMsigma_mid, sigma(SIMreg))

        SIMvar_mid <- SIMvar_mid + vcov(SIMreg)/B

      }

      SIMcoef <- rbind(SIMcoef, colMeans(SIMcoef_mid))

      if (identical(class(reg), "lm")) SIMsigma <- c(SIMsigma, mean(SIMsigma_mid))

      SIMvar[[i]]<- SIMvar_mid - cov(SIMcoef_mid)

    }

    extrap_coef <- lm(SIMcoef ~ lambda + I(lambda^2))

    SIMEXcoef <- predict(extrap_coef, newdata = data.frame(lambda = -1))

    if (identical(class(reg), "lm")) {

      extrap_sigma <- lm(SIMsigma ~ lambda + I(lambda^2))

      SIMEXsigma <- predict(extrap_sigma, newdata = data.frame(lambda = -1))

    }

    SIMvar.element <- matrix(unlist(SIMvar), nrow = length(lambda) , byrow = TRUE)

    extrap_var <- lm(SIMvar.element ~ lambda + I(lambda^2))

    SIMEXvar <- matrix(predict(extrap_var, newdata = data.frame(lambda = -1)),
                       nrow = length(SIMEXcoef))

    dimnames(SIMEXvar) <- list(colnames(SIMEXcoef), colnames(SIMEXcoef))

  } else if (MEvariable.type == "categorical") {

    if (!is.matrix(measurement.error)) {
      stop("measurement.error should be a matrix for a categorical variable measured with error")
    }

    data[, MEvariable] <- as.factor(data[, MEvariable])

    levels(data[, MEvariable]) <- 0:(length(levels(data[, MEvariable])) - 1)

    category <- levels(data[, MEvariable])

    mcmdecomp <- eigen(measurement.error)

    for (i in 2:length(lambda)) {

      SIMcoef_mid <- c()

      if (identical(class(reg), "lm")) SIMsigma_mid <- c()

      SIMvar_mid <- matrix(rep(0, ncoef^2), nrow = ncoef)

      for (b in 1:B) {

        SIMdata <- data

        SIMmcm <- mcmdecomp$vectors %*% diag(mcmdecomp$values)^lambda[i] %*% solve(mcmdecomp$vectors)

        for (j in 1:length(category)) {

          SIMdata[, MEvariable][which(data[, MEvariable] == category[j])] <-
            sample(x = category, size = length(which(data[, MEvariable] == category[j])),
                   prob = measurement.error[, j], replace = TRUE)
        }

        SIMreg <- update(reg, data = SIMdata)

        if (identical(class(SIMreg), "polr")) {
          SIMcoef_mid <- rbind(SIMcoef_mid, c(coef(SIMreg), SIMreg$zeta))
        }  else if (identical(class(SIMreg), c("multinom", "nnet"))) {
          SIMcoef_mid <- rbind(SIMcoef_mid, as.vector(t(coef(SIMreg))))
        } else {SIMcoef_mid <- rbind(SIMcoef_mid, c(coef(SIMreg)))}

        if (identical(class(reg), "lm")) SIMsigma_mid <- c(SIMsigma_mid, sigma(SIMreg))

        SIMvar_mid <- SIMvar_mid + vcov(SIMreg)/B

      }

      SIMcoef <- rbind(SIMcoef, colMeans(SIMcoef_mid))

      if (identical(class(reg), "lm")) SIMsigma <- c(SIMsigma, mean(SIMsigma_mid))

      SIMvar[[i]]<- SIMvar_mid - cov(SIMcoef_mid)

    }

    extrap_coef <- lm(SIMcoef ~ lambda + I(lambda^2))

    SIMEXcoef <- predict(extrap_coef, newdata = data.frame(lambda = -1))

    if (identical(class(reg), "lm")) {

      extrap_sigma <- lm(SIMsigma ~ lambda + I(lambda^2))

      SIMEXsigma <- predict(extrap_sigma, newdata = data.frame(lambda = -1))

    }

    SIMvar.element <- matrix(unlist(SIMvar), nrow = length(lambda) , byrow = TRUE)

    extrap_var <- lm(SIMvar.element ~ lambda + I(lambda^2))

    SIMEXvar <- matrix(predict(extrap_var, newdata = data.frame(lambda = -1)),
                       nrow = length(SIMEXcoef))

    dimnames(SIMEXvar) <- list(colnames(SIMEXcoef), colnames(SIMEXcoef))

  } else {stop("Only support MEvariable.type = 'continuous or 'categorical'")}


  if (identical(class(reg), "lm")) out <- list(NAIVEreg = reg, SIMEXcoef = SIMEXcoef,
                                               SIMEXvar = SIMEXvar, SIMEXsigma = SIMEXsigma)

  if (!identical(class(reg), "lm")) out <- list(NAIVEreg = reg, SIMEXcoef = SIMEXcoef,
                                                SIMEXvar = SIMEXvar)

  class(out) <- "SIMEXreg"

  return(out)

}

predict.SIMEXreg <- function(SIMEXreg, newdata, ...){

  if(missing(newdata)) {stop("newdata is missing")}

  reg <- SIMEXreg$NAIVEreg

  if (identical(class(reg), c("multinom", "nnet"))) {

    if(length(reg$lev) == 2) {

      coef_index <- 1+(1:length(reg$vcoefnames))

    } else {

      coef_index <- as.vector(t(matrix(1:length(reg$wts), nrow = reg$n[3],
                                       byrow=TRUE)[, 1+(1:length(reg$vcoefnames)),
                                                   drop=FALSE][-1, , drop=FALSE]))

    }

    reg$wts[coef_index] <- SIMEXreg$SIMEXcoef

  } else if (identical(class(reg), "polr")) {

    reg$coefficients <- SIMEXreg$SIMEXcoef[1:length(coef(reg))]

    reg$zeta <- SIMEXreg$SIMEXcoef[(length(coef(reg)) + 1):length(SIMEXreg$SIMEXcoef)]

  } else reg$coefficients <- SIMEXreg$SIMEXcoef

  pred <- predict(reg, newdata = newdata, ...)

  return(pred)

}


