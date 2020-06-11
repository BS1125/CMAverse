simexreg <- function (reg = NULL, data = NULL, MEvariable = NULL, MEvariable.type = NULL,
                      measurement.error = NULL, lambda = c(0.5, 1, 1.5, 2), B = 100) {

  reg <- update(reg, data = data)

  # the vector of variable names in the regression formula
  reg_formula <- formula(reg)

  if (is.list(reg_formula)) {

    var_vec <- c()

    for (i in 1:length(reg_formula)) {

      var_vec <- c(var_vec, all.vars(reg_formula[[i]]))

    }

  } else var_vec <- all.vars(reg_formula)

  # the index vector of variables measured with error in the regression formula
  MEvar_index <- which(MEvariable %in% var_vec)

  n_MEvar <- length(MEvar_index)

  if (n_MEvar == 0) {

    warning("No variable measured with error in the regression formula and the originial regression object is returned")

    out <- reg

  } else {

    if (is.null(data)) stop("Unspecified data")

    n <- nrow(data)

    if (length(MEvariable) != length(MEvariable.type)) {

      stop("length(MEvariable) != length(MEvariable.type)")

    } else if (length(MEvariable) != length(measurement.error)) {

      stop("length(MEvariable) != length(measurement.error)")

    }

    if (!all(lambda>0)) stop("lambda should be positive")

    lambda <- c(0, lambda)

    MEvariable <- MEvariable[MEvar_index]

    MEvariable.type <- MEvariable.type[MEvar_index]

    measurement.error <- measurement.error[MEvar_index]

    if (identical(class(reg), "polr")) {

      SIMcoef <- c(coef(reg), reg$zeta)

    } else if (identical(class(reg), c("multinom", "nnet"))) {

      SIMcoef <- as.vector(t(coef(reg)))

      names(SIMcoef) <- colnames(vcov(reg))

    } else SIMcoef <- coef(reg)

    SIMvar <- list(vcov(reg))

    ncoef <- length(SIMcoef)

    for (i in 2:length(lambda)) {

      SIMcoef_mid <- c()

      SIMvar_mid <- matrix(rep(0, ncoef^2), nrow = ncoef)

      for (b in 1:B) {

        SIMdata <- data

        for (v in 1:n_MEvar) {

          if (MEvariable.type[v] == "continuous") {

            if (!length(measurement.error[[v]]) == 1) {

              stop("measurement.error should be a single value for a continuous variable measured with error")

            } else if (measurement.error[[v]] <= 0) {

              stop("measurement.error should be greater than 0")

            }

            # randomly generate true MEvariable[v]
            Z <- rnorm(n = n)

            SIMdata[, MEvariable[v]] <- SIMdata[, MEvariable[v]] + (sqrt(lambda[i] * measurement.error[[v]])  *  Z )

          } else if (MEvariable.type[v] == "categorical") {

            if (!is.matrix(measurement.error[[v]])) {

              stop("measurement.error should be a matrix for a categorical variable measured with error")

            }

            SIMdata[, MEvariable[v]] <- as.factor(SIMdata[, MEvariable[v]])

            levels(SIMdata[, MEvariable[v]]) <- 0:(length(levels(SIMdata[, MEvariable[v]])) - 1)

            category <- levels(SIMdata[, MEvariable[v]])

            mcmdecomp <- eigen(measurement.error[[v]])

            SIMmcm <- mcmdecomp$vectors %*% diag(mcmdecomp$values)^lambda[i] %*% solve(mcmdecomp$vectors)

            for (j in 1:length(category)) {

              SIMdata[, MEvariable[v]][which(SIMdata[, MEvariable[v]] == category[j])] <-
                sample(x = category, size = length(which(SIMdata[, MEvariable[v]] == category[j])),
                       prob = SIMmcm[, j], replace = TRUE)
            }

          } else stop("Unsupported MEvariable.type[i]")

        }

        SIMreg <- update(reg, data = SIMdata)

        if (identical(class(SIMreg), "polr")) {
          SIMcoef_mid <- rbind(SIMcoef_mid, c(coef(SIMreg), SIMreg$zeta))
        }  else if (identical(class(SIMreg), c("multinom", "nnet"))) {
          SIMcoef_mid <- rbind(SIMcoef_mid, as.vector(t(coef(SIMreg))))
        } else {SIMcoef_mid <- rbind(SIMcoef_mid, c(coef(SIMreg)))}

        SIMvar_mid <- SIMvar_mid + vcov(SIMreg)/B

      }

      SIMcoef <- rbind(SIMcoef, colMeans(SIMcoef_mid))

      SIMvar[[i]]<- SIMvar_mid - cov(SIMcoef_mid)

    }


    extrap_coef <- lm(SIMcoef ~ lambda + I(lambda^2))

    SIMEXcoef <- predict(extrap_coef, newdata = data.frame(lambda = -1))

    SIMvar.element <- matrix(unlist(SIMvar), nrow = length(lambda) , byrow = TRUE)

    extrap_var <- lm(SIMvar.element ~ lambda + I(lambda^2))

    SIMEXvar <- matrix(predict(extrap_var, newdata = data.frame(lambda = -1)),
                       nrow = length(SIMEXcoef))

    dimnames(SIMEXvar) <- list(colnames(SIMEXcoef), colnames(SIMEXcoef))

    ME <- list(MEvariable = MEvariable, MEvariable.type = MEvariable.type,
               measurement.error = measurement.error)

    if (identical(class(reg), "lm")) {

      reg_fit <- reg

      reg_fit$coefficients <- SIMEXcoef

      SIMEXsigma <- sqrt(sum((predict(reg_fit, newdata = reg$model) - reg$model[, 1])^2)/
                           (n - length(SIMEXcoef))) -
        sum(sqrt(unlist(measurement.error[which(MEvariable.type == "continuous")])))

      out <- list(SIMEXcoef = SIMEXcoef, SIMEXvar = SIMEXvar, SIMEXsigma = SIMEXsigma,
                  NAIVEreg = reg, data = data, ME = ME, lambda = lambda[-1], B = B)

    } else  out <- list(SIMEXcoef = SIMEXcoef, SIMEXvar = SIMEXvar,
                        NAIVEreg = reg, data = data, ME = ME, lambda = lambda[-1], B = B)

    class(out) <- "simexreg"

  }

  return(out)

}


coef.simexreg <- function(simexreg) {

  return(simexreg$SIMEXcoef)

}


vcov.simexreg <- function(simexreg) {

  return(simexreg$SIMEXvar)

}


sigma.simexreg <- function(simexreg) {

  return(simexreg$SIMEXsigma)

}


formula.simexreg <- function(simexreg) {

  return(formula(simexreg$NAIVEreg))

}


print.simexreg <- function(simexreg, ...) {

  cat(paste("Naive regression object: \n"))

  print(simexreg$NAIVEreg, ...)

  cat("\nVariables measured with error:\n")

  print(simexreg$ME$MEvariable)

  cat("\nMeasurement errors:\n")

  print(simexreg$ME$measurement.error)

  cat("\nError-corrected coefficient estimates:\n")

  print(simexreg$SIMEXcoef)

}


summary.simexreg <- function(simexreg, ...) {

  n <- length(simexreg$data)

  est <- as.vector(simexreg$SIMEXcoef)

  p <- length(est)

  df <- n - p

  se <- sqrt(diag(simexreg$SIMEXvar))

  tc <- qt(0.975, df)

  t <- est/se

  pval <- 2 * pt(-abs(t), df)

  out <- data.frame(est = est, se = se, lo = est - tc * se, hi = est + tc * se,
                    t = t, pval = pval)

  colnames(out) <- c("Estimate", "Std.error", "95% CIL", "95% CIU", "t", "Pr(>|t|)")

  rownames(out) <- names(se)

  class(out) <- c("summary.simexreg", "data.frame")

  return(out)

}


print.summary.simexreg <- function(summary.simexreg, digits = 4) {

  printCoefmat(summary.simexreg, digits = digits, has.Pvalue = TRUE)

}


predict.simexreg <- function(simexreg, ...){

  reg <- simexreg$NAIVEreg

  if (identical(class(reg), c("multinom", "nnet"))) {

    if(length(reg$lev) == 2) {

      coef_index <- 1+(1:length(reg$vcoefnames))

    } else {

      coef_index <- as.vector(t(matrix(1:length(reg$wts), nrow = reg$n[3],
                                       byrow=TRUE)[, 1+(1:length(reg$vcoefnames)),
                                                   drop=FALSE][-1, , drop=FALSE]))

    }

    reg$wts[coef_index] <- simexreg$SIMEXcoef

  } else if (identical(class(reg), "polr")) {

    reg$coefficients <- simexreg$SIMEXcoef[1:length(coef(reg))]

    reg$zeta <- simexreg$SIMEXcoef[(length(coef(reg)) + 1):length(simexreg$SIMEXcoef)]

  } else reg$coefficients <- simexreg$SIMEXcoef

  out <- predict(reg, ...)

  return(out)

}


update.simexreg <- function(simexreg, data, ...) {

  reg <- simexreg$NAIVEreg

  if (missing(data)) data <- simexreg$data

  args <- list(object = reg, data = data, ...)

  reg <- do.call(update, args)

  MEvariable <- simexreg$ME$MEvariable

  MEvariable.type <- simexreg$ME$MEvariable.type

  measurement.error <- simexreg$ME$measurement.error

  lambda <- simexreg$lambda

  B <- simexreg$B

  out <- simexreg(reg = reg, data = data, MEvariable = MEvariable,
                  MEvariable.type = MEvariable.type,
                  measurement.error = measurement.error,
                  lambda = lambda, B = B)

  return(out)

}
