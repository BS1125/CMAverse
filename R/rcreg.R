rcreg <- function (reg = NULL, data = NULL, MEvariable = NULL, measurement.error = NULL) {

  reg <- update(reg, data = data)

  if (!is.numeric(data[, MEvariable])) stop("MEvariable of unsupported data type")

  # the vector of names of all variables in the regression formula
  # the vector of names of all independent variables in the regression formula
  reg_formula <- formula(reg)

  if (is.list(reg_formula)) {

    var_vec <- ind_var <- c()

    for (i in 1:length(reg_formula)) {

      var_vec <- c(var_vec, all.vars(reg_formula[[i]]))

      ind_var <- c(ind_var, all.vars(reg_formula[[i]][[3]]))

    }

  } else {

    var_vec <- all.vars(reg_formula)

    ind_var <- all.vars(reg_formula[[3]])

  }

  # the vector of variables measured with error in the regression formula
  MEvar_index <- which(MEvariable %in% var_vec)

  MEvariable <- MEvariable[MEvar_index]

  n_MEvar <- length(MEvariable)

  if (!all(MEvariable %in% ind_var)) {

    stop("Only supports independent variables in the regression formula measured with error")

  } else if (n_MEvar == 0) {

    warning("No variable measured with error in the regression formula and the originial regression object is returned")

    out <- reg

  } else if (n_MEvar == length(ind_var)) {

    stop("No independent variable measured without error in the regression formula")

  } else {

    if (is.null(data)) stop("Unspecified data")

    if (is.null(measurement.error)) {

      stop("Unspecified measurement.error")

    } else {

      measurement.error <- as.matrix(measurement.error)

      if (dim(measurement.error)[1] != n_MEvar | dim(measurement.error)[2] != n_MEvar) {

        stop("Incorrect dimension of measurement.error")

      }

      if (!all(diag(measurement.error) >= 0)) stop("Invalid measurement.error")

    }

    # measurement error matrix for independent variables in the regression formula measured with error
      measurement.error <- measurement.error[MEvar_index, MEvar_index]

      MEvar_index <- which(ind_var %in% MEvariable)

      data_rc <- model.matrix(as.formula(paste0("~", paste(c(MEvariable,
                                                             ind_var[-MEvar_index]), collapse = "+"))),
                              data = data[, ind_var])[, -1]

      n <- nrow(data_rc)

      # mean values of independent variables in the regression formula measured with error
      mean <- colMeans(data_rc)

      W <- data_rc - rbind(mean)[rep(1, n),]

      # var-cov matrix of variables measured with error with all variables
      covmat <- cov(data_rc[, MEvariable], data_rc)

      # true var-cov matrix of variables measured with error
      covmat[1:n_MEvar,
             1:n_MEvar] <- covmat[1:n_MEvar, 1:n_MEvar] - measurement.error

      # true values of variables measured with error
      EX <- rbind(mean[1:n_MEvar])[rep(1, n),] +
        W %*% t(covmat %*% solve(cov(data_rc)))

      rcdata <- data

      rcdata[, MEvariable] <- EX

      # error-corrected regression object
      rcreg <- update(reg, data = rcdata)

      ME <- list(MEvariable = MEvariable, measurement.error = measurement.error)

      out <- list(rcreg = rcreg,
                  rcdata = rcdata,
                  NAIVEreg = reg,
                  NAIVEdata = data,
                  ME = ME)

      class(out) <- c("rcreg")

    }

  return(out)

}


coef.rcreg <- function(rcreg, ...) {

  return(coef(rcreg$rcreg, ...))

}


vcov.rcreg <- function(rcreg, ...) {

  return(vcov(rcreg$rcreg, ...))

}


sigma.rcreg <- function(rcreg, ...) {

  return(sigma(rcreg$rcreg, ...))

}


formula.tcreg <- function(rcreg, ...) {

  return(formula(rcreg$NAIVEreg, ...))

}


print.rcreg <- function(rcreg, ...) {

  print(rcreg$rcreg, ...)

}


summary.rcreg <- function(rcreg, ...) {

  summary(rcreg$rcreg, ...)

}


predict.rcreg <- function(rcreg, ...) {

  predict(rcreg$rcreg, ...)

}


update.rcreg <- function(rcreg, data, ...) {

  reg <- rcreg$NAIVEreg

  if (missing(data)) data <- rcreg$NAIVEdata

  args <- list(object = object, data = data, ...)

  reg <- do.call(update, args)

  MEvariable <- rcreg$ME$MEvariable

  measurement.error <- rcreg$ME$measurement.error

  rcreg <- rcreg(reg = reg, data = data, MEvariable = MEvariable,
                  measurement.error = measurement.error)

  return(rcreg)

}
