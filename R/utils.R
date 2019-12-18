
z_p <- function(s, n, nway, yreg) {

  z <- NULL

  pval <- NULL

  if (nway == 3) {

    if (yreg == "linear") {

      z <- s$estimate / s$std.error

      pval <- 2 * pt(-abs(z), n - 1)

      } else {

      z[1:6] <- (s$estimate[1:6]-1) / s$std.error[1:6]

      z[7] <- s$estimate[7] / s$std.error[7]

      pval <- 2 * pt(-abs(z), n - 1)

      }
    } else if (nway == 4) {

        z <- s$estimate / s$std.error

        pval <- 2 * pt(-abs(z), n - 1)

        }

  return(data.frame(z = z, pval = pval))

}

format_row_boot <- function(boot.out, index = 1, conf = 0.95) {
  ci <- boot::boot.ci(boot.out, index = index, type = "perc", conf = conf)
  d <- cbind(as.data.frame(t(tail(ci$percent[1, ], 2))))
  return(d)
}


format_df_boot <- function(boot.out, conf = 0.95, n, nway, yreg) {

  d_all <- NULL

  for (i in 1:length(boot.out$t0)) {

    d <- format_row_boot(boot.out, i, conf = conf)

    d_all <- rbind(d_all, d)

  }

  estimate <- boot.out$t0

  bias <- apply(boot.out$t, 2, mean) - estimate

  std.error <- apply(boot.out$t, 2, sd)

  d_all <- cbind(estimate, bias, std.error, d_all)

  rownames(d_all) <- names(boot.out$t0)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "bias", "std.error", label_CI[1], label_CI[2])


  return(cbind(d_all, z_p(d_all, n = n, nway = nway, yreg = yreg)))

}

format_df_delta <- function(delta.out, conf = 0.95, n, nway, yreg) {

  d_all <- data.frame(matrix(NA, length(delta.out)/2, 4))

  alpha <- (1 - conf) / 2

  z <- qnorm(1 - alpha)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "std.error", label_CI[1], label_CI[2])

  rownames(d_all) <- names(delta.out)[2*(1:(length(delta.out)/2))-1]

  for (name in rownames(d_all)) {

    d_all[name, ] <- c(unname(delta.out[name]), unname(delta.out[stringr::str_c(name,"se_delta",sep="_")]),
                       unname(delta.out[name]-z*delta.out[stringr::str_c(name,"se_delta",sep="_")]),
                       unname(delta.out[name]+z*delta.out[stringr::str_c(name,"se_delta",sep="_")]))

  }

  return(cbind(d_all, z_p(d_all, n = n, nway = nway, yreg = yreg)))

}

print.delta_out <- function(delta_out, digits = 2) {

  printCoefmat(delta_out, digits = digits, has.Pvalue = TRUE)

}

print.boot_out <- function(boot_out, digits = 2) {

  printCoefmat(boot_out, digits = digits, has.Pvalue = TRUE)

}
