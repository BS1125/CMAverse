
z_p <- function(s, n, yreg, model) {

  z <- NULL

  pval <- NULL

  z <- s$Estimate / s$Std.error

  pval <- 2 * pt(-abs(z), n - 1)

  return(data.frame(z = z, pval = pval))

}


format_df <- function(results, conf, n, yreg, model) {

  d_all <- data.frame(matrix(NA, length(results$effect_estimate), 4))

  alpha <- (1 - conf) / 2

  z <- qnorm(1 - alpha)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("Estimate", "Std.error", label_CI[1], label_CI[2])

  rownames(d_all) <- names(results$effect_estimate)

  for (name in rownames(d_all)) {

    d_all[name, ] <- c(unname(results$effect_estimate[name]), unname(results$effect_se[which(rownames(d_all)==name)]),
                       unname(results$effect_estimate[name]-z*results$effect_se[which(rownames(d_all)==name)]),
                       unname(results$effect_estimate[name]+z*results$effect_se[which(rownames(d_all)==name)]))

  }

  return(cbind(d_all, z_p(d_all, n = n, yreg = yreg, model = model)))

}

print.cmest <- function(results, digits = 4) {

  printCoefmat(results, digits = digits, has.Pvalue = TRUE)

}

