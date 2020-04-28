
z_p <- function(cmest_out) {

  n <- nrow(cmest_out$data)

  z <- cmest_out$effect_estimate / cmest_out$effect_se

  pval <- 2 * pt(-abs(z), n - 1)

  return(data.frame(z = z, pval = pval))

}


format_df <- function(cmest_out, conf = 0.95) {

  d_all <- data.frame(matrix(NA, length(cmest_out$effect_estimate), 4))

  alpha <- (1 - conf) / 2

  z <- qnorm(1 - alpha)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("Estimate", "Std.error", label_CI[1], label_CI[2])

  rownames(d_all) <- names(cmest_out$effect_estimate)

  for (name in rownames(d_all)) {

    d_all[name, ] <- c(unname(cmest_out$effect_estimate[name]), unname(cmest_out$effect_se[which(rownames(d_all)==name)]),
                       unname(cmest_out$effect_estimate[name]-z*cmest_out$effect_se[which(rownames(d_all)==name)]),
                       unname(cmest_out$effect_estimate[name]+z*cmest_out$effect_se[which(rownames(d_all)==name)]))

  }

  return(cbind(d_all, z_p(cmest_out)))

}

summary.cmest <- function(cmest_out, digits = 4) {

  format_df(cmest_out = cmest_out)

}

print.summary.cmest <- function(cmest_out, digits = 4) {

  out <- format_df(cmest_out = cmest_out)

  printCoefmat(out, digits = digits, has.Pvalue = TRUE)

}

print.cmsens.me <- function(cmsens_out, digits = 4) {

  for (i in 1:length(cmsens_out)) {

    cat(paste("\n", names(cmsens_out)[i],"\n"))
    printCoefmat(cmsens_out[[i]], digits = digits, has.Pvalue = TRUE)
    cat("--------------------------------------------------------------\n")
  }
}

