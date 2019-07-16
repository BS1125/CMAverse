function(s, n) {
  s$z <- s$estimate / (s$std.error / sqrt(n))
  s$pval <- 2 * pt(-abs(s$z), n - 1)
  return(data.frame(z = s$z, pval = s$pval))
}


format_df_boot <- function(boot.out, conf = 0.95) {
  CI_lower <- apply(boot.out$boot_results, 2, function(x) quantile(x, (1-conf)/2))

  CI_upper <- apply(boot.out$boot_results, 2, function(x) quantile(x, conf+(1-conf)/2))

  estimate <- boot.out$estimate

  bias <- apply(boot.out$boot_results, 2, mean) - estimate

  std.error <- apply(boot.out$boot_results, 2, sd)

  d_all <- cbind(estimate, bias, std.error, CI_lower, CI_upper)

  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))

  colnames(d_all) <- c("estimate", "bias", "std.error", label_CI[1], label_CI[2])

  rownames(d_all) <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")

  return(d_all)
}

format_df_delta <- function(delta.out, conf = 0.95, n) {
  d_all <- data.frame(matrix(NA, 7, 4))
  alpha <- (1 - conf) / 2
  z <- qnorm(1 - alpha)
  label_CI <- paste0(round(conf * 100, 2), c("% CIL", "% CIU"))
  colnames(d_all) <- c("estimate", "std.error", label_CI[1], label_CI[2])
  rownames(d_all) <- c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm")
  ##----- cde
  d_all["cde", ] <- c(delta.out$cded, delta.out$se_cded,
                      delta.out$cded - z * delta.out$se_cded,
                      delta.out$cded + z * delta.out$se_cded)
  ##----- pnde
  d_all["pnde", ] <- c(delta.out$pnde, delta.out$se_pnde,
                       delta.out$pnde - z * delta.out$se_pnde,
                       delta.out$pnde + z * delta.out$se_pnde)
  ##----- tnde
  d_all["tnde", ] <- c(delta.out$tnde, delta.out$se_tnde,
                       delta.out$tnde - z * delta.out$se_tnde,
                       delta.out$tnde + z * delta.out$se_tnde)
  ##----- pnie
  d_all["pnie", ] <- c(delta.out$pnie, delta.out$se_pnie,
                       delta.out$pnie - z * delta.out$se_pnie,
                       delta.out$pnie + z * delta.out$se_pnie)
  ##----- tnie
  d_all["tnie", ] <- c(delta.out$tnie, delta.out$se_tnie,
                       delta.out$tnie - z * delta.out$se_tnie,
                       delta.out$tnie + z * delta.out$se_tnie)
  ##----- te
  d_all["te", ] <- c(delta.out$te, delta.out$se_te,
                     delta.out$te - z * delta.out$se_te,
                     delta.out$te + z * delta.out$se_te)
  ##----- pm
  d_all["pm", ] <- c(delta.out$pm, delta.out$se_pm,
                     delta.out$pm - z * delta.out$se_pm,
                     delta.out$pm + z * delta.out$se_pm)
  return(cbind(d_all, add_columns(d_all, n = n)))
}

print_test <- function(r, sas, type = "marginal") {
  r_string <- paste("r", type, sep = "_")
  sas_string <- paste("sas", type, sep = "_")
  cat("\n")
  cat(r_string)
  cat("\n")
  print(r)
  cat("\n")
  cat(sas_string)
  cat("\n")
  print(sas)
  cat("\nDifference\n")
  print(r - sas)
  cat("\nSum difference\n")
  print(sum(r - sas))
}

run_test <- function(cm,
                     filename = "Mbin_int/linear_delta.txt",
                     sas = read.csv("../../inst/sasoutput/Mbin_int_linear_delta.csv")) {
  cm$delta_marginal()
  cm$delta_conditional()
  cm$print_output(type = "full")

  cm$summary_coef_conditional
  r_marginal <- cm$summary_coef_marginal[1:6, c(1:4, 6)]
  r_conditional <- cm$summary_coef_marginal[1:6, c(1:4, 6)]

  sas <- sas[!is.na(sas$Obs), ]
  sas_marginal <- sas[1:6, c(3, 4, 6, 7, 5)]
  sas_conditional <- sas[7:12, c(3, 4, 6, 7, 5)]

  sink(filename)
  print_test(r_marginal, sas_marginal, type = "marginal")
  sink()
}

