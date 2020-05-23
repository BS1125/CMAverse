
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

  out <- format_df(cmest_out = cmest_out)

  class(out) <- c("summary.cmest", "data.frame")

  out

}

print.summary.cmest <- function(summary.cmest, digits = 4) {

  printCoefmat(summary.cmest, digits = digits, has.Pvalue = TRUE)

}

print.cmsens.me <- function(cmsens_out, digits = 4) {

  for (i in 1:length(cmsens_out$cmsens)) {

    cat(paste("\n", names(cmsens_out$cmsens)[[i]],"\n"))
    printCoefmat(cmsens_out$cmsens[[i]], digits = digits, has.Pvalue = TRUE)
    cat("--------------------------------------------------------------\n")

  }

}

plot.cmest <- function(cmest_out) {

  require(ggplot2)

  if (cmest_out$model %in% c("rb", "wb", "msm", "g-formula")) {

    index_out <- 1:6

  } else if (cmest_out$model == "iorw") {

    index_out <- 1:3

  }

  effect_df <- data.frame(Effect = factor(rownames(summary(cmest_out))[index_out]),
                          Point = summary(cmest_out)[index_out, 1],
                          CIlower = summary(cmest_out)[index_out, 3],
                          CIupper = summary(cmest_out)[index_out, 4])

  ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper), width = 0.3,
                    data = effect_df)+
      geom_point(aes(x = Effect, y = Point),
                 colour = "blue", data = effect_df) +
      ylab("Point Estimate and 95% CI")+
      geom_hline(yintercept = 0, color = "red")

}



plot.cmsens.me <- function(cmsens_out) {

  require(ggplot2)

  if (cmsens_out$cmest$model %in% c("rb", "wb", "msm", "g-formula")) {

    index_out <- 1:6

  } else if (cmsens_out$cmest$model == "iorw") {

    index_out <- 1:3

  }

  effect_df <- data.frame(Effect = factor(rownames(summary(cmsens_out$cmest))[index_out]),
                          Point = summary(cmsens_out$cmest)[index_out, 1],
                          CIlower = summary(cmsens_out$cmest)[index_out, 3],
                          CIupper = summary(cmsens_out$cmest)[index_out, 4],
                          ReliabilityRatio = rep(1, length(index_out)))

  if (cmsens_out$ME$MEvariable.type == "continuous") {

    for (i in 1:length(cmsens_out$cmsens)) {

      effect_df <- rbind(effect_df,
                         data.frame(Effect = rownames(cmsens_out$cmsens[[i]])[index_out],
                                    Point = cmsens_out$cmsens[[i]][index_out, 1],
                                    CIlower = cmsens_out$cmsens[[i]][index_out, 3],
                                    CIupper = cmsens_out$cmsens[[i]][index_out, 4],
                                    ReliabilityRatio = factor(round((1 - cmsens_out$ME$measurement.error[i]/
                                                                sd(cmsens_out$cmest$data[, cmsens_out$ME$MEvariable],
                                                                   na.rm = TRUE)), 2))))

    }

    ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper,
                        colour = ReliabilityRatio), width = 0.3,
                    data = effect_df,
                    position = position_dodge2(width=0.5))+
      geom_point(aes(x = Effect, y = Point, colour = ReliabilityRatio),
                 data = effect_df,
                 position = position_dodge2(width=0.3)) +
      ylab("Point Estimate and 95% CI")+
      scale_colour_hue()+
      geom_hline(yintercept = 0, color = "red")+
      theme(legend.position = "bottom")

  } else if (cmsens_out$ME$MEvariable.type == "categorical") {

    for (i in 1:length(cmsens_out$cmsens)) {

      effect_df <- rbind(effect_df,
                         data.frame(Effect = rownames(cmsens_out$cmsens[[i]])[index_out],
                                    Point = cmsens_out$cmsens[[i]][index_out, 1],
                                    CIlower = cmsens_out$cmsens[[i]][index_out, 3],
                                    CIupper = cmsens_out$cmsens[[i]][index_out, 4],
                                    MisspecificationMAtrix = factor(paste0("Matrix", i))))

    }

    ggplot() +
      geom_errorbar(aes(x = Effect, ymin = CIlower, ymax = CIupper,
                        colour = MisspecificationMAtrix), width = 0.3,
                    data = effect_df,
                    position = position_dodge2(width=0.5))+
      geom_point(aes(x = Effect, y = Point, colour = MisspecificationMAtrix),
                 data = effect_df,
                 position = position_dodge2(width=0.3)) +
      ylab("Point Estimate and 95% CI")+
      scale_colour_hue()+
      geom_hline(yintercept = 0, color = "red")+
      theme(legend.position = "bottom")

  }



}

