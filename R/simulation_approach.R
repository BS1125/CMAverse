simulation_approach <- function(data, formulas, exposure, mediator, outcome, covariates.pre,
                                covariates.post=NULL, mreg, yreg, model, nsims) {

  EY0m_sim <- EY1m_sim <- EY00_sim <- EY01_sim <- EY10_sim <- EY11_sim <-c()

  for (j in 1:nsims) {

    data_sim <- data[sample(1:nrow(data), replace = TRUE),]

    regressions_sim <- run_regressions(formulas = formulas, mreg = mreg,
                                       yreg = yreg, data = data_sim, model = model)

    if (model %in% c("rb", "msm")) {

      # counterfactual design matrices for the mediator model
      mdesign_a <- as.data.frame(cbind(c(rep(a,nrow(data))),data[,covariates.pre]))

      mdesign_a_star <- as.data.frame(cbind(c(rep(a_star,nrow(data))),data[,covariates.pre]))

      colnames(mdesign_a) <- colnames(mdesign_a_star) <- c(exposure, covariates.pre)

      mediator_regression <- regressions_sim$mediator_regression

      # simulate counterfactual mediators
      m_a <- do.call(cbind, lapply(1:length(mediator_regression),
                                   FUN = function(x) predict(mediator_regression[[x]],
                                                             newdata = mdesign_a, type = "response")))

      m_a_star <- do.call(cbind, lapply(1:length(mediator_regression),
                                        FUN = function(x) predict(mediator_regression[[x]],
                                                                  newdata = mdesign_a_star, type = "response")))

      # counterfactual design matrices for the outcome model
      ydesign0m <- as.data.frame(cbind(rep(a_star, nrow(data)), t(m_star)[rep(1,nrow(data)),],
                                       data[,covariates.pre]))

      ydesign1m <- as.data.frame(cbind(rep(a, nrow(data)), t(m_star)[rep(1,nrow(data)),],
                                       data[,covariates.pre]))

      ydesign00 <- as.data.frame(cbind(rep(a_star, nrow(data)), m_a_star, data[,covariates.pre]))

      ydesign01 <- as.data.frame(cbind(rep(a_star, nrow(data)), m_a, data[,covariates.pre]))

      ydesign10 <- as.data.frame(cbind(rep(a, nrow(data)), m_a_star, data[,covariates.pre]))

      ydesign11 <- as.data.frame(cbind(rep(a, nrow(data)), m_a, data[,covariates.pre]))

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, mediator, covariates.pre)

    }

    if (model == "g-formula") {

      mediator_regression_sim <- regressions_sim$mediator_regression

      post_covar_regression_sim <- regressions_sim$post_covar_regression

      mediator_a <- post_covar_a <- as.data.frame(cbind(c(rep(a,nrow(data))),
                                                        data[,covariates.pre]))

      mediator_a_star <- post_covar_a_star <- as.data.frame(cbind(c(rep(a_star,nrow(data))),
                                                                  data[,covariates.pre]))

      colnames(mediator_a) <- colnames(mediator_a_star) <-
        colnames(post_covar_a) <- colnames(post_covar_a_star) <-
        c(exposure, covariates.pre)

      for (i in 1:length(mediator)) {
        mediator_sim_mid <- predict(mediator_regression_sim[[i]],
                                    newdata = mediator_a, type = "response")
        mediator_a <- cbind(mediator_a, mediator_sim_mid)
        mediator_a_star <- cbind(mediator_a_star, mediator_sim_mid)
        colnames(mediator_a) <- colnames(mediator_a_star) <-
          c(exposure, covariates.pre, mediator[1:i])

      }

      for (i in 1:length(covariates.post)) {
        L_sim_mid <- predict(post_covar_regression_sim[[i]],
                             newdata = post_covar_a, type = "response")
        post_covar_a <- cbind(post_covar_a, L_sim_mid)
        post_covar_a_star <- cbind(post_covar_a_star, L_sim_mid)
        colnames(post_covar_a) <- colnames(post_covar_a_star) <-
          c(exposure, covariates.pre, covariates.post[1:i])

      }

      ydesign0m <- as.data.frame(cbind(post_covar_a_star[, c(exposure, covariates.pre, covariates.post)],
                                       t(m_star)[rep(1,nrow(data)),]))

      ydesign1m <- as.data.frame(cbind(post_covar_a[, c(exposure, covariates.pre, covariates.post)],
                                       t(m_star)[rep(1,nrow(data)),]))

      ydesign00 <- as.data.frame(cbind(post_covar_a_star[, c(exposure, covariates.pre, covariates.post)],
                                       mediator_a_star[, mediator]))

      ydesign01 <- as.data.frame(cbind(post_covar_a_star[, c(exposure, covariates.pre, covariates.post)],
                                       mediator_a[, mediator]))

      ydesign10 <- as.data.frame(cbind(post_covar_a[, c(exposure, covariates.pre, covariates.post)],
                                       mediator_a_star[, mediator]))

      ydesign11 <- as.data.frame(cbind(post_covar_a[, c(exposure, covariates.pre, covariates.post)],
                                       mediator_a[, mediator]))

      colnames(ydesign0m) <- colnames(ydesign1m) <- colnames(ydesign00) <- colnames(ydesign01) <-
        colnames(ydesign10) <- colnames(ydesign11) <- c(exposure, covariates.pre,
                                                        covariates.post, mediator)


    }

    type <- ifelse(yreg != "coxph", "response", "risk")

    y0m <- predict(regressions_sim$outcome_regression, newdata =  ydesign0m, type = type)

    y1m <- predict(regressions_sim$outcome_regression, newdata =  ydesign1m, type = type)

    y00 <- predict(regressions_sim$outcome_regression, newdata =  ydesign00, type = type)

    y01 <- predict(regressions_sim$outcome_regression, newdata =  ydesign01, type = type)

    y10 <- predict(regressions_sim$outcome_regression, newdata =  ydesign10, type = type)

    y11 <- predict(regressions_sim$outcome_regression, newdata =  ydesign11, type = type)

    # calculate causal effect components
    EY0m_sim <- c(EY0m_sim, mean(y0m)) #E(Ya0m*)

    EY1m_sim <- c(EY1m_sim, mean(y1m)) #E(Ya1m*)

    EY00_sim <- c(EY00_sim, mean(y00)) #E(Ya0Ma0)

    EY10_sim <- c(EY10_sim, mean(y10)) #E(Ya1Ma0)

    EY01_sim <- c(EY01_sim, mean(y01)) #E(Ya0Ma1)

    EY11_sim <- c(EY11_sim, mean(y11)) #E(Ya1Ma1)

  }

  # linear outcome
  if (yreg == "linear") {

    cde_sim <- EY1m_sim - EY0m_sim # cde in the 3 way or 4 way decomposition

    pnde_sim <- EY10_sim - EY00_sim # pnde in the 3 way decomposition

    tnie_sim <- EY11_sim - EY10_sim # tnie in the 3 way decomposition

    tnde_sim <- EY11_sim - EY01_sim # tnde in the 3 way decomposition

    pnie_sim <- EY01_sim - EY00_sim # pnie in the 3 way decomposition

    te_sim <- tnie_sim + pnde_sim # te in the 3 way or 4 way decomposition

    pm_sim <- tnie_sim / (pnde_sim + te_sim) # pm in the 3 way decomposition

    intref_sim <- pnde_sim - cde_sim # intref in the 4 way decomposition

    intmed_sim <- tnie_sim - pnie_sim # intmed in the 4 way decomposition

    pie_sim <- pnie_sim # pie in the 4 way decomposition

    cde_prop_sim <- cde_sim/te_sim # proportion cde in the 4 way decomposition

    intref_prop_sim <- intref_sim/te_sim # proportion intref in the 4 way decomposition

    intmed_prop_sim <- intmed_sim/te_sim # proportion intmed in the 4 way decomposition

    pie_prop_sim <- pie_sim/te_sim # proportion pie in the 4 way decomposition

    overall_pm_sim <- (pie_sim + intmed_sim)/te_sim # overall pm in the 4 way decomposition

    overall_int_sim <- (intref_sim + intmed_sim)/te_sim # overall int in the 4 way decomposition

    overall_pe_sim <- (intref_sim + intmed_sim + pie_sim)/te_sim # overall pe in the 4 way decomposition

    # calculate the mean of each simulated effects
    for (effect in c("cde", "pnde", "tnde", "pnie", "tnie", "te", "pm", "intref", "intmed",
                     "pie", "cde_prop", "intref_prop", "intmed_prop", "pie_prop",
                     "overall_pm", "overall_int", "overall_pe")) {

      mid <- mean(get(paste0(effect, "_sim")))

      assign(effect, mid)

    }

    # calculate the se of each simulated effects
    for (effect_se in c("cde_se", "pnde_se", "tnde_se", "pnie_se", "tnie_se", "te_se", "pm_se",
                        "intref_se", "intmed_se", "pie_se", "cde_prop_se", "intref_prop_se",
                        "intmed_prop_se", "pie_prop_se", "overall_pm_se",
                        "overall_int_se", "overall_pe_se")) {

      mid <- sd(get(stringr::str_replace(effect_se, "_se", "_sim")))

      assign(effect_se, mid)

    }

    # 3 way decomposition effects and se's
    decomp3way <-  c(cde = cde, cde_se = cde_se,
                     pnde = pnde, pnde_se = pnde_se,
                     tnde = tnde, tnde_se = tnde_se,
                     pnie = pnie, pnie_se = pnie_se,
                     tnie = tnie, tnie_se = tnie_se,
                     te = te, te_se = te_se,
                     pm = pm, pm_se = pm_se)

    # 4 way decomposition effects and se's
    decomp4way <-  c(cde = cde, cde_se = cde_se, intref = intref, intref_se = intref_se,
                     intmed = intmed, intmed_se = intmed_se, pie = pie, pie_se = pie_se,
                     te = te, te_se = te_se, cde_prop = cde_prop, cde_prop_se = cde_prop_se,
                     intref_prop = intref_prop, intref_prop_se = intref_prop_se,
                     intmed_prop = intmed_prop, intmed_prop_se = intmed_prop_se,
                     pie_prop = pie_prop, pie_prop_se = pie_prop_se,
                     overall_pm = overall_pm, overall_pm_se = overall_pm_se,
                     overall_int = overall_int, overall_int_se = overall_int_se,
                     overall_pe = overall_pe, overall_pe_se = overall_pe_se)

    out <- list(decomp3way = decomp3way, decomp4way = decomp4way,
                sims = data.frame(cde = cde_sim, pnde = pnde_sim, tnde = tnde_sim,
                                  pnie = pnie_sim, tnie = tnie_sim, intref = intref_sim,
                                  intmed = intmed_sim, pie_sim = pnie_sim))

  }

  # nonlinear outcome
  if (yreg != "linear") {

    cde_rr_sim <- EY1m_sim/EY0m_sim # cde in the 3 way decomposition

    pnde_rr_sim <- EY10_sim/EY00_sim # pnde in the 3 way decomposition

    tnie_rr_sim <- EY11_sim/EY10_sim # tnie in the 3 way decomposition

    tnde_rr_sim <- EY11_sim/EY01_sim # tnde in the 3 way decomposition

    pnie_rr_sim <- EY01_sim/EY00_sim # pnie in the 3 way decomposition

    te_rr_sim <- tnie_rr_sim * pnde_rr_sim # te in the 3 way or 4 way decomposition

    pm_sim <- (pnde_rr_sim * (tnie_rr_sim - 1)) / (pnde_rr_sim * tnie_rr_sim - 1) # pm in the 3 way decomposition

    cde_err_sim <- (EY1m_sim-EY0m_sim)/EY00_sim # excess RR due to cde in the 4 way decomposition

    intref_err_sim <- pnde_rr_sim - 1 - cde_err_sim # excess RR due to intref in the 4 way decomposition

    intmed_err_sim <- tnie_rr_sim * pnde_rr_sim - pnde_rr_sim - pnie_rr_sim + 1 # excess RR due to intmed in the 4 way decomposition

    pie_err_sim <- pnie_rr_sim - 1 # excess RR due to pie in the 4 way decomposition

    total_err_sim <- EY11_sim / EY00_sim - 1 # total excess RR in the 4 way decomposition

    cde_err_prop_sim <- cde_err_sim/total_err_sim # proportion cde

    intmed_err_prop_sim <- intmed_err_sim/total_err_sim # proportion intmed

    intref_err_prop_sim <- intref_err_sim/total_err_sim # proportion intref

    pie_err_prop_sim <- pie_err_sim/total_err_sim # proportion pie

    overall_pm_sim <- (pie_err_sim+intmed_err_sim)/total_err_sim # overall pm in the 4 way decomposition

    overall_int_sim <- (intref_err_sim+intmed_err_sim)/total_err_sim # overall int in the 4 way decomposition

    overall_pe_sim <- (intref_err_sim+intmed_err_sim+pie_err_sim)/total_err_sim # overall pe in the 4 way decomposition

    # calculate the mean of each simulated effects
    for (effect in c("cde_rr", "pnde_rr", "tnde_rr", "pnie_rr", "tnie_rr", "te_rr", "pm",
                     "cde_err", "intref_err", "intmed_err", "pie_err", "total_err",
                     "cde_err_prop", "intref_err_prop", "intmed_err_prop", "pie_err_prop",
                     "overall_pm", "overall_int", "overall_pe")) {

      mid <- mean(get(paste0(effect, "_sim")))

      assign(effect, mid)

    }

    # calculate the se of each simulated effects
    for (effect_se in c("cde_rr_se", "pnde_rr_se", "tnde_rr_se", "pnie_rr_se", "tnie_rr_se",
                        "te_rr_se", "pm_se", "cde_err_se", "intref_err_se", "intmed_err_se",
                        "pie_err_se", "total_err_se", "cde_err_prop_se", "intref_err_prop_se",
                        "intmed_err_prop_se", "pie_err_prop_se",
                        "overall_pm_se", "overall_int_se", "overall_pe_se")) {

      mid <- sd(get(stringr::str_replace(effect_se, "_se", "_sim")))

      assign(effect_se, mid)

    }

    # 3 way decomposition effects and se's
    decomp3way <-  c(cde_rr = cde_rr, cde_rr_se = cde_rr_se,
                     pnde_rr = pnde_rr, pnde_rr_se = pnde_rr_se,
                     tnde_rr = tnde_rr, tnde_rr_se = tnde_rr_se,
                     pnie_rr = pnie_rr, pnie_rr_se = pnie_rr_se,
                     tnie_rr = tnie_rr, tnie_rr_se = tnie_rr_se,
                     te_rr = te_rr, te_rr_se = te_rr_se,
                     pm = pm, pm_se = pm_se)

    # 4 way decomposition effects and se's
    decomp4way <-  c(cde_err = cde_err, cde_err_se = cde_err_se, intref_err = intref_err,
                     intref_err_se = intref_err_se, intmed_err = intmed_err,
                     intmed_err_se = intmed_err_se, pie_err = pie_err, pie_err_se = pie_err_se,
                     total_err = total_err, total_err_se = total_err_se,
                     cde_err_prop = cde_err_prop, cde_err_prop_se = cde_err_prop_se,
                     intref_err_prop = intref_err_prop, intref_err_prop_se = intref_err_prop_se,
                     intmed_err_prop = intmed_err_prop, intmed_err_prop_se = intmed_err_prop_se,
                     pie_err_prop = pie_err_prop, pie_err_prop_se = pie_err_prop_se,
                     overall_pm = overall_pm, overall_pm_se = overall_pm_se,
                     overall_int = overall_int, overall_int_se = overall_int_se,
                     overall_pe = overall_pe, overall_pe_se = overall_pe_se)

    out <- list(decomp3way = decomp3way, decomp4way = decomp4way,
                sims = data.frame(cde_rr = cde_rr_sim, pnde_rr = pnde_rr_sim, tnde_rr = tnde_rr_sim,
                                  pnie_rr = pnie_rr_sim, tnie_rr = tnie_rr_sim,
                                  cde_err_sim = cde_err_sim, intref_err_sim = intref_err_sim,
                                  intmed_err_sim = intmed_err_sim, pie_err_sim = pie_err_sim))



    }



  return(out)

}
