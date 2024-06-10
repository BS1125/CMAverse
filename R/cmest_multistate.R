#' Causal Mediation Analysis with the Multistate Approach
#'
#' \code{cmest_multistate} is used to implement \emph{the multistate approach} by Valeri et al. (2023) 
#' for causal mediation analysis.
#'
#' @param data a data frame (or object coercible by \link{as.data.frame} to a data frame) 
#' containing the variables in the model. 
#' @param outcome variable name of the outcome.
#' @param yevent variable name of the event for the outcome.
#' @param mediator variable name of the mediator.
#' @param mevent variable name of the event for the mediator
#' @param exposure variable name of the exposure.
#' @param EMint a logical value. \code{TRUE} indicates there is 
#' exposure-mediator interaction in \code{yreg}. 
#' @param basec a vector of variable names of the confounders.
#' @param astar the control value of the exposure. 
#' @param a the treatment value of the exposure. 
#' @param basecval (required when \code{estimation} is \code{paramfunc} and \code{EMint} is \code{TRUE}) 
#' a list of values at which each confounder is conditioned on, following the order in \code{basec}. 
#' If \code{NULL}, the mean of each confounder is used.
#' @param ymreg type of multistate survival model to be used. Currently supporting coxph only.
#' @param multistate_seed The seed to be used when generating bootstrap datasets for multistate modeling.
#' @param s The time point(s) beyond which survival probability is interested in multistate modeling.
#' @param bh_method Method for estimating baseline hazards in multistate modeling. Currently supporting "breslow" only.
#' @param mediator_event Event indicator for the mediator in multistate modeling.
#' @param nboot (used when \code{inference} is \code{bootstrap}) the number of bootstraps applied. 
#' Default is 200.
#' 
#' @details
#' \strong{Assumptions of the multistate method}
#' \itemize{
#' \item{\emph{Consistency of potential outcomes:} }{For each i and each t, the survival in a world where we intervene, i.e., setting the time to treatment to a specific value t (via a fixed or stochastic intervention) is the same as the survival in the real world where we observe a time to treatment equal to t.}
#' \item{\emph{There is no unmeasured mediator-outcome confounding:} }{Given \code{exposure} and 
#' \code{basec}, \code{mediator} is independent of \code{outcome}.}
#' \item{\emph{Non-informative censoring of event times:} }{The observed censoring time is conditionally independent of all potential event times.}
#' \item{\emph{Positivity:} }{Each exposure-covariate combination has a non-zero probability of occurring.}
#' }
#' 
#' @references
#' Valeri L, Proust-Lima C, Fan W, Chen JT, Jacqmin-Gadda H. A multistate approach for the study of interventions on an intermediate time-to-event in health disparities research. Statistical Methods in Medical Research. 2023;32(8):1445-1460.
#' 
#' @examples
#' \dontrun{
#' library(CMAverse)
#' multistate_out = cmest_multistate(data = sc_data, 
#' s = s_vec,
#' multistate_seed = 1,
#' exposure = 'A', mediator = 'M', outcome = 'S',
#' yevent = "ind_S", mevent = "ind_M",
#' basec = c("C1", "C2"),
#' basecval = c("C1" = "1", "C2" = as.character(mean(sc_data$C2))),
#' astar="0", a="1", 
#' nboot=1, EMint=F, 
#' bh_method = "breslow")
#' multistate_out
#' }
#' where in the example, \code{sc_data} is a data frame that contains a binary exposure \code{A} that has
#' the reference level \code{astar=0} and active level \code{a=1}, a time-to event \code{mediator} with event indicator \code{ind_M},
#' a time_to_event \code{outcome} with event indicator \code{ind_S}, a factor baseline covariate \code{C1} conditioned on level 1, 
#' a and continuous baseline covariate \code{C2} conditioned on its mean. 
#' 
#' @return 
#' The output is a list that consists of 4 elements:
#' \itemize{
#' \item{the model summary of the joint multistate Cox proportional hazards model fitted on the
#' original dataset}
#' \item{the point estimates of RD and SD for each of the
#' user-specified time points of interest on the original dataset}
#' \item{the summary of the bootstrapped RD, SD, and TE estimates for each of the
#' user-specified time point of interest, including the 2.5, 50, and
#' 97.5th percentiles}
#' \item{the estimated RD, SD, TD for each of the
#' user-specified time point of interest for each bootstrap dataset}}
#' 
#' @importFrom mstate msprep expand.covs msfit
#' @importFrom stats as.formula getCall
#' @importFrom survival coxph
#' @importFrom dplyr arrange filter pull group_by summarize
#' @importFrom utils txtProgressBar
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' 
#' @export

cmest_multistate <- function(data = NULL, outcome = NULL, yevent = NULL, 
                     mediator = NULL, mevent = NULL, 
                     exposure = NULL, EMint = NULL, 
                     basec = NULL, basecval = NULL,
                     ymreg = "coxph",
                     astar = NULL, a = NULL, 
                     nboot = 200, 
                     bh_method = "breslow",
                     s = NULL, 
                     multistate_seed = 123) {
  # function call
  cl <- match.call()
  # output list
  out <- list(call = cl)
  
  ###################################################################################################
  #################################Argument Restrictions And Warnings################################
  ###################################################################################################
  # data
  if (is.null(data)) stop("Unspecified data")
  data <- as.data.frame(data)
  if (!all(c(outcome, yevent, mediator, mevent, exposure, basec) %in% colnames(data))) stop("Variables specified in outcome, yevent, mediator, mevent, exposure and basec not found in data")
  #if (sum(is.na(data[, c(outcome, yevent, mediator, mevent, exposure, basec)])) > 0 && !multimp) stop("NAs in outcome, yevent, mediator, mevent, exposure or basec data; delete rows with NAs in these variables from the data or set multimp = TRUE") 
  if (sum(is.na(data[, c(outcome, yevent, mediator, mevent, exposure, basec)])) > 0) stop("NAs in outcome, yevent, mediator, mevent, exposure or basec data; delete rows with NAs in these variables from the data.") 
  out$data <- data
  n <- nrow(data)
  
  # nboot
  if (!is.numeric(nboot)) stop("nboot should be numeric")
  if (nboot <= 0) stop("nboot must be greater than 0") 
  #if (!boot.ci.type %in% c("per", "bca")) stop("Select boot.ci.type from 'per', 'bca'")
  out$methods$nboot <- nboot
  #out$methods$boot.ci.type <- boot.ci.type
  
  # outcome
  if (length(outcome) == 0) stop("Unspecified outcome")
  if (length(outcome) > 1) stop("length(outcome) > 1")
  if (length(yevent) == 0) stop("Unspecified outcome event")
  out$variables$outcome <- outcome
  out$variables$yevent <- yevent
  
  # exposure
  if (length(exposure) == 0) stop("Unspecified exposure")
  if (length(exposure) > 1) stop("length(exposure) > 1")
  out$variables$exposure <- exposure
  
  # mediator
  if (length(mediator) == 0) stop("Unspecified mediator")
  if (length(mevent) == 0) stop("Unspecified mediator event")
  out$variables$mediator <- mediator
  out$variables$mevent <- mevent
  
  # EMint
  if (!is.logical(EMint)) stop("EMint must be TRUE or FALSE")
  out$variables$EMint <- EMint
  
  # s
  if(is.null(s)) stop("Unspecified s")
  if(!is.numeric(s)) stop("Specified s must be numeric")
  
  # basec
  if (!is.null(basec)) out$variables$basec <- basec
  
  # regs
  out$reg.input <- list(ymreg = ymreg)
  
  # a, astar
  if (is.null(a) | is.null(astar)) stop("Unspecified a or astar")
  
  # basecval
  if(length(basec) != 0){
    if(!all(names(basecval) %in% basec) == T){
      stop("basec and basecval variable names mismatch.")
    }else if(length(basec) != length(names(basecval))){
      stop("There must be the same number of variables in basec and basecval. Please double check and rerun.")
    }
    for (i in length(basec)){
      curr_basec = basec[i] # current basec column name
      if (!is.numeric(data[,curr_basec]) && !is.factor(data[,curr_basec])){
        stop("Baseline covariates must be of type numeric or factor.")
      }
    }
  }
  
  # multimp
  #out$multimp <- list(multimp = multimp)
  #if (!is.logical(multimp)) stop("multimp must be TRUE or FALSE")
  # args_mice
  #if (multimp) {
  #  if (!is.null(args_mice) && !is.list(args_mice)) stop("args_mice must be a list")
  #  if (is.null(args_mice)) args_mice <- list()
  #  args_mice$print <- FALSE
  #  if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
  #  args_mice$data <- data
  #  out$multimp$args_mice <- args_mice
  #}
  
  # run multistate method
  set.seed(multistate_seed)
  data = data.frame(data)
  s_grid = s
  ## set up transition matrix
  trans = transMat(x=list(c(2, 3), c(3), c()), names=c(exposure, mediator, outcome)) 
  ## transition-dependent covariates
  covs_df = c(exposure, mediator, basec) 
  # create necessary objects
  ## create a list of mstate data that corresponds to the bootstrap samples
  boot_ind_df = make_boot_ind(data=data, nboot=nboot)
  boot_ind_list = boot_ind_df$boot_ind
  mstate_bootlist <- lapply(boot_ind_list, function(ind){
    boot_dat = data[ind,]
    mstate_boot_dat <- make_mstate_dat(dat=boot_dat, mediator, outcome, mevent, yevent, trans, covs_df)
    return(mstate_boot_dat)
  })
  ## create the formula object for joint mstate modeling
  mstate_form = mstate_formula(data=data, exposure=exposure, mediator=mediator, basec=basec, EMint=EMint)
  mstate_form <- as.formula(mstate_form)
  ## create fixed newdata dfs for msfit()
  mstate_data_orig = make_mstate_dat(dat=data, mediator=mediator, outcome=outcome, mevent=mevent, yevent=yevent, trans=trans, covs_df=covs_df)
  fixed_newd = fixed_newd(mstate_data=mstate_data_orig, trans=trans, a=a, astar=astar, exposure=exposure, mediator=mediator, basec=basec, basecval=basecval)
  
  # fit the joint multistate model on the original data and save the summary table
  mstate_fit_orig = coxph(mstate_form, data = mstate_data_orig, method = bh_method) 
  mstate_fit_orig_summ = summary(mstate_fit_orig)
  
  # compute the point estimate of RD and SD and TD on original data
  pt_est = s_point_est(i=0, mstate_bootlist, mstate_data_orig,
                       s_grid, newd_A0M0=fixed_newd[[1]], newd_A1M0=fixed_newd[[2]], mstate_form,
                       a, astar, exposure, mediator,
                       trans, bh_method)
  
  # PARALLEL computing
  ## nboot grid
  cat("\n")
  cat("Started bootstrapping...")
  i_grid = seq(1, nboot, 1)
  ## run in parallel
  no_cores <- parallel::detectCores() 
  cl <- parallel::makeCluster(no_cores-1)
  #registerDoParallel(cl)
  doSNOW::registerDoSNOW(cl)
  
  pb <- txtProgressBar(max = nboot, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  system.time({
    i_grid = seq(1, nboot, 1)
    results <- foreach::foreach(index = i_grid,
                                .options.snow = opts,
                                .verbose=F,
                                .combine = rbind,
                                #.export = c("mstate_formula", "make_mstate_dat", "s_point_est", "dynamic_newd", "fixed_newd", "make_boot_ind"),
                                .export = c("mstate_formula", "make_mstate_dat", "s_point_est", "dynamic_newd", "make_boot_ind"),
                                .packages = c("mstate", "tidyverse", "parallel")) %dopar% {
                                  s_point_est(i=index, mstate_bootlist, mstate_data_orig,
                                              s_grid, newd_A0M0=fixed_newd[[1]], newd_A1M0=fixed_newd[[2]], mstate_form,
                                              a, astar, exposure, mediator,
                                              trans, bh_method)
                                }
  })
  close(pb)
  stopCluster(cl)
  
  # clean up the data frame
  results_summ = results %>% 
    group_by(s) %>%
    summarize(RD_median = median(RD),
              RD_lower = quantile(RD, 0.025),
              RD_higher = quantile(RD, 0.975),
              SD_median = median(SD),
              SD_lower = quantile(SD, 0.025),
              SD_higher = quantile(SD, 0.975),
              TE_median = median(TE),
              TE_lower = quantile(TE, 0.025),
              TE_higher = quantile(TE, 0.975))
  
  return(list(model_summary = mstate_fit_orig_summ,
              pt_est = pt_est,
              bootstrap_summary = results_summ,
              raw_output = results))
}