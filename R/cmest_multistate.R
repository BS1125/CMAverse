#' Causal Mediation Analysis with the Multistate Approach
#'
#' \code{cmest_multistate} is used to implement \emph{the multistate approach} by Valeri et al. (2023) 
#' for causal mediation analysis.
#'
#' @param data a data frame (or object coercible by \link{as.data.frame} to a data frame) 
#' containing the variables in the model. 
#' @param outcome variable name of the outcome.
#' @param yevent variable name of the event for the outcome.
#' @param mediator a vector of variable name(s) of the mediator(s).
#' @param mevent variable name of the event for the mediator
#' @param exposure variable name of the exposure.
#' @param EMint a logical value. \code{TRUE} indicates there is 
#' exposure-mediator interaction in \code{yreg}. 
#' @param basec a vector of variable names of the confounders.
#' @param yreg outcome regression model. See \code{Details}.
#' @param mreg a list of mediator regression models following the order in \code{mediator}. See \code{Details}.
#' @param astar the control value of the exposure. 
#' @param a the treatment value of the exposure. 
#' @param mval a list of values at which each mediator is controlled to calculate the \code{cde}, following the order in \code{mediator}.
#' @param basecval (required when \code{estimation} is \code{paramfunc} and \code{EMint} is \code{TRUE}) 
#' a list of values at which each confounder is conditioned on, following the order in \code{basec}. 
#' If \code{NULL}, the mean of each confounder is used.
#' @param ymreg type of multistate survival model to be used. Currently supporting coxph only.
#' @param multistate_seed The seed to be used when generating bootstrap datasets for multistate modeling.
#' @param s The time beyond which survival probability is interested in multistate modeling.
#' @param bh_method Method for estimating baseline hazards in multistate modeling. Currently supporting "breslow" only.
#' @param mediator_event Event indicator for the mediator in multistate modeling.
#' @param nboot (used when \code{inference} is \code{bootstrap}) the number of bootstraps applied. 
#' Default is 200.
#' @param boot.ci.type (used when \code{inference} is \code{bootstrap}) the type of bootstrap confidence interval. If \code{per}, percentile bootstrap
#' confidence intervals are estimated; if \code{bca}, bias-corrected and accelerated (BCa) bootstrap 
#' confidence intervals are estimated. Default is \code{per}.
#' @param multimp a logical value. If 
#' \code{TRUE}, conduct multiple imputations using the \link[mice]{mice} function. Default is 
#' \code{FALSE}.
#' @param args_mice a list of additional arguments passed to the \link[mice]{mice} function. See \link[mice]{mice}
#' for details.
#' @param x an object of class \code{cmest}.
#' @param object an object of class \code{cmest}.
#' @param digits minimal number of significant digits. See \link{print.default}.
#' 
#' 
#' 
#' @details
#' # add modeling details/supplementary description of how to set the arguments above
#' 
#' @references
#' # references of papers you use for this code
#' 
#' @examples
#' # give some code examples
#' 
#' @importFrom
#' # add functions you use from other packages here
#' 
#' @export

cmest_multistate <- function(data = NULL, outcome = NULL, yevent = NULL, 
                     mediator = NULL, mevent = NULL, exposure = NULL, EMint = NULL, basec = NULL, 
                     yreg = NULL, mreg = NULL, 
                     astar = NULL, a = NULL, basecval = NULL, 
                     nboot = 200, boot.ci.type = "per", 
                     multimp = FALSE, args_mice = NULL) {
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
  if (sum(is.na(data[, c(outcome, yevent, mediator, mevent, exposure, basec)])) > 0 && !multimp) stop("NAs in outcome, yevent, mediator, mevent, exposure or basec data; delete rows with NAs in these variables from the data or set multimp = TRUE") 
  out$data <- data
  n <- nrow(data)
  
  # nboot
  if (!is.numeric(nboot)) stop("nboot should be numeric")
  if (nboot <= 0) stop("nboot must be greater than 0") 
  if (!boot.ci.type %in% c("per", "bca")) stop("Select boot.ci.type from 'per', 'bca'")
  out$methods$nboot <- nboot
  out$methods$boot.ci.type <- boot.ci.type
  
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
  
  # basec
  if (!is.null(basec)) out$variables$basec <- basec
  
  #regs
  out$reg.input <- list(yreg = yreg, mreg = mreg)
  
  # a, astar
  if (is.null(a) | is.null(astar)) stop("Unspecified a or astar")
  
  # basecval
  
  # multimp
  out$multimp <- list(multimp = multimp)
  if (!is.logical(multimp)) stop("multimp must be TRUE or FALSE")
  # args_mice
  if (multimp) {
    if (!is.null(args_mice) && !is.list(args_mice)) stop("args_mice must be a list")
    if (is.null(args_mice)) args_mice <- list()
    args_mice$print <- FALSE
    if (!is.null(args_mice$data)) warning("args_mice$data is overwritten by data")
    args_mice$data <- data
    out$multimp$args_mice <- args_mice
  }
  
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
    mstate_boot_dat <- make_mstate_dat(dat=boot_dat, mediator, outcome, mediator_event, event, trans, covs_df)
    return(mstate_boot_dat)
  })
  ## create the formula object for joint mstate modeling
  mstate_form = mstate_formula(data=data, exposure=exposure, mediator=mediator, basec=basec, EMint=EMint)
  mstate_form <- as.formula(mstate_form)
  ## create fixed newdata dfs for msfit()
  mstate_data_orig = make_mstate_dat(dat=data, mediator=mediator, outcome=outcome, mediator_event=mediator_event, event=event, trans=trans, covs_df=covs_df)
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