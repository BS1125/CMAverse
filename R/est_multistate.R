# updated semicompete script
#library(tidyverse)
#library(dplyr)
#library(mstate)
#library(foreach)
#library(doParallel)
#library(doSNOW)
#library(progressr) ## use progressr for procession updates
#library(doFuture)  ## attaches also foreach and future

# semicompete function code
## function to create the ymreg formula
mstate_formula <- function(data, exposure, mediator, basec, EMint) {
  # create the formula for multistate modeling if not specified
  formula_terms = c(exposure, mediator, basec)
  for (i in 1:length(formula_terms)){
    if (formula_terms[i] == exposure){
      if (class(data[,exposure]) %in% c("numeric", "integer")){
        terms = paste(exposure, c(".1", ".2", ".3"), sep = "", collapse="+")
        # if EMint == T, put the interaction term in the formula as well
        if (EMint){
          terms = paste(terms, 
                        paste(paste(exposure, ".3", sep=""), "*", paste(mediator, ".3", sep=""), sep=""),
                        sep="+")
        }
      }else{ # if the exposure is a factor variable
        nlevel = length(levels(data[,exposure]))
        if (nlevel > 2){
          terms = paste(exposure, rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="", collapse="+")
          exposure_terms = paste(exposure, rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="")
        }else {
          terms = paste(exposure, ".", seq(1,3), sep="", collapse="+")
          exposure_terms = paste(exposure, ".", seq(1,3), sep="")
        }
        
        if (EMint){
          exp_trans3_terms = exposure_terms[startsWith(exposure_terms, exposure) & endsWith(exposure_terms, ".3")]
          int_terms = paste(exp_trans3_terms, "*", paste(mediator, ".3", sep=""), sep="")
          terms = paste(terms, paste(int_terms, collapse="+"), sep="+")
        }
      }
    }else if (formula_terms[i] == mediator){
      terms = paste(terms, paste(mediator, ".3", sep = ""), sep="+")
    }else if (formula_terms[i] %in% basec){
      j = which(basec == formula_terms[i])
      if (class(data[,formula_terms[i]]) %in% c("numeric", "integer")){
        terms = paste(terms, paste(basec[j], c(".1", ".2", ".3"), sep="", collapse="+"), sep="+")
      }else{
        nlevel = length(levels(data[,basec[j]]))
        if (nlevel > 2){
          terms = paste(terms, paste(basec[j], rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="", collapse="+"), sep="+")
        }else {
          terms = paste(terms, paste(basec[j], ".", seq(1,3), sep="", collapse="+"), sep="+")
        }
      }
    }
  }
  ## create the formula for multistate modeling
  mstate_formula = as.formula(paste("Surv(Tstart, Tstop, status) ~ ", terms, "+ strata(trans)"))
  #mstate_formula = as.formula(paste("Surv(Tstart, Tstop, status) ~ ", terms))
  return(mstate_formula)
}

## function to convert the original survival data to mstate format
make_mstate_dat <- function(dat, mediator, outcome, mevent, yevent, trans, covs_df) {
  ## convert the bootstrap data to multistate format
  mstate_dat <- msprep(time = c(NA, mediator, outcome), status = c(NA, mevent, yevent),
                       data = dat, trans = trans, keep = covs_df)
  #print(str(mstate_data))  # Add this line for debugging
  mstate_dat <- expand.covs(mstate_dat, covs_df, append = TRUE, longnames = FALSE)
  return(mstate_dat)
}

## function to make data frames for msfit
fixed_newd <- function(mstate_data, trans, a, astar, exposure, mediator, basec, basecval) {
  # convert to data frame
  mstate_data = data.frame(mstate_data)
  # create individual components
  ## exposure block
  #print(class(mstate_data[,exposure]))
  if (class(mstate_data[,exposure]) %in% c("numeric", "integer")){
    A_blk = matrix(rep(as.numeric(astar),9), nrow=3)
    colnames(A_blk) = paste(exposure, ".", seq(1,3), sep="")
  } else { # if the exposure is a factor
    nlevel = length(levels(mstate_data[,exposure])) # extract the number of levels (Including ref) in the factor variable
    level_selected = which(levels(mstate_data[,exposure]) == astar) # a = astar (ie, a inactive
    #print(paste("nlevel is ", nlevel))
    A_blk = matrix(rep(0,3*3*(nlevel-1)), nrow=3) # (nlevel-1) dummy variables created in mstate formatted data
    if (level_selected > 1){ # level selected is not reference level
      start = 3 * level_selected - 5
      end = start + 2
      A_blk[,start:end] = diag(1,3)
    } 
    if (nlevel > 2){
      colnames(A_blk) = paste(exposure, rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="")
    }else {
      colnames(A_blk) = paste(exposure, ".", seq(1,3), sep="")
    }
  }
  ## mediator block
  M3_blk = matrix(rep(0,3), nrow=3) # time to mediator must be numeric
  colnames(M3_blk) = paste(mediator, ".", 3, sep="")
  ## covariate block
  C_blocks = list()
  if (length(basec) > 0) {
    for (i in 1:length(basec)){
      if (class(mstate_data[,basec[i]]) %in% c("numeric", "integer")){
        C_blocks[[i]] = diag(x=as.numeric(basecval[i]), 3)
        colnames(C_blocks[[i]]) = paste(basec[i], ".", seq(1,3), sep="")
      } else { # if variable type is factor
        nlevel = length(levels(mstate_data[,basec[i]])) # extract the number of levels (Including ref) in the factor variable
        level_selected = which(levels(mstate_data[,basec[i]]) == basecval[i])
        #print(level_selected)
        temp_mat = matrix(rep(0, 3*3*(nlevel-1)), nrow=3)
        if (nlevel > 2){
          if (level_selected > 1){
            start = 3*level_selected-5
            end = start + 2
            temp_mat[,start:end] = diag(1, 3)
          } 
          C_blocks[[i]] = temp_mat
          colnames(C_blocks[[i]]) = paste(basec[i], rep(seq(1,(nlevel-1)), each=3), ".", seq(1,3), sep="")
        } else {
          if (level_selected > 1){
            temp_mat = diag(1,3)
          }
          C_blocks[[i]] = temp_mat
          colnames(C_blocks[[i]]) = paste(basec[i], ".", seq(1,3), sep="")
        }
      }
    }
    C_blk = do.call(cbind, C_blocks)
  } 
  
  ## combine components
  if (length(basec) > 0) {
    newd_A0M0 = data.frame(cbind(A_blk, M3_blk, C_blk))
  } else {
    newd_A0M0 = data.frame(cbind(A_blk, M3_blk))
  }
  newd_A0M0$trans = c(1,2,3)
  newd_A0M0$strata = c(1,2,3)
  attr(newd_A0M0,"trans") <- trans
  class(newd_A0M0) <- c("msdata", "data.frame")
  
  # create newd010, newd100
  newd_A1M0 = newd_A0M0
  if (class(mstate_data[,exposure]) %in% c("numeric", "integer")){
    newd_A1M0[1,paste(exposure,".",1, sep="")] = as.numeric(a)
    newd_A1M0[2,paste(exposure,".",2, sep="")] = as.numeric(a)
    newd_A1M0[3,paste(exposure,".",3, sep="")] = as.numeric(a)
  } else{
    nlevel = length(levels(mstate_data[,exposure]))
    level_active = which(levels(mstate_data[,exposure]) == a) 
    if (level_active > 1){ # if active value is not the reference level
      if (nlevel > 2){
        newd_A1M0[1,paste(exposure, (level_active-1), ".", 1, sep="")] = 1
        newd_A1M0[2,paste(exposure, (level_active-1), ".", 2, sep="")] = 1
        newd_A1M0[3,paste(exposure, (level_active-1), ".", 3, sep="")] = 1
      }else {
        newd_A1M0[1,paste(exposure,".",1, sep="")] = 1
        newd_A1M0[2,paste(exposure,".",2, sep="")] = 1
        newd_A1M0[3,paste(exposure,".",3, sep="")] = 1
      }
    } 
  }
  return(list(newd_A0M0=newd_A0M0, newd_A1M0=newd_A1M0))
}

### dynamic newd (newd001): slightly different for each bootstrap sample
#### in cmest(), run this once for the biggest s. For all smaller s, only need to slice
dynamic_newd = function(mstate_data, exposure, mediator, newd_M0, time_vec, max_s, trans){
  # returns a list of newd data frames for msfit
  up_to_ind = which.max(abs(time_vec - max_s) == min(abs(time_vec - max_s)))
  #up_to_ind = min(which(time_vec >= max_s))
  #up_to_ind = max(which(time_vec <= max_s))
  ## set up newdata as needed (transition 3 now has M)
  newd_Mt_list = lapply(time_vec[1:up_to_ind], function(time){
    newd_Mt = newd_M0
    newd_Mt[3,paste(mediator, ".", 3, sep="")] = time
    attr(newd_Mt, "trans") <- trans
    class(newd_Mt) <- c("msdata", "data.frame")
    return(newd_Mt)
  })
  return(newd_Mt_list)
}

## function to make bootstrap indices
make_boot_ind = function(data, nboot){
  ## create nboot bootstrapping indices
  boot_ind = lapply(1:nboot, function(i){
    ind = sample(1:nrow(data), nrow(data), replace=T)
    data.frame(boot_iter = i, boot_ind = I(list(ind)))
  })
  boot_ind_df = do.call(rbind, boot_ind)
  return(boot_ind_df)
}

# updated function calculates the RD on all s elements in the s_grid for each bootstrap sample i 
## passed in as an argument for boot()
s_point_est <- function(i, mstate_bootlist, mstate_orig,
                        s_grid, newd_A0M0, newd_A1M0, mstate_form,
                        a, astar,
                        exposure, mediator,
                        #exposure, mediator, outcome, basec, mediator_event, event, covs_df, a
                        trans, bh_method, ymreg="coxph"){
  
  if (i > 0){
    mstate_df <- mstate_bootlist[[i]]
  } else {
    mstate_df = mstate_orig
  }
  
  if (ymreg=="coxph"){
    joint_mod <- coxph(mstate_form, data = mstate_df, method = bh_method)
    joint_mod_call <- getCall(joint_mod)
    joint_mod_call$data <- mstate_df
    joint_mod_call$formula <- mstate_form
    joint_mod_call$method <- bh_method
    joint_mod2 <- eval.parent(joint_mod_call)
  }
  cat("Fitted model.\n")
  
  # use msfit() to get predicted cumulative hazards data frames
  cumhaz_A0M0_msfit = msfit(joint_mod2, newd_A0M0, trans=trans) # for RD and SD
  cumhaz_A1M0_msfit = msfit(joint_mod2, newd_A1M0, trans=trans) # for RD and SD
  # extract cumulative hazards from the msfit objects
  cumhaz_A0M0 = cumhaz_A0M0_msfit$Haz
  cumhaz_A1M0 = cumhaz_A1M0_msfit$Haz
  # extract transitions that are needed
  ## time var the same for all data frames below; all below dfs have the same # of rows
  cumhaz_A0M0_trans1 = subset(cumhaz_A0M0, trans==1)
  cumhaz_A0M0_trans2 = subset(cumhaz_A0M0, trans==2)
  cumhaz_A1M0_trans1 = subset(cumhaz_A1M0, trans==1)
  cumhaz_A1M0_trans2 = subset(cumhaz_A1M0, trans==2)
  # populate the hazard grid
  time_vec = cumhaz_A0M0_trans1$time
  #print(paste("The range of time vec is ", range(time_vec)))
  n_grid = nrow(cumhaz_A0M0_trans1)
  haz_A0M0_trans1 = rep(NA, n_grid) # grids for hazard \alpha_{01}
  haz_A1M0_trans1 = rep(NA, n_grid)
  for (i in 1:n_grid){
    haz_A0M0_trans1[i]=(cumhaz_A0M0_trans1$Haz[i+1]-cumhaz_A0M0_trans1$Haz[i])/(time_vec[i+1]-time_vec[i])
    haz_A1M0_trans1[i]=(cumhaz_A1M0_trans1$Haz[i+1]-cumhaz_A1M0_trans1$Haz[i])/(time_vec[i+1]-time_vec[i])
  }
  
  # get the index to slice the newd_Mt_list list
  max_s = max(s_grid)
  newdA1Mt_list = dynamic_newd(mstate_df, exposure, mediator, newd_A1M0, time_vec, max_s, trans) # active - a=1
  newdA0Mt_list = dynamic_newd(mstate_df, exposure, mediator, newd_A0M0, time_vec, max_s, trans) # reference - a=0=astar
  # numerical integration
  ## use lapply to compute each integrand, corresponding to max_s
  ## for integral list for lower values of s, slice the list
  # create the msfit object, extract cumulative hazards data frame
  pb1 <- txtProgressBar(min = 1, max = length(newdA1Mt_list), style = 3)
  cat("\n")
  cat("Preparing for numerical integration (1)...\n")
  newdA1Mt_list_cumhaz = lapply(1:length(newdA1Mt_list), function(i){
    Sys.sleep(0.1)
    setTxtProgressBar(pb1, i)
    newdA1Mt_msfit = msfit(joint_mod2, newdA1Mt_list[[i]], trans=trans)
    cumhaz_A1Mt = newdA1Mt_msfit$Haz
    cumhaz_A1Mt_trans3 = subset(cumhaz_A1Mt, trans==3) 
  })
  pb2 <- txtProgressBar(min = 1, max = length(newdA0Mt_list), style = 3)
  cat("\n")
  cat("Preparing for numerical integration (2)...\n")
  newdA0Mt_list_cumhaz = lapply(1:length(newdA0Mt_list), function(i){
    Sys.sleep(0.1)
    setTxtProgressBar(pb2, i)
    newdA0Mt_msfit = msfit(joint_mod2, newdA0Mt_list[[i]], trans=trans)
    cumhaz_A0Mt = newdA0Mt_msfit$Haz
    cumhaz_A0Mt_trans3 = subset(cumhaz_A0Mt, trans==3)
  })
  
  cat("\n")
  cat("Start creating integrand list.\n")
  integrand_list = list()
  for (j in 1:length(s_grid)){
    curr_s = s_grid[j]
    up_to_ind = which.max(abs(time_vec - curr_s) == min(abs(time_vec - curr_s)))
    #up_to_ind = max(which(time_vec <= curr_s))
    #print(up_to_ind)
    # lapply(1:length(newd001_trans3_list), function(i)
    #print(paste("Creating integrand list for time point s =", curr_s))
    integrand_list[[j]] = lapply(1:up_to_ind, function(i){
      # compute cumelautive hazards
      cumhaz_A1Mt = newdA1Mt_list_cumhaz[[i]]
      cumhaz_A0Mt = newdA0Mt_list_cumhaz[[i]]
      cumhaz_A1Mt_s = cumhaz_A1Mt %>% arrange(abs(curr_s-time)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=1, C=1)
      cumhaz_A0Mt_s = cumhaz_A0Mt %>% arrange(abs(curr_s-time)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{12}(s|A=0, C=1)
      # get the row indices to evaluate the integrand over
      newdA1Mt = newdA1Mt_list[[i]]
      newdA0Mt = newdA1Mt_list[[i]]
      time_ind = which(time_vec == newdA1Mt[3, paste(mediator, ".", 3, sep="")])
      #time_ind = min(time_ind, up_to_ind)
      if (time_ind <= up_to_ind){
        # For RD
        P_01_integrand = (exp(-cumhaz_A0M0_trans1$Haz[time_ind]-cumhaz_A0M0_trans2$Haz[time_ind]) * haz_A0M0_trans1[time_ind] * exp(-cumhaz_A0Mt_s+cumhaz_A0Mt$Haz[time_ind])) * (min(curr_s,time_vec[time_ind+1])-time_vec[time_ind])
        P_g_01_integrand = (exp(-cumhaz_A0M0_trans1$Haz[time_ind]-cumhaz_A1M0_trans2$Haz[time_ind]) * haz_A0M0_trans1[time_ind] * exp(-cumhaz_A1Mt_s+cumhaz_A1Mt$Haz[time_ind])) * (min(curr_s,time_vec[time_ind+1])-time_vec[time_ind])
        # for SD
        sd_integrand1 = (exp(-cumhaz_A1M0_trans1$Haz[time_ind]-cumhaz_A1M0_trans2$Haz[time_ind]) * haz_A1M0_trans1[time_ind] * exp(-cumhaz_A1Mt_s+cumhaz_A1Mt$Haz[time_ind])) * (min(curr_s,time_vec[time_ind+1])-time_vec[time_ind])
        #time_sum = min(max_s,time_vec[time_ind+1])-time_vec[time_ind]
      }else{
        P_01_integrand = 0
        P_g_01_integrand = 0
        sd_integrand1 = 0
      }
      return(c(P_01 = P_01_integrand, P_g_01 = P_g_01_integrand, sd_integrand1 = sd_integrand1))
    })
  }
  cat("Finished creating integrand list.\n")
  
  ## sum up the individual integrands as the estimate for P_01
  RD_vec = rep(NA, length(s_grid))
  SD_vec = rep(NA, length(s_grid))
  for (i in 1:length(s_grid)){ # loop through each value in s_grid
    curr_s = s_grid[i]
    up_to_ind = which.max(abs(time_vec - curr_s) == min(abs(time_vec - curr_s)))
    sums = colSums(do.call(rbind, integrand_list[[i]][1:up_to_ind]))
    #print(sums)
    P_01 = sums[1]
    P_g_01 = sums[2]
    sd_integral = sums[3]
    # estimate cumulative hazard up to the end time s
    cumhaz_A0M0_trans1_s = cumhaz_A0M0_trans1 %>% arrange(abs(curr_s-time)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{01}(s|A=0,C=1)
    cumhaz_A0M0_trans2_s = cumhaz_A0M0_trans2 %>% arrange(abs(curr_s-time)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=0,C=1)
    cumhaz_A1M0_trans2_s = cumhaz_A1M0_trans2 %>% arrange(abs(curr_s-time)) %>% filter(row_number()==1) %>% pull(Haz) # \Lambda_{02}(s|A=1,C=1)
    cumhaz_A1M0_trans1_s = cumhaz_A1M0_trans1 %>% arrange(abs(curr_s-time)) %>% filter(row_number()==1) %>% pull(Haz)
    # P_00 and P_g_00
    P_00 = exp(-cumhaz_A0M0_trans1_s-cumhaz_A0M0_trans2_s)
    P_g_00 = exp(-cumhaz_A0M0_trans1_s-cumhaz_A1M0_trans2_s)
    P_10 = exp(-cumhaz_A1M0_trans1_s-cumhaz_A1M0_trans2_s)
    #print(paste("P_00 is ", P_00))
    #print(paste("P_g_00 is ", P_g_00))
    #print(paste("P_10 is ", P_10))
    # compute RD
    RD_vec[i] = (P_g_00 + P_g_01) - (P_00 + P_01)
    SD_vec[i] = (P_10 + sd_integral) - (P_g_00 + P_g_01)
  }
  
  out_df = data.frame(RD = RD_vec, SD = SD_vec, TE = RD_vec+SD_vec)
  out_df$s = s_grid
  
  cat("\n")
  cat("Started bootstrapping...")
  
  #out = list()
  #out$out_df = out_df
  #out$cumhaz_A0M0 = cumhaz_A0M0
  #out$cumhaz_A1M0 = cumhaz_A1M0
  
  return(out_df)
  #return(out)
}