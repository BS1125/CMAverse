context("multistate works correctly")

test_that("The multistate method correctly estimates RD and SD", {
  # generate data
  # set up coefficients 
  # M (trans 1)
  a1 = -1.9 
  a2 = 0.2 
  a3 = 0.5 
  
  # S (trans 2)
  b1 = 1 
  b3 = -0.5
  b4 = 0.3
  
  # M when generating the semi-competing observations in resample() (trans 3)
  # S (resample) (trans 3)
  c1 = 0.55
  c2 = -0.15 
  c3 = -0.1 
  c4 = -0.2 
  
  # generate dataset
  set.seed(8) 
  # build a function to generate time-to-event data
  gen_srv <- function(n, lambda, beta, X){
    X = as.matrix(X)
    beta = as.matrix(beta, ncol=1)
    time = -log(runif(n)) / (lambda * exp(X %*% beta)) # exponential distribution
    return(time)
  }
  
  n <- 1000 
  A = sample(c(0,1),replace=TRUE, size=n, c(0.5,0.5)) #binary exposure
  C1 = sample(c(0,1),replace=TRUE, size=n,c(0.6, 0.4)) #binary confounder
  C2 = rnorm(n, mean = 1, sd = 1) #continuous confounder
  id=c(1:n)
  full = data.frame(id,A,C1,C2)
  M = gen_srv(n=n, lambda = 1.5, beta = c(a1,a2,a3), X=data.frame(A,C1,C2)) #time to event mediator
  S = gen_srv(n=n, lambda = 1, beta = c(b1,b3,b4), X=data.frame(A,C1,C2)) #time to event outcome
  data = data.frame(id = c(1:n), M = M, S = S)
  # indicator for event
  data$ind_M = ifelse(data$M <= data$S, 1, 0)
  data$ind_S = 1
  data <- merge(data,full , by = "id")
  #modify Y distribution
  trans_matrix = mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("A", "M", "S"))
  covs = c("A","M", "C1","C2")
  pre_data = mstate::msprep(time = c(NA, "M", "S"), status = c(NA, "ind_M", "ind_S"),
                            data = data, trans = trans_matrix, keep = covs)
  pre_data = mstate::expand.covs(pre_data, covs, append = TRUE, longnames = FALSE)
  #pre_data$A_M.3 = pre_data$A.3*pre_data$M.3
  # resample for T < S
  data_23 = pre_data[which(pre_data$trans == 3),]
  data_23_tem = data.frame(id = rep(NA,dim(data_23)[1]),
                           new_y = rep(NA,dim(data_23)[1]))
  paste("# to resample is ", nrow(data_23))
  
  for(i in 1:dim(data_23)[1]){
    data_23_tem$id[i] = data_23$id[i]
    repeat {
      time_test = gen_srv(n = 1, 
                          lambda = 0.8,
                          beta = c(as.numeric(c1),
                                   c2,
                                   as.numeric(c3),
                                   as.numeric(c4)), 
                          X = data_23[i, c("A.3", "M.3", "C1.3","C2.3")])
      # exit if the condition is met
      if (time_test > data_23[i,"M.3"]) break
    }
    data_23_tem$new_y[i] = time_test
  }
  data_temp = merge(data, data_23_tem, by = "id", all = T)
  # modify Y and M
  data_temp$S[which(data_temp$ind_M == 1)] = data_temp$new_y[which(data_temp$ind_M == 1)]
  data_temp$M[which(data_temp$ind_M == 0)] = data_temp$S[which(data_temp$ind_M == 0)]
  data_final = data_temp
  data_final$A = as.factor(data_final$A) #generate a factor exposure
  
  sc_data = data_final %>% dplyr::select(id,A,M,S,ind_M,ind_S,C1,C2)
  
  # generate time to censoring C; compare C to S; update event indicator
  time_to_censor = runif(n, 0, 2*max(sc_data$S))
  sc_data$ind_S = ifelse(sc_data$S > time_to_censor, 0, 1)
  sc_data$A = factor(sc_data$A)
  sc_data$C1 = factor(sc_data$C1)
  
  # get point estimates
  covs = c("A","M", "C1","C2")
  sc_data_mstate = make_mstate_dat(dat=sc_data, mediator="M", outcome="S", mediator_event="ind_M", event="ind_S", trans=trans_matrix, covs_df=covs)
  sc_data_mstate_formula = mstate_formula(data=sc_data, exposure="A", mediator="M", basec=c("C1", "C2"), EMint=F)
  sc_data_fixed_newd = fixed_newd(sc_data_mstate, trans=trans_matrix, a="1", astar="0", exposure="A", mediator="M", basec=c("C1", "C2"), basecval=c("C1" = "1", "C2" = as.character(mean(sc_data$C2))))
  sc_data_boot_ind_df = make_boot_ind(data=sc_data, nboot=2)
  sc_data_boot_ind_list = sc_data_boot_ind_df$boot_ind
  sc_data_mstate_bootlist <- lapply(sc_data_boot_ind_list, function(ind){
    boot_dat = sc_data[ind,]
    mstate_boot_dat <- make_mstate_dat(dat=boot_dat, 
                                       mediator="M", outcome="S", 
                                       mediator_event="ind_M", event="ind_S", 
                                       trans=trans_matrix, covs_df=covs)
    return(mstate_boot_dat)
  })
  time_to_predict_sc = c(0.1, 0.5, 1, 2, 3, 4, 5, 6)
  sc_data_pt_est = s_point_est(i=0, mstate_bootlist=sc_data_mstate_bootlist, mstate_orig=sc_data_mstate,
                               s_grid=time_to_predict_sc, 
                               newd_A0M0=sc_data_fixed_newd[[1]], newd_A1M0=sc_data_fixed_newd[[2]],
                               mstate_form=sc_data_mstate_formula,
                               a="1", astar="0", exposure="A", mediator="M",
                               trans=trans_matrix, bh_method="breslow", ymreg="coxph")
  
  # get true values
  get_theo_RD_SD = function(s=1, C1=1, C2=1,
                            a_vec=c(a1, a2, a3),
                            b_vec=c(b1, b3, b4),
                            c_vec=c(c1, c2, c3, c4),
                            lambda_vec=c(1.5,1,0.8)){
    a1 = a_vec[1]
    a2 = a_vec[2]
    a3 = a_vec[3]
    b1 = b_vec[1]
    b3 = b_vec[2]
    b4 = b_vec[3]
    c1 = c_vec[1]
    c2 = c_vec[2]
    c3 = c_vec[3]
    c4 = c_vec[4]
    d = C2
    lambda1 = lambda_vec[1]
    lambda2 = lambda_vec[2]
    lambda3 = lambda_vec[3]
    
    P_00 = exp(-lambda1*exp(a2*C1+a3*d)*s - lambda2*exp(b3*C1+b4*d)*s)
    P_g_00 = exp(-lambda1*exp(a2*C1+a3*d)*s - lambda2*exp(b1+b3*C1+b4*d)*s)
    P_g_01_func = function(t){
      exp(-lambda1*exp(a2*C1+a3*d)*t - lambda2*exp(b1+b3*C1+b4*d)*t) * lambda1*exp(a2*C1+a3*d) * exp(lambda3*exp(c1+c2*t+c3*C1+c4*d) * (t-s)) 
    }
    P_01_func = function(t){
      exp(-lambda1*exp(a2*C1+a3*d)*t - lambda2*exp(b3*C1+b4*d)*t) * lambda1*exp(a2*C1+a3*d) * exp(lambda3*exp(c2*t+c3*C1+c4*d) * (t-s)) 
    }
    
    P_g_01 = integrate(P_g_01_func, lower = 0, upper = s)
    P_01 = integrate(P_01_func, lower=0, upper = s)
    
    # theoretical RD
    theoretical_RD = (P_g_00 + P_g_01$value) - (P_00 + P_01$value)
    
    SD_1 = exp(-lambda1*exp(a1+a2*C1+a3*d)*s - lambda2*exp(b1+b3*C1+b4*d)*s)
    SD_2_func = function(t){
      exp(-lambda1*exp(a1+a2*C1+a3*d)*t - lambda2*exp(b1+b3*C1+b4*d)*t) * lambda1*exp(a1+a2*C1+a3*d) * exp(lambda3*exp(c1+c2*t+c3*C1+c4*d) * (t-s)) 
    }
    SD_2 = integrate(SD_2_func, lower = 0, upper = s)
    # theoretical SD
    theoretical_SD = (SD_1 + SD_2$value) - (P_g_00 + P_g_01$value)
    
    return(c(RD = theoretical_RD, SD = theoretical_SD))
  }
  
  theo_list = list()
  for (i in 1:length(time_to_predict_sc)){
    theo_list[[i]] = get_theo_RD_SD(s=time_to_predict_sc[i], C1=1, C2=1,
                                    a_vec=c(a1, a2, a3),
                                    b_vec=c(b1, b3, b4),
                                    c_vec=c(c1, c2, c3, c4),
                                    lambda_vec=c(1.5,1.0,0.8))
  }
  
  theo_mat = round(do.call(rbind, theo_list), 5)
  rownames(theo_mat) = time_to_predict_sc
  
  # compare estimates with truth under tolerance criteria
  ## either error margin less than 0.02
  ## or abs((estimate-true)/true) < 20%
  cond1 = abs(sc_data_pt_est[,c("RD", "SD")] - theo_mat) <= 0.025
  cond2 = abs((sc_data_pt_est[,c("RD", "SD")] - theo_mat)/theo_mat) <= 0.25
  combined_cond = cond1 | cond2
  expect_true(all(combined_cond))
  
})