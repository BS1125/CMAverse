rm(list=ls())

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

# add sc_data as an internal dataset
usethis::use_data(sc_data, internal = TRUE, overwrite = TRUE)