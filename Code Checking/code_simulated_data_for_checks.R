rm(list=ls())

library(survival)

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2
beta31 = 0.3
beta41 = 0.1

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3
beta32 = 0.5
beta42 = 0.1

theta0 = -5
theta11 = 0.8
theta12 = 0.6
theta21 = 1.8
theta22 = 1.2
theta23 = 1.5
theta24 = 1
theta31 = 0.2
theta32 = 0.4
theta33 = 0.1
theta34 = 0.3
theta35 = 0.2
theta36 = 0.1
theta41 = 0.1
theta42 = 0.3
theta43 = 0.4
theta44 = 0.2

meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}

n=10000

set.seed(100)

A_cont = rnorm(n,5,2)
A_bin = rbinom(n,1,0.4)
A_cat = sample(x=0:2, size=n, prob = c(0.2, 0.3, 0.5), replace = TRUE)

C_cont = rnorm(n,mean=1,sd=sqrt(1))
C_bin = rbinom(n,1,0.6)

###############################################################################################
#################################################bin A#########################################
###############################################################################################

#generate Bin M, Cont M, and Cat M

#contA
linpred1 = (beta01+beta11*A_cont+beta21*C_cont+beta31*C_bin)
linpred2 = (beta02+beta12*A_cont+beta22*C_cont+beta32*C_bin)

probm1_bin = expit(linpred1)
probm0_cat = 1/(1+exp(linpred1)+exp(linpred2))
probm1_cat = exp(linpred1)/(1+exp(linpred1)+exp(linpred2))
probm2_cat = exp(linpred2)/(1+exp(linpred1)+exp(linpred2))

M_bin_contA = rbinom(n,1,probm1_bin)
M_cont_contA = rnorm(n,linpred1,0.1)
M_cat_contA = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                                   prob=c(probm0_cat[x],
                                                          probm1_cat[x],
                                                          probm2_cat[x])))

#binA
linpred1 = (beta01+beta11*A_bin+beta21*C_cont+beta31*C_bin)
linpred2 = (beta02+beta12*A_bin+beta22*C_cont+beta32*C_bin)

probm1_bin = expit(linpred1)
probm0_cat = 1/(1+exp(linpred1)+exp(linpred2))
probm1_cat = exp(linpred1)/(1+exp(linpred1)+exp(linpred2))
probm2_cat = exp(linpred2)/(1+exp(linpred1)+exp(linpred2))

M_bin_binA = rbinom(n,1,probm1_bin)
M_cont_binA = rnorm(n,linpred1,0.1)
M_cat_binA = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                                  prob=c(probm0_cat[x],
                                                         probm1_cat[x],
                                                         probm2_cat[x])))

#catA
linpred1 = (beta01+beta11*(A_cat==1)+beta21*(A_cat==2)+beta31*C_cont+beta41*C_bin)
linpred2 = (beta02+beta12*(A_cat==1)+beta22*(A_cat==2)+beta32*C_cont+beta42*C_bin)

probm1_bin = expit(linpred1)
probm0_cat = 1/(1+exp(linpred1)+exp(linpred2))
probm1_cat = exp(linpred1)/(1+exp(linpred1)+exp(linpred2))
probm2_cat = exp(linpred2)/(1+exp(linpred1)+exp(linpred2))

M_bin_catA = rbinom(n,1,probm1_bin)
M_cont_catA = rnorm(n,linpred1,0.1)
M_cat_catA = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                                  prob=c(probm0_cat[x],
                                                         probm1_cat[x],
                                                         probm2_cat[x])))

##################################################no int###################################################

#contA

linpred1 = theta0+theta11*A_cont+theta21*M_bin_contA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_cont+theta21*M_cont_contA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_cont+theta21*(M_cat_contA==1)+theta22*(M_cat_contA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binM_contA_noint = rexp(n,exp(-linpred1))
cen_binM_contA<-quantile(Y_surv_binM_contA_noint,probs = 0.75)

Y_surv_contM_contA_noint = rexp(n,exp(-linpred2))
cen_contM_contA<-quantile(Y_surv_contM_contA_noint,probs = 0.75)

Y_surv_catM_contA_noint = rexp(n,exp(-linpred3))
cen_catM_contA<-quantile(Y_surv_catM_contA_noint,probs = 0.75)

Ycen_binM_contA_noint<-pmin(Y_surv_binM_contA_noint,cen_binM_contA)
Ycen_contM_contA_noint<-pmin(Y_surv_contM_contA_noint,cen_contM_contA)
Ycen_catM_contA_noint<-pmin(Y_surv_catM_contA_noint,cen_catM_contA)

cens_binM_contA<-as.numeric(Ycen_binM_contA_noint<cen_binM_contA)
delta_binM_contA<-as.numeric(Y_surv_binM_contA_noint>cen_binM_contA)
cens_contM_contA<-as.numeric(Ycen_contM_contA_noint<cen_contM_contA)
delta_contM_contA<-as.numeric(Y_surv_contM_contA_noint>cen_contM_contA)
cens_catM_contA<-as.numeric(Ycen_catM_contA_noint<cen_catM_contA)
delta_catM_contA<-as.numeric(Y_surv_catM_contA_noint>cen_catM_contA)

contY_binM_contA_noint = rnorm(n,linpred1,0.2)
binY_binM_contA_noint = rbinom(n,1,expit(linpred1))
countY_binM_contA_noint = rpois(n,abs(linpred1))

contY_contM_contA_noint = rnorm(n,linpred2,0.2)
binY_contM_contA_noint = rbinom(n,1,expit(linpred2))
countY_contM_contA_noint = rpois(n,abs(linpred2))

contY_catM_contA_noint = rnorm(n,linpred3,0.2)
binY_catM_contA_noint = rbinom(n,1,expit(linpred3))
countY_catM_contA_noint = rpois(n,abs(linpred3))

#binA

linpred1 = theta0+theta11*A_bin+theta21*M_bin_binA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_bin+theta21*M_cont_binA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_bin+theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binM_binA_noint = rexp(n,exp(-linpred1))
cen_binM_binA<-quantile(Y_surv_binM_binA_noint,probs = 0.75)

Y_surv_contM_binA_noint = rexp(n,exp(-linpred2))
cen_contM_binA<-quantile(Y_surv_contM_binA_noint,probs = 0.75)

Y_surv_catM_binA_noint = rexp(n,exp(-linpred3))
cen_catM_binA<-quantile(Y_surv_catM_binA_noint,probs = 0.75)

Ycen_binM_binA_noint<-pmin(Y_surv_binM_binA_noint,cen_binM_binA)
Ycen_contM_binA_noint<-pmin(Y_surv_contM_binA_noint,cen_contM_binA)
Ycen_catM_binA_noint<-pmin(Y_surv_catM_binA_noint,cen_catM_binA)

cens_binM_binA<-as.numeric(Ycen_binM_binA_noint<cen_binM_binA)
delta_binM_binA<-as.numeric(Y_surv_binM_binA_noint>cen_binM_binA)
cens_contM_binA<-as.numeric(Ycen_contM_binA_noint<cen_contM_binA)
delta_contM_binA<-as.numeric(Y_surv_contM_binA_noint>cen_contM_binA)
cens_catM_binA<-as.numeric(Ycen_catM_binA_noint<cen_catM_binA)
delta_catM_binA<-as.numeric(Y_surv_catM_binA_noint>cen_catM_binA)

contY_binM_binA_noint = rnorm(n,linpred1,0.2)
binY_binM_binA_noint = rbinom(n,1,expit(linpred1))
countY_binM_binA_noint = rpois(n,abs(linpred1))

contY_contM_binA_noint = rnorm(n,linpred2,0.2)
binY_contM_binA_noint = rbinom(n,1,expit(linpred2))
countY_contM_binA_noint = rpois(n,abs(linpred2))

contY_catM_binA_noint = rnorm(n,linpred3,0.2)
binY_catM_binA_noint = rbinom(n,1,expit(linpred3))
countY_catM_binA_noint = rpois(n,abs(linpred3))

#catA

linpred1 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_bin_catA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_cont_catA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binM_catA_noint = rexp(n,exp(-linpred1))
cen_binM_catA<-quantile(Y_surv_binM_catA_noint,probs = 0.75)

Y_surv_contM_catA_noint = rexp(n,exp(-linpred2))
cen_contM_catA<-quantile(Y_surv_contM_catA_noint,probs = 0.75)

Y_surv_catM_catA_noint = rexp(n,exp(-linpred3))
cen_catM_catA<-quantile(Y_surv_catM_catA_noint,probs = 0.75)

Ycen_binM_catA_noint<-pmin(Y_surv_binM_catA_noint,cen_binM_catA)
Ycen_contM_catA_noint<-pmin(Y_surv_contM_catA_noint,cen_contM_catA)
Ycen_catM_catA_noint<-pmin(Y_surv_catM_catA_noint,cen_catM_catA)

cens_binM_catA<-as.numeric(Ycen_binM_catA_noint<cen_binM_catA)
delta_binM_catA<-as.numeric(Y_surv_binM_catA_noint>cen_binM_catA)
cens_contM_catA<-as.numeric(Ycen_contM_catA_noint<cen_contM_catA)
delta_contM_catA<-as.numeric(Y_surv_contM_catA_noint>cen_contM_catA)
cens_catM_catA<-as.numeric(Ycen_catM_catA_noint<cen_catM_catA)
delta_catM_catA<-as.numeric(Y_surv_catM_catA_noint>cen_catM_catA)

contY_binM_catA_noint = rnorm(n,linpred1,0.2)
binY_binM_catA_noint = rbinom(n,1,expit(linpred1))
countY_binM_catA_noint = rpois(n,abs(linpred1))

contY_contM_catA_noint = rnorm(n,linpred2,0.2)
binY_contM_catA_noint = rbinom(n,1,expit(linpred2))
countY_contM_catA_noint = rpois(n,abs(linpred2))

contY_catM_catA_noint = rnorm(n,linpred3,0.2)
binY_catM_catA_noint = rbinom(n,1,expit(linpred3))
countY_catM_catA_noint = rpois(n,abs(linpred3))

noint_data = cbind(A_cont,A_bin,A_cat,
                   M_cont_contA,M_bin_contA,M_cat_contA,
                   M_cont_binA,M_bin_binA,M_cat_binA,
                   M_cont_catA,M_bin_catA,M_cat_catA,
                   contY_binM_contA_noint,contY_contM_contA_noint,contY_catM_contA_noint,
                   binY_binM_contA_noint,binY_contM_contA_noint,binY_catM_contA_noint,
                   countY_binM_contA_noint,countY_contM_contA_noint,countY_catM_contA_noint,
                   contY_binM_binA_noint,contY_contM_binA_noint,contY_catM_binA_noint,
                   binY_binM_binA_noint,binY_contM_binA_noint,binY_catM_binA_noint,
                   countY_binM_binA_noint,countY_contM_binA_noint,countY_catM_binA_noint,
                   contY_binM_catA_noint,contY_contM_catA_noint,contY_catM_catA_noint,
                   binY_binM_catA_noint,binY_contM_catA_noint,binY_catM_catA_noint,
                   countY_binM_catA_noint,countY_contM_catA_noint,countY_catM_catA_noint,
                   Ycen_binM_contA_noint,Ycen_contM_contA_noint,Ycen_catM_contA_noint,
                   cens_binM_contA,cens_contM_contA,cens_catM_contA,
                   delta_binM_contA,delta_contM_contA,delta_catM_contA,
                   Ycen_binM_binA_noint,Ycen_contM_binA_noint,Ycen_catM_binA_noint,
                   cens_binM_binA,cens_contM_binA,cens_catM_binA,
                   delta_binM_binA,delta_contM_binA,delta_catM_binA,
                   Ycen_binM_catA_noint,Ycen_contM_catA_noint,Ycen_catM_catA_noint,
                   cens_binM_catA,cens_contM_catA,cens_catM_catA,
                   delta_binM_catA,delta_contM_catA,delta_catM_catA,
                   C_cont,C_bin)

write.table(noint_data,file="./Data Simulation and Code Checking/noint_data_10000.txt")

######################################################int################################################

#contA

linpred1 = theta0+theta11*A_cont+theta21*M_bin_contA+theta31*A_cont*M_bin_contA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_cont+theta21*M_cont_contA+theta31*A_cont*M_cont_contA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_cont+theta21*(M_cat_contA==1)+theta31*A_cont*(M_cat_contA==1)+
  theta32*A_cont*(M_cat_contA==2)+theta22*(M_cat_contA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binM_contA_int = rexp(n,exp(-linpred1))
cen_binM_contA<-quantile(Y_surv_binM_contA_int,probs = 0.75)

Y_surv_contM_contA_int = rexp(n,exp(-linpred2))
cen_contM_contA<-quantile(Y_surv_contM_contA_int,probs = 0.75)

Y_surv_catM_contA_int = rexp(n,exp(-linpred3))
cen_catM_contA<-quantile(Y_surv_catM_contA_int,probs = 0.75)

Ycen_binM_contA_int<-pmin(Y_surv_binM_contA_int,cen_binM_contA)
Ycen_contM_contA_int<-pmin(Y_surv_contM_contA_int,cen_contM_contA)
Ycen_catM_contA_int<-pmin(Y_surv_catM_contA_int,cen_catM_contA)

cens_binM_contA<-as.numeric(Ycen_binM_contA_int<cen_binM_contA)
delta_binM_contA<-as.numeric(Y_surv_binM_contA_int>cen_binM_contA)
cens_contM_contA<-as.numeric(Ycen_contM_contA_int<cen_contM_contA)
delta_contM_contA<-as.numeric(Y_surv_contM_contA_int>cen_contM_contA)
cens_catM_contA<-as.numeric(Ycen_catM_contA_int<cen_catM_contA)
delta_catM_contA<-as.numeric(Y_surv_catM_contA_int>cen_catM_contA)

contY_binM_contA_int = rnorm(n,linpred1,0.2)
binY_binM_contA_int = rbinom(n,1,expit(linpred1))
countY_binM_contA_int = rpois(n,abs(linpred1))

contY_contM_contA_int = rnorm(n,linpred2,0.2)
binY_contM_contA_int = rbinom(n,1,expit(linpred2))
countY_contM_contA_int = rpois(n,abs(linpred2))

contY_catM_contA_int = rnorm(n,linpred3,0.2)
binY_catM_contA_int = rbinom(n,1,expit(linpred3))
countY_catM_contA_int = rpois(n,abs(linpred3))

#binA

linpred1 = theta0+theta11*A_bin+theta21*M_bin_binA+theta31*A_bin*M_bin_binA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_bin+theta21*M_cont_binA+theta31*A_bin*M_cont_binA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_bin+theta21*(M_cat_binA==1)+theta31*A_bin*(M_cat_binA==1)+
  theta32*A_bin*(M_cat_binA==2)+theta22*(M_cat_binA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binM_binA_int = rexp(n,exp(-linpred1))
cen_binM_binA<-quantile(Y_surv_binM_binA_int,probs = 0.75)

Y_surv_contM_binA_int = rexp(n,exp(-linpred2))
cen_contM_binA<-quantile(Y_surv_contM_binA_int,probs = 0.75)

Y_surv_catM_binA_int = rexp(n,exp(-linpred3))
cen_catM_binA<-quantile(Y_surv_catM_binA_int,probs = 0.75)

Ycen_binM_binA_int<-pmin(Y_surv_binM_binA_int,cen_binM_binA)
Ycen_contM_binA_int<-pmin(Y_surv_contM_binA_int,cen_contM_binA)
Ycen_catM_binA_int<-pmin(Y_surv_catM_binA_int,cen_catM_binA)

cens_binM_binA<-as.numeric(Ycen_binM_binA_int<cen_binM_binA)
delta_binM_binA<-as.numeric(Y_surv_binM_binA_int>cen_binM_binA)
cens_contM_binA<-as.numeric(Ycen_contM_binA_int<cen_contM_binA)
delta_contM_binA<-as.numeric(Y_surv_contM_binA_int>cen_contM_binA)
cens_catM_binA<-as.numeric(Ycen_catM_binA_int<cen_catM_binA)
delta_catM_binA<-as.numeric(Y_surv_catM_binA_int>cen_catM_binA)

contY_binM_binA_int = rnorm(n,linpred1,0.2)
binY_binM_binA_int = rbinom(n,1,expit(linpred1))
countY_binM_binA_int = rpois(n,abs(linpred1))

contY_contM_binA_int = rnorm(n,linpred2,0.2)
binY_contM_binA_int = rbinom(n,1,expit(linpred2))
countY_contM_binA_int = rpois(n,abs(linpred2))

contY_catM_binA_int = rnorm(n,linpred3,0.2)
binY_catM_binA_int = rbinom(n,1,expit(linpred3))
countY_catM_binA_int = rpois(n,abs(linpred3))

#catA

linpred1 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_bin_catA+
  theta31*(A_cat==1)*M_bin_catA+theta32*(A_cat==2)*M_bin_catA+
  theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_cont_catA+
  theta31*(A_cat==1)*M_cont_catA+theta32*(A_cat==2)*M_cont_catA+
  theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*(M_cat_catA==1)+
  theta22*(M_cat_catA==2)+theta31*(A_cat==1)*(M_cat_catA==1)+
  theta32*(A_cat==1)*(M_cat_catA==2)+theta33*(A_cat==2)*(M_cat_catA==1)+
  theta34*(A_cat==2)*(M_cat_catA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binM_catA_int = rexp(n,exp(-linpred1))
cen_binM_catA<-quantile(Y_surv_binM_catA_int,probs = 0.75)

Y_surv_contM_catA_int = rexp(n,exp(-linpred2))
cen_contM_catA<-quantile(Y_surv_contM_catA_int,probs = 0.75)

Y_surv_catM_catA_int = rexp(n,exp(-linpred3))
cen_catM_catA<-quantile(Y_surv_catM_catA_int,probs = 0.75)

Ycen_binM_catA_int<-pmin(Y_surv_binM_catA_int,cen_binM_catA)
Ycen_contM_catA_int<-pmin(Y_surv_contM_catA_int,cen_contM_catA)
Ycen_catM_catA_int<-pmin(Y_surv_catM_catA_int,cen_catM_catA)

cens_binM_catA<-as.numeric(Ycen_binM_catA_int<cen_binM_catA)
delta_binM_catA<-as.numeric(Y_surv_binM_catA_int>cen_binM_catA)
cens_contM_catA<-as.numeric(Ycen_contM_catA_int<cen_contM_catA)
delta_contM_catA<-as.numeric(Y_surv_contM_catA_int>cen_contM_catA)
cens_catM_catA<-as.numeric(Ycen_catM_catA_int<cen_catM_catA)
delta_catM_catA<-as.numeric(Y_surv_catM_catA_int>cen_catM_catA)

contY_binM_catA_int = rnorm(n,linpred1,0.2)
binY_binM_catA_int = rbinom(n,1,expit(linpred1))
countY_binM_catA_int = rpois(n,abs(linpred1))

contY_contM_catA_int = rnorm(n,linpred2,0.2)
binY_contM_catA_int = rbinom(n,1,expit(linpred2))
countY_contM_catA_int = rpois(n,abs(linpred2))

contY_catM_catA_int = rnorm(n,linpred3,0.2)
binY_catM_catA_int = rbinom(n,1,expit(linpred3))
countY_catM_catA_int = rpois(n,abs(linpred3))

int_data = cbind(A_cont,A_bin,A_cat,
                 M_cont_contA,M_bin_contA,M_cat_contA,
                 M_cont_binA,M_bin_binA,M_cat_binA,
                 M_cont_catA,M_bin_catA,M_cat_catA,
                 contY_binM_contA_int,contY_contM_contA_int,contY_catM_contA_int,
                 binY_binM_contA_int,binY_contM_contA_int,binY_catM_contA_int,
                 countY_binM_contA_int,countY_contM_contA_int,countY_catM_contA_int,
                 contY_binM_binA_int,contY_contM_binA_int,contY_catM_binA_int,
                 binY_binM_binA_int,binY_contM_binA_int,binY_catM_binA_int,
                 countY_binM_binA_int,countY_contM_binA_int,countY_catM_binA_int,
                 contY_binM_catA_int,contY_contM_catA_int,contY_catM_catA_int,
                 binY_binM_catA_int,binY_contM_catA_int,binY_catM_catA_int,
                 countY_binM_catA_int,countY_contM_catA_int,countY_catM_catA_int,
                 Ycen_binM_contA_int,Ycen_contM_contA_int,Ycen_catM_contA_int,
                 cens_binM_contA,cens_contM_contA,cens_catM_contA,
                 delta_binM_contA,delta_contM_contA,delta_catM_contA,
                 Ycen_binM_binA_int,Ycen_contM_binA_int,Ycen_catM_binA_int,
                 cens_binM_binA,cens_contM_binA,cens_catM_binA,
                 delta_binM_binA,delta_contM_binA,delta_catM_binA,
                 Ycen_binM_catA_int,Ycen_contM_catA_int,Ycen_catM_catA_int,
                 cens_binM_catA,cens_contM_catA,cens_catM_catA,
                 delta_binM_catA,delta_contM_catA,delta_catM_catA,
                 C_cont,C_bin)

write.table(int_data,file="./Data Simulation and Code Checking/int_data_10000.txt")

##############################################multiple M without int####################################################

# contA

linpred1 = theta0+theta11*A_cont+theta21*M_cont_contA+theta22*M_bin_contA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_cont+theta21*(M_cat_contA==1)+theta22*(M_cat_contA==2)+theta23*M_bin_contA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_cont+theta21*M_cont_contA+theta22*(M_cat_contA==1)+theta23*(M_cat_contA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcontA_noint = rexp(n,exp(-(linpred1)))
cen_binMcontMcontA<-quantile(Y_surv_binMcontMcontA_noint,probs = 0.75)
Y_surv_catMbinMcontA_noint = rexp(n,exp(-(linpred2)))
cen_catMbinMcontA<-quantile(Y_surv_catMbinMcontA_noint,probs = 0.75)
Y_surv_catMcontMcontA_noint = rexp(n,exp(-(linpred3)))
cen_catMcontMcontA<-quantile(Y_surv_catMcontMcontA_noint,probs = 0.75)

Ycen_binMcontMcontA_noint<-pmin(Y_surv_binMcontMcontA_noint,cen_binMcontMcontA)
Ycen_catMbinMcontA_noint<-pmin(Y_surv_catMbinMcontA_noint,cen_catMbinMcontA)
Ycen_catMcontMcontA_noint<-pmin(Y_surv_catMcontMcontA_noint,cen_catMcontMcontA)

cens_binMcontMcontA<-as.numeric(Ycen_binMcontMcontA_noint<cen_binMcontMcontA)
delta_binMcontMcontA<-as.numeric(Y_surv_binMcontMcontA_noint>cen_binMcontMcontA)
cens_catMbinMcontA<-as.numeric(Ycen_catMbinMcontA_noint<cen_catMbinMcontA)
delta_catMbinMcontA<-as.numeric(Y_surv_catMbinMcontA_noint>cen_catMbinMcontA)
cens_catMcontMcontA<-as.numeric(Ycen_catMcontMcontA_noint<cen_catMcontMcontA)
delta_catMcontMcontA<-as.numeric(Y_surv_catMcontMcontA_noint>cen_catMcontMcontA)

contY_binMcontMcontA_noint = rnorm(n,linpred1,0.2)
contY_catMbinMcontA_noint = rnorm(n,linpred2,0.2)
contY_catMcontMcontA_noint = rnorm(n,linpred3,0.2)
binY_binMcontMcontA_noint = rbinom(n,1,expit(linpred1))
binY_catMbinMcontA_noint = rbinom(n,1,expit(linpred2))
binY_catMcontMcontA_noint = rbinom(n,1,expit(linpred3))
countY_binMcontMcontA_noint = rpois(n,abs(linpred1))
countY_catMbinMcontA_noint = rpois(n,abs(linpred2))
countY_catMcontMcontA_noint = rpois(n,abs(linpred3))

# binA

linpred1 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*M_bin_binA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_bin+theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta23*M_bin_binA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*(M_cat_binA==1)+theta23*(M_cat_binA==2)+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMbinA_noint = rexp(n,exp(-(linpred1)))
cen_binMcontMbinA<-quantile(Y_surv_binMcontMbinA_noint,probs = 0.75)
Y_surv_catMbinMbinA_noint = rexp(n,exp(-(linpred2)))
cen_catMbinMbinA<-quantile(Y_surv_catMbinMbinA_noint,probs = 0.75)
Y_surv_catMcontMbinA_noint = rexp(n,exp(-(linpred3)))
cen_catMcontMbinA<-quantile(Y_surv_catMcontMbinA_noint,probs = 0.75)

Ycen_binMcontMbinA_noint<-pmin(Y_surv_binMcontMbinA_noint,cen_binMcontMbinA)
Ycen_catMbinMbinA_noint<-pmin(Y_surv_catMbinMbinA_noint,cen_catMbinMbinA)
Ycen_catMcontMbinA_noint<-pmin(Y_surv_catMcontMbinA_noint,cen_catMcontMbinA)

cens_binMcontMbinA<-as.numeric(Ycen_binMcontMbinA_noint<cen_binMcontMbinA)
delta_binMcontMbinA<-as.numeric(Y_surv_binMcontMbinA_noint>cen_binMcontMbinA)
cens_catMbinMbinA<-as.numeric(Ycen_catMbinMbinA_noint<cen_catMbinMbinA)
delta_catMbinMbinA<-as.numeric(Y_surv_catMbinMbinA_noint>cen_catMbinMbinA)
cens_catMcontMbinA<-as.numeric(Ycen_catMcontMbinA_noint<cen_catMcontMbinA)
delta_catMcontMbinA<-as.numeric(Y_surv_catMcontMbinA_noint>cen_catMcontMbinA)

contY_binMcontMbinA_noint = rnorm(n,linpred1,0.2)
contY_catMbinMbinA_noint = rnorm(n,linpred2,0.2)
contY_catMcontMbinA_noint = rnorm(n,linpred3,0.2)
binY_binMcontMbinA_noint = rbinom(n,1,expit(linpred1))
binY_catMbinMbinA_noint = rbinom(n,1,expit(linpred2))
binY_catMcontMbinA_noint = rbinom(n,1,expit(linpred3))
countY_binMcontMbinA_noint = rpois(n,abs(linpred1))
countY_catMbinMbinA_noint = rpois(n,abs(linpred2))
countY_catMcontMbinA_noint = rpois(n,abs(linpred3))

# catA

linpred1 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+
  theta21*M_cont_catA+theta22*M_bin_catA+
  theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+
  theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta23*M_bin_catA+
  theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+
  theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta23*M_cont_catA+
  theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcatA_noint = rexp(n,exp(-(linpred1)))
cen_binMcontMcatA<-quantile(Y_surv_binMcontMcatA_noint,probs = 0.75)
Y_surv_catMbinMcatA_noint = rexp(n,exp(-(linpred2)))
cen_catMbinMcatA<-quantile(Y_surv_catMbinMcatA_noint,probs = 0.75)
Y_surv_catMcontMcatA_noint = rexp(n,exp(-(linpred3)))
cen_catMcontMcatA<-quantile(Y_surv_catMcontMcatA_noint,probs = 0.75)

Ycen_binMcontMcatA_noint<-pmin(Y_surv_binMcontMcatA_noint,cen_binMcontMcatA)
Ycen_catMbinMcatA_noint<-pmin(Y_surv_catMbinMcatA_noint,cen_catMbinMcatA)
Ycen_catMcontMcatA_noint<-pmin(Y_surv_catMcontMcatA_noint,cen_catMcontMcatA)

cens_binMcontMcatA<-as.numeric(Ycen_binMcontMcatA_noint<cen_binMcontMcatA)
delta_binMcontMcatA<-as.numeric(Y_surv_binMcontMcatA_noint>cen_binMcontMcatA)
cens_catMbinMcatA<-as.numeric(Ycen_catMbinMcatA_noint<cen_catMbinMcatA)
delta_catMbinMcatA<-as.numeric(Y_surv_catMbinMcatA_noint>cen_catMbinMcatA)
cens_catMcontMcatA<-as.numeric(Ycen_catMcontMcatA_noint<cen_catMcontMcatA)
delta_catMcontMcatA<-as.numeric(Y_surv_catMcontMcatA_noint>cen_catMcontMcatA)

contY_binMcontMcatA_noint = rnorm(n,linpred1,0.2)
contY_catMbinMcatA_noint = rnorm(n,linpred2,0.2)
contY_catMcontMcatA_noint = rnorm(n,linpred3,0.2)
binY_binMcontMcatA_noint = rbinom(n,1,expit(linpred1))
binY_catMbinMcatA_noint = rbinom(n,1,expit(linpred2))
binY_catMcontMcatA_noint = rbinom(n,1,expit(linpred3))
countY_binMcontMcatA_noint = rpois(n,abs(linpred1))
countY_catMbinMcatA_noint = rpois(n,abs(linpred2))
countY_catMcontMcatA_noint = rpois(n,abs(linpred3))

multipleM_noint_data = cbind(A_cont,A_bin,A_cat,
                             M_cont_contA,M_bin_contA,M_cat_contA,
                             M_cont_binA,M_bin_binA,M_cat_binA,
                             M_cont_catA,M_bin_catA,M_cat_catA,
                             contY_binMcontMcontA_noint,contY_catMbinMcontA_noint,contY_catMcontMcontA_noint,
                             binY_binMcontMcontA_noint,binY_catMbinMcontA_noint,binY_catMcontMcontA_noint,
                             countY_binMcontMcontA_noint,countY_catMbinMcontA_noint,countY_catMcontMcontA_noint,
                             Ycen_binMcontMcontA_noint,cens_binMcontMcontA,delta_binMcontMcontA,
                             contY_binMcontMbinA_noint,contY_catMbinMbinA_noint,contY_catMcontMbinA_noint,
                             binY_binMcontMbinA_noint,binY_catMbinMbinA_noint,binY_catMcontMbinA_noint,
                             countY_binMcontMbinA_noint,countY_catMbinMbinA_noint,countY_catMcontMbinA_noint,
                             Ycen_binMcontMbinA_noint,cens_binMcontMbinA,delta_binMcontMbinA,
                             contY_binMcontMcatA_noint,contY_catMbinMcatA_noint,contY_catMcontMcatA_noint,
                             binY_binMcontMcatA_noint,binY_catMbinMcatA_noint,binY_catMcontMcatA_noint,
                             countY_binMcontMcatA_noint,countY_catMbinMcatA_noint,countY_catMcontMcatA_noint,
                             Ycen_binMcontMcatA_noint,cens_binMcontMcatA,delta_binMcontMcatA,
                             contY_catMbinMcontA_noint,contY_catMbinMcontA_noint,contY_catMcontMcontA_noint,
                             binY_catMbinMcontA_noint,binY_catMbinMcontA_noint,binY_catMcontMcontA_noint,
                             countY_catMbinMcontA_noint,countY_catMbinMcontA_noint,countY_catMcontMcontA_noint,
                             Ycen_catMbinMcontA_noint,cens_catMbinMcontA,delta_catMbinMcontA,
                             contY_catMbinMbinA_noint,contY_catMbinMbinA_noint,contY_catMcontMbinA_noint,
                             binY_catMbinMbinA_noint,binY_catMbinMbinA_noint,binY_catMcontMbinA_noint,
                             countY_catMbinMbinA_noint,countY_catMbinMbinA_noint,countY_catMcontMbinA_noint,
                             Ycen_catMbinMbinA_noint,cens_catMbinMbinA,delta_catMbinMbinA,
                             contY_catMbinMcatA_noint,contY_catMbinMcatA_noint,contY_catMcontMcatA_noint,
                             binY_catMbinMcatA_noint,binY_catMbinMcatA_noint,binY_catMcontMcatA_noint,
                             countY_catMbinMcatA_noint,countY_catMbinMcatA_noint,countY_catMcontMcatA_noint,
                             Ycen_catMbinMcatA_noint,cens_catMbinMcatA,delta_catMbinMcatA,
                             contY_catMcontMcontA_noint,contY_catMbinMcontA_noint,contY_catMcontMcontA_noint,
                             binY_catMcontMcontA_noint,binY_catMbinMcontA_noint,binY_catMcontMcontA_noint,
                             countY_catMcontMcontA_noint,countY_catMbinMcontA_noint,countY_catMcontMcontA_noint,
                             Ycen_catMcontMcontA_noint,cens_catMcontMcontA,delta_catMcontMcontA,
                             contY_catMcontMbinA_noint,contY_catMbinMbinA_noint,contY_catMcontMbinA_noint,
                             binY_catMcontMbinA_noint,binY_catMbinMbinA_noint,binY_catMcontMbinA_noint,
                             countY_catMcontMbinA_noint,countY_catMbinMbinA_noint,countY_catMcontMbinA_noint,
                             Ycen_catMcontMbinA_noint,cens_catMcontMbinA,delta_catMcontMbinA,
                             contY_catMcontMcatA_noint,contY_catMbinMcatA_noint,contY_catMcontMcatA_noint,
                             binY_catMcontMcatA_noint,binY_catMbinMcatA_noint,binY_catMcontMcatA_noint,
                             countY_catMcontMcatA_noint,countY_catMbinMcatA_noint,countY_catMcontMcatA_noint,
                             Ycen_catMcontMcatA_noint,cens_catMcontMcatA,delta_catMcontMcatA,
                             C_cont,C_bin)

write.table(multipleM_noint_data,file="./Data Simulation and Code Checking/multipleM_noint_data_10000.txt")

##############################################multiple M with exposure-mediator int##########################

# contA

linpred1 = theta0+theta11*A_cont+
  theta21*M_cont_contA+theta22*M_bin_contA+
  theta31*M_cont_contA*A_cont+theta32*M_bin_contA*A_cont+
  theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_cont+
  theta21*(M_cat_contA==1)+theta22*(M_cat_contA==2)+theta23*M_bin_contA+
  theta31*(M_cat_contA==1)*A_cont+theta32*(M_cat_contA==2)*A_cont+
  theta33*M_bin_contA*A_cont+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_cont+
  theta21*(M_cat_contA==1)+theta22*(M_cat_contA==2)+theta23*M_cont_contA+
  theta31*(M_cat_contA==1)*A_cont+theta32*(M_cat_contA==2)*A_cont+
  theta33*M_cont_contA*A_cont+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcontA_EMint = rexp(n,exp(-(linpred1)))
cen_binMcontMcontA<-quantile(Y_surv_binMcontMcontA_EMint,probs = 0.75)
Y_surv_catMbinMcontA_EMint = rexp(n,exp(-(linpred2)))
cen_catMbinMcontA<-quantile(Y_surv_catMbinMcontA_EMint,probs = 0.75)
Y_surv_catMcontMcontA_EMint = rexp(n,exp(-(linpred3)))
cen_catMcontMcontA<-quantile(Y_surv_catMcontMcontA_EMint,probs = 0.75)

Ycen_binMcontMcontA_EMint<-pmin(Y_surv_binMcontMcontA_EMint,cen_binMcontMcontA)
Ycen_catMbinMcontA_EMint<-pmin(Y_surv_catMbinMcontA_EMint,cen_catMbinMcontA)
Ycen_catMcontMcontA_EMint<-pmin(Y_surv_catMcontMcontA_EMint,cen_catMcontMcontA)

cens_binMcontMcontA<-as.numeric(Ycen_binMcontMcontA_EMint<cen_binMcontMcontA)
delta_binMcontMcontA<-as.numeric(Y_surv_binMcontMcontA_EMint>cen_binMcontMcontA)
cens_catMbinMcontA<-as.numeric(Ycen_catMbinMcontA_EMint<cen_catMbinMcontA)
delta_catMbinMcontA<-as.numeric(Y_surv_catMbinMcontA_EMint>cen_catMbinMcontA)
cens_catMcontMcontA<-as.numeric(Ycen_catMcontMcontA_EMint<cen_catMcontMcontA)
delta_catMcontMcontA<-as.numeric(Y_surv_catMcontMcontA_EMint>cen_catMcontMcontA)

contY_binMcontMcontA_EMint = rnorm(n,linpred1,0.2)
contY_catMbinMcontA_EMint = rnorm(n,linpred2,0.2)
contY_catMcontMcontA_EMint = rnorm(n,linpred3,0.2)
binY_binMcontMcontA_EMint = rbinom(n,1,expit(linpred1))
binY_catMbinMcontA_EMint = rbinom(n,1,expit(linpred2))
binY_catMcontMcontA_EMint = rbinom(n,1,expit(linpred3))
countY_binMcontMcontA_EMint = rpois(n,abs(linpred1))
countY_catMbinMcontA_EMint = rpois(n,abs(linpred2))
countY_catMcontMcontA_EMint = rpois(n,abs(linpred3))

# binA

linpred1 = theta0+theta11*A_bin+
  theta21*M_cont_binA+theta22*M_bin_binA+
  theta31*M_cont_binA*A_bin+theta32*M_bin_binA*A_bin+
  theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_bin+
  theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta23*M_bin_binA+
  theta31*(M_cat_binA==1)*A_bin+theta32*(M_cat_binA==2)*A_bin+
  theta33*M_bin_binA*A_bin+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_bin+
  theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta23*M_cont_binA+
  theta31*(M_cat_binA==1)*A_bin+theta32*(M_cat_binA==2)*A_bin+
  theta33*M_cont_binA*A_bin+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMbinA_EMint = rexp(n,exp(-(linpred1)))
cen_binMcontMbinA<-quantile(Y_surv_binMcontMbinA_EMint,probs = 0.75)
Y_surv_catMbinMbinA_EMint = rexp(n,exp(-(linpred2)))
cen_catMbinMbinA<-quantile(Y_surv_catMbinMbinA_EMint,probs = 0.75)
Y_surv_catMcontMbinA_EMint = rexp(n,exp(-(linpred3)))
cen_catMcontMbinA<-quantile(Y_surv_catMcontMbinA_EMint,probs = 0.75)

Ycen_binMcontMbinA_EMint<-pmin(Y_surv_binMcontMbinA_EMint,cen_binMcontMbinA)
Ycen_catMbinMbinA_EMint<-pmin(Y_surv_catMbinMbinA_EMint,cen_catMbinMbinA)
Ycen_catMcontMbinA_EMint<-pmin(Y_surv_catMcontMbinA_EMint,cen_catMcontMbinA)

cens_binMcontMbinA<-as.numeric(Ycen_binMcontMbinA_EMint<cen_binMcontMbinA)
delta_binMcontMbinA<-as.numeric(Y_surv_binMcontMbinA_EMint>cen_binMcontMbinA)
cens_catMbinMbinA<-as.numeric(Ycen_catMbinMbinA_EMint<cen_catMbinMbinA)
delta_catMbinMbinA<-as.numeric(Y_surv_catMbinMbinA_EMint>cen_catMbinMbinA)
cens_catMcontMbinA<-as.numeric(Ycen_catMcontMbinA_EMint<cen_catMcontMbinA)
delta_catMcontMbinA<-as.numeric(Y_surv_catMcontMbinA_EMint>cen_catMcontMbinA)

contY_binMcontMbinA_EMint = rnorm(n,linpred1,0.2)
contY_catMbinMbinA_EMint = rnorm(n,linpred2,0.2)
contY_catMcontMbinA_EMint = rnorm(n,linpred3,0.2)
binY_binMcontMbinA_EMint = rbinom(n,1,expit(linpred1))
binY_catMbinMbinA_EMint = rbinom(n,1,expit(linpred2))
binY_catMcontMbinA_EMint = rbinom(n,1,expit(linpred3))
countY_binMcontMbinA_EMint = rpois(n,abs(linpred1))
countY_catMbinMbinA_EMint = rpois(n,abs(linpred2))
countY_catMcontMbinA_EMint = rpois(n,abs(linpred3))

# catA

linpred1 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+
  theta21*M_cont_catA+theta22*M_bin_catA+
  theta31*M_cont_catA*(A_cat==1)+theta32*M_cont_catA*(A_cat==2)+
  theta33*M_bin_catA*(A_cat==1)+theta34*M_bin_catA*(A_cat==2)+
  theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+
  theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta23*M_bin_catA+
  theta31*(M_cat_catA==1)*(A_cat==1)+theta32*(M_cat_catA==1)*(A_cat==2)+
  theta33*(M_cat_catA==2)*(A_cat==1)+theta34*(M_cat_catA==2)*(A_cat==2)+
  theta35*M_bin_catA*(A_cat==1)+theta36*M_bin_catA*(A_cat==2)+
  theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+
  theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta23*M_cont_catA+
  theta31*(M_cat_catA==1)*(A_cat==1)+theta32*(M_cat_catA==1)*(A_cat==2)+
  theta33*(M_cat_catA==2)*(A_cat==1)+theta34*(M_cat_catA==2)*(A_cat==2)+
  theta35*M_cont_catA*(A_cat==1)+theta36*M_cont_catA*(A_cat==2)+
  theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcatA_EMint = rexp(n,exp(-(linpred1)))
cen_binMcontMcatA<-quantile(Y_surv_binMcontMcatA_EMint,probs = 0.75)
Y_surv_catMbinMcatA_EMint = rexp(n,exp(-(linpred2)))
cen_catMbinMcatA<-quantile(Y_surv_catMbinMcatA_EMint,probs = 0.75)
Y_surv_catMcontMcatA_EMint = rexp(n,exp(-(linpred3)))
cen_catMcontMcatA<-quantile(Y_surv_catMcontMcatA_EMint,probs = 0.75)

Ycen_binMcontMcatA_EMint<-pmin(Y_surv_binMcontMcatA_EMint,cen_binMcontMcatA)
Ycen_catMbinMcatA_EMint<-pmin(Y_surv_catMbinMcatA_EMint,cen_catMbinMcatA)
Ycen_catMcontMcatA_EMint<-pmin(Y_surv_catMcontMcatA_EMint,cen_catMcontMcatA)

cens_binMcontMcatA<-as.numeric(Ycen_binMcontMcatA_EMint<cen_binMcontMcatA)
delta_binMcontMcatA<-as.numeric(Y_surv_binMcontMcatA_EMint>cen_binMcontMcatA)
cens_catMbinMcatA<-as.numeric(Ycen_catMbinMcatA_EMint<cen_catMbinMcatA)
delta_catMbinMcatA<-as.numeric(Y_surv_catMbinMcatA_EMint>cen_catMbinMcatA)
cens_catMcontMcatA<-as.numeric(Ycen_catMcontMcatA_EMint<cen_catMcontMcatA)
delta_catMcontMcatA<-as.numeric(Y_surv_catMcontMcatA_EMint>cen_catMcontMcatA)

contY_binMcontMcatA_EMint = rnorm(n,linpred1,0.2)
contY_catMbinMcatA_EMint = rnorm(n,linpred2,0.2)
contY_catMcontMcatA_EMint = rnorm(n,linpred3,0.2)
binY_binMcontMcatA_EMint = rbinom(n,1,expit(linpred1))
binY_catMbinMcatA_EMint = rbinom(n,1,expit(linpred2))
binY_catMcontMcatA_EMint = rbinom(n,1,expit(linpred3))
countY_binMcontMcatA_EMint = rpois(n,abs(linpred1))
countY_catMbinMcatA_EMint = rpois(n,abs(linpred2))
countY_catMcontMcatA_EMint = rpois(n,abs(linpred3))

multipleM_EMint_data = cbind(A_cont,A_bin,A_cat,
                             M_cont_contA,M_bin_contA,M_cat_contA,
                             M_cont_binA,M_bin_binA,M_cat_binA,
                             M_cont_catA,M_bin_catA,M_cat_catA,
                             contY_binMcontMcontA_EMint,contY_catMbinMcontA_EMint,contY_catMcontMcontA_EMint,
                             binY_binMcontMcontA_EMint,binY_catMbinMcontA_EMint,binY_catMcontMcontA_EMint,
                             countY_binMcontMcontA_EMint,countY_catMbinMcontA_EMint,countY_catMcontMcontA_EMint,
                             Ycen_binMcontMcontA_EMint,cens_binMcontMcontA,delta_binMcontMcontA,
                             contY_binMcontMbinA_EMint,contY_catMbinMbinA_EMint,contY_catMcontMbinA_EMint,
                             binY_binMcontMbinA_EMint,binY_catMbinMbinA_EMint,binY_catMcontMbinA_EMint,
                             countY_binMcontMbinA_EMint,countY_catMbinMbinA_EMint,countY_catMcontMbinA_EMint,
                             Ycen_binMcontMbinA_EMint,cens_binMcontMbinA,delta_binMcontMbinA,
                             contY_binMcontMcatA_EMint,contY_catMbinMcatA_EMint,contY_catMcontMcatA_EMint,
                             binY_binMcontMcatA_EMint,binY_catMbinMcatA_EMint,binY_catMcontMcatA_EMint,
                             countY_binMcontMcatA_EMint,countY_catMbinMcatA_EMint,countY_catMcontMcatA_EMint,
                             Ycen_binMcontMcatA_EMint,cens_binMcontMcatA,delta_binMcontMcatA,
                             contY_catMbinMcontA_EMint,contY_catMbinMcontA_EMint,contY_catMcontMcontA_EMint,
                             binY_catMbinMcontA_EMint,binY_catMbinMcontA_EMint,binY_catMcontMcontA_EMint,
                             countY_catMbinMcontA_EMint,countY_catMbinMcontA_EMint,countY_catMcontMcontA_EMint,
                             Ycen_catMbinMcontA_EMint,cens_catMbinMcontA,delta_catMbinMcontA,
                             contY_catMbinMbinA_EMint,contY_catMbinMbinA_EMint,contY_catMcontMbinA_EMint,
                             binY_catMbinMbinA_EMint,binY_catMbinMbinA_EMint,binY_catMcontMbinA_EMint,
                             countY_catMbinMbinA_EMint,countY_catMbinMbinA_EMint,countY_catMcontMbinA_EMint,
                             Ycen_catMbinMbinA_EMint,cens_catMbinMbinA,delta_catMbinMbinA,
                             contY_catMbinMcatA_EMint,contY_catMbinMcatA_EMint,contY_catMcontMcatA_EMint,
                             binY_catMbinMcatA_EMint,binY_catMbinMcatA_EMint,binY_catMcontMcatA_EMint,
                             countY_catMbinMcatA_EMint,countY_catMbinMcatA_EMint,countY_catMcontMcatA_EMint,
                             Ycen_catMbinMcatA_EMint,cens_catMbinMcatA,delta_catMbinMcatA,
                             contY_catMcontMcontA_EMint,contY_catMbinMcontA_EMint,contY_catMcontMcontA_EMint,
                             binY_catMcontMcontA_EMint,binY_catMbinMcontA_EMint,binY_catMcontMcontA_EMint,
                             countY_catMcontMcontA_EMint,countY_catMbinMcontA_EMint,countY_catMcontMcontA_EMint,
                             Ycen_catMcontMcontA_EMint,cens_catMcontMcontA,delta_catMcontMcontA,
                             contY_catMcontMbinA_EMint,contY_catMbinMbinA_EMint,contY_catMcontMbinA_EMint,
                             binY_catMcontMbinA_EMint,binY_catMbinMbinA_EMint,binY_catMcontMbinA_EMint,
                             countY_catMcontMbinA_EMint,countY_catMbinMbinA_EMint,countY_catMcontMbinA_EMint,
                             Ycen_catMcontMbinA_EMint,cens_catMcontMbinA,delta_catMcontMbinA,
                             contY_catMcontMcatA_EMint,contY_catMbinMcatA_EMint,contY_catMcontMcatA_EMint,
                             binY_catMcontMcatA_EMint,binY_catMbinMcatA_EMint,binY_catMcontMcatA_EMint,
                             countY_catMcontMcatA_EMint,countY_catMbinMcatA_EMint,countY_catMcontMcatA_EMint,
                             Ycen_catMcontMcatA_EMint,cens_catMcontMcatA,delta_catMcontMcatA,
                             C_cont,C_bin)

write.table(multipleM_EMint_data,file="./Data Simulation and Code Checking/multipleM_EMint_data_10000.txt")


##############################multiple M with mediator-mediator int##########################

# contA

linpred1 = theta0+theta11*A_cont+theta21*M_cont_contA+theta22*M_bin_contA+theta31*M_cont_contA*M_bin_contA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_cont+theta21*(M_cat_contA==1)+theta22*(M_cat_contA==2)+theta23*M_bin_contA+theta31*(M_cat_contA==1)*M_bin_contA+theta32*(M_cat_contA==2)*M_bin_contA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_cont+theta21*M_cont_contA+theta22*(M_cat_contA==1)+theta23*(M_cat_contA==2)+theta31*(M_cat_contA==1)*M_cont_contA+theta32*(M_cat_contA==2)*M_cont_contA+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcontA_MMint = rexp(n,exp(-(linpred1)))
cen_binMcontMcontA<-quantile(Y_surv_binMcontMcontA_MMint,probs = 0.75)
Y_surv_catMbinMcontA_MMint = rexp(n,exp(-(linpred2)))
cen_catMbinMcontA<-quantile(Y_surv_catMbinMcontA_MMint,probs = 0.75)
Y_surv_catMcontMcontA_MMint = rexp(n,exp(-(linpred3)))
cen_catMcontMcontA<-quantile(Y_surv_catMcontMcontA_MMint,probs = 0.75)

Ycen_binMcontMcontA_MMint<-pmin(Y_surv_binMcontMcontA_MMint,cen_binMcontMcontA)
Ycen_catMbinMcontA_MMint<-pmin(Y_surv_catMbinMcontA_MMint,cen_catMbinMcontA)
Ycen_catMcontMcontA_MMint<-pmin(Y_surv_catMcontMcontA_MMint,cen_catMcontMcontA)

cens_binMcontMcontA<-as.numeric(Ycen_binMcontMcontA_MMint<cen_binMcontMcontA)
delta_binMcontMcontA<-as.numeric(Y_surv_binMcontMcontA_MMint>cen_binMcontMcontA)
cens_catMbinMcontA<-as.numeric(Ycen_catMbinMcontA_MMint<cen_catMbinMcontA)
delta_catMbinMcontA<-as.numeric(Y_surv_catMbinMcontA_MMint>cen_catMbinMcontA)
cens_catMcontMcontA<-as.numeric(Ycen_catMcontMcontA_MMint<cen_catMcontMcontA)
delta_catMcontMcontA<-as.numeric(Y_surv_catMcontMcontA_MMint>cen_catMcontMcontA)

contY_binMcontMcontA_MMint = rnorm(n,linpred1,0.2)
contY_catMbinMcontA_MMint = rnorm(n,linpred2,0.2)
contY_catMcontMcontA_MMint = rnorm(n,linpred3,0.2)
binY_binMcontMcontA_MMint = rbinom(n,1,expit(linpred1))
binY_catMbinMcontA_MMint = rbinom(n,1,expit(linpred2))
binY_catMcontMcontA_MMint = rbinom(n,1,expit(linpred3))
countY_binMcontMcontA_MMint = rpois(n,abs(linpred1))
countY_catMbinMcontA_MMint = rpois(n,abs(linpred2))
countY_catMcontMcontA_MMint = rpois(n,abs(linpred3))

# binA

linpred1 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*M_bin_binA+theta31*M_cont_binA*M_bin_binA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_bin+theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta23*M_bin_binA+theta31*(M_cat_binA==1)*M_bin_binA+theta32*(M_cat_binA==2)*M_bin_binA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*(M_cat_binA==1)+theta23*(M_cat_binA==2)+theta31*(M_cat_binA==1)*M_cont_binA+theta32*(M_cat_binA==2)*M_cont_binA+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMbinA_MMint = rexp(n,exp(-(linpred1)))
cen_binMcontMbinA<-quantile(Y_surv_binMcontMbinA_MMint,probs = 0.75)
Y_surv_catMbinMbinA_MMint = rexp(n,exp(-(linpred2)))
cen_catMbinMbinA<-quantile(Y_surv_catMbinMbinA_MMint,probs = 0.75)
Y_surv_catMcontMbinA_MMint = rexp(n,exp(-(linpred3)))
cen_catMcontMbinA<-quantile(Y_surv_catMcontMbinA_MMint,probs = 0.75)

Ycen_binMcontMbinA_MMint<-pmin(Y_surv_binMcontMbinA_MMint,cen_binMcontMbinA)
Ycen_catMbinMbinA_MMint<-pmin(Y_surv_catMbinMbinA_MMint,cen_catMbinMbinA)
Ycen_catMcontMbinA_MMint<-pmin(Y_surv_catMcontMbinA_MMint,cen_catMcontMbinA)

cens_binMcontMbinA<-as.numeric(Ycen_binMcontMbinA_MMint<cen_binMcontMbinA)
delta_binMcontMbinA<-as.numeric(Y_surv_binMcontMbinA_MMint>cen_binMcontMbinA)
cens_catMbinMbinA<-as.numeric(Ycen_catMbinMbinA_MMint<cen_catMbinMbinA)
delta_catMbinMbinA<-as.numeric(Y_surv_catMbinMbinA_MMint>cen_catMbinMbinA)
cens_catMcontMbinA<-as.numeric(Ycen_catMcontMbinA_MMint<cen_catMcontMbinA)
delta_catMcontMbinA<-as.numeric(Y_surv_catMcontMbinA_MMint>cen_catMcontMbinA)

contY_binMcontMbinA_MMint = rnorm(n,linpred1,0.2)
contY_catMbinMbinA_MMint = rnorm(n,linpred2,0.2)
contY_catMcontMbinA_MMint = rnorm(n,linpred3,0.2)
binY_binMcontMbinA_MMint = rbinom(n,1,expit(linpred1))
binY_catMbinMbinA_MMint = rbinom(n,1,expit(linpred2))
binY_catMcontMbinA_MMint = rbinom(n,1,expit(linpred3))
countY_binMcontMbinA_MMint = rpois(n,abs(linpred1))
countY_catMbinMbinA_MMint = rpois(n,abs(linpred2))
countY_catMcontMbinA_MMint = rpois(n,abs(linpred3))

# catA

linpred1 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_cont_catA+theta22*M_bin_catA+theta31*M_cont_catA*M_bin_catA+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta23*M_bin_catA+theta31*(M_cat_catA==1)*M_bin_catA+theta32*(M_cat_catA==2)*M_bin_catA+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_cont_catA+theta22*(M_cat_catA==1)+theta23*(M_cat_catA==2)+theta31*(M_cat_catA==1)*M_cont_catA+theta32*(M_cat_catA==2)*M_cont_catA+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcatA_MMint = rexp(n,exp(-(linpred1)))
cen_binMcontMcatA<-quantile(Y_surv_binMcontMcatA_MMint,probs = 0.75)
Y_surv_catMbinMcatA_MMint = rexp(n,exp(-(linpred2)))
cen_catMbinMcatA<-quantile(Y_surv_catMbinMcatA_MMint,probs = 0.75)
Y_surv_catMcontMcatA_MMint = rexp(n,exp(-(linpred3)))
cen_catMcontMcatA<-quantile(Y_surv_catMcontMcatA_MMint,probs = 0.75)

Ycen_binMcontMcatA_MMint<-pmin(Y_surv_binMcontMcatA_MMint,cen_binMcontMcatA)
Ycen_catMbinMcatA_MMint<-pmin(Y_surv_catMbinMcatA_MMint,cen_catMbinMcatA)
Ycen_catMcontMcatA_MMint<-pmin(Y_surv_catMcontMcatA_MMint,cen_catMcontMcatA)

cens_binMcontMcatA<-as.numeric(Ycen_binMcontMcatA_MMint<cen_binMcontMcatA)
delta_binMcontMcatA<-as.numeric(Y_surv_binMcontMcatA_MMint>cen_binMcontMcatA)
cens_catMbinMcatA<-as.numeric(Ycen_catMbinMcatA_MMint<cen_catMbinMcatA)
delta_catMbinMcatA<-as.numeric(Y_surv_catMbinMcatA_MMint>cen_catMbinMcatA)
cens_catMcontMcatA<-as.numeric(Ycen_catMcontMcatA_MMint<cen_catMcontMcatA)
delta_catMcontMcatA<-as.numeric(Y_surv_catMcontMcatA_MMint>cen_catMcontMcatA)

contY_binMcontMcatA_MMint = rnorm(n,linpred1,0.2)
contY_catMbinMcatA_MMint = rnorm(n,linpred2,0.2)
contY_catMcontMcatA_MMint = rnorm(n,linpred3,0.2)
binY_binMcontMcatA_MMint = rbinom(n,1,expit(linpred1))
binY_catMbinMcatA_MMint = rbinom(n,1,expit(linpred2))
binY_catMcontMcatA_MMint = rbinom(n,1,expit(linpred3))
countY_binMcontMcatA_MMint = rpois(n,abs(linpred1))
countY_catMbinMcatA_MMint = rpois(n,abs(linpred2))
countY_catMcontMcatA_MMint = rpois(n,abs(linpred3))

multipleM_MMint_data = cbind(A_cont,A_bin,A_cat,
                             M_cont_contA,M_bin_contA,M_cat_contA,
                             M_cont_binA,M_bin_binA,M_cat_binA,
                             M_cont_catA,M_bin_catA,M_cat_catA,
                             contY_binMcontMcontA_MMint,contY_catMbinMcontA_MMint,contY_catMcontMcontA_MMint,
                             binY_binMcontMcontA_MMint,binY_catMbinMcontA_MMint,binY_catMcontMcontA_MMint,
                             countY_binMcontMcontA_MMint,countY_catMbinMcontA_MMint,countY_catMcontMcontA_MMint,
                             Ycen_binMcontMcontA_MMint,cens_binMcontMcontA,delta_binMcontMcontA,
                             contY_binMcontMbinA_MMint,contY_catMbinMbinA_MMint,contY_catMcontMbinA_MMint,
                             binY_binMcontMbinA_MMint,binY_catMbinMbinA_MMint,binY_catMcontMbinA_MMint,
                             countY_binMcontMbinA_MMint,countY_catMbinMbinA_MMint,countY_catMcontMbinA_MMint,
                             Ycen_binMcontMbinA_MMint,cens_binMcontMbinA,delta_binMcontMbinA,
                             contY_binMcontMcatA_MMint,contY_catMbinMcatA_MMint,contY_catMcontMcatA_MMint,
                             binY_binMcontMcatA_MMint,binY_catMbinMcatA_MMint,binY_catMcontMcatA_MMint,
                             countY_binMcontMcatA_MMint,countY_catMbinMcatA_MMint,countY_catMcontMcatA_MMint,
                             Ycen_binMcontMcatA_MMint,cens_binMcontMcatA,delta_binMcontMcatA,
                             contY_catMbinMcontA_MMint,contY_catMbinMcontA_MMint,contY_catMcontMcontA_MMint,
                             binY_catMbinMcontA_MMint,binY_catMbinMcontA_MMint,binY_catMcontMcontA_MMint,
                             countY_catMbinMcontA_MMint,countY_catMbinMcontA_MMint,countY_catMcontMcontA_MMint,
                             Ycen_catMbinMcontA_MMint,cens_catMbinMcontA,delta_catMbinMcontA,
                             contY_catMbinMbinA_MMint,contY_catMbinMbinA_MMint,contY_catMcontMbinA_MMint,
                             binY_catMbinMbinA_MMint,binY_catMbinMbinA_MMint,binY_catMcontMbinA_MMint,
                             countY_catMbinMbinA_MMint,countY_catMbinMbinA_MMint,countY_catMcontMbinA_MMint,
                             Ycen_catMbinMbinA_MMint,cens_catMbinMbinA,delta_catMbinMbinA,
                             contY_catMbinMcatA_MMint,contY_catMbinMcatA_MMint,contY_catMcontMcatA_MMint,
                             binY_catMbinMcatA_MMint,binY_catMbinMcatA_MMint,binY_catMcontMcatA_MMint,
                             countY_catMbinMcatA_MMint,countY_catMbinMcatA_MMint,countY_catMcontMcatA_MMint,
                             Ycen_catMbinMcatA_MMint,cens_catMbinMcatA,delta_catMbinMcatA,
                             contY_catMcontMcontA_MMint,contY_catMbinMcontA_MMint,contY_catMcontMcontA_MMint,
                             binY_catMcontMcontA_MMint,binY_catMbinMcontA_MMint,binY_catMcontMcontA_MMint,
                             countY_catMcontMcontA_MMint,countY_catMbinMcontA_MMint,countY_catMcontMcontA_MMint,
                             Ycen_catMcontMcontA_MMint,cens_catMcontMcontA,delta_catMcontMcontA,
                             contY_catMcontMbinA_MMint,contY_catMbinMbinA_MMint,contY_catMcontMbinA_MMint,
                             binY_catMcontMbinA_MMint,binY_catMbinMbinA_MMint,binY_catMcontMbinA_MMint,
                             countY_catMcontMbinA_MMint,countY_catMbinMbinA_MMint,countY_catMcontMbinA_MMint,
                             Ycen_catMcontMbinA_MMint,cens_catMcontMbinA,delta_catMcontMbinA,
                             contY_catMcontMcatA_MMint,contY_catMbinMcatA_MMint,contY_catMcontMcatA_MMint,
                             binY_catMcontMcatA_MMint,binY_catMbinMcatA_MMint,binY_catMcontMcatA_MMint,
                             countY_catMcontMcatA_MMint,countY_catMbinMcatA_MMint,countY_catMcontMcatA_MMint,
                             Ycen_catMcontMcatA_MMint,cens_catMcontMcatA,delta_catMcontMcatA,
                             C_cont,C_bin)

write.table(multipleM_MMint_data,file="./Data Simulation and Code Checking/multipleM_MMint_data_10000.txt")

##############################################multiple M with multiway int##########################

# contA

linpred1 = theta0+theta11*A_cont+theta21*M_cont_contA+theta22*M_bin_contA+theta31*M_cont_contA*M_bin_contA*A_cont+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_cont+theta21*(M_cat_contA==1)+theta22*(M_cat_contA==2)+theta23*M_bin_contA+theta31*(M_cat_contA==1)*M_bin_contA*A_cont+theta32*(M_cat_contA==2)*M_bin_contA*A_cont+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_cont+theta21*M_cont_contA+theta22*(M_cat_contA==1)+theta23*(M_cat_contA==2)+theta31*(M_cat_contA==1)*M_cont_contA*A_cont+theta32*(M_cat_contA==2)*M_cont_contA*A_cont+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcontA_EMMint = rexp(n,exp(-(linpred1)))
cen_binMcontMcontA<-quantile(Y_surv_binMcontMcontA_EMMint,probs = 0.75)
Y_surv_catMbinMcontA_EMMint = rexp(n,exp(-(linpred2)))
cen_catMbinMcontA<-quantile(Y_surv_catMbinMcontA_EMMint,probs = 0.75)
Y_surv_catMcontMcontA_EMMint = rexp(n,exp(-(linpred3)))
cen_catMcontMcontA<-quantile(Y_surv_catMcontMcontA_EMMint,probs = 0.75)

Ycen_binMcontMcontA_EMMint<-pmin(Y_surv_binMcontMcontA_EMMint,cen_binMcontMcontA)
Ycen_catMbinMcontA_EMMint<-pmin(Y_surv_catMbinMcontA_EMMint,cen_catMbinMcontA)
Ycen_catMcontMcontA_EMMint<-pmin(Y_surv_catMcontMcontA_EMMint,cen_catMcontMcontA)

cens_binMcontMcontA<-as.numeric(Ycen_binMcontMcontA_EMMint<cen_binMcontMcontA)
delta_binMcontMcontA<-as.numeric(Y_surv_binMcontMcontA_EMMint>cen_binMcontMcontA)
cens_catMbinMcontA<-as.numeric(Ycen_catMbinMcontA_EMMint<cen_catMbinMcontA)
delta_catMbinMcontA<-as.numeric(Y_surv_catMbinMcontA_EMMint>cen_catMbinMcontA)
cens_catMcontMcontA<-as.numeric(Ycen_catMcontMcontA_EMMint<cen_catMcontMcontA)
delta_catMcontMcontA<-as.numeric(Y_surv_catMcontMcontA_EMMint>cen_catMcontMcontA)

contY_binMcontMcontA_EMMint = rnorm(n,linpred1,0.2)
contY_catMbinMcontA_EMMint = rnorm(n,linpred2,0.2)
contY_catMcontMcontA_EMMint = rnorm(n,linpred3,0.2)
binY_binMcontMcontA_EMMint = rbinom(n,1,expit(linpred1))
binY_catMbinMcontA_EMMint = rbinom(n,1,expit(linpred2))
binY_catMcontMcontA_EMMint = rbinom(n,1,expit(linpred3))
countY_binMcontMcontA_EMMint = rpois(n,abs(linpred1))
countY_catMbinMcontA_EMMint = rpois(n,abs(linpred2))
countY_catMcontMcontA_EMMint = rpois(n,abs(linpred3))

# binA

linpred1 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*M_bin_binA+theta31*M_cont_binA*M_bin_binA*A_bin+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*A_bin+theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta23*M_bin_binA+theta31*(M_cat_binA==1)*M_bin_binA*A_bin+theta32*(M_cat_binA==2)*M_bin_binA*A_bin+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*(M_cat_binA==1)+theta23*(M_cat_binA==2)+theta31*(M_cat_binA==1)*M_cont_binA*A_bin+theta32*(M_cat_binA==2)*M_cont_binA*A_bin+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMbinA_EMMint = rexp(n,exp(-(linpred1)))
cen_binMcontMbinA<-quantile(Y_surv_binMcontMbinA_EMMint,probs = 0.75)
Y_surv_catMbinMbinA_EMMint = rexp(n,exp(-(linpred2)))
cen_catMbinMbinA<-quantile(Y_surv_catMbinMbinA_EMMint,probs = 0.75)
Y_surv_catMcontMbinA_EMMint = rexp(n,exp(-(linpred3)))
cen_catMcontMbinA<-quantile(Y_surv_catMcontMbinA_EMMint,probs = 0.75)

Ycen_binMcontMbinA_EMMint<-pmin(Y_surv_binMcontMbinA_EMMint,cen_binMcontMbinA)
Ycen_catMbinMbinA_EMMint<-pmin(Y_surv_catMbinMbinA_EMMint,cen_catMbinMbinA)
Ycen_catMcontMbinA_EMMint<-pmin(Y_surv_catMcontMbinA_EMMint,cen_catMcontMbinA)

cens_binMcontMbinA<-as.numeric(Ycen_binMcontMbinA_EMMint<cen_binMcontMbinA)
delta_binMcontMbinA<-as.numeric(Y_surv_binMcontMbinA_EMMint>cen_binMcontMbinA)
cens_catMbinMbinA<-as.numeric(Ycen_catMbinMbinA_EMMint<cen_catMbinMbinA)
delta_catMbinMbinA<-as.numeric(Y_surv_catMbinMbinA_EMMint>cen_catMbinMbinA)
cens_catMcontMbinA<-as.numeric(Ycen_catMcontMbinA_EMMint<cen_catMcontMbinA)
delta_catMcontMbinA<-as.numeric(Y_surv_catMcontMbinA_EMMint>cen_catMcontMbinA)

contY_binMcontMbinA_EMMint = rnorm(n,linpred1,0.2)
contY_catMbinMbinA_EMMint = rnorm(n,linpred2,0.2)
contY_catMcontMbinA_EMMint = rnorm(n,linpred3,0.2)
binY_binMcontMbinA_EMMint = rbinom(n,1,expit(linpred1))
binY_catMbinMbinA_EMMint = rbinom(n,1,expit(linpred2))
binY_catMcontMbinA_EMMint = rbinom(n,1,expit(linpred3))
countY_binMcontMbinA_EMMint = rpois(n,abs(linpred1))
countY_catMbinMbinA_EMMint = rpois(n,abs(linpred2))
countY_catMcontMbinA_EMMint = rpois(n,abs(linpred3))

# catA

linpred1 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_cont_catA+theta22*M_bin_catA+theta31*M_cont_catA*M_bin_catA*(A_cat==1)+theta32*M_cont_catA*M_bin_catA*(A_cat==2)+theta41*C_cont+theta42*C_bin
linpred2 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*(M_cat_catA==1)+theta22*(M_cat_catA==2)+theta23*M_bin_catA+theta31*(M_cat_catA==1)*M_bin_catA*(A_cat==1)+theta32*(M_cat_catA==2)*M_bin_catA*(A_cat==1)+theta33*(M_cat_catA==1)*M_bin_catA*(A_cat==2)+theta34*(M_cat_catA==2)*M_bin_catA*(A_cat==2)+theta41*C_cont+theta42*C_bin
linpred3 = theta0+theta11*(A_cat==1)+theta12*(A_cat==2)+theta21*M_cont_catA+theta22*(M_cat_catA==1)+theta23*(M_cat_catA==2)+theta31*(M_cat_catA==1)*M_cont_catA*(A_cat==1)+theta32*(M_cat_catA==2)*M_cont_catA*(A_cat==1)+theta33*(M_cat_catA==1)*M_cont_catA*(A_cat==2)+theta34*(M_cat_catA==2)*M_cont_catA*(A_cat==2)+theta41*C_cont+theta42*C_bin

Y_surv_binMcontMcatA_EMMint = rexp(n,exp(-(linpred1)))
cen_binMcontMcatA<-quantile(Y_surv_binMcontMcatA_EMMint,probs = 0.75)
Y_surv_catMbinMcatA_EMMint = rexp(n,exp(-(linpred2)))
cen_catMbinMcatA<-quantile(Y_surv_catMbinMcatA_EMMint,probs = 0.75)
Y_surv_catMcontMcatA_EMMint = rexp(n,exp(-(linpred3)))
cen_catMcontMcatA<-quantile(Y_surv_catMcontMcatA_EMMint,probs = 0.75)

Ycen_binMcontMcatA_EMMint<-pmin(Y_surv_binMcontMcatA_EMMint,cen_binMcontMcatA)
Ycen_catMbinMcatA_EMMint<-pmin(Y_surv_catMbinMcatA_EMMint,cen_catMbinMcatA)
Ycen_catMcontMcatA_EMMint<-pmin(Y_surv_catMcontMcatA_EMMint,cen_catMcontMcatA)

cens_binMcontMcatA<-as.numeric(Ycen_binMcontMcatA_EMMint<cen_binMcontMcatA)
delta_binMcontMcatA<-as.numeric(Y_surv_binMcontMcatA_EMMint>cen_binMcontMcatA)
cens_catMbinMcatA<-as.numeric(Ycen_catMbinMcatA_EMMint<cen_catMbinMcatA)
delta_catMbinMcatA<-as.numeric(Y_surv_catMbinMcatA_EMMint>cen_catMbinMcatA)
cens_catMcontMcatA<-as.numeric(Ycen_catMcontMcatA_EMMint<cen_catMcontMcatA)
delta_catMcontMcatA<-as.numeric(Y_surv_catMcontMcatA_EMMint>cen_catMcontMcatA)

contY_binMcontMcatA_EMMint = rnorm(n,linpred1,0.2)
contY_catMbinMcatA_EMMint = rnorm(n,linpred2,0.2)
contY_catMcontMcatA_EMMint = rnorm(n,linpred3,0.2)
binY_binMcontMcatA_EMMint = rbinom(n,1,expit(linpred1))
binY_catMbinMcatA_EMMint = rbinom(n,1,expit(linpred2))
binY_catMcontMcatA_EMMint = rbinom(n,1,expit(linpred3))
countY_binMcontMcatA_EMMint = rpois(n,abs(linpred1))
countY_catMbinMcatA_EMMint = rpois(n,abs(linpred2))
countY_catMcontMcatA_EMMint = rpois(n,abs(linpred3))


multipleM_EMMint_data = cbind(A_cont,A_bin,A_cat,
                             M_cont_contA,M_bin_contA,M_cat_contA,
                             M_cont_binA,M_bin_binA,M_cat_binA,
                             M_cont_catA,M_bin_catA,M_cat_catA,
                             contY_binMcontMcontA_EMMint,contY_catMbinMcontA_EMMint,contY_catMcontMcontA_EMMint,
                             binY_binMcontMcontA_EMMint,binY_catMbinMcontA_EMMint,binY_catMcontMcontA_EMMint,
                             countY_binMcontMcontA_EMMint,countY_catMbinMcontA_EMMint,countY_catMcontMcontA_EMMint,
                             Ycen_binMcontMcontA_EMMint,cens_binMcontMcontA,delta_binMcontMcontA,
                             contY_binMcontMbinA_EMMint,contY_catMbinMbinA_EMMint,contY_catMcontMbinA_EMMint,
                             binY_binMcontMbinA_EMMint,binY_catMbinMbinA_EMMint,binY_catMcontMbinA_EMMint,
                             countY_binMcontMbinA_EMMint,countY_catMbinMbinA_EMMint,countY_catMcontMbinA_EMMint,
                             Ycen_binMcontMbinA_EMMint,cens_binMcontMbinA,delta_binMcontMbinA,
                             contY_binMcontMcatA_EMMint,contY_catMbinMcatA_EMMint,contY_catMcontMcatA_EMMint,
                             binY_binMcontMcatA_EMMint,binY_catMbinMcatA_EMMint,binY_catMcontMcatA_EMMint,
                             countY_binMcontMcatA_EMMint,countY_catMbinMcatA_EMMint,countY_catMcontMcatA_EMMint,
                             Ycen_binMcontMcatA_EMMint,cens_binMcontMcatA,delta_binMcontMcatA,
                             contY_catMbinMcontA_EMMint,contY_catMbinMcontA_EMMint,contY_catMcontMcontA_EMMint,
                             binY_catMbinMcontA_EMMint,binY_catMbinMcontA_EMMint,binY_catMcontMcontA_EMMint,
                             countY_catMbinMcontA_EMMint,countY_catMbinMcontA_EMMint,countY_catMcontMcontA_EMMint,
                             Ycen_catMbinMcontA_EMMint,cens_catMbinMcontA,delta_catMbinMcontA,
                             contY_catMbinMbinA_EMMint,contY_catMbinMbinA_EMMint,contY_catMcontMbinA_EMMint,
                             binY_catMbinMbinA_EMMint,binY_catMbinMbinA_EMMint,binY_catMcontMbinA_EMMint,
                             countY_catMbinMbinA_EMMint,countY_catMbinMbinA_EMMint,countY_catMcontMbinA_EMMint,
                             Ycen_catMbinMbinA_EMMint,cens_catMbinMbinA,delta_catMbinMbinA,
                             contY_catMbinMcatA_EMMint,contY_catMbinMcatA_EMMint,contY_catMcontMcatA_EMMint,
                             binY_catMbinMcatA_EMMint,binY_catMbinMcatA_EMMint,binY_catMcontMcatA_EMMint,
                             countY_catMbinMcatA_EMMint,countY_catMbinMcatA_EMMint,countY_catMcontMcatA_EMMint,
                             Ycen_catMbinMcatA_EMMint,cens_catMbinMcatA,delta_catMbinMcatA,
                             contY_catMcontMcontA_EMMint,contY_catMbinMcontA_EMMint,contY_catMcontMcontA_EMMint,
                             binY_catMcontMcontA_EMMint,binY_catMbinMcontA_EMMint,binY_catMcontMcontA_EMMint,
                             countY_catMcontMcontA_EMMint,countY_catMbinMcontA_EMMint,countY_catMcontMcontA_EMMint,
                             Ycen_catMcontMcontA_EMMint,cens_catMcontMcontA,delta_catMcontMcontA,
                             contY_catMcontMbinA_EMMint,contY_catMbinMbinA_EMMint,contY_catMcontMbinA_EMMint,
                             binY_catMcontMbinA_EMMint,binY_catMbinMbinA_EMMint,binY_catMcontMbinA_EMMint,
                             countY_catMcontMbinA_EMMint,countY_catMbinMbinA_EMMint,countY_catMcontMbinA_EMMint,
                             Ycen_catMcontMbinA_EMMint,cens_catMcontMbinA,delta_catMcontMbinA,
                             contY_catMcontMcatA_EMMint,contY_catMbinMcatA_EMMint,contY_catMcontMcatA_EMMint,
                             binY_catMcontMcatA_EMMint,binY_catMbinMcatA_EMMint,binY_catMcontMcatA_EMMint,
                             countY_catMcontMcatA_EMMint,countY_catMbinMcatA_EMMint,countY_catMcontMcatA_EMMint,
                             Ycen_catMcontMcatA_EMMint,cens_catMcontMcatA,delta_catMcontMcatA,
                             C_cont,C_bin)

write.table(multipleM_EMMint_data,file="./Data Simulation and Code Checking/multipleM_EMMint_data_10000.txt")

#########################################post-exposure confounding##############################

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2
beta31 = 0.1
beta41 = 0.3
beta51 = 0.25

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3
beta32 = 0.5
beta42 = 0.1
beta52 = 0.2

#generate Bin M, Cont M, and Cat M
linpred1 = (beta01+beta11*A_bin+beta21*C_cont)
linpred2 = (beta02+beta12*A_bin+beta22*C_cont)

probl1_bin = expit(linpred1)
probl0_cat = 1/(1+exp(linpred1)+exp(linpred2))
probl1_cat = exp(linpred1)/(1+exp(linpred1)+exp(linpred2))
probl2_cat = exp(linpred2)/(1+exp(linpred1)+exp(linpred2))

L_bin_binA = rbinom(n,1,probl1_bin)
L_cat_binA = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                             prob=c(probl0_cat[x],
                                                    probl1_cat[x],
                                                    probl2_cat[x])))

#generate Bin M, Cont M, and Cat M
linpred3 = (beta01+beta11*A_bin+beta21*C_cont+beta31*L_bin_binA+beta41*(L_cat_binA==1)+beta51*(L_cat_binA==2))
linpred4 = (beta02+beta12*A_bin+beta22*C_cont+beta32*L_bin_binA+beta42*(L_cat_binA==1)+beta52*(L_cat_binA==2))

probm1_bin = expit(linpred3)
probm0_cat = 1/(1+exp(linpred3)+exp(linpred4))
probm1_cat = exp(linpred3)/(1+exp(linpred3)+exp(linpred4))
probm2_cat = exp(linpred4)/(1+exp(linpred3)+exp(linpred4))

M_bin_binA = rbinom(n,1,probm1_bin)
M_cont_binA = rnorm(n,linpred3,0.1)
M_cat_binA = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                             prob=c(probm0_cat[x],
                                                    probm1_cat[x],
                                                    probm2_cat[x])))

#generate Cont Y, Bin Y, Count Y with int
linpred1 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*M_bin_binA+theta31*A_bin*M_bin_binA+theta41*C_cont+theta42*L_bin_binA+theta43*(L_cat_binA==1)+theta44*(L_cat_binA==2)
linpred2 = theta0+theta11*A_bin+theta21*(M_cat_binA==1)+theta22*(M_cat_binA==2)+theta23*M_bin_binA+theta31*A_bin*M_bin_binA+theta41*C_cont+theta42*L_bin_binA+theta43*(L_cat_binA==1)+theta44*(L_cat_binA==2)
linpred3 = theta0+theta11*A_bin+theta21*M_cont_binA+theta22*(M_cat_binA==1)+theta23*(M_cat_binA==2)+theta31*A_bin*M_cont_binA+theta41*C_cont+theta42*L_bin_binA+theta43*(L_cat_binA==1)+theta44*(L_cat_binA==2)

contY_binMcontMbinA_EMint = rnorm(n,linpred1,0.2)
contY_catMbinMbinA_EMint = rnorm(n,linpred2,0.2)
contY_catMcontMbinA_EMint = rnorm(n,linpred3,0.2)
binY_binMcontMbinA_EMint = rbinom(n,1,expit(linpred1))
binY_catMbinMbinA_EMint = rbinom(n,1,expit(linpred2))
binY_catMcontMbinA_EMint = rbinom(n,1,expit(linpred3))
countY_binMcontMbinA_EMint = rpois(n,abs(linpred1))
countY_catMbinMbinA_EMint = rpois(n,abs(linpred2))
countY_catMcontMbinA_EMint = rpois(n,abs(linpred3))

multipleM_EMint_postcovar_data = cbind(A_bin, L_bin_binA, L_cat_binA,
                                       M_cont_binA,M_bin_binA,M_cat_binA,
                                       contY_binMcontMbinA_EMint,contY_catMbinMbinA_EMint,contY_catMcontMbinA_EMint,
                                       binY_binMcontMbinA_EMint,binY_catMbinMbinA_EMint,binY_catMcontMbinA_EMint,
                                       countY_binMcontMbinA_EMint,countY_catMbinMbinA_EMint,countY_catMcontMbinA_EMint,
                                       C_cont)

write.table(multipleM_EMint_postcovar_data,file="./Data Simulation and Code Checking/multipleM_EMint_postcovar_data_10000.txt")
