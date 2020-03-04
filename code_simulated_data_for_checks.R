rm(list=ls())

library(survival)

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3

beta03 = -0.1
beta13 = 0.2
beta23 = 0.3

beta04 = -0.4
beta14 = 0.5
beta24 = 0.2

theta0 = -5
theta1 =0.8
theta21 = 1.8
theta22 = 1.2
theta23 = 1.5
theta24 = 1
theta31 = 0.2
theta32 = 0.4
theta33 = 0.1
theta34 = 0.3
theta4 = 0.1

meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M, Cont M, and Cat M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)

probm1_bin = expit(linpred1)
probm0_cat = 1/(1+exp(linpred1)+exp(linpred2))
probm1_cat = exp(linpred1)/(1+exp(linpred1)+exp(linpred2))
probm2_cat = exp(linpred2)/(1+exp(linpred1)+exp(linpred2))

M_bin = rbinom(n,1,probm1_bin)
M_cont = rnorm(n,linpred1,0.1)
M_cat = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
               prob=c(probm0_cat[x],
                      probm1_cat[x],
                      probm2_cat[x])))

##################################################no int###################################################

#generate survival data
linpred = exp(-(theta0+theta1*A+theta21*M_bin+theta4*C))
Y_surv_binM_noint = rexp(n,linpred)
cen_binM<-quantile(Y_surv_binM_noint,probs = 0.75)

linpred = exp(-(theta0+theta1*A+theta21*M_cont+theta4*C))
Y_surv_contM_noint = rexp(n,linpred)
cen_contM<-quantile(Y_surv_contM_noint,probs = 0.75)

linpred = exp(-(theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta4*C))
Y_surv_catM_noint = rexp(n,linpred)
cen_catM<-quantile(Y_surv_catM_noint,probs = 0.75)

#censor observations after censoring time
Ycen_binM_noint<-pmin(Y_surv_binM_noint,cen_binM)
Ycen_contM_noint<-pmin(Y_surv_contM_noint,cen_contM)
Ycen_catM_noint<-pmin(Y_surv_catM_noint,cen_catM)
#generate event indicator
cens_binM<-as.numeric(Ycen_binM_noint<cen_binM)
delta_binM<-as.numeric(Y_surv_binM_noint>cen_binM)
cens_contM<-as.numeric(Ycen_contM_noint<cen_contM)
delta_contM<-as.numeric(Y_surv_contM_noint>cen_contM)
cens_catM<-as.numeric(Ycen_catM_noint<cen_catM)
delta_catM<-as.numeric(Y_surv_catM_noint>cen_catM)

#generate Cont Y, Bin Y, Count Y without int
linpred = theta0+theta1*A+theta21*M_bin+theta4*C
contY_binM_noint = rnorm(n,linpred,0.2)
binY_binM_noint = rbinom(n,1,expit(linpred))
countY_binM_noint = rpois(n,abs(linpred))

linpred = theta0+theta1*A+theta21*M_cont+theta4*C
contY_contM_noint = rnorm(n,linpred,0.2)
binY_contM_noint = rbinom(n,1,expit(linpred))
countY_contM_noint = rpois(n,abs(linpred))

linpred = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta4*C
contY_catM_noint = rnorm(n,linpred,0.2)
binY_catM_noint = rbinom(n,1,expit(linpred))
countY_catM_noint = rpois(n,abs(linpred))

noint_data = cbind(A,
                   M_cont,M_bin,M_cat,
                   contY_binM_noint,contY_contM_noint,contY_catM_noint,
                   binY_binM_noint,binY_contM_noint,binY_catM_noint,
                   countY_binM_noint,countY_contM_noint,countY_catM_noint,
                   Ycen_binM_noint,Ycen_contM_noint,Ycen_catM_noint,
                   cens_binM,cens_contM,cens_catM,
                   delta_binM,delta_contM,delta_catM,
                   C)

write.table(noint_data,file="noint_data_10000.txt")

######################################################int################################################

#generate survival data
linpred = exp(-(theta0+theta1*A+theta21*M_bin+theta31*A*M_bin+theta4*C))
Y_surv_binM_int = rexp(n,linpred)
cen_binM<-quantile(Y_surv_binM_int,probs = 0.75)

linpred = exp(-(theta0+theta1*A+theta21*M_cont+theta31*A*M_cont+theta4*C))
Y_surv_contM_int = rexp(n,linpred)
cen_contM<-quantile(Y_surv_contM_int,probs = 0.75)

linpred = exp(-(theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+
                  theta31*A*(M_cat==1)+theta32*A*(M_cat==2)+theta4*C))
Y_surv_catM_int = rexp(n,linpred)
cen_catM<-quantile(Y_surv_catM_int,probs = 0.75)

#censor observations after censoring time
Ycen_binM_int<-pmin(Y_surv_binM_int,cen_binM)
Ycen_contM_int<-pmin(Y_surv_contM_int,cen_contM)
Ycen_catM_int<-pmin(Y_surv_catM_int,cen_catM)
#generate event indicator
cens_binM<-as.numeric(Ycen_binM_int<cen_binM)
delta_binM<-as.numeric(Y_surv_binM_int>cen_binM)
cens_contM<-as.numeric(Ycen_contM_int<cen_contM)
delta_contM<-as.numeric(Y_surv_contM_int>cen_contM)
cens_catM<-as.numeric(Ycen_catM_int<cen_catM)
delta_catM<-as.numeric(Y_surv_catM_int>cen_catM)

#generate Cont Y, Bin Y, Count Y without int
linpred = theta0+theta1*A+theta21*M_bin+theta31*A*M_bin+theta4*C
contY_binM_int = rnorm(n,linpred,0.2)
binY_binM_int = rbinom(n,1,expit(linpred))
countY_binM_int = rpois(n,abs(linpred))

linpred = theta0+theta1*A+theta21*M_cont+theta31*A*M_cont+theta4*C
contY_contM_int = rnorm(n,linpred,0.2)
binY_contM_int = rbinom(n,1,expit(linpred))
countY_contM_int = rpois(n,abs(linpred))

linpred = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+
  theta31*A*(M_cat==1)+theta32*A*(M_cat==2)+theta4*C
contY_catM_int = rnorm(n,linpred,0.2)
binY_catM_int = rbinom(n,1,expit(linpred))
countY_catM_int = rpois(n,abs(linpred))

int_data = cbind(A,
                   M_cont,M_bin,M_cat,
                   contY_binM_int,contY_contM_int,contY_catM_int,
                   binY_binM_int,binY_contM_int,binY_catM_int,
                   countY_binM_int,countY_contM_int,countY_catM_int,
                   Ycen_binM_int,Ycen_contM_int,Ycen_catM_int,
                   cens_binM,cens_contM,cens_catM,
                   delta_binM,delta_contM,delta_catM,
                   C)

write.table(int_data,file="int_data_10000.txt")


##############################################multiple M without int####################################################

#generate Cont Y, Bin Y, Count Y without int
linpred1 = theta0+theta1*A+theta21*M_cont+theta22*M_bin+theta4*C
linpred2 = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta23*M_bin+theta4*C
linpred3 = theta0+theta1*A+theta21*M_cont+theta22*(M_cat==1)+theta23*(M_cat==2)+theta4*C

contY_binMcontM_noint = rnorm(n,linpred1,0.2)
contY_catMbinM_noint = rnorm(n,linpred2,0.2)
contY_catMcontM_noint = rnorm(n,linpred3,0.2)
binY_binMcontM_noint = rbinom(n,1,expit(linpred1))
binY_catMbinM_noint = rbinom(n,1,expit(linpred2))
binY_catMcontM_noint = rbinom(n,1,expit(linpred3))
countY_binMcontM_noint = rpois(n,abs(linpred1))
countY_catMbinM_noint = rpois(n,abs(linpred2))
countY_catMcontM_noint = rpois(n,abs(linpred3))

multipleM_noint_data = cbind(A,
                             M_cont,M_bin,M_cat,
                             contY_binMcontM_noint,contY_catMbinM_noint,contY_catMcontM_noint,
                             binY_binMcontM_noint,binY_catMbinM_noint,binY_catMcontM_noint,
                             countY_binMcontM_noint,countY_catMbinM_noint,countY_catMcontM_noint,
                             C)

write.table(multipleM_noint_data,file="multipleM_noint_data_10000.txt")


##############################################multiple M with exposure-mediator int##########################

#generate Cont Y, Bin Y, Count Y with int
linpred1 = theta0+theta1*A+theta21*M_cont+theta22*M_bin+theta31*A*M_bin+theta4*C
linpred2 = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta23*M_bin+theta31*A*M_bin+theta4*C
linpred3 = theta0+theta1*A+theta21*M_cont+theta22*(M_cat==1)+theta23*(M_cat==2)+theta31*A*M_cont+theta4*C

contY_binMcontM_EMint = rnorm(n,linpred1,0.2)
contY_catMbinM_EMint = rnorm(n,linpred2,0.2)
contY_catMcontM_EMint = rnorm(n,linpred3,0.2)
binY_binMcontM_EMint = rbinom(n,1,expit(linpred1))
binY_catMbinM_EMint = rbinom(n,1,expit(linpred2))
binY_catMcontM_EMint = rbinom(n,1,expit(linpred3))
countY_binMcontM_EMint = rpois(n,abs(linpred1))
countY_catMbinM_EMint = rpois(n,abs(linpred2))
countY_catMcontM_EMint = rpois(n,abs(linpred3))

multipleM_EMint_data = cbind(A,
                             M_cont,M_bin,M_cat,
                             contY_binMcontM_EMint,contY_catMbinM_EMint,contY_catMcontM_EMint,
                             binY_binMcontM_EMint,binY_catMbinM_EMint,binY_catMcontM_EMint,
                             countY_binMcontM_EMint,countY_catMbinM_EMint,countY_catMcontM_EMint,
                             C)

write.table(multipleM_EMint_data,file="multipleM_EMint_data_10000.txt")


##############################multiple M with mediator-mediator int##########################

linpred1 = theta0+theta1*A+theta21*M_cont+theta22*M_bin+theta31*M_cont*M_bin+theta4*C
linpred2 = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta23*M_bin+theta31*(M_cat==1)*M_bin+theta32*(M_cat==2)*M_bin+theta4*C
linpred3 = theta0+theta1*A+theta21*M_cont+theta22*(M_cat==1)+theta23*(M_cat==2)+theta31*(M_cat==1)*M_cont+theta32*(M_cat==2)*M_cont+theta4*C

contY_binMcontM_MMint = rnorm(n,linpred1,0.2)
contY_catMbinM_MMint = rnorm(n,linpred2,0.2)
contY_catMcontM_MMint = rnorm(n,linpred3,0.2)
binY_binMcontM_MMint = rbinom(n,1,expit(linpred1))
binY_catMbinM_MMint = rbinom(n,1,expit(linpred2))
binY_catMcontM_MMint = rbinom(n,1,expit(linpred3))
countY_binMcontM_MMint = rpois(n,abs(linpred1))
countY_catMbinM_MMint = rpois(n,abs(linpred2))
countY_catMcontM_MMint = rpois(n,abs(linpred3))

multipleM_MMint_data = cbind(A,
                             M_cont,M_bin,M_cat,
                             contY_binMcontM_MMint,contY_catMbinM_MMint,contY_catMcontM_MMint,
                             binY_binMcontM_MMint,binY_catMbinM_MMint,binY_catMcontM_MMint,
                             countY_binMcontM_MMint,countY_catMbinM_MMint,countY_catMcontM_MMint,
                             C)

write.table(multipleM_MMint_data,file="multipleM_MMint_data_10000.txt")

##############################################multiple M with multiway int##########################

linpred1 = theta0+theta1*A+theta21*M_cont+theta22*M_bin+theta31*A*M_cont*M_bin+theta4*C
linpred2 = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta23*M_bin+theta31*A*M_cat*M_bin+theta4*C
linpred3 = theta0+theta1*A+theta21*M_cont+theta22*(M_cat==1)+theta23*(M_cat==2)+theta31*A*M_cat*M_cont+theta4*C

contY_binMcontM_EMMint = rnorm(n,linpred1,0.2)
contY_catMbinM_EMMint = rnorm(n,linpred2,0.2)
contY_catMcontM_EMMint = rnorm(n,linpred3,0.2)
binY_binMcontM_EMMint = rbinom(n,1,expit(linpred1))
binY_catMbinM_EMMint = rbinom(n,1,expit(linpred2))
binY_catMcontM_EMMint = rbinom(n,1,expit(linpred3))
countY_binMcontM_EMMint = rpois(n,abs(linpred1))
countY_catMbinM_EMMint = rpois(n,abs(linpred2))
countY_catMcontM_EMMint = rpois(n,abs(linpred3))

multipleM_EMMint_data = cbind(A,
                             M_cont,M_bin,M_cat,
                             contY_binMcontM_EMMint,contY_catMbinM_EMMint,contY_catMcontM_EMMint,
                             binY_binMcontM_EMMint,binY_catMbinM_EMMint,binY_catMcontM_EMMint,
                             countY_binMcontM_EMMint,countY_catMbinM_EMMint,countY_catMcontM_EMMint,
                             C)

write.table(multipleM_EMMint_data,file="multipleM_EMMint_data_10000.txt")

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

theta0 = -5
theta1 =0.8
theta21 = 1.8
theta22 = 1.2
theta23 = 1.5
theta24 = 1
theta31 = 0.2
theta32 = 0.4
theta33 = 0.1
theta34 = 0.3
theta4 = 0.1
theta5 = 0.3
theta6 = 0.4
theta7 = 0.2

meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}
n=10000

set.seed(100)
#generate binary exposure and confounder
A = rbinom(n,1,0.4)
C = rnorm(n,mean=1,sd=sqrt(1))

#generate Bin M, Cont M, and Cat M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)

probl1_bin = expit(linpred1)
probl0_cat = 1/(1+exp(linpred1)+exp(linpred2))
probl1_cat = exp(linpred1)/(1+exp(linpred1)+exp(linpred2))
probl2_cat = exp(linpred2)/(1+exp(linpred1)+exp(linpred2))

L_bin = rbinom(n,1,probl1_bin)
L_cat = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                          prob=c(probl0_cat[x],
                                                 probl1_cat[x],
                                                 probl2_cat[x])))

#generate Bin M, Cont M, and Cat M
linpred3 = (beta01+beta11*A+beta21*C+beta31*L_bin+beta41*(L_cat==1)+beta51*(L_cat==2))
linpred4 = (beta02+beta12*A+beta22*C+beta32*L_bin+beta42*(L_cat==1)+beta52*(L_cat==2))

probm1_bin = expit(linpred3)
probm0_cat = 1/(1+exp(linpred3)+exp(linpred4))
probm1_cat = exp(linpred3)/(1+exp(linpred3)+exp(linpred4))
probm2_cat = exp(linpred4)/(1+exp(linpred3)+exp(linpred4))

M_bin = rbinom(n,1,probm1_bin)
M_cont = rnorm(n,linpred3,0.1)
M_cat = sapply(1:n, FUN = function(x) sample(c(0,1,2),size=1,replace=TRUE,
                                             prob=c(probm0_cat[x],
                                                    probm1_cat[x],
                                                    probm2_cat[x])))

#generate Cont Y, Bin Y, Count Y with int
linpred1 = theta0+theta1*A+theta21*M_cont+theta22*M_bin+theta31*A*M_bin+theta4*C+theta5*L_bin+theta6*(L_cat==1)+theta7*(L_cat==2)
linpred2 = theta0+theta1*A+theta21*(M_cat==1)+theta22*(M_cat==2)+theta23*M_bin+theta31*A*M_bin+theta4*C+theta5*L_bin+theta6*(L_cat==1)+theta7*(L_cat==2)
linpred3 = theta0+theta1*A+theta21*M_cont+theta22*(M_cat==1)+theta23*(M_cat==2)+theta31*A*M_cont+theta4*C+theta5*L_bin+theta6*(L_cat==1)+theta7*(L_cat==2)

contY_binMcontM_EMint = rnorm(n,linpred1,0.2)
contY_catMbinM_EMint = rnorm(n,linpred2,0.2)
contY_catMcontM_EMint = rnorm(n,linpred3,0.2)
binY_binMcontM_EMint = rbinom(n,1,expit(linpred1))
binY_catMbinM_EMint = rbinom(n,1,expit(linpred2))
binY_catMcontM_EMint = rbinom(n,1,expit(linpred3))
countY_binMcontM_EMint = rpois(n,abs(linpred1))
countY_catMbinM_EMint = rpois(n,abs(linpred2))
countY_catMcontM_EMint = rpois(n,abs(linpred3))

multipleM_EMint_postcovar_data = cbind(A, L_bin, L_cat,
                             M_cont,M_bin,M_cat,
                             contY_binMcontM_EMint,contY_catMbinM_EMint,contY_catMcontM_EMint,
                             binY_binMcontM_EMint,binY_catMbinM_EMint,binY_catMcontM_EMint,
                             countY_binMcontM_EMint,countY_catMbinM_EMint,countY_catMcontM_EMint,
                             C)

write.table(multipleM_EMint_postcovar_data,file="multipleM_EMint_postcovar_data_10000.txt")
