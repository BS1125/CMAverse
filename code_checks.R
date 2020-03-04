rm(list=ls())

beta01 = -0.25
beta11 = 0.5
beta21 = 0.2

beta02 = -0.3
beta12 = 0.4
beta22 = 0.3

theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 1.2
theta4 = 0.6
theta5 = 0.4
theta6 = 0.3
theta7 = 0.5
theta8 = 0.1

meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}

n=10000

set.seed(100)

#generate binary exposure and confounder
A = rbinom(n,1,0.4)

C = rnorm(n,mean=1,sd=sqrt(1))
C_post = rnorm(n,mean=0,sd=sqrt(1))

#generate Bin M and Cont M
linpred1 = (beta01+beta11*A+beta21*C)
linpred2 = (beta02+beta12*A+beta22*C)
probm1 = expit(linpred1)
probm2 = expit(linpred2)
M_bin1 = rbinom(n,1,probm1)
M_bin2 = rbinom(n,1,probm2)
M_cont1 = rnorm(n,linpred1,0.1)
M_cont2 = rnorm(n,linpred2,0.1)

#generate Cont Y, Bin Y, Count Y without int
linpred1 = theta0+theta1*A+theta2*M_cont1+theta3*M_cont2+theta4*A*M_cont1+theta5*A*M_cont2+
  theta6*M_cont1*M_cont2+theta7*A*M_cont1*M_cont2+theta8*C
linpred2 = theta0+theta1*A+theta2*M_bin1+theta3*M_cont1+theta4*A*M_bin1+theta5*A*M_cont1+
  theta6*M_bin1*M_cont1+theta7*A*M_bin1*M_cont1+theta8*C
linpred3 = theta0+theta1*A+theta2*M_bin1+theta3*M_bin2+theta4*A*M_bin1+theta5*A*M_bin2+
  theta6*M_bin1*M_bin2+theta7*A*M_bin1*M_bin2+theta8*C
contY_2contM_mint = rnorm(n,linpred1,0.2)
contY_binMcontM_mint = rnorm(n,linpred2,0.2)
contY_2binM_mint = rnorm(n,linpred3,0.2)
binY_2contM_mint = rbinom(n,1,expit(linpred1))
binY_binMcontM_mint = rbinom(n,1,expit(linpred2))
binY_2binM_mint = rbinom(n,1,expit(linpred3))

data = as.data.frame(cbind(A,M_cont1,M_cont2,M_bin1,M_bin2,contY_2contM_mint,
                            contY_binMcontM_mint,contY_2binM_mint,binY_2contM_mint,
                            binY_binMcontM_mint,binY_2binM_mint,C,C_post))



exposure="A"
exposure.type="binary"
outcome="contY_binMcontM_mint"
yreg="linear"
mediator="M_cont1"
mreg="linear"
covariates.pre="C"
covariates.post="C_post"
EMint=TRUE
model="msm"
a=1
a_star=0
m_star=1
nsims=100
