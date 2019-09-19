beta0 = -0.25
beta1 = 0.5
beta2 = 0.2
theta0 = -5
theta1 =0.8
theta2 = 1.8
theta3 = 0.2
theta4 = 0.1
m<-0
a<-1
a_star<-0
meanc<-1
expit<-function(x){
  exp(x)/(1+exp(x))
}

#MCONT YCONT

cde_contcont<-(theta1+theta3*m)*(a-a_star)
pnde_contcont<-(theta1+theta3*(beta0+beta1*a_star+beta2*meanc))*(a-a_star)
tnde_contcont<-(theta1+theta3*(beta0+beta1*a+beta2*meanc))*(a-a_star)
pnie_contcont<-(theta2*beta1+theta3*beta1*a_star)*(a-a_star)
tnie_contcont<-(theta2*beta1+theta3*beta1*a)*(a-a_star)
te_contcont<-pnde_contcont+tnie_contcont
pm_contcont <- tnie_contcont/(pnde_contcont+te_contcont)
intref_contcont <- pnde_contcont - cde_contcont
intmed_contcont <- tnie_contcont - pnie_contcont
pie_contcont <- pnie_contcont
cde_prop_contcont <- cde_contcont/te_contcont
intref_prop_contcont <- intref_contcont/te_contcont
intmed_prop_contcont <- intmed_contcont/te_contcont
pie_prop_contcont <- pie_contcont/te_contcont
overall_pm_contcont<- (pie_contcont+intmed_contcont)/te_contcont
overall_int_contcont <- (intref_contcont+intmed_contcont)/te_contcont
overall_pe_contcont <- (intref_contcont+intmed_contcont+pie_contcont)/te_contcont
cde_contcont
pnde_contcont
tnde_contcont
pnie_contcont
tnie_contcont
te_contcont
pm_contcont
intref_contcont
intmed_contcont
pie_contcont
cde_prop_contcont
intref_prop_contcont
intmed_prop_contcont
pie_prop_contcont
overall_pm_contcont
overall_int_contcont
overall_pe_contcont

#MBIN YCONT

cde_bincont<-(theta1+theta3*m)*(a-a_star)
pnde_bincont<-theta1*(a-a_star) + theta3*(a-a_star)*expit(beta0+beta1*a_star+beta2*meanc)
tnde_bincont<-theta1*(a-a_star) + theta3*(a-a_star)*expit(beta0+beta1*a+beta2*meanc)
pnie_bincont<-(theta2+theta3*a_star)*(expit(beta0+beta1*a+beta2*meanc)-expit(beta0+beta1*a_star+beta2*meanc))
tnie_bincont<-(theta2+theta3*a)*(expit(beta0+beta1*a+beta2*meanc)-expit(beta0+beta1*a_star+beta2*meanc))
te_bincont<-pnde_bincont+tnie_bincont
pm_bincont <- tnie_bincont/(pnde_bincont+te_bincont)
intref_bincont <- pnde_bincont - cde_bincont
intmed_bincont <- tnie_bincont - pnie_bincont
pie_bincont <- pnie_bincont
cde_prop_bincont <- cde_bincont/te_bincont
intref_prop_bincont <- intref_bincont/te_bincont
intmed_prop_bincont <- intmed_bincont/te_bincont
pie_prop_bincont <- pie_bincont/te_bincont
overall_pm_bincont<- (pie_bincont+intmed_bincont)/te_bincont
overall_int_bincont <- (intref_bincont+intmed_bincont)/te_bincont
overall_pe_bincont <- (intref_bincont+intmed_bincont+pie_bincont)/te_bincont
cde_bincont
pnde_bincont
tnde_bincont
pnie_bincont
tnie_bincont
te_bincont
pm_bincont
intref_bincont
intmed_bincont
pie_bincont
cde_prop_bincont
intref_prop_bincont
intmed_prop_bincont
pie_prop_bincont
overall_pm_bincont
overall_int_bincont
overall_pe_bincont

#MCONT YBIN
cde_contbin<-exp((theta1+theta3*m)*(a-a_star))
pnde_contbin<-exp((theta1+theta3*(beta0+beta1*a_star+beta2*meanc+theta2))*(a-a_star)+0.5*theta3^2*(a-a_star)^2)
tnde_contbin<-exp((theta1+theta3*(beta0+beta1*a+beta2*meanc+theta2))*(a-a_star)+0.5*theta3^2*(a-a_star)^2)
pnie_contbin<-exp((theta2*beta1+theta3*beta1*a_star)*(a-a_star))
tnie_contbin<-exp((theta2*beta1+theta3*beta1*a)*(a-a_star))
cde_contbin
pnde_contbin
tnde_contbin
pnie_contbin
tnie_contbin
te_contbin<-pnde_contbin*tnie_contbin
pm_contbin <- (pnde_contbin * (tnie_contbin - 1)) / (pnde_contbin * tnie_contbin - 1)
cde_err_contbin <- exp(theta1*(a-a_star)+theta2*m_star+
             theta3*a*m_star- (theta2+theta3*a_star)*
             (beta0+beta1*a_star+beta2)-
             0.5*(theta2+theta3*a_star)^2)-
         exp(theta2*m_star+theta3*a_star*m_star-
               (theta2+theta3*a_star)*(beta0+beta1*a_star+beta2)-
               0.5*(theta2+theta3*a_star)^2)
intref_contbin <- pnde_contbin - 1 - cde_err_contbin
intmed_contbin <- tnie_contbin * pnde_contbin - pnde_contbin - pnie_contbin + 1
pie_contbin <- pnie_contbin - 1
total_err_contbin <- te_contbin - 1
cde_err_prop_contbin <- cde_err_contbin/total_err_contbin
intmed_err_prop_contbin <- intmed_contbin/total_err_contbin
intref_err_prop_contbin <- intref_contbin/total_err_contbin
pie_err_prop_contbin <- pie_contbin/total_err_contbin
overall_pm_contbin <- (pie_contbin+intmed_contbin)/total_err_contbin
overall_int_contbin <- (intref_contbin+intmed_contbin)/total_err_contbin
overall_pe_contbin <- (intref_contbin+intmed_contbin+pie_contbin)/total_err_contbin
cde_contbin
pnde_contbin
tnde_contbin
pnie_contbin
tnie_contbin
te_contbin
pm_contbin
cde_err_contbin
intref_contbin
intmed_contbin
pie_contbin
total_err_contbin
cde_err_prop_contbin
intmed_err_prop_contbin
intref_err_prop_contbin
pie_err_prop_contbin
overall_pm_contbin
overall_int_contbin
overall_pe_contbin

#MBIN YBIN

cde_binbin<-exp((theta1+theta3*m)* (a - a_star))
pnde_binbin<-(exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a_star+beta2*meanc)))/
  (exp(theta1*a_star)*(1+exp(theta2+theta3*a_star+beta0+beta1*a_star+beta2*meanc)))
tnde_binbin<-(exp(theta1*a)*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc)))/
  (exp(theta1*a_star)*(1+exp(theta2+theta3*a_star+beta0+beta1*a+beta2*meanc)))
pnie_binbin<-((1+exp(beta0+beta1*a_star+beta2*meanc))*(1+exp(theta2+theta3*a_star+beta0+beta1*a+beta2*meanc)))/
  ((1+exp(beta0+beta1*a+beta2*meanc))*(1+exp(theta2+theta3*a_star+beta0+beta1*a_star+beta2*meanc)))
tnie_binbin<-((1+exp(beta0+beta1*a_star+beta2*meanc))*(1+exp(theta2+theta3*a+beta0+beta1*a+beta2*meanc)))/
  ((1+exp(beta0+beta1*a+beta2*meanc))*(1+exp(theta2+theta3*a+beta0+beta1*a_star+beta2*meanc)))
te_binbin<-pnde_binbin*tnie_binbin
pm_binbin <- (pnde_binbin * (tnie_binbin - 1)) / (pnde_binbin * tnie_binbin - 1)
cde_err_binbin <- (exp(theta1*(a-a_star)+theta2*m_star+
                  theta3*a*m_star)*(1+exp(beta0+
                  beta1*a_star+sum(beta2)))/(1+exp(beta0+
                  beta1*a_star+sum(beta2)+theta2+
                 theta3*a_star))-exp(theta2*m_star+
                 theta3*a_star*m_star)*(1+exp(beta0+
                 beta1*a_star+sum(beta2)))/(1+exp(beta0+
                 beta1*a_star+sum(beta2)+theta2+
                 theta3*a_star)))
intref_binbin <- pnde_binbin - 1 - cde_err_binbin
intmed_binbin <- tnie_binbin * pnde_binbin - pnde_binbin - pnie_binbin + 1
pie_binbin <- pnie_binbin - 1
total_err_binbin <- te_binbin - 1
cde_err_prop_binbin <- cde_err_binbin/total_err_binbin
intmed_err_prop_binbin <- intmed_binbin/total_err_binbin
intref_err_prop_binbin <- intref_binbin/total_err_binbin
pie_err_prop_binbin <- pie_binbin/total_err_binbin
overall_pm_binbin <- (pie_binbin+intmed_binbin)/total_err_binbin
overall_int_binbin <- (intref_binbin+intmed_binbin)/total_err_binbin
overall_pe_binbin <- (intref_binbin+intmed_binbin+pie_binbin)/total_err_binbin
cde_binbin
pnde_binbin
tnde_binbin
pnie_binbin
tnie_binbin
te_binbin
pm_binbin
cde_err_binbin
intref_binbin
intmed_binbin
pie_binbin
total_err_binbin
cde_err_prop_binbin
intmed_err_prop_binbin
intref_err_prop_binbin
pie_err_prop_binbin
overall_pm_binbin
overall_int_binbin
overall_pe_binbin
