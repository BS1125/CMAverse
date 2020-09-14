## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=F,warning=F,message=F,fig.width=6,fig.height=3----------------------
library(dagitty)
g1 <- dagitty( "dag {
    A -> Y
    A -> M -> Y
    C -> A
    C -> M
    C -> Y
}")
coordinates(g1) <- list(x=c(A=0,Y=4,M=2,C=0.2),y=c(A=0,Y=0,M=-0.8,C=-1))
plot(g1)

## ----echo=F,warning=F,message=F-----------------------------------------------
library(knitr)
library(dplyr)
library(kableExtra)
tb1 <- data.frame(c("Controlled Direct Effect", "Pure Natural Direct Effect",
                    "Total Natural Direct Effect", "Pure Natural Indirect Effect",
                    "Total Natural Indirect Effect", "Total Effect", 
                    "Reference Interaction", "Mediated Interaction", "Proportion $CDE$",
                    "Proportion $INT_{ref}$", "Proportion $INT_{med}$", "Proportion $PNIE$", "Proportion Mediated", "Proportion Attributable to Interaction", 
                    "Proportion Eliminated"),
           c("$CDE$", "$PNDE$", "$TNDE$", "$PNIE$", "$TNIE$", "$TE$",  "$INT_{ref}$", 
             "$INT_{med}$", "$prop^{CDE}$", "$prop^{INT_{ref}}$", "$prop^{INT_{med}}$", 
             "$prop^{PNIE}$", "$PM$", "$INT$", "$PE$"),
           c("$E[Y_{am}-Y_{a^*m}]$", "$E[Y_{aM_a^*}-Y_{a^*M_a^*}]$", "$E[Y_{aM_a}-Y_{a^*M_a}]$", 
             "$E[Y_{a^*M_a}-Y_{a^*M_a^*}]$", "$E[Y_{aM_a}-Y_{aM_a^*}]$", "$PNDE+TNIE$", 
             "$PNDE-CDE$", "$TNIE-PNIE$", "$CDE/TE$", 
             "$INT_{ref}/TE$", "$INT_{med}/TE$", "$PNIE/TE$", "$TNIE/TE$", 
             "$(INT_{ref}+INT_{med})/TE$", 
             "$(INT_{ref}+INT_{med}+PNIE)/TE$"))
colnames(tb1) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb1, escape = FALSE, caption = "Table 1: Causal Effects on the Difference Scale") %>%
  footnote(general = "$a$ and $a^*$ are two reference values for $A$. $m$ is the reference value for $M$. $M_a$ denotes the counterfactual value of $M$ that would have been observed had $A$ been set to be $a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aMa*}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $M_{a*}$.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive", "condensed"))

## ----echo=F,warning=F,message=F-----------------------------------------------
tb2 <- data.frame(c("Controlled Direct Effect", "Pure Natural Direct Effect",
                    "Total Natural Direct Effect", "Pure Natural Indirect Effect",
                    "Total Natural Indirect Effect", "Total Effect", 
                    "Excess RR due to Controlled Direct Effect", 
                    "Excess RR due to Reference Interaction", 
                    "Excess RR due to Mediated Interaction", 
                    "Excess RR due to Pure Natural Indirect Effect", 
                    "Proportion $err^{CDE}$",
                    "Proportion $err^{INT_{ref}}$", "Proportion $err^{INT_{med}}$", 
                    "Proportion $err^{PNIE}$", "Proportion Mediated",
                    "Proportion Attributable to Interaction", 
                    "Proportion Eliminated"),
           c("$rr^{CDE}$", "$rr^{PNDE}$", "$rr^{TNDE}$", "$rr^{PNIE}$", "$rr^{TNIE}$", "$rr^{TE}$", 
             "$err^{CDE}$", "$err^{INT_{ref}}$", "$err^{INT_{med}}$", "$err^{PNIE}$", 
             "$prop^{err^{CDE}}$", 
             "$prop^{err^{INT_{ref}}}$", "$prop^{err^{INT_{med}}}$", 
             "$prop^{err^{PNIE}}$", "$PM$", "$INT$", "$PE$"),
           c("$E[Y_{am}]/E[Y_{a^*m}]$", "$E[Y_{aM_a^*}]/E[Y_{a^*M_a^*}]$", 
             "$E[Y_{aM_a}]/E[Y_{a^*M_a}]$", "$E[Y_{a^*M_a}]/E[Y_{a^*M_a^*}]$", 
             "$E[Y_{aM_a}]/E[Y_{aM_a^*}]$", "$rr^{PNDE}\\times rr^{TNIE}$", 
              "$(E[Y_{am}-Y_{a^*m}])/E[Y_{a^*M_a^*}]$",
             "$rr^{PNDE}-1-err^{CDE}$", "$rr^{TNIE}*rr^{PNDE}-rr^{PNDE}-rr^{PNIE}+1$", 
             "$rr^{PNIE}-1$", "$err^{CDE}/(rr^{TE}-1)$", 
             "$err^{INT_{ref}}/(rr^{TE}-1)$", 
             "$err^{INT_{med}}/(rr^{TE}-1)$", "$err^{PNIE}/(rr^{TE}-1)$", 
             "$(rr^{PNDE}*(rr^{TNIE}-1))/(rr^{TE}-1)$",
             "$(err^{INT_{ref}}+err^{INT_{med}})/(rr^{TE}-1)$", 
             "$(err^{INT_{ref}}+err^{INT_{med}}+err^{PNIE})/(rr^{TE}-1)$"))
colnames(tb2) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb2, escape = FALSE, caption = "Table 2: Causal Effects on the Ratio Scale") %>%
  footnote(general = "$a$ and $a^*$ are two reference values for $A$. $m$ is the reference value for $M$. $M_a$ denotes the counterfactual value of $M$ that would have been observed had $A$ been set to be $a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aMa*}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $M_{a*}$. $RR$ represents risk or rate ratio. If $Y$ is categorical, a reference value of $Y$, denoted as $y$, is needed and E[Y] represents the probability of $Y=y$.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive"))

## ----echo=F,warning=F,message=F,fig.width=6,fig.height=3----------------------
g2 <- dagitty::dagitty( "dag {
    A -> Y
    A -> M -> Y
    C -> A
    C -> M
    C -> Y
    A -> L
    L -> Y
    L -> M
    C -> L
}")

dagitty::coordinates(g2) <- list(x=c(A=0,Y=4,M=2,C=0.2,L=2),y=c(A=0,Y=0,M=-0.8,C=-1,L=0.5))

plot(g2)

## ----echo=F,warning=F,message=F-----------------------------------------------
library(knitr)
library(dplyr)
library(kableExtra)
tb3 <- data.frame(c("Controlled Direct Effect", "Randomized Analogue of $PNDE$",
                    "Randomized Analogue of $TNDE$", "Randomized Analogue of $PNIE$",
                    "Randomized Analogue of $TNIE$", "Total Effect", 
                    "Randomized Analogue of $INT_{ref}$", "Randomized Analogue of $INT_{med}$", "Proportion $CDE$",
                    "Proportion $rINT_{ref}$", "Proportion $rINT_{med}$", "Proportion $rPNIE$", "Randomized Analogue of $PM$", "Randomized Analogue of $INT$", 
                    "Randomized Analogue of $PE$"),
           c("$CDE$", "$rPNDE$", "$rTNDE$", "$rPNIE$", "$rTNIE$", "$TE$",  "$rINT_{ref}$", 
             "$rINT_{med}$", "$prop^{CDE}$", "$prop^{rINT_{ref}}$", "$prop^{rINT_{med}}$", 
             "$prop^{rPNIE}$", "$rPM$", "$rINT$", "$rPE$"),
           c("$E[Y_{am}-Y_{a^*m}]$", "$E[Y_{aG_{a^*}}-Y_{a^*G_{a^*}}]$", "$E[Y_{aG_{a}}-Y_{a^*G_{a}}]$", 
             "$E[Y_{a^*G_{a}}-Y_{a^*G_{a^*}}]$", "$E[Y_{aG_{a}}-Y_{aG_{a^*}}]$", "$rPNDE+rTNIE$", 
             "$rPNDE-CDE$", "$rTNIE-rPNIE$", "$CDE/TE$", 
             "$rINT_{ref}/TE$", "$rINT_{med}/TE$", "$rPNIE/TE$", "$rTNIE/TE$", 
             "$(rINT_{ref}+rINT_{med})/TE$", 
             "$(rINT_{ref}+rINT_{med}+rPNIE)/TE$"))
colnames(tb3) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb3, escape = FALSE, caption = "Table 3: Causal Effects on the Difference Scale") %>%
  footnote(general = "$a$ and $a^*$ are two reference values for $A$. $m$ is the reference value for $M$. $G_{a}$ denotes a random draw from the distribution of $M$ among those with $A=a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aG_{a*}}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $G_{a*}$.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive", "condensed"))

## ----echo=F,warning=F,message=F-----------------------------------------------
tb4 <- data.frame(c("Controlled Direct Effect", "Randomized Analogue of $PNDE$",
                    "Randomized Analogue of $TNDE$", "Randomized Analogue of $PNIE$",
                    "Randomized Analogue of $TNIE$", "Total Effect", 
                    "Excess RR due to Controlled Direct Effect", 
                    "Randomized Analogue of $err^{INT_{ref}}$", 
                    "Randomized Analogue of $err^{INT_{med}}$", 
                    "Randomized Analogue of $err^{PNIE}$", 
                    "Proportion $err^{CDE}$",
                    "Proportion $rerr^{INT_{ref}}$", "Proportion $rerr^{INT_{med}}$", 
                    "Proportion $rerr^{PNIE}$", "Randomized Analogue of $PM$", "Randomized Analogue of $INT$", 
                    "Randomized Analogue of $PE$"),
           c("$rr^{CDE}$", "$rr^{rPNDE}$", "$rr^{rTNDE}$", "$rr^{rPNIE}$", "$rr^{rTNIE}$", "$rr^{TE}$", 
             "$err^{CDE}$", "$rerr^{INT_{ref}}$", "$rerr^{INT_{med}}$", "$rerr^{PNIE}$", 
             "$prop^{err^{CDE}}$", 
             "$prop^{rerr^{INT_{ref}}}$", "$prop^{rerr^{INT_{med}}}$", 
             "$prop^{rerr^{PNIE}}$", "$rPM$", "$rINT$", "$rPE$"),
           c("$E[Y_{am}]/E[Y_{a^*m}]$", "$E[Y_{aG_{a^*}}]/E[Y_{a^*G_{a^*}}]$", 
             "$E[Y_{aG_{a}}]/E[Y_{a^*G_{a}}]$", "$E[Y_{a^*G_{a}}]/E[Y_{a^*G_{a^*}}]$", 
             "$E[Y_{aG_{a}}]/E[Y_{aG_{a^*}}]$", "$rr^{rPNDE}\\times rr^{rTNIE}$", 
              "$(E[Y_{am}-Y_{a^*m}])/E[Y_{a^*M_a^*}]$",
             "$rr^{rPNDE}-1-err^{CDE}$", "$rr^{rTNIE}*rr^{rPNDE}-rr^{rPNDE}-rr^{rPNIE}+1$", 
             "$rr^{rPNIE}-1$", "$err^{CDE}/(rr^{TE}-1)$", 
             "$rerr^{INT_{ref}}/(rr^{TE}-1)$", 
             "$rerr^{INT_{med}}/(rr^{TE}-1)$", "$rerr^{PNIE}/(rr^{TE}-1)$", 
             "$(rr^{rPNDE}*(rr^{rTNIE}-1))/(rr^{TE}-1)$",
             "$(rerr^{INT_{ref}}+rerr^{INT_{med}})/(rr^{TE}-1)$", 
             "$(rerr^{INT_{ref}}+rerr^{INT_{med}}+rerr^{PNIE})/(rr^{TE}-1)$"))
colnames(tb4) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb4, escape = FALSE, caption = "Table 4: Causal Effects on the Ratio Scale") %>%
  footnote(general = "$a$ and $a^*$ are two reference values for $A$. $m$ is the reference value for $M$. $G_{a}$ denotes a random draw from the distribution of $M$ among those with $A=a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aG_{a*}}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $G_{a*}$. $RR$ represents risk or rate ratio.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive"))

