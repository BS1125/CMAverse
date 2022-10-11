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
             "$E[Y_{a^*M_a}-Y_{a^*M_a^*}]$", "$E[Y_{aM_a}-Y_{aM_a^*}]$", "$PNDE+TNIE$ or $TNDE+PNIE$", 
             "$PNDE-CDE$", "$TNIE-PNIE$", "$CDE/TE$", 
             "$INT_{ref}/TE$", "$INT_{med}/TE$", "$PNIE/TE$", "$TNIE/TE$", 
             "$(INT_{ref}+INT_{med})/TE$", 
             "$(INT_{ref}+INT_{med}+PNIE)/TE$"))
colnames(tb1) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb1, escape = FALSE, caption = "Table 1: Causal Effects on the Difference Scale") %>%
  footnote(general = "$a$ and $a^*$ are the active and control values for $A$ respectively. $m$ is the value at which $M$ is controlled. $M_a$ denotes the counterfactual value of $M$ that would have been observed had $A$ been set to be $a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aMa*}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be the counterfactual value $M_{a*}$.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive", "condensed"))

## ----echo=F,warning=F,message=F-----------------------------------------------
tb2 <- data.frame(c("Controlled Direct Effect", "Pure Natural Direct Effect",
                    "Total Natural Direct Effect", "Pure Natural Indirect Effect",
                    "Total Natural Indirect Effect", "Total Effect", 
                    "Excess Ratio due to Controlled Direct Effect", 
                    "Excess Ratio due to Reference Interaction", 
                    "Excess Ratio due to Mediated Interaction", 
                    "Excess Ratio due to Pure Natural Indirect Effect", 
                    "Proportion $ER^{CDE}$",
                    "Proportion $ER^{INT_{ref}}$", "Proportion $ER^{INT_{med}}$", 
                    "Proportion $ER^{PNIE}$", "Proportion Mediated",
                    "Proportion Attributable to Interaction", 
                    "Proportion Eliminated"),
           c("$R^{CDE}$", "$R^{PNDE}$", "$R^{TNDE}$", "$R^{PNIE}$", "$R^{TNIE}$", "$R^{TE}$", 
             "$ER^{CDE}$", "$ER^{INT_{ref}}$", "$ER^{INT_{med}}$", "$ER^{PNIE}$", 
             "$prop^{ER^{CDE}}$", 
             "$prop^{ER^{INT_{ref}}}$", "$prop^{ER^{INT_{med}}}$", 
             "$prop^{ER^{PNIE}}$", "$PM$", "$INT$", "$PE$"),
           c("$E[Y_{am}]/E[Y_{a^*m}]$", "$E[Y_{aM_a^*}]/E[Y_{a^*M_a^*}]$", 
             "$E[Y_{aM_a}]/E[Y_{a^*M_a}]$", "$E[Y_{a^*M_a}]/E[Y_{a^*M_a^*}]$", 
             "$E[Y_{aM_a}]/E[Y_{aM_a^*}]$", "$R^{PNDE}\\times R^{TNIE}$ or $R^{TNDE}\\times R^{PNIE}$", 
              "$(E[Y_{am}-Y_{a^*m}])/E[Y_{a^*M_a^*}]$",
             "$R^{PNDE}-1-ER^{CDE}$", "$R^{TNIE}*R^{PNDE}-R^{PNDE}-R^{PNIE}+1$", 
             "$R^{PNIE}-1$", "$ER^{CDE}/(R^{TE}-1)$", 
             "$ER^{INT_{ref}}/(R^{TE}-1)$", 
             "$ER^{INT_{med}}/(R^{TE}-1)$", "$ER^{PNIE}/(R^{TE}-1)$", 
             "$(R^{PNDE}*(R^{TNIE}-1))/(R^{TE}-1)$",
             "$(ER^{INT_{ref}}+ER^{INT_{med}})/(R^{TE}-1)$", 
             "$(ER^{INT_{ref}}+ER^{INT_{med}}+ER^{PNIE})/(R^{TE}-1)$"))
colnames(tb2) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb2, escape = FALSE, caption = "Table 2: Causal Effects on the Ratio Scale") %>%
  footnote(general = "$a$ and $a^*$ are the active and control values for $A$ respectively. $m$ is the value at which $M$ is controlled. $M_a$ denotes the counterfactual value of $M$ that would have been observed had $A$ been set to be $a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aMa*}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be the counterfactual value $M_{a*}$. If $Y$ is categorical, $E[Y]$ represents the probability of $Y=y$ where $y$ is a pre-specified value of $Y$.") %>%
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
             "$E[Y_{a^*G_{a}}-Y_{a^*G_{a^*}}]$", "$E[Y_{aG_{a}}-Y_{aG_{a^*}}]$", "$rPNDE+rTNIE$ or $rTNDE+rPNIE$", 
             "$rPNDE-CDE$", "$rTNIE-rPNIE$", "$CDE/TE$", 
             "$rINT_{ref}/TE$", "$rINT_{med}/TE$", "$rPNIE/TE$", "$rTNIE/TE$", 
             "$(rINT_{ref}+rINT_{med})/TE$", 
             "$(rINT_{ref}+rINT_{med}+rPNIE)/TE$"))
colnames(tb3) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb3, escape = FALSE, caption = "Table 3: Causal Effects on the Difference Scale") %>%
  footnote(general = "$a$ and $a^*$ are the active and control values for $A$. $m$ is the value at which $M$ is controlled. $G_{a}$ denotes a random draw from the distribution of $M$ had $A=a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aG_{a*}}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be the counterfactual value $G_{a*}$.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive", "condensed"))

## ----echo=F,warning=F,message=F-----------------------------------------------
tb4 <- data.frame(c("Controlled Direct Effect", "Randomized Analogue of $PNDE$",
                    "Randomized Analogue of $TNDE$", "Randomized Analogue of $PNIE$",
                    "Randomized Analogue of $TNIE$", "Total Effect", 
                    "Excess Ratio due to Controlled Direct Effect", 
                    "Randomized Analogue of $ER^{INT_{ref}}$", 
                    "Randomized Analogue of $ER^{INT_{med}}$", 
                    "Randomized Analogue of $ER^{PNIE}$", 
                    "Proportion $ER^{CDE}$",
                    "Proportion $rER^{INT_{ref}}$", "Proportion $rER^{INT_{med}}$", 
                    "Proportion $rER^{PNIE}$", "Randomized Analogue of $PM$", "Randomized Analogue of $INT$", 
                    "Randomized Analogue of $PE$"),
           c("$R^{CDE}$", "$rR^{PNDE}$", "$rR^{TNDE}$", "$rR^{PNIE}$", "$rR^{TNIE}$", "$R^{TE}$", 
             "$ER^{CDE}$", "$rER^{INT_{ref}}$", "$rER^{INT_{med}}$", "$rER^{PNIE}$", 
             "$prop^{ER^{CDE}}$", 
             "$prop^{rER^{INT_{ref}}}$", "$prop^{rER^{INT_{med}}}$", 
             "$prop^{rER^{PNIE}}$", "$rPM$", "$rINT$", "$rPE$"),
           c("$E[Y_{am}]/E[Y_{a^*m}]$", "$E[Y_{aG_{a^*}}]/E[Y_{a^*G_{a^*}}]$", 
             "$E[Y_{aG_{a}}]/E[Y_{a^*G_{a}}]$", "$E[Y_{a^*G_{a}}]/E[Y_{a^*G_{a^*}}]$", 
             "$E[Y_{aG_{a}}]/E[Y_{aG_{a^*}}]$", "$rR^{PNDE}\\times rR^{TNIE}$ or $rR^{TNDE}\\times rR^{PNIE}$", 
              "$(E[Y_{am}-Y_{a^*m}])/E[Y_{a^*M_a^*}]$",
             "$rR^{PNDE}-1-ER^{CDE}$", "$rR^{TNIE}*rR^{PNDE}-rR^{PNDE}-rR^{PNIE}+1$", 
             "$rR^{PNIE}-1$", "$ER^{CDE}/(rR^{TE}-1)$", 
             "$rER^{INT_{ref}}/(R^{TE}-1)$", 
             "$rER^{INT_{med}}/(R^{TE}-1)$", "$rER^{PNIE}/(R^{TE}-1)$", 
             "$(rR^{PNDE}*(rR^{TNIE}-1))/(R^{TE}-1)$",
             "$(rER^{INT_{ref}}+rER^{INT_{med}})/(R^{TE}-1)$", 
             "$(rER^{INT_{ref}}+rER^{INT_{med}}+rER^{PNIE})/(R^{TE}-1)$"))
colnames(tb4) <- c("Full Name", "Abbreviation", "Formula")
knitr::kable(tb4, escape = FALSE, caption = "Table 4: Causal Effects on the Ratio Scale") %>%
  footnote(general = "$a$ and $a^*$ are the active and control values for $A$. $m$ is the value at which $M$ is controlled. $G_{a}$ denotes a random draw from the distribution of $M$ among those with $A=a$. $Y_{am}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be $m$. $Y_{aG_{a*}}$ denotes the counterfactual value of $Y$ that would have been observed had $A$ been set to be $a$, and $M$ to be the counterfactual value $G_{a*}$. If $Y$ is categorical, $E[Y]$ represents the probability of $Y=y$ where $y$ is a pre-specified value of $Y$.") %>%
  column_spec(1, width = "20em") %>%
  column_spec(2, width = "8em") %>%
  column_spec(3, width = "10em") %>%
  kable_styling(bootstrap_options = c("bordered", "striped", "hover", "responsive"))

