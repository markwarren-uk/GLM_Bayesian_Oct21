
#=======================================

# R code for Chapter 8 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Load packages
library(ggplot2)
library(GGally)
library(tidyverse)
library(mgcv)
library(lme4)
library(car)
library(plyr)
library(devtools)
library(ggpubr)
library(qqplotr)
library(geiger)
library(INLA)
library(brinla)
library(inlatools)
library(gridExtra) 
library(rlang)

#=======================================

# Import the data

# To import the data, change the working directory to the one
# in which the .csv file is saved, and use the read.csv function

seal <- read_csv(file = "seal.csv")

str(seal)

# 'data.frame':	29 obs. of  5 variables:
# $ id        : Factor w/ 29 levels "A20","A31","A33"
# $ sex       : int  1 1 2 1 2 2 2 2 1 2 ...
# $ age       : int  8 10 10 10 7 10 10  ...
# $ lean_mass : int  77 105 100 85 86    ...
# $ dive_dur  : num  816 1128 1026 870   ...

# We have 29 observations, each an individual seal

#=======================================

# The 9 steps to fitting a Bayesian GLM are:

# 1. State the question
# 2. Perform data exploration
# 3. Select a statistical model
# 4. Specify and justify priors
# 5. Fit the model
# 6. Obtain the posterior distribution
# 7. Conduct model checks
# 8. Interpret and present model output
# 9. Visualise the results

#=======================================

# 1. STATE THE QUESTION

#=======================================

# The aim of this study is to understand whether the dive duration of
# common seals:
# 1. Increases as a function of lean mass 
# 2. Decreases as a function of age
# 3. Is greater for males

#=======================================

# 2. PERFORM DATA EXPLORATION

#=======================================

# MISSING VALUES?
colSums(is.na(seal))

# id   sex   age   lean_mass  dive_dur 
# 0    0     0     0          0 
#No missing data


# OUTLIERS

# A multi-panel Cleveland dotplot to examine continuous variables
# First make a vector of variables of interest
# Var <- c("age", "lean_mass", "dive_dur")
# 
# And plot
# Fig. 8.1
# dotplot(as.matrix(as.matrix(seal[,Var])),
#         groups=FALSE,
#         strip = strip.custom(bg = 'white',
#                              par.strip.text = list(cex = 1.2)),
#         scales = list(x = list(relation = "free", draw = TRUE),
#                       y = list(relation = "free", draw = FALSE)),
#         col = 1, cex  = 1, pch = 16,
#         xlab = list(label = "Data range", cex = 1.2),
#         ylab = list(label = "Data order", cex = 1.2))
# No obvious outliers in continuous variables

seal <- seal %>%
  mutate(order = seq(1:nrow(seal)))

# Set preferred format
My_theme <- theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, size = 1),
                  strip.background = element_rect(fill = "white", 
                                                  color = "white", size = 1),
                  text = element_text(size = 14),
                  panel.grid.major = element_line(colour = "white", size = 0.1),
                  panel.grid.minor = element_line(colour = "white", size = 0.1))

#Write a function
multi_dotplot <- function(filename, Xvar, Yvar){
  filename %>% 
    ggplot(aes(x = {{Xvar}})) +
    #Note the double curly brackets from rlang package - it allows us to pass unquoted col names
    geom_point(aes(y = {{Yvar}})) +
    theme_bw() +
    My_theme +
    coord_flip() +
    labs(x = "Order of Data")
  
}

#Choose variables

p1 <- multi_dotplot(seal, order, age) 
p2 <- multi_dotplot(seal, order, lean_mass) 
p3 <- multi_dotplot(seal, order, dive_dur)

#Create grip and plot
grid.arrange(p1, p2, p3, nrow = 1)


#NORMALITY AND HOMOGENEITY OF DEPENDENT VARIABLE

#Fig. 8.2 Density plot of dive duration

seal %>% 
  ggplot(aes(dive_dur)) +
  geom_density() +
  xlab("Dive duration (sec.)") + ylab("Density") +
  xlim(500,1500) +
  My_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))

# A positive skew


# BALANCE OF CATEGORICAL VARIABLES

# Sex is coded numerically
# and should be a factor not an integer
seal$fSex <- factor(seal$sex)

# #Rename levels
# levels(seal$fSex)[levels(seal$fSex)=="1"] <- "Male"
# levels(seal$fSex)[levels(seal$fSex)=="2"] <- "Female"

seal$fSex <- fct_recode(seal$fSex, Male = "1", Female = "2")

#Check balance
table(seal$fSex)
# Male Female 
#   15     14
#Good balance


# COLLINEARITY

# Coll <- c("fSex", "age", "lean_mass")
# 
# # Fig. 8.3
# # Obtain summary using the ggpairs command from the GGally library
# ggpairs(seal[,Coll], ggplot2::aes(colour=fSex, alpha = 0.8))
# # Lean mass and age are correlated

seal %>% 
    ggpairs(columns = c("fSex", "age", "lean_mass"), aes(colour=fSex, alpha = 0.8), lower = list(combo = wrap("facethist", binwidth = 15))) + My_theme


#Calculate Variance Inflation Factor (VIF) for a frequentist model
round(vif(gamma1 <- glm(dive_dur ~ lean_mass + age + fSex,
                                   family = Gamma(link = log),
                                   data = seal)),2)

# lean_mass  age      fSex 
# 1.90      1.88      1.10 
# No serious variance inflation


# ZEROS IN THE RESPONSE VARIABLE

sum(seal$dive_dur == 0)
#no zeros


# RELATIONSHIPS

# Plot lean mass
mass_plot <- seal %>% 
  ggplot() +
  geom_point(aes(y = dive_dur, x = lean_mass, size = 1,
                 alpha = 0.8)) +
  geom_smooth(method = "lm", se = FALSE, 
              aes(y = dive_dur, x = lean_mass)) + 
  xlab("Lean mass (kg)") + ylab("Dive duration (sec.)") +
  xlim(58,122) + ylim(700,1400) +
  theme(text = element_text(size=15))  +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  facet_grid(. ~ fSex,
             scales = "fixed", space = "fixed") +
  theme(legend.position = "none")
mass_plot

# 1. A positive effect of lean mass on dive duration
# 2. Possible effect of sex on dive duration

# Plot age
age_plot <- seal %>% 
  ggplot() +
  geom_point(aes(y = dive_dur, x = age, size = 1,
                 alpha = 0.8)) +
  geom_smooth(method = "lm", se = FALSE, 
              aes(y = dive_dur, x = age)) + 
  xlab("Age (years)") + ylab("") +
  xlim(4,12) + ylim(700,1400) +
  theme(text = element_text(size=15))  +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  facet_grid(. ~ fSex, 
             scales = "fixed", space = "fixed") +
  theme(legend.position = "none")
age_plot

# 1. A positive effect of age on dive duration 
# (though mass correlates with age)
# 2. Possible effect of sex on dive duration

# Combine plots (Fig 8.4)
ggarrange(mass_plot, age_plot,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

#=======================================

# 3. SELECT A STATISTICAL MODEL

#=======================================

# An option is to fit a Gaussian GLM
# However, the response variable is strictly positive, positively skewed
# and with no zeros, and an alternative is a gamma distribution

#=======================================

# 4. SPECIFY PRIORS

#=======================================

# Fit one model with default priors and second with weakly informative
# priors based on previous studies on elephant and Weddell seals 

# Previous research in other seal species has shown:
# 1. A weakly positive effect of lean mass on dive duration
# 2. A weakly negative effect of age on dive duration
# 3. A marginally greater dive duration in males than females

# Formulate priors:
# 1. (B1~N(5,25)) tau = 0.04
# 2. (B2~N(0.01,0.0025)) tau = 400
# 3. (B3~N(-0.01,0.0025)) tau = 400
# 3. (B4~N(-0.1,0.25)) tau = 4

# ======================================

# 5. FIT MODELS

#=======================================

# Model M01 with default priors
M01 <- inla(dive_dur ~ lean_mass + age + fSex,
                       family = "gamma", data = seal,
                       control.compute = list(dic = TRUE))

# Model I01 with informative priors

# Model I01 is fitted with weakly informative priors, including on hyperparameter

# The hyperparameter is the precision parameter phi, which is represented 
# as phi = exp(theta). The prior is defined on theta in INLA. So:

# Obtain estimate of phi from M01
phi <- M01$summary.hyperpar["Precision parameter for the Gamma observations", "mean"]
round(log(phi),2)
# 5.43

# Weakly informative priors on fixed effects and hyperparameter
I01 <- inla(dive_dur ~ lean_mass + age + fSex, 
                       data = seal, family = "gamma",
                       control.compute = list(dic = TRUE),
                       control.family = list(hyper =
                            list(prec = list(prior = "loggamma", 
                                             param = c(5.43, 31.6^(-2))))),
              control.fixed = list(mean.intercept = 5,
                                   prec.intercept = 5^(-2),
                            mean = list(lean_mass = 0.01, 
                                              age = -0.01,
                                       fSexFemale = -0.1), 
                            prec = list(lean_mass = 0.05^(-2), 
                                              age = 0.05^(-2),
                                       fSexFemale = 0.5^(-2))))

#=======================================

# 6. OBTAIN THE POSTERIOR DISTRIBUTION

#=======================================

# Plot posterior distributions for fixed parameters 

# Output for the fixed effects of M01

M01Betas <- I01$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(M01Betas, digits = 3)

#              mean    sd  0.025quant 0.975quant
# (Intercept)  6.161 0.060      6.043      6.279
# lean_mass    0.008 0.001      0.007      0.010
# age         -0.002 0.008     -0.018      0.013
# fSexFemale  -0.069 0.021     -0.111     -0.027

# Model intercept (Beta1)
PosteriorBeta1.M01 <- as.data.frame(M01$marginals.fixed$`(Intercept)`)
PriorBeta1.M01     <- data.frame(x = PosteriorBeta1.M01[,"x"],
                                 y = dnorm(PosteriorBeta1.M01[,"x"],0,0))
Beta1mean.M01 <- M01Betas["(Intercept)", "mean"]
Beta1lo.M01   <- M01Betas["(Intercept)", "0.025quant"]
Beta1up.M01   <- M01Betas["(Intercept)", "0.975quant"]

beta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.M01, xmax = Beta1up.M01,
           ymin = 0, ymax = 6, fill = "gray88") +
  geom_line(data = PosteriorBeta1.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta1.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") + ylab("Density") +
  xlim(5.8,6.5) + ylim(0,6) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# beta1

# lean mass (Beta2)
PosteriorBeta2.M01 <- as.data.frame(M01$marginals.fixed$`lean_mass`)
PriorBeta2.M01 <- data.frame(x = PosteriorBeta2.M01[,"x"], 
                  y = dnorm(PosteriorBeta2.M01[,"x"],0,31.6))
Beta2mean.M01 <- M01Betas["lean_mass", "mean"]
Beta2lo.M01   <- M01Betas["lean_mass", "0.025quant"]
Beta2up.M01   <- M01Betas["lean_mass", "0.975quant"]

beta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.M01, xmax = Beta2up.M01,
           ymin = 0, ymax = 400, fill = "gray88") +
  geom_line(data = PosteriorBeta2.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for lean mass") + ylab("Density") +
  xlim(0.0025,0.015) + ylim(0,400) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta2mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# beta2

# age (Beta3)
PosteriorBeta3.M01 <- as.data.frame(M01$marginals.fixed$`age`)
PriorBeta3.M01     <- data.frame(x = PosteriorBeta3.M01[,"x"],
                                 y = dnorm(PosteriorBeta3.M01[,"x"],0,31.6))
Beta3mean.M01 <- M01Betas["age", "mean"]
Beta3lo.M01   <- M01Betas["age", "0.025quant"]
Beta3up.M01   <- M01Betas["age", "0.975quant"]

beta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.M01, xmax = Beta3up.M01,
           ymin = 0, ymax = 45, fill = "gray88") +
  geom_line(data = PosteriorBeta3.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for age") + ylab("Density") +
  xlim(-0.06,0.06) + ylim(0,45) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# beta3

# fSexFemale (beta4)
PosteriorBeta4.M01 <- as.data.frame(M01$marginals.fixed$`fSexFemale`)
PriorBeta4.M01     <- data.frame(x = PosteriorBeta4.M01[,"x"],
                                 y = dnorm(PosteriorBeta4.M01[,"x"],0,31.6))
Beta4mean.M01 <- M01Betas["fSexFemale", "mean"]
Beta4lo.M01   <- M01Betas["fSexFemale", "0.025quant"]
Beta4up.M01   <- M01Betas["fSexFemale", "0.975quant"]

beta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.M01, xmax = Beta4up.M01,
           ymin = 0, ymax = 16, fill = "gray88") +
  geom_line(data = PosteriorBeta4.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Sex (female)") + ylab("Density") +
  xlim(-0.2,0.075) + ylim(0,16) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta4mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# beta4

# Combine plots (Fig 8.5)
ggarrange(beta1, beta2, beta3, beta4,
                        labels = c("A", "B", "C", "D"),
                        ncol = 2, nrow = 2)


#=======================================
# Model hyperparameter

# The model contains a hyperparameter for the precision parameter

# Obtain posterior distribution of precision with
M01hyp <- M01$summary.hyper[,c("mean", "sd", "mode", "0.025quant", "0.975quant")] 
#round(M01hyp, 1)

#                         mean  sd   mode  0.025quant 0.975quant
# Precision for Gamma obs 228.2 61.9 211.2      123.3      365.1

# Plot posterior distribution of precision
PosteriorHyp.M01 <- as.data.frame(M01$marginals.hyperpar$
                                    `Precision parameter for the Gamma observations`)

PriorHyp.M01 <- data.frame(x = PosteriorHyp.M01[,"x"], 
                           y = dgamma(PosteriorHyp.M01[,"x"],1,0.01))

Hypmean.M01 <- M01hyp["Precision parameter for the Gamma observations", "mode"]
Hyplo.M01   <- M01hyp["Precision parameter for the Gamma observations", "0.025quant"]
Hypup.M01   <- M01hyp["Precision parameter for the Gamma observations", "0.975quant"]

prec.dhyp <- ggplot() +
  annotate("rect", xmin = Hyplo.M01, xmax = Hypup.M01,
           ymin = 0, ymax = 0.007, fill = "gray88") +
  geom_line(data = PosteriorHyp.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorHyp.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  ylab("Density") + xlab("Precision") +
  xlim(-10,500) + ylim(0,0.007) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Hypmean.M01, linetype = "dashed") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# prec.dhyp

# We can obtain the posterior for dispersion using `bri.hyperpar.summary`
M01var <- bri.hyperpar.summary(M01)[,c("mean", "mode", "q0.025", "q0.975")]
#round(M01var,4)

# mean   mode q0.025 q0.975 
# 0.0681 0.0650 0.0524 0.0899 

# And plot the posterior distribution
Hypvmean.M01 <- M01var["mode"]
Hypvlo.M01   <- M01var["q0.025"]
Hypvup.M01   <- M01var["q0.975"]

TauM01   <- M01$marginals.hyperpar$`Precision parameter for the Gamma observations`
dispM01 <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), TauM01))
PriorVar.M01 <- data.frame(x = dispM01[,"x"], 
                           y = dgamma(dispM01[,"x"],1,0.01))

dis.dhyp <- ggplot() +
  annotate("rect", xmin = Hypvlo.M01, xmax = Hypvup.M01,
           ymin = 0, ymax = 46, fill = "gray88") +
  geom_line(data = dispM01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorVar.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  ylab("Density") + xlab("Dispersion") +
  xlim(0.04,0.12) + ylim(0,46) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Hypvmean.M01, linetype = "dashed") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# dis.dhyp

# Combine plots
ggarrange(prec.dhyp, dis.dhyp,
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1)

#=======================================

# B. Informative model (I01)

# Posterior mean values and 95% CI for fixed effects (informative)
I01Betas <- I01$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(I01Betas, digits = 3)
#               mean    sd 0.025quant 0.975quant
# (Intercept)   6.161 0.060      6.043      6.279
# lean_mass     0.008 0.001      0.007      0.010
# age          -0.002 0.008     -0.018      0.013
# fSexFemale   -0.069 0.021     -0.111     -0.027

# The first column shows the posterior mean, the second the
# standard deviation and the third and fourth columns show the 
# 2.5% and 97.5% quantiles of the posterior distribution.

# Plot posterior distributions for fixed parameters 

# Model intercept (Beta1)
PosteriorBeta1.I01 <- as.data.frame(I01$marginals.fixed$`(Intercept)`)
PriorBeta1.I01     <- data.frame(x = PosteriorBeta1.I01[,"x"],
                                 y = dnorm(PosteriorBeta1.I01[,"x"],5,5))
Beta1mean.I01 <- I01Betas["(Intercept)", "mean"]
Beta1lo.I01   <- I01Betas["(Intercept)", "0.025quant"]
Beta1up.I01   <- I01Betas["(Intercept)", "0.975quant"]

Ibeta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.I01, xmax = Beta1up.I01,
           ymin = 0, ymax = 7, fill = "gray88") +
  geom_line(data = PosteriorBeta1.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta1.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") + ylab("Density") +
  xlim(5.8,6.5) + ylim(0,7) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta1

# mass (Beta2)
PosteriorBeta2.I01 <- as.data.frame(I01$marginals.fixed$`lean_mass`)
PriorBeta2.I01 <- data.frame(x = PosteriorBeta2.I01[,"x"], 
                             y = dnorm(PosteriorBeta2.I01[,"x"],0.01,0.05))
Beta2mean.I01 <- I01Betas["lean_mass", "mean"]
Beta2lo.I01   <- I01Betas["lean_mass", "0.025quant"]
Beta2up.I01   <- I01Betas["lean_mass", "0.975quant"]

Ibeta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.I01, xmax = Beta2up.I01,
           ymin = 0, ymax = 480, fill = "gray88") +
  geom_line(data = PosteriorBeta2.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for lean mass") + ylab("Density") +
  xlim(0.0025,0.014) + ylim(0,480) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta2mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta2

# age (Beta3)
PosteriorBeta3.I01 <- as.data.frame(I01$marginals.fixed$`age`)
PriorBeta3.I01     <- data.frame(x = PosteriorBeta3.I01[,"x"],
                                 y = dnorm(PosteriorBeta3.I01[,"x"],-0.01,0.05))

Beta3mean.I01 <- I01Betas["age", "mean"]
Beta3lo.I01   <- I01Betas["age", "0.025quant"]
Beta3up.I01   <- I01Betas["age", "0.975quant"]

Ibeta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.I01, xmax = Beta3up.I01,
           ymin = 0, ymax = 55, fill = "gray88") +
  geom_line(data = PosteriorBeta3.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for age") + ylab("Density") +
  xlim(-0.05,0.05) + ylim(0,58) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta3

# fSexFemale (beta4)
PosteriorBeta4.I01 <- as.data.frame(I01$marginals.fixed$`fSexFemale`)
PriorBeta4.I01     <- data.frame(x = PosteriorBeta4.I01[,"x"],
                                 y = dnorm(PosteriorBeta4.I01[,"x"],-0.1,0.5))
Beta4mean.I01 <- I01Betas["fSexFemale", "mean"]
Beta4lo.I01   <- I01Betas["fSexFemale", "0.025quant"]
Beta4up.I01   <- I01Betas["fSexFemale", "0.975quant"]

Ibeta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.I01, xmax = Beta4up.I01,
           ymin = 0, ymax = 21, fill = "gray88") +
  geom_line(data = PosteriorBeta4.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Sex (female)") + ylab("Density") +
  xlim(-0.18,0.05) + ylim(0,21) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta4mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta4

# Combine plots (Fig 8.6)
ggarrange(Ibeta1, Ibeta2, Ibeta3, Ibeta4,
                        labels = c("A", "B", "C", "D"),
                        ncol = 2, nrow = 2)


#=======================================
# Model hyperparameter

# The model contains a hyperparameter for the precision parameter

# Obtain posterior distribution of precision with
I01hyp <- I01$summary.hyper[,c("mean", "sd", "mode", "0.025quant", "0.975quant")] 
#round(I01hyp, 1)

#                         mean  sd   mode  0.025quant 0.975quant
# Precision for Gamma obs 358.1 84.6  338        212      542.2

# Plot posterior distribution of precision
PosteriorHyp.I01 <- as.data.frame(I01$marginals.hyperpar$
                                    `Precision parameter for the Gamma observations`)

PriorHyp.I01 <- data.frame(x = PosteriorHyp.I01[,"x"], 
                           y = dgamma(PosteriorHyp.I01[,"x"],5.4,0.001))

Hypmean.I01 <- I01hyp["Precision parameter for the Gamma observations", "mode"]
Hyplo.I01   <- I01hyp["Precision parameter for the Gamma observations", "0.025quant"]
Hypup.I01   <- I01hyp["Precision parameter for the Gamma observations", "0.975quant"]

prec.Ihyp <- ggplot() +
  annotate("rect", xmin = Hyplo.I01, xmax = Hypup.I01,
           ymin = 0, ymax = 0.005, fill = "gray88") +
  geom_line(data = PosteriorHyp.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorHyp.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  ylab("Density") + xlab("Precision") +
  xlim(100,700) + ylim(0,0.005) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Hypmean.I01, linetype = "dashed") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# prec.Ihyp

# We can obtain the posterior for dispersion using `bri.hyperpar.summary`
I01var <- bri.hyperpar.summary(I01)[,c("mean", "mode", "q0.025", "q0.975")]
#round(I01var,4)

# mean   mode q0.025 q0.975 
# 0.0540 0.0521 0.0430 0.0686 

# And plot the posterior distribution
Hypvmean.I01 <- I01var["mode"]
Hypvlo.I01   <- I01var["q0.025"]
Hypvup.I01   <- I01var["q0.975"]

TauI01   <- I01$marginals.hyperpar$`Precision parameter for the Gamma observations`
dispI01 <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), TauI01))
PriorVar.I01 <- data.frame(x = dispI01[,"x"], 
                           y = dgamma(dispI01[,"x"],5.4,0.001))

dis.Ihyp <- ggplot() +
  annotate("rect", xmin = Hypvlo.I01, xmax = Hypvup.I01,
           ymin = 0, ymax = 65, fill = "gray88") +
  geom_line(data = dispI01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorVar.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  ylab("Density") + xlab("Dispersion") +
  xlim(0.035,0.08) + ylim(0,65) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Hypvmean.I01, linetype = "dashed") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# dis.Ihyp

# Combine plots 
ggarrange(prec.Ihyp, dis.Ihyp,
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1)

#=======================================

# Compare M01 and I01

# Extract DICs
InfDIC <- c(M01$dic$dic, I01$dic$dic)

# Add weighting
InfDIC.weights <- aicw(InfDIC)

# Add names
rownames(InfDIC.weights) <- c("default","informative")

# Print DICs
dprint.inf <- print (InfDIC.weights,
                     abbrev.names = FALSE)

# Order DICs by fit
round(dprint.inf[order(dprint.inf$fit),],2)

#               fit    delta    w
# informative 323.84  0.00 0.63
# default     324.87  1.03 0.37

# Informative model more probable
 
#=======================================
# Compare with frequentist GLM

Freq <- glm(dive_dur ~ lean_mass + age + fSex,
                       family = Gamma(link = log),
                       data = seal)

round(summary(Freq)$coef[,1:4],3)

#               Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)    6.160      0.070  88.307    0.000
# lean_mass      0.008      0.001   8.138    0.000
# age           -0.002      0.009  -0.223    0.825
# fSexFemale    -0.069      0.024  -2.845    0.009

# Dispersion
round(summary(Freq)$dispersion,3)
# 0.004

# Posterior mean values and 95% CI for fixed effects (default priors)
round(M01Betas, digits = 3)

#              mean    sd  0.025quant 0.975quant
# (Intercept)  6.161 0.060      6.043      6.279
# lean_mass    0.008 0.001      0.007      0.010
# age         -0.002 0.008     -0.018      0.013
# fSexFemale  -0.069 0.021     -0.111     -0.027

# Posterior mean values and 95% CI for fixed effects (informative priors)
round(I01Betas, digits = 3)

#              mean    sd   0.025quant 0.975quant
# (Intercept)  6.161 0.060      6.043      6.279
# lean_mass    0.008 0.001      0.007      0.010
# age         -0.002 0.008     -0.018      0.013
# fSexFemale  -0.069 0.021     -0.111     -0.027

#=======================================

# 7. MODEL CHECKS

#=======================================

# A. Model selection

# For Bayesian GLMs we can choose among models using the 
# Deviance Information Criterion (DIC)- equivalent to AIC
# Compare alternative formulations of the model

f01 <- dive_dur ~ lean_mass + age + fSex
f02 <- dive_dur ~ lean_mass + age
f03 <- dive_dur ~ lean_mass + fSex
f04 <- dive_dur ~ age + fSex

I01.full <- inla(f01,
            data = seal, family = "gamma",
            control.compute = list(dic = TRUE),
            control.family = list(hyper =
            list(prec = list(prior = "loggamma", 
            param = c(5.43, 31.6^(-2))))),
            control.fixed = list(mean.intercept = 5,
                                 prec.intercept = 5^(-2),
                                 mean = list(lean_mass = 0.01, 
                                             age = -0.01,
                                             fSexFemale = -0.1), 
                                 prec = list(lean_mass = 0.05^(-2), 
                                             age = 0.05^(-2),
                                             fSexFemale = 0.5^(-2))))

I01.1 <- inla(f02,
                 data = seal, family = "gamma",
                 control.compute = list(dic = TRUE),
                 control.family = list(hyper =
                 list(prec = list(prior = "loggamma", 
                 param = c(5.43, 31.6^(-2))))),
                 control.fixed = list(mean.intercept = 5,
                                      prec.intercept = 5^(-2),
                                      mean = list(lean_mass = 0.01, 
                                                  age = -0.01,
                                                  fSexFemale = -0.1), 
                                      prec = list(lean_mass = 0.05^(-2), 
                                                  age = 0.05^(-2),
                                                  fSexFemale = 0.5^(-2))))

I01.2 <- inla(f03,
                 data = seal, family = "gamma",
                 control.compute = list(dic = TRUE),
                 control.family = list(hyper =
                 list(prec = list(prior = "loggamma", 
                 param = c(5.43, 31.6^(-2))))),
                 control.fixed = list(mean.intercept = 5,
                                      prec.intercept = 5^(-2),
                                      mean = list(lean_mass = 0.01, 
                                                  age = -0.01,
                                                  fSexFemale = -0.1), 
                                      prec = list(lean_mass = 0.05^(-2), 
                                                  age = 0.05^(-2),
                                                  fSexFemale = 0.5^(-2))))

I01.3 <- inla(f04,
                 data = seal, family = "gamma",
                 control.compute = list(dic = TRUE),
                 control.family = list(hyper =
                 list(prec = list(prior = "loggamma", 
                 param = c(5.43, 31.6^(-2))))),
                 control.fixed = list(mean.intercept = 5,
                                      prec.intercept = 5^(-2),
                                      mean = list(lean_mass = 0.01, 
                                                  age = -0.01,
                                                  fSexFemale = -0.1), 
                                      prec = list(lean_mass = 0.05^(-2), 
                                                  age = 0.05^(-2),
                                                  fSexFemale = 0.5^(-2))))


# Compare models with the DIC
I01dic <- c(I01.full$dic$dic, I01.1$dic$dic, 
            I01.2$dic$dic,    I01.3$dic$dic)
DIC <- cbind(I01dic)
rownames(DIC) <- c("lean_mass + age + fSex", 
                   "lean_mass + age",
                   "lean_mass + fSex", 
                   "age + fSex")
round(DIC,1)

# lean_mass + age + fSex  323.8
# lean_mass + age         330.0
# lean_mass + fSex        322.1 <-
# age + fSex              359.3

#=======================================

# B. Posterior predictive check

I01.pred <- inla(f03,
                 data = seal, family = "gamma",
                 control.predictor = list(link = 1,
                                          compute = TRUE),
                 control.compute = list(dic = TRUE, 
                                        cpo = TRUE),
                 control.family = list(hyper =
                 list(prec = list(prior = "loggamma", 
                 param = c(5.43, 31.6^(-2))))),
                 control.fixed = list(mean.intercept = 5,
                                      prec.intercept = 5^(-2),
                                      mean = list(lean_mass = 0.01, 
                                                  age = -0.01,
                                                  fSexFemale = -0.1), 
                                      prec = list(lean_mass = 0.05^(-2), 
                                                  age = 0.05^(-2),
                                                  fSexFemale = 0.5^(-2))))


ppp <- vector(mode = "numeric", length = nrow(seal))
for(i in (1:nrow(seal))) {
  ppp[i] <- inla.pmarginal(q = seal$dive_dur[i],
                           marginal = I01.pred$marginals.fitted.values[[i]])
}
                      
# Fig. 8.7
ggplot() +
  geom_histogram(aes(ppp), binwidth = 0.1, 
                 colour = "black", fill = "gray88") +
  xlab("Posterior predictive p-values") +
  ylab("Frequency") +
  geom_vline(xintercept = 0.5, linetype = "dotted") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))


# Problem.

#=======================================

# C. Cross-Validation Model Checking

# Use CPO and PIT
sum(I01.pred$cpo$failure)
# [1] 0
# 0 indicates CPO is reliable

#Extract pit values
PIT <- (I01.pred$cpo$pit)

#And plot
Pit1 <- ggplot() +
  geom_histogram(aes(PIT), binwidth = 0.1, 
                 colour = "black", fill = "gray88") +
  xlab("PIT") + ylab("Frequency") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Pit1

Pit2 <- ggplot(mapping = aes(sample = I01.pred$cpo$pit)) +
  stat_qq_band(distribution = "unif", alpha = 0.5) +
  stat_qq_line(distribution = "unif", qprobs = c(0.1, 0.9)) +
  stat_qq_point(distribution = "unif", size = 2.5, alpha = 0.7) +
  xlab("Theoretical quantiles") + ylab("Sample quantiles") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Pit2

# Combine plots
# Fig. 8.8
ggarrange(Pit1, Pit2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
#PIT values are relatively uniform


#=======================================

# D. Bayesian residual analysis

# Plot residuals versus fitted values & each covariate in the model

# Obtain fitted values
# (Ensure  compute = TRUE in `control.predictor`)

# Plot residuals versus fitted values & each covariate in the model

# Obtain fitted values
# (Ensure  compute = TRUE in `control.predictor`)

Fit <- I01.pred$summary.fitted.values[, "mean"]

# Calculate residuals
Res <- seal$dive_dur - Fit
ResPlot <- cbind.data.frame(Fit,Res,seal$lean_mass,seal$fSex)

# Plot residuals against fitted
Res1 <- ggplot(ResPlot, aes(x=Fit, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("Bayesian residuals") + xlab("Fitted values") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# And plot residuals against variables in the model
Res2 <- ggplot(ResPlot, aes(x=seal$lean_mass, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Lean mass (kg)") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

Res3 <- ggplot(ResPlot, aes(x=seal$fSex, y=Res)) + 
  geom_boxplot(fill='gray88', color="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Sex") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
# Fig. 8.9
ggarrange(Res1, Res2, Res3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)


# Residuals analysis indicates no violation of model assumptions


# =======================================

# F. Prior sensitivity analysis

# A prior sensitivity analysis is advised but is not presented here.
# See Chapter 5 for details.

#=======================================

# 8. INTERPRET AND PRESENT MODEL OUTPUT

#=======================================

# The final model is: 
Final <- inla(f03,
              data = seal, family = "gamma",
              control.predictor = list(link = 1,
                                       compute = TRUE),
              control.compute = list(dic = TRUE, 
                                     cpo = TRUE),
              control.family = list(hyper =
              list(prec = list(prior = "loggamma", 
              param = c(5.43, 31.6^(-2))))),
              control.fixed = list(mean.intercept = 5,
                                   prec.intercept = 5^(-2),
                                   mean = list(lean_mass = 0.01, 
                                               age = -0.01,
                                               fSexFemale = -0.1), 
                                   prec = list(lean_mass = 0.05^(-2), 
                                               age = 0.05^(-2),
                                               fSexFemale = 0.5^(-2))))
                       
# Posterior mean values and 95% CI for fixed effects
BetasFinal <- Final$summary.fixed[,c("mean", "sd", 
                                     "0.025quant", 
                                     "0.975quant")] 
round(BetasFinal, digits = 3)

#              mean  sd   0.025quant 0.975quant
# (Intercept)  6.158 0.058      6.043      6.273
# lean_mass    0.008 0.001      0.007      0.009
# fSexFemale  -0.071 0.020     -0.111     -0.031


# Model interpretation

# Positive effect of lean mass on dive duration
# An effect of sex on dive duration (lower in females)

#=======================================

# 9. VISUALISE THE RESULTS

#=======================================

# Start by defining a data frame that contains SL and fSupp using 'ddply'

# MyData <- ddply(seal, 
#                 .(fSex), summarize,
#                    lean_mass = seq(
#                         from = min(lean_mass), 
#                           to = max(lean_mass), 
#                       length = 50))

MyData <- expand.grid(
  fSex = levels(seal$fSex),
  lean_mass = seq(
    from = min(seal$lean_mass), 
    to = max(seal$lean_mass), 
    length = 50))

# Make a design matrix
Xmat <- model.matrix(~ lean_mass + fSex, data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)

# Re-run the model in R-INLA using the combined data set, ensuring
# that `compute = TRUE` is selected in the `control.predictor` argument

# The final model is: 
Final.Pred <- inla(f03,
              data = seal, family = "gamma",
              lincomb = lcb,
              control.inla = list(lincomb.derived.only = TRUE),
              control.predictor = list(link = 1,
                                       compute = TRUE),
              control.family = list(hyper =
              list(prec = list(prior = "loggamma", 
              param = c(5.43, 31.6^(-2))))),
              control.fixed = list(mean.intercept = 5,
                                   prec.intercept = 5^(-2),
                                   mean = list(lean_mass = 0.01, 
                                               age = -0.01,
                                               fSexFemale = -0.1), 
                                   prec = list(lean_mass = 0.05^(-2), 
                                               age = 0.05^(-2),
                                               fSexFemale = 0.5^(-2))))

# Get the marginal distributions:
Pred.marg <- Final.Pred$marginals.lincomb.derived

# The output above is on the log-scale. We need to convert the 
# marginal distributions of x into the distribution of exp(x). 
# Starting point is the $marginals.lincomb.derived, which 
# contains the marginal posterior distribution for 
# each predicted value, and we have nrow(MyData) of them.
# The process is for 99% identical to that in the Poisson GLMM.

# This is a function that converts x into exp(x)
MyLog <- function(x) {exp(x)}

# Get mu, selo and seup
MyData$mu <- unlist(lapply(Pred.marg,function(x)inla.emarginal(MyLog,x)))

MyData$selo <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.025), 
                      inla.tmarginal(MyLog, x))))

MyData$seup <- unlist(lapply(Pred.marg,function(x)inla.qmarginal(c(0.975), 
                             inla.tmarginal(MyLog, x))))

# Define labels
label_sex <- c("Female" = "Female seals", 
               "Male"   = "Male seals")

# Plot
ggplot() + 
  geom_jitter(data = seal, aes(y = dive_dur, x = lean_mass),
              shape = 19, size = 2.5, height = 0.1, 
              width = 0.1, alpha = 0.7) +
  xlab("Seal lean mass (kg)") + 
  ylab("Posterior mean dive duration (sec.)") +
  # ylim(680,1300) + 
  theme(text = element_text(size = 12)) + 
  theme(panel.background = element_blank()) + 
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) + 
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  geom_line(data = MyData, aes(x = lean_mass, y = mu), size = 1) +
  geom_ribbon(data = MyData, aes(x = lean_mass, 
                                 ymax = seup, ymin = selo), alpha = 0.5) +
  facet_grid(. ~ fSex, scales = "fixed", space = "fixed",
             labeller = labeller (fSex = label_sex)) +
  theme(strip.text = element_text(size = 12, face="italic"))

#===================END===================#
