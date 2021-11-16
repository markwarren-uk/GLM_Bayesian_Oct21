
#=======================================

# R code for Chapter 4 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Load packages
library(lattice)  
library(ggplot2)
library(kableExtra)
library(GGally)
library(tidyverse)
library(mgcv)
library(lme4)
library(car)
library(devtools)
library(plyr)
library(ggpubr)
library(qqplotr)
library(gridExtra) 
library(rlang)

# Install the latest stable version of INLA:
#install.packages("INLA", repos=c(getOption("repos"), 
#                INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

# Also install brinla
#install_github("julianfaraway/brinla")
library(brinla)

#=======================================

# Import the data

# To import the data, change the working directory to the one
# in which the .csv file is saved. Here we use readr in Tidyverse.
# It is faster for bigger datasets and doesn't require a header

bitt <- read_csv("bitterling.csv") 

#Tibble output
# ── Column specification ──
# cols(
#   male = col_double(),
#   sl = col_double(),
#   supp_feed = col_double(),
#   resp_dist = col_double()
# )

str(bitt)

# We have 48 observations of 4 variables.
# Each row is an individual male bitterling
# that was observed and measured.

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

# The aim of this study is to understand whether the response distance of
# territorial male European bitterling:
# 1. Is positively related to body size
# 2. Increases with supplemental feeding
# 3. Whether the association with body sizes changes with supplemental feeding
# (i.e. there is a body size x supplemental feeding interaction)

#=======================================

# 2. PERFORM DATA EXPLORATION

#=======================================

# MISSING VALUES?
colSums(is.na(bitt))

# male  sl   supp_feed resp_dist 
# 0     0         0         0 
#No missing data


# OUTLIERS

# A multi-panel Cleveland dotplot to examine continuous variables
# Fig. 4.1

bitt <- bitt %>%
  mutate(order = seq(1:nrow(bitt)))

# Define preferred figure format
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

#Function for dotplots
multi_dotplot <- function(filename, Xvar, Yvar){
  filename %>% 
    ggplot(aes(x = {{Xvar}})) +
    geom_point(aes(y = {{Yvar}})) +
    theme_bw() +
    My_theme +
    coord_flip() +
    labs(x = "Order of Data")
}

#Choose continuous variables to plot
p1 <- multi_dotplot(bitt, male, sl) 
p2 <- multi_dotplot(bitt, male, resp_dist)

#Plot as a grid
grid.arrange(p1, p2, nrow = 1)

# No obvious outliers in the continuous variables

#NORMALITY AND HOMOGENEITY OF DEPENDENT VARIABLE

# Frequency polygon plot
# Fig 4.2
bitt %>% 
  ggplot(aes(resp_dist)) +
  geom_freqpoly( bins = 6) +
  labs(x = "Response distance (cm)", y = "Frequency") +
  My_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))

# A broadly normal distribution


# BALANCE OF CATEGORICAL VARIABLES

# The supplementary feeding covariate (supp_feed) is categorical
# but currently coded numerically. Designate as a factor:

bitt$fSupp <- factor(bitt$supp_feed)

#And examine balance
table(bitt$fSupp)

# No Yes 
# 25  23 
# A little unbalanced but acceptable


# COLLINEARITY

# Obtain summary using the ggpairs command from the GGally library
# Fig. 4.3

bitt %>% 
    ggpairs(columns = c("sl", "fSupp"),
            aes(colour=fSupp, alpha = 0.8),
            lower = list(continuous = "smooth_loess", combo = wrap("facethist", binwidth = 5))) +
    My_theme

#Calculate Variance Inflation Factor (VIF)
round(vif(lm(resp_dist ~ sl + fSupp, data = bitt)),2)
# No collinearity

# ZEROS IN THE RESPONSE VARIABLE

sum(bitt$resp_dist == 0)
#no zeros

# RELATIONSHIPS
# Fig. 4.4

# Define labels
label_supp <- c("0" = "No supplement", 
                "1" = "Food supplement")

bitt %>% 
  ggplot(aes(x = sl, y = resp_dist, size = 1)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE, size = 1.2, colour = "black") +
  facet_grid(.~ fSupp, 
             scales = "fixed", space = "fixed", 
             labeller=labeller (fSupp = label_supp)) +
  xlab("Male standard length (mm)") + 
  ylab("Male response distance (cm)") +
  theme(panel.background = element_blank()) +
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") +
  theme(text = element_text(size=14)) +
  theme(strip.text = element_text(size = 13, face="italic")) +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))

# 1. A positive effect of male size on response distance
# 2. A positive effect of feeding on response distance
# 3. Possible interaction of male size and feeding on response distance?

#=======================================

# 3. SELECT A STATISTICAL MODEL

#=======================================

# Response distance is a continuous variable and a  
# Gaussian distribution is an appropriate starting point

# The model will include main terms and a 2-way interaction

#=======================================

# 4. SPECIFY PRIORS

#=======================================

# Fit one model with default priors and second with informative
# priors based on pilot observations conducted 2 days before starting 
# data collection proper

# Import pilot data with readr:

pilot <- read_tsv("pilot.txt")

str(pilot)
# The data comprise 8 rows of 4 variables

#Obtain description of the data using (frequentist) lm
p1 <- lm(distance ~ length + supplement, data = pilot)

broom::tidy(p1)

# term          estimate std.error statistic p.value
# <chr>            <dbl>     <dbl>     <dbl>   <dbl>
# 1 (Intercept)   19.9     36.1       0.552  0.605 
# 2 length         1.27     0.684     1.86   0.123 
# 3 supplement1   34.0     12.2       2.78   0.0390

# intercept ~ 20  (sd = 40)
# sl        ~ 1.5 (sd = 1)
# supp1     ~ 35  (sd = 15)

# Pilot observations show:
# 1. Minimum response distance was about 20 cm (sd = 40) 
#    (B1~N(20,1600)) tau = 0.000625
# 2. A 1 mm increase in male length increased response
#    distance by about 1.3 cm (sd = 0.7)
#    (B2~N(1.3,0.49)) tau = 2.04082 
# 3. Supplementary feeding increased response distance 
#    by about 35 cm (sd = 15)
#    (B3~N(35,225)) tau = 0.004444
# 4 no clear prediction for an interaction - use default
#    (B4~N(0,1000) tau = 0.001

# ======================================

# 5. FIT MODELS

#=======================================

# Model M0 with default priors
M0 <- inla(resp_dist ~ sl * fSupp,
                       data = bitt)

inla.priors.used(M0)

# Model M1 with informative priors based on pilot study
M1 <- inla(resp_dist ~ sl * fSupp,  data = bitt,
             control.family = list(hyper =
                               list(prec = list(prior = "gaussian", 
                                                param = c(0, 1)))),
              control.fixed = list(mean.intercept = 20,
                                   prec.intercept = 40^(-2),
                                   mean = list(sl = 1.3, 
                                           fSupp1 = 35,
                                          default = 0), 
                                   prec = list(sl = 0.7^(-2), 
                                           fSupp1 = 15^(-2),
                                          default = 31.62^(-2))))

inla.priors.used(M1)

#=======================================

# 6. OBTAIN THE POSTERIOR DISTRIBUTION

#=======================================

# A. Default model

#Model output can be obtained with
summary(M0)

#But this gives a lot of information...


# Posterior mean values and 95% CI for fixed effects (default)
M0Betas <- M0$summary.fixed[,c("mean", "sd", 
                               "0.025quant", 
                               "0.975quant")] 
round(M0Betas, digits = 2)

#              mean    sd 0.025quant 0.975quant
# (Intercept) 59.07 16.99      25.93      92.81
# sl           1.77  0.32       1.13       2.40
# fSupp1      32.37 21.48     -10.33      74.12
# sl:fSupp1   -0.07  0.42      -0.88       0.75

# The first column shows the posterior mean, the second the
# standard deviation and the third and fourth columns show the 
# 2.5% and 97.5% quantiles of the posterior distribution.

# These two latter values are the 95% credible intervals; there is 
# a 95% probability that the regression parameters fall within these 
# intervals. 

# If 0 lies outside this interval we conclude that the regression parameter 
# is  statistically 'important’. This is the Bayesian equivalent of 
# ‘significant’ in a frequentist model.
 
# Instead of summarising the posterior distribution of the fixed effects 
# with a posterior mean and a 95% credible interval, we can plot
# the posterior distributions - available in the object `M0$marginals.fixed`. 

# Plot posterior distributions for fixed parameters 
# Fig. 4.5

# Model intercept (Beta1)
PosteriorBeta1.M0 <- data.frame(M0$marginals.fixed$`(Intercept)`)
PriorBeta1.M0     <- data.frame(x = PosteriorBeta1.M0[,"x"], 
                          y = dnorm(PosteriorBeta1.M0[,"x"],0,0))
Beta1mean.M0 <- M0Betas["(Intercept)", "mean"]
Beta1lo.M0   <- M0Betas["(Intercept)", "0.025quant"]
Beta1up.M0   <- M0Betas["(Intercept)", "0.975quant"]

#Create plot object
dbeta1 <- PosteriorBeta1.M0 %>% 
  ggplot(aes(y = y, x = x)) + 
  annotate("rect", xmin = Beta1lo.M0, xmax = Beta1up.M0,
           ymin = 0, ymax = 0.027, fill = "gray88") +
  geom_line(lwd = 1.2) +
  geom_line(data = PriorBeta1.M0,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") +
  ylab("Density") +
  xlim(-30,140) + 
  ylim(0,0.027) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.M0, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Male sl (Beta2)
PosteriorBeta2.M0 <- as.data.frame(M0$marginals.fixed$`sl`)
PriorBeta2.M0 <- data.frame(x = PosteriorBeta2.M0[,"x"], 
                         y = dnorm(PosteriorBeta2.M0[,"x"],0,0))
Beta2mean.M0 <- M0Betas["sl", "mean"]
Beta2lo.M0   <- M0Betas["sl", "0.025quant"]
Beta2up.M0   <- M0Betas["sl", "0.975quant"]

dbeta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.M0, xmax = Beta2up.M0,
           ymin = 0, ymax = 1.5, fill = "gray88") +
  geom_line(data = PosteriorBeta2.M0,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.M0,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for standard length") +
  ylab("Density") +
  xlim(-0.5,3.5) + 
  ylim(0,1.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta2mean.M0, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Supplementary feeding (Beta3)
PosteriorBeta3.M0 <- as.data.frame(M0$marginals.fixed$`fSupp`)
PriorBeta3.M0     <- data.frame(x = PosteriorBeta3.M0[,"x"],
                          y = dnorm(PosteriorBeta3.M0[,"x"],0,0))
inla.priors.used(M0)
Beta3mean.M0 <- M0Betas["fSupp", "mean"]
Beta3lo.M0   <- M0Betas["fSupp", "0.025quant"]
Beta3up.M0   <- M0Betas["fSupp", "0.975quant"]

dbeta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.M0, xmax = Beta3up.M0,
           ymin = 0, ymax = 0.022, fill = "gray88") +
  geom_line(data = PosteriorBeta3.M0,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.M0,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for suppl. feeding") +
  ylab("Density") +
  xlim(-50,120) + 
  ylim(0,0.022) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.M0, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# 2-way interaction - `sl:fSupp1`
PosteriorBeta4.M0 <- as.data.frame(M0$marginals.fixed$`sl:fSupp1`)
PriorBeta4.M0     <- data.frame(x = PosteriorBeta4.M0[,"x"],
                                y = dnorm(PosteriorBeta4.M0[,"x"],0,0))
Beta4mean.M0 <- M0Betas["sl:fSupp1", "mean"]
Beta4lo.M0   <- M0Betas["sl:fSupp1", "0.025quant"]
Beta4up.M0   <- M0Betas["sl:fSupp1", "0.975quant"]

dbeta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.M0, xmax = Beta4up.M0,
           ymin = 0, ymax = 1.18, fill = "gray88") +
  geom_line(data = PosteriorBeta4.M0,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.M0,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Interaction") +
  ylab("Density") +
  xlim(-2.1,2.1) + ylim(0,1.18) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta4mean.M0, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
ggarrange(dbeta1, dbeta2, dbeta3, dbeta4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)
# These density functions look normally distributed because INLA assumes
# that the posterior distributions of the betas are Gaussian.

#=======================================
# Model hyperparameter

# The model contains a parameter sigma that is used for the variance (sigma^2) 
# of the normal distribution for male response distance  The sigma is called a 
# 'hyperparameter'. Whereas the `lm` function gives the estimated value of 
# sigma, getting the posterior mean value of sigma from R-INLA is slightly 
# complicated. R-INLA uses  precision (tau), which is defined as:
# tau = sigma^(-2). 

# Obtain posterior distribution of precision (tau)
M0hyp <- M0$summary.hyper[,c("mean","mode","0.025quant","0.975quant")] 
round(M0hyp, 4)

#                                         mean  mode   0.025quant 0.975quant
# Precision for the Gaussian observations 0.004 0.0038 0.0025     0.0057

# Plot posterior distribution of precision (tau)
# Fig. 4.6

PosteriorHyp.M0 <- as.data.frame(M0$marginals.hyperpar$
                   `Precision for the Gaussian observations`)
PriorHyp.M0 <- data.frame(x = PosteriorHyp.M0[,"x"], 
                   y = dgamma(PosteriorHyp.M0[,"x"],1,2^5, log = TRUE))

Hypmean.M0 <- M0hyp["Precision for the Gaussian observations", "mode"]
Hyplo.M0   <- M0hyp["Precision for the Gaussian observations", "0.025quant"]
Hypup.M0   <- M0hyp["Precision for the Gaussian observations", "0.975quant"]

prec.M0 <- ggplot() +
  annotate("rect", xmin = Hyplo.M0, xmax = Hypup.M0,
           ymin = 0, ymax = 550, fill = "gray88") +
  geom_line(data = PosteriorHyp.M0,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorHyp.M0,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  ylab("Density") +
  xlab(expression(paste("Tau (", tau ,")"))) +
  xlim(0,0.009) + ylim(0,550) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Hypmean.M0, linetype = "dashed") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
prec.M0

# We can obtain sigma using `bri.hyperpar.summary`
round(bri.hyperpar.summary(M0), 3)

#                                  mean   sd    q0.025 q0.5   q0.975 mode
# SD for the Gaussian observations 16.168 1.712 13.236 16.008 19.966 15.687

# And plot the posterior distribution of sigma (given
# the skewed distribution, the mode is more appropriate)
# Fig. 4.7

M0var <- bri.hyperpar.summary(M0)[,c("mode", "q0.025", "q0.975")]

Hypvmean.M0 <- M0var["mode"]
Hypvlo.M0   <- M0var["q0.025"]
Hypvup.M0   <- M0var["q0.975"]

TauM0   <- M0$marginals.hyperpar$`Precision for the Gaussian observations`
SigmaM0 <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), TauM0))
PriorVar.M0 <- data.frame(x = SigmaM0[,"x"], 
                          y = dgamma(SigmaM0[,"x"],1,2^(-5)))

sigma.M0 <- ggplot()  +
  annotate("rect", xmin = Hypvlo.M0, xmax = Hypvup.M0,
           ymin = 0, ymax = 0.27, fill = "gray88") +
  geom_line(data = SigmaM0,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorVar.M0,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  ylab("") +
  xlab(expression(paste("SD (", sigma ,")"))) +
  xlim(10,25) + ylim(0,0.27) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Hypvmean.M0, linetype = "dashed") +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
sigma.M0

#=======================================
# B. Informative model

# Posterior mean values and 95% CI for fixed effects (informative)
M1Betas <- M1$summary.fixed[,c("mean", "sd", 
                               "0.025quant", 
                               "0.975quant")] 
round(M1Betas, digits = 2)

#              mean    sd 0.025quant 0.975quant
# (Intercept) 55.08 12.53      30.48      79.75
# sl           1.84  0.24       1.38       2.31
# fSupp1      41.15 12.81      15.95      66.22
# sl:fSupp1   -0.24  0.26      -0.74       0.27

# The first column shows the posterior mean, the second the
# standard deviation and the third and fourth columns show the 
# 2.5% and 97.5% quantiles of the posterior distribution.

# Plot posterior distributions for fixed parameters 
# Fig. 4.8

# Model intercept (Beta1)
PosteriorBeta1.M1 <- as.data.frame(M1$marginals.fixed$`(Intercept)`)
PriorBeta1.M1     <- data.frame(x = PosteriorBeta1.M1[,"x"], 
                          y = dnorm(PosteriorBeta1.M1[,"x"],20,40))
Beta1mean.M1 <- M1Betas["(Intercept)", "mean"]
Beta1lo.M1   <- M1Betas["(Intercept)", "0.025quant"]
Beta1up.M1   <- M1Betas["(Intercept)", "0.975quant"]

ibeta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.M1, xmax = Beta1up.M1,
           ymin = 0, ymax = 0.035, fill = "gray88")+
  geom_line(data = PosteriorBeta1.M1,
            aes(y = y, x = x), lwd = 1.2)+
  geom_line(data = PriorBeta1.M1,
            aes(y = y, x = x), color = "gray55", lwd = 1.2)+
  xlab("Intercept")+
  ylab("Density")+
  xlim(-30,140) + ylim(0,0.035)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = Beta1mean.M1, linetype = "dashed")+
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1))+
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Male sl (Beta2)
PosteriorBeta2.M1 <- as.data.frame(M1$marginals.fixed$`sl`)
PriorBeta2.M1     <- data.frame(x = PosteriorBeta2.M1[,"x"],
                          y = dnorm(PosteriorBeta2.M1[,"x"],1.3,0.7))
Beta2mean.M1 <- M1Betas["sl", "mean"]
Beta2lo.M1   <- M1Betas["sl", "0.025quant"]
Beta2up.M1   <- M1Betas["sl", "0.975quant"]

ibeta2 <- ggplot()+
  annotate("rect", xmin = Beta2lo.M1, xmax = Beta2up.M1,
           ymin = 0, ymax = 1.8, fill = "gray88")+
  geom_line(data = PosteriorBeta2.M1,
            aes(y = y, x = x), lwd = 1.2)+
  geom_line(data = PriorBeta2.M1,
            aes(y = y, x = x), color = "gray55", lwd = 1.2)+
  xlab("Slope for Standard Length")+
  ylab("Density")+
  xlim(-0.5,3.5) + ylim(0,1.8)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = Beta2mean.M1, linetype = "dashed")+
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1))+
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))


# Supplementary feeding (Beta3)
PosteriorBeta3.M1 <- as.data.frame(M1$marginals.fixed$`fSupp`)
PriorBeta3.M1     <- data.frame(x = PosteriorBeta3.M1[,"x"],
                                y = dnorm(PosteriorBeta3.M1[,"x"],35,15))

Beta3mean.M1 <- M1Betas["fSupp", "mean"]
Beta3lo.M1   <- M1Betas["fSupp", "0.025quant"]
Beta3up.M1   <- M1Betas["fSupp", "0.975quant"]

ibeta3 <- ggplot()+
  annotate("rect", xmin = Beta3lo.M1, xmax = Beta3up.M1,
           ymin = 0, ymax = 0.035, fill = "gray88")+
  geom_line(data = PosteriorBeta3.M1,
            aes(y = y, x = x), lwd = 1.2)+
  geom_line(data = PriorBeta3.M1,
            aes(y = y, x = x), color = "gray55", lwd = 1.2)+
  xlab("Slope for suppl. feeding")+
  ylab("Density")+
  xlim(-50,120) + ylim(0,0.035)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = Beta3mean.M1, linetype = "dashed")+
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1))+
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# 2-way interaction - `sl:fSupp1`
PosteriorBeta4.M1 <- as.data.frame(M1$marginals.fixed$`sl:fSupp1`)
PriorBeta4.M1     <- data.frame(x = PosteriorBeta4.M1[,"x"],
                                y = dnorm(PosteriorBeta4.M1[,"x"],0,31.62))
Beta4mean.M1 <- M1Betas["sl:fSupp1", "mean"]
Beta4lo.M1   <- M1Betas["sl:fSupp1", "0.025quant"]
Beta4up.M1   <- M1Betas["sl:fSupp1", "0.975quant"]

ibeta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.M1, xmax = Beta4up.M1,
           ymin = 0, ymax = 1.8, fill = "gray88")+
  geom_line(data = PosteriorBeta4.M1,
            aes(y = y, x = x), lwd = 1.2)+
  geom_line(data = PriorBeta4.M1,
            aes(y = y, x = x), color = "gray55", lwd = 1.2)+
  xlab("Interaction")+
  ylab("Density")+
  xlim(-2.1,2.1) + ylim(0,1.8)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = Beta4mean.M1, linetype = "dashed")+
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1))+
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
ggarrange(ibeta1, ibeta2, ibeta3, ibeta4,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)

# Model hyperparameter

# Obtain posterior distribution of precision (tau)
M1hyp <- M1$summary.hyper[,c("mean","mode","0.025quant","0.975quant")] 
round(M1hyp, digits = 4)

#                                         mean   mode    0.025quant 0.975quant
# Precision for the Gaussian observations 0.0048 0.0046  0.0032     0.0067

# Plot posterior distribution of precision (tau)
# Fig. 4.9

PosteriorHyp.M1 <- as.data.frame(M1$marginals.hyperpar$
                   `Precision for the Gaussian observations`)
PriorHyp.M1 <- data.frame(x = PosteriorHyp.M1[,"x"], 
                        y = dnorm(PosteriorHyp.M1[,"x"],0,1))
Hypmean.M1 <- M1hyp["Precision for the Gaussian observations", "mode"]
Hyplo.M1   <- M1hyp["Precision for the Gaussian observations", "0.025quant"]
Hypup.M1   <- M1hyp["Precision for the Gaussian observations", "0.975quant"]

prec.M1 <- ggplot() +
  annotate("rect", xmin = Hyplo.M1, xmax = Hypup.M1,
           ymin = 0, ymax = 510, fill = "gray88")+
  geom_line(data = PosteriorHyp.M1,
            aes(y = y, x = x), lwd = 1.2)+
  geom_line(data = PriorHyp.M1,
            aes(y = y, x = x), color = "gray55", lwd = 1.2)+
  ylab("Density")+
  xlab(expression(paste("Tau (", tau ,")")))+
  xlim(0,0.009) + ylim(0,510)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = Hypmean.M1, linetype = "dashed")+
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1))+
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
prec.M1

# We can obtain sigma using `bri.hyperpar.summary`
round(bri.hyperpar.summary(M1), 3)

# And plot the posterior distribution of sigma
# Fig. 4.10

M1var <- bri.hyperpar.summary(M1)[,c("mean","mode", "q0.025", "q0.975")]
Hypvmean.M1 <- M1var["mode"]
Hypvlo.M1   <- M1var["q0.025"]
Hypvup.M1   <- M1var["q0.975"]

TauM1 <- M1$marginals.hyperpar$`Precision for the Gaussian observations`
SigmaM1 <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), TauM1))
PriorVar.M1 <- data.frame(x = SigmaM1[,"x"], 
                          y = dnorm(SigmaM1[,"x"],0,1))

sigma.M1 <- ggplot() +
  annotate("rect", xmin = Hypvlo.M1, xmax = Hypvup.M1,
           ymin = 0, ymax = 0.33, fill = "gray88")+
  geom_line(data = SigmaM1,
            aes(y = y, x = x), lwd = 1.2)+
  geom_line(data = PriorVar.M1,
            aes(y = y, x = x), color = "gray55", lwd = 1.2)+
  ylab("")+
  xlab(expression(paste("SD (", sigma ,")")))+
  xlim(10,20) + ylim(0,0.33)+
  geom_vline(xintercept = 0, linetype = "dotted")+
  geom_vline(xintercept = Hypvmean.M1, linetype = "dashed")+
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank())+
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1))+
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
sigma.M1

#=======================================
# Compare with frequentist GLM

Freq <- lm(resp_dist ~ sl * fSupp, 
           data = bitt)

broom::tidy(Freq)%>%
  mutate_if(is.numeric, round, 4)

# A tibble: 4 x 5
# term         estimate  std.error statistic  p.value
# <chr>          <dbl>     <dbl>     <dbl>   <dbl>
# 1 (Intercept)   47.7      18.8        2.54  0.0148
# 2 sl             1.98      0.353      5.61  0     
# 3 fSupp1        59.1      28.9        2.05  0.0467
# 4 sl:fSupp1     -0.580     0.554     -1.05  0.301 

# Posterior mean values and 95% CI for fixed effects (default priors)
round(M0Betas, digits = 2)

#              mean    sd  0.025quant 0.975quant
# (Intercept) 59.07 16.99      25.93      92.81
# sl           1.77  0.32       1.13       2.40
# fSupp1      32.37 21.48     -10.33      74.12
# sl:fSupp1   -0.07  0.42      -0.88       0.75

# Posterior mean values and 95% CI for fixed effects (informative priors)
round(M1Betas, digits = 2)

#              mean    sd   0.025quant 0.975quant
# (Intercept) 55.29 12.47      30.79      79.84
# sl           1.84  0.24       1.37       2.30
# fSupp1      40.49 12.44      16.01      64.85
# sl:fSupp1   -0.22  0.25      -0.71       0.27


#Also obtain estimates of sigma

# For frequentist model
round(summary(Freq)$sigma,1)
# 16.2

# For Bayesian model with default priors
round(bri.hyperpar.summary(M0)[,c("mean")],1)
# 16.2

# For Bayesian model with informative priors
round(bri.hyperpar.summary(M1)[,c("mean")],1)
# 14.7

#=======================================

# 7. MODEL CHECKS

#=======================================

# A. Model selection

# For Bayesian GLMs we can choose among models using the 
# Deviance Information Criterion (DIC)  - equivalent to AIC

# To use DIC we must re-run the model and specify its calculation using 
# 'control.compute'

M0 <- inla(resp_dist ~ sl * fSupp,
                       control.compute = list(dic = TRUE),
                       data = bitt)

# Model selection
# Full model with default priors
M0.full <- inla(resp_dist ~ sl * fSupp,
                            control.compute = list(dic = TRUE),
                            data = bitt)
# Model with default priors with interaction dropped
M0.1 <- inla(resp_dist ~ sl + fSupp,
                         control.compute = list(dic = TRUE),
                         data = bitt)
# Model with default priors with supplementary feeding dropped
M0.2 <- inla(resp_dist ~ sl,
                         control.compute = list(dic = TRUE),
                         data = bitt)
# Model with default priors with standard length dropped
M0.3 <- inla(resp_dist ~ fSupp,
                         control.compute = list(dic = TRUE),
                         data = bitt)

# Compare models with DIC
DIC <- cbind(c(M0.full$dic$dic, M0.1$dic$dic, 
               M0.2$dic$dic,    M0.3$dic$dic))

#DIC <- cbind(M0dic)
rownames(DIC) <- c("full","no inter","no suppl","no sl")
colnames(DIC) <- c("DIC")
round(DIC,1)

# full     409.6
# no inter 408.9 <- lowest DIC
# no suppl 436.7
# no sl    438.0

# Now with informative priors
M1.full <- inla(resp_dist ~ sl * fSupp,  data = bitt,
           control.compute = list(dic = TRUE),
           control.family = list(hyper =
                                   list(prec = list(prior = "gaussian", 
                                                    param = c(0, 1)))),
           control.fixed = list(mean.intercept = 20,
                                prec.intercept = 40^(-2),
                                mean = list(sl = 1.3, 
                                            fSupp1 = 35,
                                            default = 0), 
                                prec = list(sl = 0.7^(-2), 
                                            fSupp1 = 15^(-2),
                                            default = 1000)))

M1.1 <- inla(resp_dist ~ sl + fSupp,  data = bitt,
                control.compute = list(dic = TRUE),
                control.family = list(hyper =
                                        list(prec = list(prior = "gaussian", 
                                                         param = c(0, 1)))),
                control.fixed = list(mean.intercept = 20,
                                     prec.intercept = 40^(-2),
                                     mean = list(sl = 1.3, 
                                                 fSupp1 = 35,
                                                 default = 0), 
                                     prec = list(sl = 0.7^(-2), 
                                                 fSupp1 = 15^(-2),
                                                 default = 1000)))

M1.2 <- inla(resp_dist ~ sl,  data = bitt,
                control.compute = list(dic = TRUE),
                control.family = list(hyper =
                                        list(prec = list(prior = "gaussian", 
                                                         param = c(0, 1)))),
                control.fixed = list(mean.intercept = 20,
                                     prec.intercept = 40^(-2),
                                     mean = list(sl = 1.3, 
                                                 fSupp1 = 35,
                                                 default = 0), 
                                     prec = list(sl = 0.7^(-2), 
                                                 fSupp1 = 15^(-2),
                                                 default = 1000)))

M1.3 <- inla(resp_dist ~ fSupp,  data = bitt,
                control.compute = list(dic = TRUE),
                control.family = list(hyper =
                                        list(prec = list(prior = "gaussian", 
                                                         param = c(0, 1)))),
                control.fixed = list(mean.intercept = 20,
                                     prec.intercept = 40^(-2),
                                     mean = list(sl = 1.3, 
                                                 fSupp1 = 35,
                                                 default = 0), 
                                     prec = list(sl = 0.7^(-2), 
                                                 fSupp1 = 15^(-2),
                                                 default = 1000)))

DIC1 <- cbind(c(M1.full$dic$dic, M1.1$dic$dic, 
           M1.2$dic$dic,    M1.3$dic$dic))
rownames(DIC1) <- c("full","no inter","no suppl","no sl")
round(DIC1,1)

# full     408.0
# no inter 408.6
# no suppl 436.7
# no sl    438.2

# Compare best-fitting models with non-informative and informative priors
DIC2 <- cbind(c(M0.1$dic$dic, M1.1$dic$dic))
rownames(DIC2) <- c("default priors","informative priors")
colnames(DIC2) <- "DIC"
round(DIC2,2)

# default priors     408.92
# informative priors 408.60

# These are essentially the same

#=======================================

# B. Posterior predictive check
# (Just for model with informative priors)
M1.pred <- inla(resp_dist ~ sl + fSupp,  data = bitt,
                control.predictor = list(link = 1,
                                      compute = TRUE),
                control.compute = list(dic = TRUE, 
                                       cpo = TRUE),
               control.family = list(hyper =
                                 list(prec = list(prior = "gaussian",
                                                  param = c(0,1)))),
                control.fixed = list(mean.intercept = 20,
                                     prec.intercept = 40^(-2),
                                     mean = list(sl = 1.3, 
                                             fSupp1 = 35,
                                            default = 0), 
                                     prec = list(sl = 0.7^(-2), 
                                             fSupp1 = 15^(-2),
                                            default = 1000)))

ppp <- vector(mode = "numeric", length = nrow(bitt))
for(i in (1:nrow(bitt))) {
  ppp[i] <- inla.pmarginal(q = bitt$resp_dist[i],
                    marginal = M1.pred$marginals.fitted.values[[i]])
}

# Fig. 4.11 
ggplot() +
geom_histogram(aes(ppp), binwidth = 0.1, 
                    colour = "black", fill = "gray88") +
xlab("Posterior predictive p-values") +
ylab("Frequency") +
geom_vline(xintercept = 0.5, linetype = "dotted") +
theme(text = element_text(size=15))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                         colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1)) 


# Problem. This outcome may be affected by the 
# strong effect of supplemental feeding, which is binomial.
# As a result we need additional model checks...

#=======================================

# C. Cross-Validation Model Checking

# Use CPO and PIT
sum(M1.pred$cpo$failure)
# [1] 0
# 0 indicates CPO is reliable

#Extract pit values
PIT <- (M1.pred$cpo$pit)

#And plot
# Fig. 4.12
Pit1 <- ggplot() +
geom_histogram(aes(PIT), binwidth = 0.11, 
                        colour = "black", fill = "gray88") +
xlab("PIT") +
ylab("Frequency") +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1)) 


Pit2 <- ggplot(mapping = aes(sample = M1.pred$cpo$pit)) +
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


# Combine plots
ggarrange(Pit1, Pit2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)
#PIT values are uniform

#=======================================

# D. Bayesian residual analysis

# Plot residuals versus fitted values & each covariate in the model

# Obtain fitted values
# (Ensure  compute = TRUE in `control.predictor`)

Fit <- M1.pred$summary.fitted.values[, "mean"]

# Calculate residuals
Res <- bitt$resp_dist - Fit
ResPlot <- cbind.data.frame(Fit,Res,bitt$sl,bitt$fSupp)

# Plot residuals against fitted
Fig.A <- ggplot(ResPlot, aes(x=Fit, y=Res)) + 
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
Fig.B <- ggplot(ResPlot, aes(x=bitt$sl, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Male SL (mm)") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

Fig.C <- ggplot(ResPlot, aes(x=bitt$fSupp, y=Res)) + 
  geom_boxplot(fill='gray88', color="black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Suppl. feeding") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
# Fig. 4.13
ggarrange(Fig.A, Fig.B, Fig.C,
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)

# Residuals analysis indicates no violation of model assumptions

# =======================================
  
# E. Prior sensitivity analysis

# Model with priors unchanged
M1.inform <- inla(resp_dist ~ sl + fSupp,  data = bitt,
                control.predictor = list(link = 1,
                                         compute = TRUE),
                control.compute = list(dic = TRUE, 
                                       cpo = TRUE),
                control.family = list(hyper =
                                        list(prec = list(prior = "gaussian",
                                                         param = c(0,1)))),
                control.fixed = list(mean.intercept = 20,
                                     prec.intercept = 40^(-2),
                                     mean = list(sl = 1.3, 
                                             fSupp1 = 35,
                                            default = 0), 
                                     prec = list(sl = 0.7^(-2), 
                                             fSupp1 = 15^(-2),
                                            default = 1000)))

Betas.Inform <- tibble(M1.inform$summary.fixed[,c("mean", "sd", 
                                           "0.025quant", 
                                           "0.975quant")]%>% 
  mutate(term = c("Intercept", "sl", "fSupp"),
         change_prior = 0)) 

Betas.Inform

# A tibble: 3 x 6
#   mean     sd  `0.025quant` `0.975quant` term      change_prior
#   <dbl>  <dbl>        <dbl>        <dbl> <chr>            <dbl>
# 1 58.6  12.0          34.9         82.1  Intercept            0
# 2  1.77  0.223         1.33         2.21 sl                   0
# 3 30.0   4.11         21.9         38.1  fSupp                0

#== 1. Increase priors by 20% ==
#   Intercept from 20 to 24 (var 1600 to 1920)
#   sl from 1.3 to 1.56 (var 0.49 to 0.59)
#   fSupp1 from 35 to 42 (var 225 to 270)

M1.plus20 <- inla(resp_dist ~ sl + fSupp,  data = bitt,
           control.predictor = list(compute = TRUE),
             control.compute = list(dic = TRUE, 
                                    cpo = TRUE),
            control.family = list(hyper =
                 list(prec = list(prior = "gaussian",
                                  param = c(0,1)))),
           control.fixed = list(mean.intercept = 24,
                                prec.intercept = 43.8^(-2),
                                mean = list(sl = 1.56, 
                                        fSupp1 = 42,
                                       default = 0), 
                                prec = list(sl = 0.77^(-2), 
                                        fSupp1 = 16.4^(-2),
                                       default = 1000)))

# Obtain estimates of betas
Betas.plus20 <- tibble(M1.plus20$summary.fixed[,c("mean", "sd", 
                                           "0.025quant", 
                                           "0.975quant")]%>% 
  mutate(term = c("Intercept", "sl", "fSupp"),
         change_prior = 20)) 

Betas.plus20

# A tibble: 3 x 6
#   mean     sd  `0.025quant` `0.975quant` term      change_prior
#   <dbl>  <dbl>        <dbl>        <dbl> <chr>            <dbl>
# 1 57.6  12.2          33.5         81.5  Intercept           20
# 2  1.78  0.227         1.34         2.23 sl                  20
# 3 30.3   4.15         22.2         38.5  fSupp               20

#== 2. Decrease priors by 20% ==
#   Intercept from 20 to 16 (var 1600 to 1280)
#   sl from 1.3 to 1.04 (var 0.49 to 0.39)
#   fSupp1 from 35 to 28 (var 225 to 180)

M1.minus20 <- inla(resp_dist ~ sl + fSupp,  data = bitt,
                   control.predictor = list(compute = TRUE),
                   control.compute = list(dic = TRUE, 
                                          cpo = TRUE),
                   control.family = list(hyper =
                                           list(prec = list(prior = "gaussian",
                                                            param = c(0,1)))),
                   control.fixed = list(mean.intercept = 16,
                                        prec.intercept = 35.8^(-2),
                                        mean = list(sl = 1.04, 
                                                    fSupp1 = 28,
                                                    default = 0), 
                                        prec = list(sl = 0.62^(-2), 
                                                    fSupp1 = 13.4^(-2),
                                                    default = 1000)))

# Obtain estimates of betas
Betas.minus20 <- tibble(M1.minus20$summary.fixed[,c("mean", "sd", 
                                             "0.025quant", 
                                             "0.975quant")] %>% 
  mutate(term = c("Intercept", "sl", "fSupp"),
         change_prior = -20)) 

Betas.minus20

# A tibble: 3 x 6
#   mean     sd  `0.025quant` `0.975quant` term      change_prior
#   <dbl>  <dbl>        <dbl>        <dbl> <chr>            <dbl>
# 1 60.0  11.7          37.0         83.0  Intercept          -20
# 2  1.74  0.218         1.31         2.17 sl                 -20
# 3 29.4   4.08         21.4         37.5  fSupp              -20

#=========================
# Plot posterior distributions of the alternative models

# Model intercept (Beta1)
PostBeta1.M1.inform  <- as.data.frame(M1.inform$marginals.fixed$`(Intercept)`)
PostBeta1.M1.plus20  <- as.data.frame(M1.plus20$marginals.fixed$`(Intercept)`)
PostBeta1.M1.minus20 <- as.data.frame(M1.minus20$marginals.fixed$`(Intercept)`)

beta1.sens <- ggplot() +
geom_line(data = PostBeta1.M1.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid")+
geom_line(data = PostBeta1.M1.plus20, 
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44") +
geom_line(data = PostBeta1.M1.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44") +
xlab("Intercept") +
ylab("Density") +
xlim(20,100) +
theme(text = element_text(size=13)) +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                  colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# Male sl (Beta2)
PostBeta2.M1.inform  <- as.data.frame(M1.inform$marginals.fixed$`sl`)
PostBeta2.M1.plus20  <- as.data.frame(M1.plus20$marginals.fixed$`sl`)
PostBeta2.M1.minus20 <- as.data.frame(M1.minus20$marginals.fixed$`sl`)

beta2.sens <- ggplot() +
geom_line(data = PostBeta2.M1.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid")+
geom_line(data = PostBeta2.M1.plus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44")+
geom_line(data = PostBeta2.M1.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44")+
xlab("Slope for male SL")+
ylab("Density")+
xlim(1,2.5)+
theme(text = element_text(size=13)) +
theme(panel.background = element_blank())+
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1))+
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# Supplementary feeding (Beta3)
PostBeta3.M1.inform  <- as.data.frame(M1.inform$marginals.fixed$`fSupp`)
PostBeta3.M1.plus20  <- as.data.frame(M1.plus20$marginals.fixed$`fSupp`)
PostBeta3.M1.minus20 <- as.data.frame(M1.minus20$marginals.fixed$`fSupp`)

beta3.sens <- ggplot() +
geom_line(data = PostBeta3.M1.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid")+
geom_line(data = PostBeta3.M1.plus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44")+
geom_line(data = PostBeta3.M1.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44")+
xlab("Slope for suppl. feeding")+
ylab("Density")+
xlim(15,45)+
theme(text = element_text(size=13)) +
theme(panel.background = element_blank())+
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1))+
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# And plot the posterior distributions of sigma
PostBeta3.M1.inform  <- as.data.frame(M1.inform$marginals.fixed$`fSupp`)
PostBeta3.M1.plus20  <- as.data.frame(M1.plus20$marginals.fixed$`fSupp`)
PostBeta3.M1.minus20 <- as.data.frame(M1.minus20$marginals.fixed$`fSupp`)

Tau.M1.inform  <- M1.inform$marginals.hyperpar$`Precision for the Gaussian observations`
Tau.M1.plus20  <- M1.plus20$marginals.hyperpar$`Precision for the Gaussian observations`
Tau.M1.minus20 <- M1.minus20$marginals.hyperpar$`Precision for the Gaussian observations`

Sigma.M1.inform  <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), Tau.M1.inform))
Sigma.M1.plus20  <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), Tau.M1.plus20))
Sigma.M1.minus20 <- as.data.frame(inla.tmarginal(function(x) sqrt(1/x), Tau.M1.minus20))

sigma.sens <- ggplot() +
geom_line(data = Sigma.M1.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid")+
geom_line(data = Sigma.M1.plus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44")+
geom_line(data = Sigma.M1.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44")+
ylab("Density")+
xlab(expression(paste(sigma)))+
theme(text = element_text(size=13)) +
theme(panel.background = element_blank())+
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1))+
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# Combine plots
# Fig. 4.14
ggarrange(beta1.sens, beta2.sens, 
                    beta3.sens, sigma.sens,
                    labels = c("A", "B", 
                               "C", "D"),
                    ncol = 2, nrow = 2)


# Changes to posterior distributions are negligible

#=======================================

# 8. INTERPRET AND PRESENT MODEL OUTPUT

#=======================================

# The final model is: 
Final <- inla(resp_dist ~ sl + fSupp,  data = bitt,
              control.predictor = list(link = 1,
                                    compute = TRUE),
                 control.compute = list(dic = TRUE, 
                                        cpo = TRUE),
                control.family = list(hyper =
                     list(prec = list(prior = "gaussian",
                                      param = c(0,1)))),
                  control.fixed = list(mean.intercept = 20,
                                       prec.intercept = 40^(-2),
                                       mean = list(sl = 1.3, 
                                               fSupp1 = 35,
                                              default = 0), 
                                       prec = list(sl = 0.7^(-2), 
                                               fSupp1 = 15^(-2),
                                              default = 1000)))

# Posterior mean values and 95% CI for fixed effects
BetasFinal <- Final$summary.fixed[,c("mean", "sd", 
                                     "0.025quant", 
                                     "0.975quant")] 
round(BetasFinal, digits = 2)

#              mean  sd   0.025quant 0.975quant
# (Intercept) 58.56 11.99      34.92      82.09
# sl           1.77  0.22       1.33       2.21
# fSupp1      29.95  4.11      21.88      38.07

SigmaFinal <- bri.hyperpar.summary(Final)[,c("mean", 
                                             "q0.025", 
                                             "q0.975")]

Sigma_df <- data.frame(as.list(SigmaFinal))

Sigma_df %>% 
  mutate(term = "sigma")

## Model interpretation

# Positive effect of sl on response distance
# Positive effect of supplementary feeding on response distance
# Effect of sl does not vary with supplementary feeding

#=======================================

# 9. VISUALISE THE RESULTS

#=======================================

# Start by defining a data frame that contains SL and fSupp using 'ddply'

MyData <- ddply(bitt,
                .(fSupp), summarize,
                          sl = seq(
                        from = min(sl),
                          to = max(sl),
                         length = 50))

# This creates 100 artificial covariate values. 
# There is no `predict` function in INLA, but we can obtain fitted 
# values manually with a design matrix for the values in `MyData` and 
# then multiplying this with the posterior mean values of the model. 

# We must add an extra variable for the response variable to `MyData`
# and set it to NA. We will then combine the `bitt2` and `MyData` objects, 
# and apply R-INLA to this combined data set. R-INLA will predict the 
# response variable where an NA occurs.

MyData$resp_dist <- NA
head(MyData, 6)

# Because the data frame `bitt` has various covariates, we have to be careful how we
# combine `bitt` and `MyData`.

bitt.Pred <- bitt[, colnames(MyData)]

# It is important to ensure that the order of the columns in
# `bitt.Pred` is the same as in `MyData`, because we will combine
# them with the `rbind` function. The `colnames(MyData)` construction
# ensures that this is indeed the case.

bitt.Comb <- rbind(bitt.Pred, MyData)

# We next re-run the model in R-INLA using the combined data set, ensuring
# that `compute = TRUE` is selected in the `control.predictor` argument

Final.Pred <- inla(resp_dist ~ sl + fSupp,  data = bitt.Comb,
              control.predictor = list(compute = TRUE),
              control.family = list(hyper =
                                   list(prec = list(prior="gaussian",
                                                    param =c(0,1)))),
              control.fixed = list(mean.intercept = 20,
                                   prec.intercept = 40^(-2),
                                   mean = list(sl = 1.3, 
                                               fSupp1 = 35,
                                               default = 0), 
                                   prec = list(sl = 0.7^(-2), 
                                               fSupp1 = 15^(-2),
                                               default = 1000)))

# Rows 1 to nrow(bitt) in `Final.Pred$summary.fitted.values` are for
# the observed data and rows nrow(bitt)+1 to  nrow(bitt) + nrow(MyData)
# are for the added data with NA.

# N <- (nrow(bitt))
# K <- nrow(MyData)
Pred <- Final.Pred$summary.fitted.values[((nrow(bitt))+1):
                                           (nrow(bitt) + 
                                              nrow(MyData)),]

# Add the relevant pieces to MyData and plot the whole thing 
MyData$mu    <- Pred[,"mean"]
MyData$selow <- Pred[,"0.025quant"]
MyData$seup  <- Pred[,"0.975quant"]

# Labels
label_supp <- c("0" = "No food supplement", 
                "1" = "With food supplement")

# Plot
# Fig. 4.15
ggplot() +
geom_jitter(data = bitt, 
                     aes(y = resp_dist, x = sl),
                     shape = 19, size = 2.5,
                     height = 0.25, width = 0.25, alpha = 0.6) +
xlab("Male standard length (mm)") +
ylab("Posterior mean response distance (cm)") +
ylim(100,225)+
theme(text = element_text(size = 13)) +
theme(panel.background = element_blank())+
theme(panel.border = element_rect(fill = NA, colour = "black", size = 1))+
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))+
geom_line(data = MyData, aes(x = sl, y = mu), size = 1)+
geom_ribbon(data = MyData,
                     aes(x = sl, ymax = seup, 
                         ymin = selow), alpha = 0.5)+
facet_grid(. ~ fSupp, scales = "fixed", space = "fixed", 
                    labeller=labeller (fSupp = label_supp)) +
    theme(strip.text = element_text(size = 12, face="italic")) +
theme(legend.position = "none")

#===================END===================#
