
#=======================================

# R code for Chapter 6 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Load packages
library(lattice)  
library(ggplot2)
library(GGally)
library(tidyverse)
library(mgcv)
library(lme4)
library(car)
library(devtools)
library(ggpubr)
library(qqplotr)
library(geiger)
library(gridExtra)
library(rlang)
# Install the latest stable version of INLA:
# install.packages("INLA", repos=c(getOption("repos"),
#                 INLA="https://inla.r-inla-download.org/R/stable"),
#                 dep=TRUE)
library(INLA)

# Also install brinla
# install_github("julianfaraway/brinla")
library(brinla)

# ip <- rownames(installed.packages())
# if (!"remotes" %in% ip) {
#   install.packages("remotes")
# }
# if (!"INLA" %in% ip) {
#   install.packages(
#     "INLA", 
#     repos = c(getOption("repos"), "https://inla.r-inla-download.org/R/stable")
#   )
# }
#remotes::install_github("inbo/inlatools")
library(inlatools)

#=======================================

# Import the data

# To import the data, change the working directory to the one
# in which the .csv file is saved, and use the read.csv function

coral <- read_csv("coralsulu.csv")

str(coral)

# 'data.frame':	32 obs. of 4 variables:
# $ site    : int  1 2 3 4 5 6 7 8 9 10 ...
# $ status  : Factor w/ 2 levels "protected","unprotected"
# $ depth   : int  8 4 5 6 8 5 8 10 10  ...
# $ species : int  50 52 86 63 49 62 43 ...

# We have 32 observations, each a separate reef

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

# The aim of this study is to understand the effect of protected status
# on coral reef diversity while accommodating a strong depth effect

#=======================================

# 2. PERFORM DATA EXPLORATION

#=======================================

# MISSING VALUES?
colSums(is.na(coral))

# site  status   depth species 
# 0       0       0       0 
#No missing data

#=======================================

# OUTLIERS

# No obvious outliers in continuous variables depth and species

#Add row numbers
coral <- coral %>%
  mutate(order = seq(1:nrow(coral)))

# Set theme
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
    geom_point(aes(y = {{Yvar}})) +
    theme_bw() +
    My_theme +
    coord_flip() +
    labs(x = "Order of Data")
}

#Choose variables
p1 <- multi_dotplot(coral, order, depth) 
p2 <- multi_dotplot(coral, order, species) 

#Create grid and plot
# Fig. 6.1
grid.arrange(p1, p2, nrow = 1)

#=======================================

# DISTRIBUTION OF THE DEPENDENT VARIABLE

#Fig. 6.2 Density plot of species

coral %>% 
  ggplot(aes(species)) +
  geom_density() +
  xlab("Number of coral species") + ylab("Density") +
  xlim(25,125) +
  My_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))

# A positively skewed distribution

#=======================================

# BALANCE
# Examine the balance of categorical variables

table(coral$status)
#  protected unprotected 
#     16          16 

# The design is balanced

#=======================================

# COLLINEARITY
# Obtain summary using the ggpairs command from the GGally library
# Fig. 6.3
coral %>% 
    ggpairs(columns = c("status","depth"), aes(colour=status, alpha = 0.8), 
            lower = list(combo = wrap("facethist", binwidth = 2))) + 
    My_theme
# No obvious problems

#=======================================

# ZEROS IN THE RESPONSE VARIABLE
round((sum(coral$species == 0) / nrow(coral))*100,0)
# No zeros

#=======================================

# RELATIONSHIPS
# Fig. 6.4

label_status <- c("protected"   = "Protected reef", 
                  "unprotected" = "Unprotected reef")

coral %>% ggplot() +
  geom_point(aes(y = species, x = depth, 
                               size = 1, alpha = 0.8)) +
  geom_smooth(method = "lm", se = FALSE, 
              aes(y = species, x = depth)) + 
  xlab("Depth (m)") + ylab("Number of coral species") +
  xlim(2,10) + ylim(25,125) +
  theme(text = element_text(size=15))  +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  facet_grid(. ~ status, 
             scales = "fixed", space = "fixed", 
             labeller=labeller (status = label_status)) +
    theme(strip.text = element_text(size = 14, face="italic")) +
  theme(legend.position = "none")

# The plot strongly suggests a depth x status interaction

#=======================================

# 3. SELECT A STATISTICAL MODEL

#=======================================

# Species number is a count, and includes zeros but
# no negative values. A Poisson distribution is an  
# appropriate starting point.

#=======================================

# 4. SPECIFY PRIORS

#=======================================

# Fit each model with default priors and with informative
# priors based on data from Waheed et al. (2015)

# An intercept of ~N(5,0.5) tau = 4

# A negative effect of depth
# ~N(-0.18,0.06) tau = 277.8

# A negative effect of unprotected status
# of similar magnitude as depth
# ~N(-0.5,0.25) tau = 16

# default for positive interaction depth:statusunprotected
# (default~N(0.05,0.05)) tau = 25

# ======================================

# 5. FIT MODELS

#=======================================

# Specify the formulae for the a priori models

f01 <- species ~ depth * status

# Models M01-10 with default priors
M01 <- inla(f01, 
            control.compute = list(dic = TRUE), 
            family = "poisson", 
            data = coral)

# Models I01 with informative priors
I01 <- inla(f01, family = "poisson", data = coral,
        control.compute = list(dic = TRUE),
          control.fixed = list(mean.intercept = 5.0,
                               prec.intercept = 0.5^(-2),
                            mean = list(depth = -0.18,
                            statusunprotected = -0.5,
                                      default = 0.05),
                            prec = list(depth = 0.06^(-2),
                            statusunprotected = 0.25^(-2),
                                      default = 1^(-2))))

#=======================================

# 6. OBTAIN THE POSTERIOR DISTRIBUTION

#=======================================

# Output for the fixed effects of M01

M01Betas <- M01$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(M01Betas, digits = 2)

#                          mean   sd 0.025quant 0.975quant
# (Intercept)              4.88 0.08       4.72       5.04
# depth                   -0.10 0.01      -0.12      -0.07
# statusunprotected       -0.73 0.13      -0.99      -0.48
# depth:statusunprotected  0.06 0.02       0.02       0.10


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
# with a posterior mean and a 95% credible interval, we can also plot
# the posterior distributions - available in the object `M10$marginals.fixed`. 

# Plot posterior distributions for fixed parameters 
# Fig. 6.5

# Model intercept (Beta1)
PosteriorBeta1.M01 <- as.data.frame(M01$marginals.fixed$`(Intercept)`)
PriorBeta1.M01     <- data.frame(x = PosteriorBeta1.M01[,"x"], 
                           y = dnorm(PosteriorBeta1.M01[,"x"],0,0))
Beta1mean.M01 <- M01Betas["(Intercept)", "mean"]
Beta1lo.M01   <- M01Betas["(Intercept)", "0.025quant"]
Beta1up.M01   <- M01Betas["(Intercept)", "0.975quant"]

beta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.M01, xmax = Beta1up.M01,
           ymin = 0, ymax = 5.1, fill = "gray88") +
  geom_line(data = PosteriorBeta1.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta1.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") + ylab("Density") +
  xlim(4.5,5.3) + ylim(0,5.1) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta1

# depth (Beta2)
PosteriorBeta2.M01 <- as.data.frame(M01$marginals.fixed$`depth`)
PriorBeta2.M01 <- data.frame(x = PosteriorBeta2.M01[,"x"], 
                       y = dnorm(PosteriorBeta2.M01[,"x"],0,0.001))
Beta2mean.M01 <- M01Betas["depth", "mean"]
Beta2lo.M01   <- M01Betas["depth", "0.025quant"]
Beta2up.M01   <- M01Betas["depth", "0.975quant"]

beta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.M01, xmax = Beta2up.M01,
         ymin = 0, ymax = 33, fill = "gray88") +
  geom_line(data = PosteriorBeta2.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for depth") + ylab("Density") +
  xlim(-0.18,0.01) + ylim(0,33) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta2mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta2

# status (Beta3)
PosteriorBeta3.M01 <- as.data.frame(M01$marginals.fixed$`statusunprotected`)
PriorBeta3.M01     <- data.frame(x = PosteriorBeta3.M01[,"x"],
                           y = dnorm(PosteriorBeta3.M01[,"x"],0,0.001))

Beta3mean.M01 <- M01Betas["statusunprotected", "mean"]
Beta3lo.M01   <- M01Betas["statusunprotected", "0.025quant"]
Beta3up.M01   <- M01Betas["statusunprotected", "0.975quant"]

beta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.M01, xmax = Beta3up.M01,
           ymin = 0, ymax = 3.5, fill = "gray88") +
  geom_line(data = PosteriorBeta3.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for status") + ylab("Density") +
  xlim(-1.5,0.25) + ylim(0,3.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta3

# 2-way interaction - `depth:statusunprotected`
PosteriorBeta4.M01 <- as.data.frame(M01$marginals.fixed$`depth:statusunprotected`)
PriorBeta4.M01     <- data.frame(x = PosteriorBeta4.M01[,"x"],
                                 y = dnorm(PosteriorBeta4.M01[,"x"],0,0))
Beta4mean.M01 <- M01Betas["depth:statusunprotected", "mean"]
Beta4lo.M01   <- M01Betas["depth:statusunprotected", "0.025quant"]
Beta4up.M01   <- M01Betas["depth:statusunprotected", "0.975quant"]

beta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.M01, xmax = Beta4up.M01,
           ymin = 0, ymax = 22, fill = "gray88") +
  geom_line(data = PosteriorBeta4.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Interaction") + ylab("Density") +
  xlim(-0.05,0.2) + ylim(0,22) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta4mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta4

# Combine plots (Fig 6.5)
 ggarrange(beta1, beta2, beta3, beta4,
                        labels = c("A", "B", "C", "D"),
                        ncol = 2, nrow = 2)


# These density functions look normally distributed because INLA assumes
# that the posterior distributions of the betas are Gaussian.

#=======================================

# B. Informative model (I01)

# Posterior mean values and 95% CI for fixed effects (informative)
I01Betas <- I01$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(I01Betas, digits = 2)
#                         mean   sd  0.025quant  0.975quant
# (Intercept)              4.89 0.07       4.74       5.03
# depth                   -0.10 0.01      -0.12      -0.07
# statusunprotected       -0.70 0.11      -0.93      -0.48
# depth:statusunprotected  0.06 0.02       0.02       0.09

# The first column shows the posterior mean, the second the
# standard deviation and the third and fourth columns show the 
# 2.5% and 97.5% quantiles of the posterior distribution.

# Plot posterior distributions for fixed parameters 
# Fig. 6.6

# Model intercept (Beta1)
PosteriorBeta1.I01 <- as.data.frame(I01$marginals.fixed$`(Intercept)`)
PriorBeta1.I01     <- data.frame(x = PosteriorBeta1.I01[,"x"], 
                                 y = dnorm(PosteriorBeta1.I01[,"x"],5,0.5))
Beta1mean.I01 <- I01Betas["(Intercept)", "mean"]
Beta1lo.I01   <- I01Betas["(Intercept)", "0.025quant"]
Beta1up.I01   <- I01Betas["(Intercept)", "0.975quant"]

Ibeta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.I01, xmax = Beta1up.I01,
           ymin = 0, ymax = 5.5, fill = "gray88") +
  geom_line(data = PosteriorBeta1.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta1.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") + ylab("Density") +
  xlim(4.4,5.4) + ylim(0,5.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta1

# depth (Beta2)
PosteriorBeta2.I01 <- as.data.frame(I01$marginals.fixed$`depth`)
PriorBeta2.I01 <- data.frame(x = PosteriorBeta2.I01[,"x"], 
                             y = dnorm(PosteriorBeta2.I01[,"x"],-0.18, 0.06))
Beta2mean.I01 <- I01Betas["depth", "mean"]
Beta2lo.I01   <- I01Betas["depth", "0.025quant"]
Beta2up.I01   <- I01Betas["depth", "0.975quant"]

Ibeta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.I01, xmax = Beta2up.I01,
           ymin = 0, ymax = 35, fill = "gray88") +
  geom_line(data = PosteriorBeta2.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for depth") + ylab("Density") +
  xlim(-0.2,0.02) + ylim(0,35) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta2mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta2

# status (Beta3)
PosteriorBeta3.I01 <- as.data.frame(I01$marginals.fixed$`statusunprotected`)
PriorBeta3.I01     <- data.frame(x = PosteriorBeta3.I01[,"x"],
                                 y = dnorm(PosteriorBeta3.I01[,"x"],-0.5, 0.25))

Beta3mean.I01 <- I01Betas["statusunprotected", "mean"]
Beta3lo.I01   <- I01Betas["statusunprotected", "0.025quant"]
Beta3up.I01   <- I01Betas["statusunprotected", "0.975quant"]

Ibeta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.I01, xmax = Beta3up.I01,
           ymin = 0, ymax = 4, fill = "gray88") +
  geom_line(data = PosteriorBeta3.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for status") + ylab("Density") +
  xlim(-1.5,0.25) + ylim(0,4) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta3

# 2-way interaction - `depth:statusunprotected`
PosteriorBeta4.I01 <- as.data.frame(I01$marginals.fixed$`depth:statusunprotected`)
PriorBeta4.I01     <- data.frame(x = PosteriorBeta4.I01[,"x"],
                                 y = dnorm(PosteriorBeta4.I01[,"x"],5e-02, 1))
Beta4mean.I01 <- I01Betas["depth:statusunprotected", "mean"]
Beta4lo.I01   <- I01Betas["depth:statusunprotected", "0.025quant"]
Beta4up.I01   <- I01Betas["depth:statusunprotected", "0.975quant"]

Ibeta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.I01, xmax = Beta4up.I01,
           ymin = 0, ymax = 25, fill = "gray88") +
  geom_line(data = PosteriorBeta4.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Interaction") + ylab("Density") +
  xlim(-0.05,0.15) + ylim(0,25) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta4mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta4

# Combine plots
ggarrange(Ibeta1, Ibeta2, Ibeta3, Ibeta4,
                        labels = c("A", "B", "C", "D"),
                          ncol = 2, nrow = 2)

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

#               fit   delta    w
# informative   252    0.00   0.57
# default       253    0.59   0.43

# Almost the same

#=======================================

# Compare with frequentist Poisson GLM

Freq <- glm(species ~ depth * status,
                      family = "poisson",
                      data = coral)

round(summary(Freq)$coef[,1:4],3)

#                         Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)                4.884      0.082  59.530    0.000
# depth                     -0.097      0.013  -7.453    0.000
# statusunprotected         -0.734      0.131  -5.593    0.000
# depth:statusunprotected    0.061      0.020   3.040    0.002

# Posterior mean values and 95% CI for fixed effects (default priors)
round(M01Betas, digits = 3)

#                          mean    sd   0.025quant 0.975quant
# (Intercept)              4.885 0.082      4.722      5.045
# depth                   -0.097 0.013     -0.123     -0.072
# statusunprotected       -0.734 0.131     -0.992     -0.477
# depth:statusunprotected  0.061 0.020      0.022      0.100


# Posterior mean values and 95% CI for fixed effects (informative priors)
round(I01Betas, digits = 3)

#                          mean    sd   0.025quant 0.975quant
# (Intercept)              4.886 0.075      4.738      5.033
# depth                   -0.098 0.012     -0.121     -0.075
# statusunprotected       -0.699 0.109     -0.913     -0.485
# depth:statusunprotected  0.057 0.017      0.023      0.090

#=======================================

# 7. MODEL CHECKS

#=======================================

# # A. Model selection

# For Bayesian GLMs we can choose among models using the 
# Deviance Information Criterion (DIC)- equivalent to AIC
# Compare alternative formulations of the model

f01 <- species ~ depth * status
f02 <- species ~ depth + status
f03 <- species ~ depth
f04 <- species ~ status

# To use DIC we must re-run the model and specify its calculation using 
# 'control.compute'

# Full model with default priors
M01.full <- inla(f01, 
                 control.compute = list(dic = TRUE), 
                 family = "poisson", 
                 data = coral)

# Model with default priors with interaction dropped
M01.1 <- inla(f02, 
                 control.compute = list(dic = TRUE), 
                 family = "poisson", 
                 data = coral)

# Model with default priors with reef status dropped
M01.2 <- inla(f03, 
                 control.compute = list(dic = TRUE), 
                 family = "poisson", 
                 data = coral)

# Model with default priors with depth dropped
M01.3 <- inla(f04, 
                 control.compute = list(dic = TRUE), 
                 family = "poisson", 
                 data = coral)


# Compare models with DIC
M01dic <- c(M01.full$dic$dic, M01.1$dic$dic, 
           M01.2$dic$dic,    M01.3$dic$dic)
DIC <- cbind(M01dic)
rownames(DIC) <- c("full","no inter","no status","no depth")
round(DIC,0)

# full         253 <- lowest
# no inter     260
# no status    321
# no depth     311

# Now with informative priors
I01.full <- inla(f01, family = "poisson", data = coral,
            control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 5.0,
                                 prec.intercept = 0.5^(-2),
                                 mean = list(depth = -0.18,
                                             statusunprotected = -0.5,
                                             default = 0.05),
                                 prec = list(depth = 0.06^(-2),
                                             statusunprotected = 0.25^(-2),
                                             default = 1.0^(-2))))

I01.1 <- inla(f02, family = "poisson", data = coral,
                 control.compute = list(dic = TRUE),
                 control.fixed = list(mean.intercept = 5.0,
                                      prec.intercept = 0.5^(-2),
                                      mean = list(depth = -0.18,
                                                  statusunprotected = -0.5,
                                                  default = 0.05),
                                      prec = list(depth = 0.06^(-2),
                                                  statusunprotected = 0.25^(-2),
                                                  default = 1.0^(-2))))

I01.2 <- inla(f03, family = "poisson", data = coral,
                 control.compute = list(dic = TRUE),
                 control.fixed = list(mean.intercept = 5.0,
                                      prec.intercept = 0.5^(-2),
                                      mean = list(depth = -0.18,
                                                  statusunprotected = -0.5,
                                                  default = 0.05),
                                      prec = list(depth = 0.06^(-2),
                                                  statusunprotected = 0.25^(-2),
                                                  default = 1.0^(-2))))

I01.3 <- inla(f04, family = "poisson", data = coral,
                 control.compute = list(dic = TRUE),
                 control.fixed = list(mean.intercept = 5.0,
                                      prec.intercept = 0.5^(-2),
                                      mean = list(depth = -0.18,
                                                  statusunprotected = -0.5,
                                                  default = 0.05),
                                      prec = list(depth = 0.06^(-2),
                                                  statusunprotected = 0.25^(-2),
                                                  default = 1.0^(-2))))

# Compare models with DIC
I01dic <- c(I01.full$dic$dic, I01.1$dic$dic, 
            I01.2$dic$dic,    I01.3$dic$dic)
DIC <- cbind(I01dic)
rownames(DIC) <- c("full","no inter","no status","no depth")
round(DIC,0)

# full         252 <- lowest 
# no inter     260
# no status    321
# no depth     311

# Compare models with non-informative and informative priors
dic2 <- c(M01.full$dic$dic, I01.full$dic$dic)
DIC2 <- cbind(dic2)
rownames(DIC2) <- c("default priors","informative priors")
round(DIC2,0)

# default priors     253
# informative priors 252

# These are essentially the same - continue with informative model

#=======================================

# B. Dispersion
# (Just for model with informative priors)

# Poisson GLMs assume that the mean and variance of the response variable 
# increase at the same rate. This assumption must be confirmed. 
# Overdispersion means that a Poisson distribution does not adequately 
# model the variance and is not appropriate for the analysis.

# Dispersion can be assessed by summing the squared Pearson residuals 
# and dividing them by the number of observations minus the degrees of 
# freedom. This value should be close to one. Values above one indicate 
# overdispersion, while values below one indicate underdispersion.

# However this approach is rather arbitrary and a better comparison 
# for assessing dispersion is to simulate data from the fitted model and 
# calculate the dispersion for each of the simulated data sets.

# Simulation can be undertaken using the 'inlatools' package

# Start by refitting the model with the 'config = TRUE' option,
# which permits us to simulate regression parameters.

I01 <- inla(f01, family = "poisson", data = coral,
               control.compute = list(config = TRUE, dic = TRUE),
             control.predictor = list(compute = TRUE), 
                 control.fixed = list(mean.intercept = 5.0,
                                      prec.intercept = 0.5^(-2),
                                   mean = list(depth = -0.18,
                                   statusunprotected = -0.5,
                                             default = 0.05),
                                   prec = list(depth = 0.06^(-2),
                                   statusunprotected = 0.25^(-2),
                                             default = 1^(-2))))

# We simulate data from the model with 'dispersion_check'

dis_pois <- dispersion_check(I01)

# Which generates the dispersion value based on the data and 
# a vector of dispersion values for 1000 simulated data sets

# The dispersion value for the data is given by:
round(dis_pois$data,2)

# This should be close to 1 (<1 underdispersion, >1 overdispersion)

# Plot dispersion values for simulated data sets 
# (with dispersion value for the data added as a dashed line)
# Fig. 6.7

pois_plot <- ggplot() +
    geom_density(aes(dis_pois$model)) +
    geom_vline(xintercept = dis_pois$data, linetype = "dashed") +
    xlab("Dispersion") + ylab("Density") +
    xlim(0,2.1) + ylim(0,1.8) +
    theme(text = element_text(size=13)) +
    theme(panel.background = element_blank()) +
    theme(panel.border = element_rect(fill = NA, 
                                      colour = "black", size = 1)) +
    theme(strip.background = element_rect
          (fill = "white", color = "white", size = 1))
pois_plot

# The plot displays the density of the dispersion of the simulated data sets.
# The vertical dashed line shows the dispersion of the original data. 
# If the dispersion value for the data (dashed line) falls within
# the distribution of dispersion values for simulated data sets
# the model is not under- or overdispersed

# In this case there is clear evidence of overdispersion

#Why is the model overdispersed?
#1. Missing covariates or interactions
#2. Zero inflation
#3. Influential outliers
#4. Non-independence
#5. Wrong link function
#6. Non-linearity
#7. True overdispersion

# Fit the model with negative binomial distribution
I01.nb <- inla(f01, family = "nbinomial", data = coral,
            control.compute = list(config = TRUE, dic = TRUE),
            control.predictor = list(compute = TRUE), 
            control.fixed = list(mean.intercept = 5.0,
                                 prec.intercept = 0.5^(-2),
                                 mean = list(depth = -0.18,
                                             statusunprotected = -0.5,
                                             default = 0.05),
                                 prec = list(depth = 0.06^(-2),
                                             statusunprotected = 0.25^(-2),
                                             default = 1^(-2))))


# Then repeat simulation exercise
dis_nbin <- dispersion_check(I01.nb)

# The dispersion value for the data is given by:
round(dis_nbin$data,2)
# This should be close to 1.

# Plot
nb_plot <- ggplot() +
  geom_density(aes(dis_nbin$model)) +
  geom_vline(xintercept = dis_nbin$data, linetype = "dashed") +
  xlab("Dispersion") + ylab("") +
  xlim(0,2.5) + ylim(0,1.8) +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
# Fig. 6.8

FigDisp <- ggarrange(pois_plot, nb_plot,
                     labels = c("A", "B"),
                       ncol = 2, nrow = 1)
FigDisp

# The model with the Poisson distribution (A) has clear overdispersion.
# The plot of the negative binomial model (B) indicates no over- or 
# underdispersion. 

# Compare DICs for Poisson and negative binomial models
dic3 <- c(I01$dic$dic, I01.nb$dic$dic)
DIC3 <- cbind(dic3)
rownames(DIC3) <- c("Poisson","negative binomial")
round(DIC3,0)

# Poisson            252
# negative binomial  248

# The negative binomial model gives the better fit

#=======================================

# C. Posterior predictive check of negative binomial model

I01.pred <- inla(f01, family = "nbinomial", data = coral,
               control.predictor = list(link = 1,
                                     compute = TRUE),
               control.compute = list(dic = TRUE, 
                                      cpo = TRUE),
               control.fixed = list(mean.intercept = 5.0,
                                    prec.intercept = 0.5^(-2),
                                 mean = list(depth = -0.18,
                                 statusunprotected = -0.5,
                                           default = 0.05),
                                 prec = list(depth = 0.06^(-2),
                                 statusunprotected = 0.25^(-2),
                                           default = 1^(-2))))

ppp <- vector(mode = "numeric", length = nrow(coral))
for(i in (1:nrow(coral))) {
  ppp[i] <- inla.pmarginal(q = coral$species[i],
                    marginal = I01.pred$marginals.fitted.values[[i]])
}

# Fig. 6.9
ggplot() +
  geom_histogram(aes(ppp), binwidth = 0.09, 
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


# Problem

#=======================================

# D. Cross-Validation Model Checking

# Use CPO and PIT
sum(I01.pred$cpo$failure)
# [1] 0
# 0 indicates CPO is reliable

#Extract pit values
PIT <- (I01.pred$cpo$pit)

#And plot
Pit1 <- ggplot() +
  geom_histogram(aes(PIT), binwidth = 0.11, 
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
# Fig. 6.10
ggarrange(Pit1, Pit2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

#PIT values are uniform

#=======================================

# E. Bayesian residual analysis

# Plot residuals versus fitted values & each covariate in the model

# Obtain fitted values
# (Ensure  compute = TRUE in `control.predictor`)

Fit <- I01.pred$summary.fitted.values[, "mean"]

# Calculate residuals
Res <- coral$species - Fit
ResPlot <- cbind.data.frame(Fit,Res,coral$depth,coral$status)

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
Res2 <- ggplot(ResPlot, aes(x=coral$depth, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Depth (m)") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

Res3 <- ggplot(ResPlot, aes(x=coral$status, y=Res)) + 
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Reef status") +
  theme(text = element_text(size=13)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
# Fig. 6.11
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
Final <- inla(f01, family = "nbinomial", data = coral,
          control.compute = list(config = TRUE),
        control.predictor = list(compute=TRUE), 
            control.fixed = list(mean.intercept = 5.0,
                                 prec.intercept = 0.5^(-2),
                              mean = list(depth = -0.18,
                              statusunprotected = -0.5,
                                        default = 0.05),
                              prec = list(depth = 0.06^(-2),
                              statusunprotected = 0.25^(-2),
                                        default = 1^(-2))))

# Posterior mean values and 95% CI for fixed effects
BetasFinal <- Final$summary.fixed[,c("mean", "sd", 
                                     "0.025quant", 
                                     "0.975quant")] 
round(BetasFinal, digits = 2)

#                          mean   sd 0.025quant 0.975quant
# (Intercept)              4.90 0.10       4.69       5.11
# depth                   -0.10 0.02      -0.13      -0.07
# statusunprotected       -0.69 0.13      -0.95      -0.42
# depth:statusunprotected  0.06 0.02       0.02       0.10


## Model interpretation

# Negative effect of depth on species number
# Positive effect of protected status on species number
# Interaction of depth and status (depth relationship steeper
# in protected sites)

#=======================================

# 9. VISUALISE THE RESULTS

#=======================================

#Plot figure Fig. 6.12
# MyData <- ddply(coral, .(status), summarize,
#                  depth = seq(from = min(depth), 
#                     to = max(depth), length = 50))

MyData <- expand.grid(
  status = c("protected", "unprotected"),
  depth = seq(from = min(coral$depth), 
              to = max(coral$depth), length = 50))

# 2. Make a design matrix
Xmat <- model.matrix(~ depth * status, data = MyData)
Xmat <- as.data.frame(Xmat)

lcb <- inla.make.lincombs(Xmat)

# Re-run the model in R-INLA using the combined data set, ensuring
# that `compute = TRUE` is selected in the `control.predictor` argument
Final.Pred <- inla(f01, family = "nbinomial", data = coral,
              lincomb = lcb,
              control.inla = list(lincomb.derived.only = TRUE),
              control.predictor = list(compute = TRUE), 
              control.fixed = list(mean.intercept = 5.0,
                                   prec.intercept = 0.5^(-2),
                                mean = list(depth = -0.18,
                                statusunprotected = -0.5,
                                          default = 0.05),
                                prec = list(depth = 0.06^(-2),
                                statusunprotected = 0.25^(-2),
                                          default = 1^(-2))))

# Run loop to get mu, selo and seup
Pred.marg <- Final.Pred$marginals.lincomb.derived

for (i in 1:nrow(MyData)){
  MyData$mu[i]  <- inla.emarginal(exp, Pred.marg[[i]])
  lo.up <- inla.qmarginal(c(0.025, 0.975), 
                          inla.tmarginal(exp, Pred.marg[[i]]))
  MyData$selo[i] <- lo.up[1]
  MyData$seup[i] <- lo.up[2]    	
}               

# Labels
label_status <- c("protected" = "Protected reef", 
                  "unprotected" = "Unprotected reef")

# Plot
ggplot() + 
  geom_jitter(data = coral, aes(y = species, x = depth),
              shape = 19, size = 2.5, height = 1, 
              width = 0.1, alpha = 0.7) +
  xlab("Water depth (m)") + 
  ylab("Posterior mean number of coral species") +
  ylim(30,120) + xlim(2,10.5) + 
  theme(text = element_text(size = 14)) + 
  theme(panel.background = element_blank()) + 
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) + 
  theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1)) +
  geom_line(data = MyData, aes(x = depth, y = mu), size = 1) +
  geom_ribbon(data = MyData, aes(x = depth, 
                                 ymax = seup, ymin = selo), alpha = 0.5) +
  theme(strip.text = element_text(size = 14, face="italic")) +
  facet_grid(. ~ status, scales = "fixed", space = "fixed", 
                    labeller=labeller (status = label_status))

#===================END===================#
