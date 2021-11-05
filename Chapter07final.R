
#=======================================

# R code for Chapter 7 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Load packages
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
library(INLA)
library(brinla)
library(inlatools)
library(gridExtra)
library(rlang)
#=======================================

# Import the data

# To import the data, change the working directory to the one
# in which the .csv file is saved, and use the read.csv function

cc <- read_csv("cuckoo.csv")

str(cc)

# 'data.frame':	18 obs. of  5 variables:
# $ nest   : int  1 2 3 4 5 6 7 8 9 10 ...
# $ water  : int  168 14 114 247 184 159 100 ...
# $ tree   : int  54 73 12 55 21 48 3 89 17 66 ...
# $ aggress: Factor w/ 2 levels "high","low":  ...
# $ egg    : int  1 0 1 0 1 0 1 0 1 0 ...

# We have 18 observations, each a separate great reed warbler nest

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

# The aim of this study is to understand what ecological factors predict
# parasitism of great reed warbler nests by the common cuckoo

#=======================================

# 2. PERFORM DATA EXPLORATION

#=======================================

# MISSING VALUES?
colSums(is.na(cc))

# nest   water    tree aggress     egg 
# 0       0       0       0       0 
#No missing data

#=======================================

# OUTLIERS

cc <- cc %>%
  mutate(order = seq(1:nrow(cc)))

# Set preferred theme
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
p1 <- multi_dotplot(cc, order, water) 
p2 <- multi_dotplot(cc, order, tree) 

#Create grid
# Fig. 7.1
grid.arrange(p1, p2, nrow = 1)

# No obvious outliers in continuous variables water and tree

#=======================================

# DISTRIBUTION OF THE DEPENDENT VARIABLE

# The response variable is binomial - whether a nest is 
# parasitised ('1') or not ('0')

# Balance of response variable
table(cc$egg)

# 0  1 
# 8 10 
# 10 nests contained at least 1 cuckoo egg

#=======================================

# BALANCE

table(cc$aggress)
# high low 
# 9    9
# Perfect balance of aggressive parents

#=======================================

# COLLINEARITY
# Fig. 7.2
cc %>% 
    ggpairs(columns = c("water", "tree", "aggress"), 
            aes(colour=aggress, alpha = 0.8), 
            lower = list(combo = wrap("facethist", binwidth = 15))) + 
    My_theme

# No obvious problems

#=======================================

# ZEROS IN THE RESPONSE VARIABLE

round((sum(cc$egg == 0) / nrow(cc))*100,0)
# 44% zeros - but binomial response variable

#=======================================

# RELATIONSHIPS
# Fig. 7.3

label_agg <- c("low"   = "Low aggression", 
               "high"  = "High aggression")

water_plot <- ggplot() +
  geom_jitter(data = cc, position = position_jitter
              (width = 0.05, height = 0.05),
              aes(y = egg, x = water, size = 1, alpha = 0.8)) +
  xlab("Distance to open water (m)") + 
  ylab("Probability of cuckoo egg") +
  ylim(-0.1,1.1) +
  theme(text = element_text(size=11))  +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
        colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  facet_grid(. ~ aggress, 
             scales = "fixed", space = "fixed", 
             labeller=labeller (aggress = label_agg)) +
  theme(strip.text = element_text(size = 12, face="italic")) +
  theme(legend.position = "none")

tree_plot <- ggplot() +
  geom_jitter(data = cc, position = position_jitter
              (width = 0.05, height = 0.05),
              aes(y = egg, x = tree, size = 1, alpha = 0.8)) +
  xlab("Distance to tree (m)") + 
  ylab("") +
  ylim(-0.1,1.1) +
  theme(text = element_text(size=11))  +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  facet_grid(. ~ aggress, 
             scales = "fixed", space = "fixed", 
             labeller=labeller (aggress = label_agg)) +
  theme(strip.text = element_text(size = 12, face="italic")) +
  theme(legend.position = "none")

# Combine plots
ggarrange(water_plot, tree_plot,
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

# The plot suggests:
# 1. Greater probability of cuckoo parasitism with low aggression parents
# 2. There may also be a negative effect of distance to tree
# 3. No clear relationship with distance to open water

#=======================================

# 3. SELECT A STATISTICAL MODEL

#=======================================

# Probability of cuckoo parasitism is strictly binomial. 
# A Bernoulli distribution is an appropriate starting point.

#=======================================

# 4. SPECIFY PRIORS

#=======================================

# Fit one model with default priors and second with informative
# priors based on pilot observations conducted from previous year

# Import pilot data:
# ccpilot <- read.csv(file = "ccpilot.csv", 
#                     header = TRUE, 
#                     dec = ".", 
#                     stringsAsFactors = TRUE)

ccpilot <- read_csv("ccpilot.csv")

#Obtain structure of the data
str(ccpilot)

# 'data.frame':	7 obs. of  5 variables:
# $ nest   : Factor w/ 7 levels "a","b","c","d",
# $ water  : int  100 136 175 69 66 84 133
# $ tree   : int  42 97 61 54 43 21 31
# $ aggress: Factor w/ 2 levels "high","low"
# $ egg    : int  1 0 0 1 0 1 1

# Formulate model
f01 <- egg ~ water + tree + aggress

# Run with INLA
P01 <- inla(f01,
            family = "binomial", 
            Ntrials = 1,
            data = ccpilot)

# Obtain summary of fixed effects
P01Betas <- P01$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(P01Betas, digits = 2)

#             mean     sd 0.025quant 0.975quant
# (Intercept)  1.68  1.88      -1.99       5.39
# water       -0.01  0.02      -0.04       0.02
# tree        -0.03  0.03      -0.10       0.03
# aggresslow   1.83  1.60      -1.13       5.14


# Priors for model

# intercept  ~  1.68 (sd = 1.88) tau = 0.28
# water      ~ -0.01 (sd = 0.02) tau = 2500
# tree       ~ -0.03 (sd = 0.03) tau = 1111
# aggresslow ~  1.83 (sd = 1.60) tau = 0.39

# ======================================

# 5. FIT MODELS

#=======================================

# Specify model formula

f01 <- egg ~ water + tree + aggress

# Model M01 with default priors
M01 <- inla(f01, family = "binomial", 
                 Ntrials = 1, data = cc,
                 control.compute = list(dic = TRUE))

# Models I01 with informative priors (priors based on pilot)
I01 <- inla(f01, family = "binomial", Ntrials = 1, data = cc,
            control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 1.68,
                                 prec.intercept = 1.88^(-2),
                              mean = list(water = -0.01,
                                           tree = -0.03,
                                     aggresslow = 1.83),
                              prec = list(water = 0.02^(-2),
                                          tree  = 0.03^(-2),
                                     aggresslow = 1.6^(-2))))

#=======================================

# 6. OBTAIN THE POSTERIOR DISTRIBUTION

#=======================================

# Output for the fixed effects of M01

M01Betas <- M01$summary.fixed[,c("mean", "mode", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(M01Betas, digits = 2)

#               mean   mode    sd 0.025quant 0.975quant
#(Intercept)  195.84 164.09 45.04     131.41     301.86
# water         0.37   0.27  0.18       0.08       0.79
# tree         -5.41  -4.43  1.41      -8.74      -3.37
# aggresslow   98.63  90.90 13.35      78.45     129.66


# The first column shows the posterior mean, the second the
# mode, the third the standard deviation and the fourth and fifth columns show the 
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

# Model intercept (Beta1)
PosteriorBeta1.M01 <- as.data.frame(M01$marginals.fixed$`(Intercept)`)
PriorBeta1.M01     <- data.frame(x = PosteriorBeta1.M01[,"x"], 
                           y = dnorm(PosteriorBeta1.M01[,"x"],0,0))
Beta1mean.M01 <- M01Betas["(Intercept)", "mode"]
Beta1lo.M01   <- M01Betas["(Intercept)", "0.025quant"]
Beta1up.M01   <- M01Betas["(Intercept)", "0.975quant"]

beta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.M01, xmax = Beta1up.M01,
           ymin = 0, ymax = 0.011, fill = "gray88") +
  geom_line(data = PosteriorBeta1.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta1.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") + ylab("Density") +
  xlim(0,410) + ylim(0,0.011) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta1

# water (Beta2)
PosteriorBeta2.M01 <- as.data.frame(M01$marginals.fixed$`water`)
PriorBeta2.M01 <- data.frame(x = PosteriorBeta2.M01[,"x"], 
                       y = dnorm(PosteriorBeta2.M01[,"x"],0,0.001))
Beta2mean.M01 <- M01Betas["water", "mode"]
Beta2lo.M01   <- M01Betas["water", "0.025quant"]
Beta2up.M01   <- M01Betas["water", "0.975quant"]

beta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.M01, xmax = Beta2up.M01,
         ymin = 0, ymax = 2.5, fill = "gray88") +
  geom_line(data = PosteriorBeta2.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for water") + ylab("Density") +
  xlim(-0.5,1.5) + ylim(0,2.5) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta2mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta2

# tree (Beta3)
PosteriorBeta3.M01 <- as.data.frame(M01$marginals.fixed$`tree`)
PriorBeta3.M01     <- data.frame(x = PosteriorBeta3.M01[,"x"],
                           y = dnorm(PosteriorBeta3.M01[,"x"],0,0.001))

Beta3mean.M01 <- M01Betas["tree", "mode"]
Beta3lo.M01   <- M01Betas["tree", "0.025quant"]
Beta3up.M01   <- M01Betas["tree", "0.975quant"]

beta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.M01, xmax = Beta3up.M01,
           ymin = 0, ymax = 0.35, fill = "gray88") +
  geom_line(data = PosteriorBeta3.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for tree") + ylab("Density") +
  xlim(-12,2) + ylim(0,0.35) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta3

# aggresslow (beta4)
PosteriorBeta4.M01 <- as.data.frame(M01$marginals.fixed$`aggresslow`)
PriorBeta4.M01     <- data.frame(x = PosteriorBeta4.M01[,"x"],
                                 y = dnorm(PosteriorBeta4.M01[,"x"],0,0))
Beta4mean.M01 <- M01Betas["aggresslow", "mode"]
Beta4lo.M01   <- M01Betas["aggresslow", "0.025quant"]
Beta4up.M01   <- M01Betas["aggresslow", "0.975quant"]

beta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.M01, xmax = Beta4up.M01,
           ymin = 0, ymax = 0.038, fill = "gray88") +
  geom_line(data = PosteriorBeta4.M01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.M01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Aggression") + ylab("Density") +
  xlim(40,160) + ylim(0,0.038) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta4mean.M01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                  colour = "black", size = 1)) +
  theme(strip.background = element_rect
       (fill = "white", color = "white", size = 1))
# beta4

# Combine plots
# Fig. 7.4
ggarrange(beta1, beta2, beta3, beta4,
                        labels = c("A", "B", "C", "D"),
                        ncol = 2, nrow = 2)

#=======================================

# Informative model (I01)

# Posterior mean values and 95% CI for fixed effects (informative)
I01Betas <- I01$summary.fixed[,c("mean", "mode", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(I01Betas, digits = 2)
#              mean  mode   sd 0.025quant 0.975quant
#(Intercept)   2.90  2.85 1.23       0.54       5.37
# water       -0.01 -0.01 0.01      -0.03       0.01
# tree        -0.06 -0.06 0.02      -0.11      -0.03
# aggresslow   2.47  2.37 1.08       0.44       4.69

# The first column shows the posterior mean, the second the mode,
# the third the standard deviation and the fourth and fifth columns
# show the 2.5% and 97.5% quantiles of the posterior distribution.

# Plot posterior distributions for fixed parameters 

# Model intercept (Beta1)
PosteriorBeta1.I01 <- as.data.frame(I01$marginals.fixed$`(Intercept)`)
PriorBeta1.I01     <- data.frame(x = PosteriorBeta1.I01[,"x"], 
                                 y = dnorm(PosteriorBeta1.I01[,"x"],1.68,1.88))
Beta1mean.I01 <- I01Betas["(Intercept)", "mode"]
Beta1lo.I01   <- I01Betas["(Intercept)", "0.025quant"]
Beta1up.I01   <- I01Betas["(Intercept)", "0.975quant"]

Ibeta1 <- ggplot() +
  annotate("rect", xmin = Beta1lo.I01, xmax = Beta1up.I01,
           ymin = 0, ymax = 0.38, fill = "gray88") +
  geom_line(data = PosteriorBeta1.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta1.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Intercept") + ylab("Density") +
  xlim(-3,7.5) + ylim(0,0.38) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta1mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta1

# water (Beta2)
PosteriorBeta2.I01 <- as.data.frame(I01$marginals.fixed$`water`)
PriorBeta2.I01 <- data.frame(x = PosteriorBeta2.I01[,"x"], 
                             y = dnorm(PosteriorBeta2.I01[,"x"], -0.01, 0.02))
Beta2mean.I01 <- I01Betas["water", "mode"]
Beta2lo.I01   <- I01Betas["water", "0.025quant"]
Beta2up.I01   <- I01Betas["water", "0.975quant"]

Ibeta2 <- ggplot() +
  annotate("rect", xmin = Beta2lo.I01, xmax = Beta2up.I01,
           ymin = 0, ymax = 45, fill = "gray88") +
  geom_line(data = PosteriorBeta2.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta2.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for water") + ylab("Density") +
  xlim(-0.06,0.05) + ylim(0,45) +
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
PosteriorBeta3.I01 <- as.data.frame(I01$marginals.fixed$`tree`)
PriorBeta3.I01     <- data.frame(x = PosteriorBeta3.I01[,"x"],
                                 y = dnorm(PosteriorBeta3.I01[,"x"],-0.03, 0.03))

Beta3mean.I01 <- I01Betas["tree", "mode"]
Beta3lo.I01   <- I01Betas["tree", "0.025quant"]
Beta3up.I01   <- I01Betas["tree", "0.975quant"]

Ibeta3 <- ggplot() +
  annotate("rect", xmin = Beta3lo.I01, xmax = Beta3up.I01,
           ymin = 0, ymax = 22, fill = "gray88") +
  geom_line(data = PosteriorBeta3.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta3.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Slope for tree") + ylab("Density") +
  xlim(-0.15,0.05) + ylim(0,22) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = Beta3mean.I01, linetype = "dashed") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))
# Ibeta3

# aggresslow (beta4)

PosteriorBeta4.I01 <- as.data.frame(I01$marginals.fixed$`aggresslow`)
PriorBeta4.I01     <- data.frame(x = PosteriorBeta4.I01[,"x"],
                                 y = dnorm(PosteriorBeta4.I01[,"x"],1.83, 1.6))
Beta4mean.I01 <- I01Betas["aggresslow", "mode"]
Beta4lo.I01   <- I01Betas["aggresslow", "0.025quant"]
Beta4up.I01   <- I01Betas["aggresslow", "0.975quant"]

Ibeta4 <- ggplot() +
  annotate("rect", xmin = Beta4lo.I01, xmax = Beta4up.I01,
           ymin = 0, ymax = 0.42, fill = "gray88") +
  geom_line(data = PosteriorBeta4.I01,
            aes(y = y, x = x), lwd = 1.2) +
  geom_line(data = PriorBeta4.I01,
            aes(y = y, x = x), color = "gray55", lwd = 1.2) +
  xlab("Aggression") + ylab("Density") +
  xlim(-2.5,7.5) + ylim(0,0.42) +
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
# Fig. 7.5
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

#               fit    delta    w
# informative   13.12  0.00     1
# default       40.83 27.72     0

# Informative model much more probable

#=======================================

# Compare with frequentist Poisson GLM

Freq <- glm(egg ~ water + tree + aggress,
                  family = binomial(link = "logit"),
                  data = cc)

round(summary(Freq)$coef[,1:4],2)

#             Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)   185.18  480365.37       0        1
# water           0.39    3689.23       0        1
# tree           -5.71   13125.29       0        1
# aggresslow     82.41  346541.08       0        1

# This looks strange....

# Posterior mean values and 95% CI for fixed effects (default priors)
round(M01Betas, digits = 2)

#               mean    sd   0.025quant 0.975quant
# (Intercept)  197.29 44.98     132.92     303.15
# water          0.36  0.18       0.08       0.79
# tree          -5.45  1.41      -8.77      -3.42
# aggresslow    98.42 13.34      78.25     129.44

# Posterior mean values and 95% CI for fixed effects (informative priors)
round(I01Betas, digits = 2)

#              mean    sd   0.025quant 0.975quant
# (Intercept)  2.90 1.23       0.54       5.37
# water       -0.01 0.01      -0.03       0.01
# tree        -0.06 0.02      -0.11      -0.03
# aggresslow   2.47 1.08       0.44       4.69

#=======================================

# 7. MODEL CHECKS

#=======================================

# Model validation of binomial GLMs is difficult 

# A. Model selection

# For Bayesian GLMs we can choose among models using the 
# Deviance Information Criterion (DIC)- equivalent to AIC
# Compare alternative formulations of the model

# However, is there justification for model selection? All variables
# are considered to play a role in influencing risk of parasitism.

f01 <- egg ~ water + tree + aggress
f02 <- egg ~ water + tree
f03 <- egg ~ water + aggress
f04 <- egg ~ tree + aggress

# To use DIC we must re-run the model and specify its calculation using 
# 'control.compute'

# Compare models with informative priors
I01.full <- inla(f01, family = "binomial", Ntrials = 1, data = cc,
                     control.compute = list(dic = TRUE),
                     control.fixed = list(mean.intercept = 1.68,
                                          prec.intercept = 1.88^(-2),
                                          mean = list(water = -0.01,
                                                      tree = -0.03,
                                                      aggresslow = 1.83),
                                          prec = list(water = 0.02^(-2),
                                                      tree  = 0.03^(-2),
                                                      aggresslow = 1.6^(-2))))

I01.1 <- inla(f02, family = "binomial", Ntrials = 1, data = cc,
              control.compute = list(dic = TRUE),
              control.fixed = list(mean.intercept = 1.68,
                                   prec.intercept = 1.88^(-2),
                                   mean = list(water = -0.01,
                                               tree = -0.03,
                                               aggresslow = 1.83),
                                   prec = list(water = 0.02^(-2),
                                               tree  = 0.03^(-2),
                                               aggresslow = 1.6^(-2))))

I01.2 <- inla(f03, family = "binomial", Ntrials = 1, data = cc,
              control.compute = list(dic = TRUE),
              control.fixed = list(mean.intercept = 1.68,
                                   prec.intercept = 1.88^(-2),
                                   mean = list(water = -0.01,
                                               tree = -0.03,
                                               aggresslow = 1.83),
                                   prec = list(water = 0.02^(-2),
                                               tree  = 0.03^(-2),
                                               aggresslow = 1.6^(-2))))

I01.3 <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
              control.compute = list(dic = TRUE),
              control.fixed = list(mean.intercept = 1.68,
                                   prec.intercept = 1.88^(-2),
                                   mean = list(water = -0.01,
                                               tree = -0.03,
                                               aggresslow = 1.83),
                                   prec = list(water = 0.02^(-2),
                                               tree  = 0.03^(-2),
                                               aggresslow = 1.6^(-2))))

# Compare models with the DIC
I01dic <- c(I01.full$dic$dic, I01.1$dic$dic, 
            I01.2$dic$dic,    I01.3$dic$dic)
DIC <- cbind(I01dic)
rownames(DIC) <- c("water + tree + aggress", 
                   "water + tree",
                   "water + aggress", 
                   "tree + aggress")
round(DIC,1)

# water + tree + aggress   13.1
# water + tree             15.9
# water + aggress          23.3
# tree + aggress           12.5 <- lowest

# Full model and with water dropped are the same

# Further refine best model
f05 <- egg ~ tree
f06 <- egg ~ aggress

I01.4 <- inla(f05, family = "binomial", Ntrials = 1, data = cc,
              control.compute = list(dic = TRUE),
              control.fixed = list(mean.intercept = 1.68,
                                   prec.intercept = 1.88^(-2),
                                   mean = list(water = -0.01,
                                               tree = -0.03,
                                               aggresslow = 1.83),
                                   prec = list(water = 0.02^(-2),
                                               tree  = 0.03^(-2),
                                               aggresslow = 1.6^(-2))))

I01.5 <- inla(f06, family = "binomial", Ntrials = 1, data = cc,
              control.compute = list(dic = TRUE),
              control.fixed = list(mean.intercept = 1.68,
                                   prec.intercept = 1.88^(-2),
                                   mean = list(water = -0.01,
                                               tree = -0.03,
                                               aggresslow = 1.83),
                                   prec = list(water = 0.02^(-2),
                                               tree  = 0.03^(-2),
                                               aggresslow = 1.6^(-2))))

# Compare models with DIC
I01dic2 <- c(I01.3$dic$dic, I01.4$dic$dic, I01.5$dic$dic)
DIC2 <- cbind(I01dic2)
rownames(DIC2) <- c("tree + aggress ",
                    "tree","aggress")
round(DIC2,1)

# tree + aggress     12.5 <- lowest
# tree               14.8
# aggress            24.2

# No improvement with further model selection

#=======================================

# B. Posterior predictive check of informative model

I01.pred <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
                 control.predictor = list(link = 1,
                                          compute = TRUE),
                 control.compute = list(dic = TRUE, 
                                        cpo = TRUE),
                 control.fixed = list(mean.intercept = 1.68,
                                      prec.intercept = 1.88^(-2),
                                      mean = list(water = -0.01,
                                                  tree = -0.03,
                                                  aggresslow = 1.83),
                                      prec = list(water = 0.02^(-2),
                                                  tree  = 0.03^(-2),
                                                  aggresslow = 1.6^(-2))))

ppp <- vector(mode = "numeric", length = nrow(cc))
for(i in (1:nrow(cc))) {
  ppp[i] <- inla.pmarginal(q = cc$egg[i],
                           marginal = I01.pred$marginals.fitted.values[[i]])
}

# Do model predictions and observed correspond?
# Fig. 7.6

cc$ppp <- ppp-0.04
ggplot() + 
  geom_vline(xintercept = 1:18, linetype = "dotted") +
  geom_point(data = cc, aes(y = egg, x = nest),
              shape = 19, size = 4, colour = "black",) +
  geom_point(data = cc, aes(y = ppp, x = nest),
              shape = 19, size = 4, colour = "gray60",) +
  ylab("ppp/egg") + xlab("Nest number") +
  scale_x_continuous(breaks = c(1:18)) +
  theme(text = element_text(size=15)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))


# Yes, perfect correspondance

#=======================================

# C. Bayesian residual analysis

# Plot residuals versus fitted values & each covariate in the model

# Obtain fitted values
# (Ensure  compute = TRUE in `control.predictor`)

Fit <- I01.pred$summary.fitted.values[, "mean"]

# Calculate residuals
Res     <- cc$egg - Fit
ResPlot <- cbind.data.frame(Fit,Res,cc$tree,cc$aggress)

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
Res2 <- ggplot(ResPlot, aes(x=cc$tree, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Nearest tree (m)") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

Res3 <- ggplot(ResPlot, aes(x=cc$aggress, y=Res)) + 
  geom_boxplot(fill = "grey88", colour = "black") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("Aggression") +
  theme(text = element_text(size=13)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
               colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1)) +
  theme(legend.position = "none")

# Combine plots
# Fig. 7.7
ggarrange(Res1, Res2, Res3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)

# These are difficult to interpret, unless there are many observations

# =======================================
  
# F. Prior sensitivity analysis

# Model with priors unchanged
I01.un <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
                control.predictor = list(link = 1,
                                         compute = TRUE),
                control.compute = list(dic = TRUE, 
                                       cpo = TRUE),
                control.fixed = list(mean.intercept = 1.68,
                                     prec.intercept = 1.88^(-2),
                                     mean = list(tree = -0.03,
                                                 aggresslow = 1.83),
                                     prec = list(tree  = 0.03^(-2),
                                                 aggresslow = 1.6^(-2))))

Betas.un <- I01.un$summary.fixed[,c("mean", "sd", 
                                    "0.025quant", 
                                    "0.975quant")] 
round(Betas.un, 2)

#              mean   sd 0.025quant 0.975quant
# (Intercept)  2.27 1.03       0.32       4.34
# tree        -0.07 0.02      -0.11      -0.03
# aggresslow   2.11 1.01       0.21       4.17

#== 1. Increase priors by 20% ==
#   Intercept from 1.68 to 2.02 (sd 1.88 to 2.26)
#   tree from -0.03 to -0.04 (sd 0.03 to 0.036)
#   aggresslow from 1.83 to 2.2 (sd 1.60 to 1.92)


M1.plus20 <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
                  control.predictor = list(link = 1,
                                           compute = TRUE),
                  control.compute = list(dic = TRUE, 
                                         cpo = TRUE),
                  control.fixed = list(mean.intercept = 2.02,
                                       prec.intercept = 2.26^(-2),
                                       mean = list(tree = -0.04,
                                                   aggresslow = 2.2),
                                       prec = list(tree  = 0.036^(-2),
                                                   aggresslow = 1.92^(-2))))
  
Betas.plus20 <- M1.plus20$summary.fixed[,c("mean", "sd", 
                                           "0.025quant", 
                                           "0.975quant")] 
round(Betas.plus20, 2)

#              mean   sd 0.025quant 0.975quant
# (Intercept)  2.83 1.18       0.61       5.23
# tree        -0.08 0.02      -0.13      -0.04
# aggresslow   2.38 1.14       0.25       4.73
 

#== 2. Decrease priors by 20% ==
#   Intercept from 1.68 to 1.34 (sd 1.88 to 1.50)
#   tree from -0.03 to -0.024 (sd 0.03 to 0.024)
#   aggresslow from 1.83 to 1.46 (sd 1.60 to 1.28)

M1.minus20 <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
                  control.predictor = list(link = 1,
                                           compute = TRUE),
                  control.compute = list(dic = TRUE, 
                                         cpo = TRUE),
                  control.fixed = list(mean.intercept = 1.34,
                                       prec.intercept = 1.50^(-2),
                                       mean = list(tree = -0.024,
                                                   aggresslow = 1.46),
                                       prec = list(tree  = 0.024^(-2),
                                                   aggresslow = 1.28^(-2))))

#Re-run model and obtain estimates of betas
Betas.minus20 <- M1.minus20$summary.fixed[,c("mean", "sd", 
                                             "0.025quant", 
                                             "0.975quant")] 
round(Betas.minus20, digits = 2)

#              mean   sd 0.025quant 0.975quant
# (Intercept)  1.77 0.88       0.09       3.52
# tree        -0.05 0.02      -0.09      -0.02
# aggresslow   1.83 0.87       0.17       3.59


#=======================================

# 8. INTERPRET AND PRESENT MODEL OUTPUT

#=======================================

# The final model is: 
Final <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
                       control.compute = list(config = TRUE),
                    control.predictor = list(compute = TRUE),
                 control.fixed = list(mean.intercept = 1.68,
                                      prec.intercept = 1.88^(-2),
                                    mean = list(tree = -0.03,
                                          aggresslow = 1.83),
                                   prec = list(tree  = 0.03^(-2),
                                          aggresslow = 1.6^(-2))))

# Posterior mean values and 95% CI for fixed effects
BetasFinal <- Final$summary.fixed[,c("mean", "sd", 
                                     "0.025quant", 
                                     "0.975quant")] 
round(BetasFinal, digits = 2)

#              mean   sd 0.025quant 0.975quant
# (Intercept)  2.27 1.03       0.32       4.34
# tree        -0.07 0.02      -0.11      -0.03
# aggresslow   2.11 1.01       0.21       4.17


## Model interpretation

# Negative effect of distance to nearest tree on probability of cuckoo parasitism
# Positive effect of low parental aggression of cuckoo parasitism

#=======================================

# 9. VISUALISE THE RESULTS

#=======================================

#Plot figure Fig. 7.8
# MyData <- ddply(cc, .(aggress), summarize,
#                  tree = seq(from = min(tree), 
#                     to = max(tree), length = 50))

MyData <- expand.grid(
  aggress = c("high", "low"),
  tree = seq(from = min(cc$tree), 
              to = max(cc$tree), length = 50))

# Make a design matrix
Xmat <- model.matrix(~ tree + aggress, data = MyData)
Xmat <- as.data.frame(Xmat)
lcb <- inla.make.lincombs(Xmat)

# Re-run the model in R-INLA using the combined data set, ensuring
# that `compute = TRUE` is selected in the `control.predictor` argument

Final.Pred <- inla(f04, family = "binomial", Ntrials = 1, data = cc,
                       lincomb = lcb,
                  control.inla = list(lincomb.derived.only = TRUE),
             control.predictor = list(compute = TRUE), 
                 control.fixed = list(mean.intercept = 1.68,
                                      prec.intercept = 1.88^(-2),
                                    mean = list(tree = -0.03,
                                          aggresslow = 1.83),
                                   prec = list(tree  = 0.03^(-2),
                                          aggresslow = 1.6^(-2))))

# Get the marginal distributions:
Pred.marg <- Final.Pred$marginals.lincomb.derived

# Results are on the logit-scale and need to be converted.
# This function converts x into exp(x) / (1 + exp(x))
MyLogit <- function(x) {exp(x)/(1+exp(x))}

# Get mu, selo and seup
MyData$mu <- unlist(lapply(Pred.marg,
             function(x) inla.emarginal(MyLogit,x)))

MyData$selo <- unlist(lapply(Pred.marg,
               function(x) inla.qmarginal(c(0.025), 
                           inla.tmarginal(MyLogit, x))))

MyData$seup <- unlist(lapply(Pred.marg,
               function(x)inla.qmarginal(c(0.975), 
                          inla.tmarginal(MyLogit, x))))

# Define labels
label_agg <- c("high" = "High aggression parents", 
                "low" = "Low aggression parents")

# Plot
ggplot() + 
  geom_jitter(data = cc, aes(y = egg, x = tree),
              shape = 19, size = 2.5, height = 0.01, 
              width = 0.01, alpha = 0.7) +
  xlab("Distance to nearest tree (m)") + 
  ylab("Posterior mean probability of cuckoo parasitism") +
  xlim(-1,101) + 
  theme(text = element_text(size = 13)) + 
  theme(panel.background = element_blank()) + 
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) + 
  theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1)) +
  geom_line(data = MyData, aes(x = tree, y = mu), size = 1) +
  geom_ribbon(data = MyData, aes(x = tree, 
                                 ymax = seup, ymin = selo), alpha = 0.5) +
  facet_grid(. ~ aggress, scales = "fixed", space = "fixed",
                    labeller = labeller (aggress = label_agg)) +
  theme(strip.text = element_text(size = 12, face="italic"))

#===================END===================#
