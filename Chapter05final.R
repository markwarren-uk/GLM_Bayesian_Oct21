
#=======================================

# R code for Chapter 5 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Load packages
library(lattice)  
library(ggplot2)
library(GGally)
library(tidyverse)
library(mgcv)
library(plyr)
library(lme4)
library(car)
library(devtools)
library(ggpubr)
library(qqplotr)
library(geiger)
library(INLA)
library(brinla)
library(gridExtra)
library(rlang)
#=======================================

# Import the data
ga <- read_csv(file = "stickleback.csv")

str(ga)

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

# The aim of this study is to understand what ecological variables select
# for lateral plate number in 3-spined sticklebacks on North Uist:

# We will use an Information Theory approach to formulate a set of  
# competing hypothesis and examine which is best supported by the data 

models <- tibble(
model = c("M01", "M02", "M03", "M04", "M05", "M06", "M07", "M08", "M09"),
model_form = c("vert", "ph", "inve", "alt + dist", "alt + area", "comp", "vert + ph", "vert + comp", "ph + sl"),
ref = c("Hoogland et al. (1956)", "Giles (1983)", "Reimchen (1994)", "Raeymaekers et al. (2007)", 
         "Lucek et al. (2016)", "MacColl et al (2013)", "Spence et al. (2013)", "Magalhaes et al. (2016)", "Smith et al. (2020)")
)
models

#=======================================

# 2. PERFORM DATA EXPLORATION

#=======================================

# MISSING VALUES?
colSums(is.na(ga))

# loch  dist   alt  area  ph   comp  inve  vert  sl  plate 
# 0     0     0     0     0     0     0     0    0     0 
#No missing data

#=======================================

# OUTLIERS

# A multi-panel Cleveland dotplot to examine continuous variables
# Fig. 5.1

#Add row numbers
ga <- ga %>%
  mutate(order = seq(1:nrow(ga)))

My_theme <- theme(panel.background = element_blank(),
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
p1 <- multi_dotplot(ga, order, dist) 
p2 <- multi_dotplot(ga, order, alt) 
p3 <- multi_dotplot(ga, order, area) 
p4 <- multi_dotplot(ga, order, ph) 
p5 <- multi_dotplot(ga, order, sl) 
p6 <- multi_dotplot(ga, order, plate) 

#Create grid and plot
grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 2)

#Log transform area
ga$log_area <- log10((ga$area))

# Compare area and log_area
# Fig 5.2
p7 <- multi_dotplot(ga, order, log_area)
grid.arrange(p3, p7, nrow = 1)

# Better

#=======================================

# DISTRIBUTION OF THE DEPENDENT VARIABLE

#Fig 5.3 Frequency polygon plot of plate number
ga %>% 
  ggplot(aes(plate)) +
  geom_freqpoly( bins = 9) +
  labs(x = "Number of lateral plates", y = "Frequency") +
  My_theme +
  theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))

#=======================================

# BALANCE
# Examine the balance of categorical variables

table(ga$comp)
# no yes 
# 33  24 

table(ga$inve)
# no yes 
# 32  25 

table(ga$vert)
# no yes 
# 8  49 

#=======================================

# COLLINEARITY

ga %>% 
    ggpairs(columns = c("dist", "alt", "log_area", "ph", "comp", "inve", "vert", "sl"), 
                            aes(alpha = 0.8), lower = list(continuous = "smooth_loess", 
                                combo = wrap("facethist", binwidth = 5))) + My_theme

# A weak positive correlation between sl and pH

#Calculate Variance Inflation Factor (VIF)
round(vif(glm(plate ~ ph + sl, family = "poisson", data = ga)),2)
# ph   sl 
# 1.95 1.95 

# Some variance inflation but <3

#=======================================

# ZEROS IN THE RESPONSE VARIABLE

round((sum(ga$plate == 0) / nrow(ga))*100,0)
#11% zeros
# Not clear whether this is problematic at this stage...

#=======================================

# RELATIONSHIPS
# Fig. 5.5
My_theme <- theme(panel.background = element_blank(),
                  panel.border = element_rect(fill = NA, size = 1),
                  strip.background = element_rect(fill = "white", 
                                                  color = "white", size = 1),
                  text = element_text(size = 14),
                  panel.grid.major = element_line(colour = "white", size = 0.1),
                  panel.grid.minor = element_line(colour = "white", size = 0.1))

grid.arrange(
    ga %>% 
  ggplot() +
  geom_point(aes(x = dist, y = plate), alpha = 0.5) +
  geom_smooth(aes(x = dist, y = plate), method = 'lm', se=FALSE)+
  labs(x = "Distance (km)", y = "Plate number") +
  My_theme, 

ga %>% 
  ggplot() +
  geom_point(aes(x = alt, y = plate), alpha = 0.5) +
  geom_smooth(aes(x = alt, y = plate), method = 'lm', se=FALSE)+
  labs(x = "Altitude (m)", y = "Plate number") +
  My_theme,   

ga %>% 
  ggplot() +
  geom_point(aes(x = ph, y = plate), alpha = 0.5) +
  geom_smooth(aes(x = ph, y = plate), method = 'lm', se=FALSE)+
  labs(x = "pH", y = "Plate number") +
  My_theme,  

ga %>% 
  ggplot() +
  geom_point(aes(x = sl, y = plate), alpha = 0.5) +
  geom_smooth(aes(x = sl, y = plate), method = 'lm', se=FALSE)+
  labs(x = "SL (mm)", y = "Plate number") +
  My_theme,  
nrow = 2)


# Categorical covariates
# Fig. 5.6
label_inve <- c("no" = "No invertebrate predators", 
                "yes" = "Invertebrate predators")
label_vert <- c("no" = "No vertebrate predators", 
                "yes" = "Vertebrate predators")
ggplot() +
geom_boxplot(data = ga, aes(y = plate, x = comp), fill = "gray88") +
geom_jitter(data = ga, aes(y = plate, x = comp),
                     size = 2, width = 0.05, height = 0.1) +
xlab("Competitors present") + ylab("Plate number")+
theme(text = element_text(size=13)) +
theme(panel.background = element_blank())+
theme(panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))+
  theme(strip.text = element_text(size = 12, face="italic")) +
facet_grid(vert ~ inve, 
                    scales = "fixed", space = "fixed", 
                    labeller=labeller (vert = label_vert, 
                                       inve = label_inve))
         
#=======================================

# 3. SELECT A STATISTICAL MODEL

#=======================================

# Plate number is a count, and includes zeros but
# no negative values
# Poisson distribution is an appropriate starting point

# Use an IT approach and formulate a series of a priori models

#=======================================

# 4. SPECIFY PRIORS

#=======================================

# Fit each model with default priors and with informative
# priors based on pilot data from Smith et al. (2020) J Zool

# ======================================

# 5. FIT MODELS

#=======================================

# Specify the formulae for the a priori models

f01 <- plate ~ vert
f02 <- plate ~ ph 
f03 <- plate ~ inve
f04 <- plate ~ alt + dist 
f05 <- plate ~ alt + log_area 
f06 <- plate ~ comp 
f07 <- plate ~ vert + ph 
f08 <- plate ~ vert + comp
f09 <- plate ~ ph + sl 

# Models M01-10 with default priors
M01 <- inla(f01, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M02 <- inla(f02, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M03 <- inla(f03, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M04 <- inla(f04, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M05 <- inla(f05, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M06 <- inla(f06, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M07 <- inla(f07, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M08 <- inla(f08, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)
M09 <- inla(f09, control.compute = list(dic = TRUE), 
            family = "poisson", data = ga)

# Extract DICs
DefDIC <- c(M01$dic$dic,M02$dic$dic,M03$dic$dic,
            M04$dic$dic,M05$dic$dic,M06$dic$dic,
            M07$dic$dic,M08$dic$dic,M09$dic$dic)

# Add weighting
DefDIC.weights <- aicw(DefDIC)

# Add names
DefDIC.weights %>% 
  mutate(model = c("M01","M02","M03",
                              "M04","M05","M06",
                              "M07","M08","M09")) %>% 
  select(model, everything()) %>% 
  arrange(fit) %>% 
  mutate(across(2:4, round, 2))

# Print DICs
dprint.def <- print (DefDIC.weights,
              abbrev.names = FALSE)

# Order DICs by fit
round(dprint.def[order(dprint.def$fit),],2)

#     fit     delta   w
# M09 190.97  0.00 0.79
# M02 194.29  3.32 0.15
# M07 196.01  5.04 0.06
# M04 219.77 28.81 0.00
# M05 224.46 33.49 0.00
# M08 231.96 40.99 0.00
# M01 236.48 45.51 0.00
# M06 243.70 52.73 0.00
# M03 246.54 55.57 0.00

#=======================================

# Models I01-10 with informative priors
I01 <- inla(f01, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 2.1,
                                 prec.intercept = 0.5^(-2),
                            mean = list(vertyes = 1.4), 
                            prec = list(vertyes = 0.6^(-2))))

I02 <- inla(f02, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = -3.3,
                                 prec.intercept = 0.8^(-2),
                                 mean = list(ph = 0.7), 
                                 prec = list(ph = 0.2^(-2))))

I03 <- inla(f03, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 1.1,
                                 prec.intercept = 0.3^(-2),
                            mean = list(inveyes = 0.4),
                            prec = list(inveyes = 0.3^(-2))))

I04 <- inla(f04, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 2.5,
                                 prec.intercept = 0.3^(-2),
                                mean = list(alt = -0.01,
                                           dist = -0.1),
                                prec = list(alt = 0.1^(-2),
                                           dist = 0.1^(-2))))

I05 <- inla(f05, family = "poisson", data = ga,
            control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 2.9,
                                 prec.intercept = 0.7^(-2),
                                mean = list(alt = -0.01,
                                       log_area = -0.02),
                                prec = list(alt = 0.1^(-2),
                                       log_area = 0.1^(-2))))

I06 <- inla(f06, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 1.6,
                                 prec.intercept = 0.3^(-2),
                            mean = list(compyes = 0.5),
                            prec = list(compyes = 0.3^(-2))))

I07 <- inla(f07, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = -4.8,
                                 prec.intercept = 1.9^(-2),
                            mean = list(vertyes = 0.6,
                                             ph = 0.6),
                            prec = list(vertyes = 0.5^(-2),
                                             ph = 0.2^(-2))))

I08 <- inla(f08, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = 1.7,
                                 prec.intercept = 0.4^(-2),
                            mean = list(vertyes = 0.7,
                                        compyes = 0.9),
                            prec = list(vertyes = 0.3^(-2),
                                        compyes = 0.3^(-2))))

I09 <- inla(f09, family = "poisson", data = ga,
          control.compute = list(dic = TRUE),
            control.fixed = list(mean.intercept = -3.4,
                                 prec.intercept = 0.6^(-2),
                                 mean = list(ph = 0.5,
                                             sl = 0.1),
                                 prec = list(ph = 0.2^(-2),
                                             sl = 0.05^(-2))))

# Extract DICs
InfDIC <- c(I01$dic$dic, I02$dic$dic, I03$dic$dic, 
            I04$dic$dic, I05$dic$dic, I06$dic$dic,
            I07$dic$dic, I08$dic$dic, I09$dic$dic)

# Add weighting
InfDIC.weights <- aicw(InfDIC)

# Add names
rownames(InfDIC.weights) <- c("I01","I02","I03",
                              "I04","I05","I06",
                              "I07","I08","I09")
# Print DICs
dprint.inf <- print (InfDIC.weights,
                     abbrev.names = FALSE)

# Order DICs by fit
round(dprint.inf[order(dprint.inf$fit),],2)

#      fit    delta  w
# I09 189.39  0.00 0.84
# I02 193.34  3.95 0.12
# I07 195.43  6.03 0.04
# I04 219.46 30.07 0.00
# I05 224.38 34.98 0.00
# I08 234.29 44.90 0.00
# I01 236.59 47.19 0.00
# I06 243.45 54.06 0.00
# I03 245.97 56.58 0.00

# Compare best default and informative models
ComDIC <- c(M09$dic$dic, I09$dic$dic)
DIC <- cbind(ComDIC)
rownames(DIC) <- c("default priors","informative priors")
round(DIC,0)

#                      ComDIC
# default priors        191
# informative priors    189

# I09 has a marginally better fit

#=======================================

# 6. OBTAIN THE POSTERIOR DISTRIBUTION

#=======================================

# A. Best-fitting default model (M09)

# Posterior mean values and 95% CI for fixed effects (default)
M09Betas <- M09$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(M09Betas, digits = 3)
#              mean   sd   0.025quant 0.975quant
# (Intercept) -4.011 0.657     -5.305     -2.724
# ph           0.488 0.126      0.242      0.735
# sl           0.052 0.022      0.007      0.095


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
# Fig. 5.7

# Model intercept (Beta1)
PosteriorBeta1.M09 <- as.data.frame(M09$marginals.fixed$`(Intercept)`)
PriorBeta1.M09     <- data.frame(x = PosteriorBeta1.M09[,"x"], 
                           y = dnorm(PosteriorBeta1.M09[,"x"],0,0))
Beta1mean.M09 <- M09Betas["(Intercept)", "mean"]
Beta1lo.M09   <- M09Betas["(Intercept)", "0.025quant"]
Beta1up.M09   <- M09Betas["(Intercept)", "0.975quant"]

beta1 <- ggplot() +
annotate("rect", xmin = Beta1lo.M09, xmax = Beta1up.M09,
                  ymin = 0, ymax = 0.62, fill = "gray88") +
geom_line(data = PosteriorBeta1.M09,
                  aes(y = y, x = x), lwd = 1.2) +
geom_line(data = PriorBeta1.M09,
                   aes(y = y, x = x), color = "gray55", lwd = 1.2) +
xlab("Intercept") +
ylab("Density") +
xlim(-8,1) + ylim(0,0.62) +
geom_vline(xintercept = 0, linetype = "dotted") +
geom_vline(xintercept = Beta1mean.M09, linetype = "dashed") +
theme(text = element_text(size=13)) +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                     colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# ph (Beta2)
PosteriorBeta2.M09 <- as.data.frame(M09$marginals.fixed$`ph`)
PriorBeta2.M09 <- data.frame(x = PosteriorBeta2.M09[,"x"], 
                       y = dnorm(PosteriorBeta2.M09[,"x"],0,31.6))
Beta2mean.M09 <- M09Betas["ph", "mean"]
Beta2lo.M09   <- M09Betas["ph", "0.025quant"]
Beta2up.M09   <- M09Betas["ph", "0.975quant"]

beta2 <- ggplot() +
annotate("rect", xmin = Beta2lo.M09, xmax = Beta2up.M09,
                  ymin = 0, ymax = 3.3, fill = "gray88") +
geom_line(data = PosteriorBeta2.M09,
                   aes(y = y, x = x), lwd = 1.2) +
geom_line(data = PriorBeta2.M09,
                   aes(y = y, x = x), color = "gray55", lwd = 1.2) +
xlab("Slope for pH") +
ylab("Density") +
xlim(-0.5,1.5) + ylim(0,3.3) +
geom_vline(xintercept = 0, linetype = "dotted") +
geom_vline(xintercept = Beta2mean.M09, linetype = "dashed") +
theme(text = element_text(size=13)) +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# SL (Beta3)
PosteriorBeta3.M09 <- as.data.frame(M09$marginals.fixed$`sl`)
PriorBeta3.M09     <- data.frame(x = PosteriorBeta3.M09[,"x"],
                           y = dnorm(PosteriorBeta3.M09[,"x"],0,31.6))

Beta3mean.M09 <- M09Betas["sl", "mean"]
Beta3lo.M09   <- M09Betas["sl", "0.025quant"]
Beta3up.M09   <- M09Betas["sl", "0.975quant"]

beta3 <- ggplot() +
annotate("rect", xmin = Beta3lo.M09, xmax = Beta3up.M09,
                  ymin = 0, ymax = 20, fill = "gray88") +
geom_line(data = PosteriorBeta3.M09,
                   aes(y = y, x = x), lwd = 1.2) +
geom_line(data = PriorBeta3.M09,
                   aes(y = y, x = x), color = "gray55", lwd = 1.2) +
xlab("Slope for SL") +
ylab("Density") +
xlim(-0.05,0.15) + ylim(0,20) +
geom_vline(xintercept = 0, linetype = "dotted") +
geom_vline(xintercept = Beta3mean.M09, linetype = "dashed") +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))

# Combine plots
ggarrange(beta1, beta2, beta3,
                        labels = c("A", "B", "C"),
                        ncol = 2, nrow = 2)


# These density functions look normally distributed because INLA assumes
# that the posterior distributions of the betas are Gaussian.

#=======================================
# B. Informative model (I09)

# Posterior mean values and 95% CI for fixed effects (informative)
I09Betas <- I09$summary.fixed[,c("mean", "sd", 
                                 "0.025quant", 
                                 "0.975quant")] 
round(I09Betas, digits = 2)
#              mean   sd  0.025quant  0.975quant
# (Intercept) -3.74 0.43      -4.59      -2.89
# ph           0.44 0.10       0.25       0.62
# sl           0.05 0.02       0.02       0.09

# The first column shows the posterior mean, the second the
# standard deviation and the third and fourth columns show the 
# 2.5% and 97.5% quantiles of the posterior distribution.

# Plot posterior distributions for fixed parameters 
# Fig. 5.8

# Model intercept (Beta1)
PosteriorBeta1.I09 <- as.data.frame(I09$marginals.fixed$`(Intercept)`)
PriorBeta1.I09     <- data.frame(x = PosteriorBeta1.I09[,"x"], 
                                 y = dnorm(PosteriorBeta1.I09[,"x"],-3.4,0.6))
Beta1mean.I09 <- I09Betas["(Intercept)", "mean"]
Beta1lo.I09   <- I09Betas["(Intercept)", "0.025quant"]
Beta1up.I09   <- I09Betas["(Intercept)", "0.975quant"]

Ibeta1 <- ggplot() +
annotate("rect", xmin = Beta1lo.I09, xmax = Beta1up.I09,
                  ymin = 0, ymax = 1, fill = "gray88") +
geom_line(data = PosteriorBeta1.I09,
                   aes(y = y, x = x), lwd = 1.2) +
geom_line(data = PriorBeta1.I09,
                   aes(y = y, x = x), color = "gray55", lwd = 1.2) +
xlab("Intercept") +
ylab("Density") +
xlim(-8,1) + ylim(0,1) +
geom_vline(xintercept = 0, linetype = "dotted") +
geom_vline(xintercept = Beta1mean.I09, linetype = "dashed") +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


# ph (Beta2)
PosteriorBeta2.I09 <- as.data.frame(I09$marginals.fixed$`ph`)
PriorBeta2.I09 <- data.frame(x = PosteriorBeta2.I09[,"x"], 
                             y = dnorm(PosteriorBeta2.I09[,"x"],0.5,0.2))
Beta2mean.I09 <- I09Betas["ph", "mean"]
Beta2lo.I09   <- I09Betas["ph", "0.025quant"]
Beta2up.I09   <- I09Betas["ph", "0.975quant"]

Ibeta2 <- ggplot() +
annotate("rect", xmin = Beta2lo.I09, xmax = Beta2up.I09,
                  ymin = 0, ymax = 4.5, fill = "gray88") +
geom_line(data = PosteriorBeta2.I09,
                   aes(y = y, x = x), lwd = 1.2) +
geom_line(data = PriorBeta2.I09,
                   aes(y = y, x = x), color = "gray55", lwd = 1.2) +
xlab("Slope for pH") +
ylab("Density") +
xlim(-0.5,1.5) + ylim(0,4.5) +
geom_vline(xintercept = 0, linetype = "dotted") +
geom_vline(xintercept = Beta2mean.I09, linetype = "dashed") +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))



# SL (Beta3)
PosteriorBeta3.I09 <- as.data.frame(I09$marginals.fixed$`sl`)
PriorBeta3.I09     <- data.frame(x = PosteriorBeta3.I09[,"x"],
                                 y = dnorm(PosteriorBeta3.I09[,"x"],0.1,0.05))

Beta3mean.I09 <- I09Betas["sl", "mean"]
Beta3lo.I09   <- I09Betas["sl", "0.025quant"]
Beta3up.I09   <- I09Betas["sl", "0.975quant"]

Ibeta3 <- ggplot() +
annotate("rect", xmin = Beta3lo.I09, xmax = Beta3up.I09,
                  ymin = 0, ymax = 25, fill = "gray88") +
geom_line(data = PosteriorBeta3.I09,
                   aes(y = y, x = x), lwd = 1.2) +
geom_line(data = PriorBeta3.I09,
                   aes(y = y, x = x), color = "gray55", lwd = 1.2) +
xlab("Slope for SL") +
ylab("Density") +
xlim(-0.05,0.15) + ylim(0,25) +
geom_vline(xintercept = 0, linetype = "dotted") +
geom_vline(xintercept = Beta3mean.I09, linetype = "dashed") +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1))+
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1)) 

# Combine plots
ggarrange(Ibeta1, Ibeta2, Ibeta3,
                        labels = c("A", "B", "C"),
                        ncol = 2, nrow = 2)

#=======================================
# Compare with frequentist Poisson GLM

Freq <- glm(plate ~ ph + sl,
            family = "poisson",
            data = ga)

round(summary(Freq)$coef[,1:4],3)

#               Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept)   -4.013      0.657  -6.106    0.000
# ph             0.488      0.126   3.885    0.000
# sl             0.052      0.022   2.316    0.021

# Posterior mean values and 95% CI for fixed effects (default priors)
round(M09Betas, digits = 3)

#              mean    sd   0.025quant 0.975quant
# (Intercept) -4.011 0.657     -5.305     -2.724
# ph           0.488 0.126      0.242      0.735
# sl           0.052 0.022      0.007      0.095

# Posterior mean values and 95% CI for fixed effects (informative priors)
round(I09Betas, digits = 3)

#              mean    sd   0.025quant 0.975quant
# (Intercept) -3.737 0.432     -4.585     -2.891
# ph           0.436 0.095      0.249      0.623
# sl           0.054 0.019      0.018      0.091


pivot_wider(round(M09Betas, digits = 2) %>% 
              mutate(term = c("(Intercept)", "ph", "sl")) %>% 
              select(term, mean, sd) %>% 
              unite(param, mean, sd, sep = ", sd="),
            names_from = term,
            values_from = param)

kable(
  bind_rows(
    pivot_wider(broom::tidy(Freq)%>% 
                  mutate_if(is.numeric, round, 2) %>% 
                  select(term, estimate, std.error) %>% 
                  unite(param, estimate, std.error, sep = ", sd="),
                names_from = term, 
                values_from = c(param)),
    
    pivot_wider(round(M09Betas, digits = 2) %>% 
                  mutate(term = c("(Intercept)", "ph", "sl")) %>% 
                  select(term, mean, sd) %>% 
                  unite(param, mean, sd, sep = ", sd="),
                names_from = term,
                values_from = param),
    
    pivot_wider(round(I09Betas, digits = 2) %>% 
                  mutate(term = c("(Intercept)", "ph", "sl")) %>% 
                  select(term, mean, sd) %>% 
                  unite(param, mean, sd, sep = ", sd="),
                names_from = term,
                values_from = param)
  )%>% 
    mutate(Model = c("Frequentist", "Bayesian (default)", "Bayesian (informative)")) %>% 
    select(Model, everything()),
  align = "lccrr",
  format = "html", col.names = c("Model", "Intercept", "pH", "sl"),
  caption = 'Comparison of model parameters for frequentist, 
  Bayesian model with non-informative and informative priors of Poisson 
  GLM model to investigate the evolution of lateral plate number in three-spined sticklebacks on North Uist.'
) %>% 
  kable_styling(bootstrap_options = "striped", full_width = F)

#=======================================

# 7. MODEL CHECKS

#=======================================

# A. Model selection

# We have used an IT approach to conduct model selection, setting up
# alternative plausible models and testing among them using the DIC

# The rank order of models with default priors is:
round(dprint.def[order(dprint.def$fit),],2)

#     fit     delta  w
# M09 190.97  0.00 0.79
# M02 194.29  3.32 0.15
# M07 196.01  5.04 0.06
# M04 219.77 28.81 0.00
# M05 224.46 33.49 0.00
# M08 231.96 40.99 0.00
# M01 236.48 45.51 0.00
# M06 243.70 52.73 0.00
# M03 246.54 55.57 0.00

# And with informative priors:
round(dprint.inf[order(dprint.inf$fit),],2)

#     fit    delta   w
# I09 189.39  0.00 0.84
# I02 193.34  3.95 0.12
# I07 195.40  6.00 0.04
# I04 219.46 30.07 0.00
# I05 224.38 34.98 0.00
# I08 231.71 42.32 0.00
# I01 236.40 47.00 0.00
# I06 244.07 54.68 0.00
# I03 246.33 56.94 0.00

# In each case model 09 is the best-fitting, with the model 
# with informative priors marginally better than the model
# with default priors

# Compare best-fitting models with non-informative and informative priors
dic2 <- c(M09$dic$dic, I09$dic$dic)
DIC2 <- cbind(dic2)
rownames(DIC2) <- c("default priors","informative priors")
round(DIC2,0)

# default priors      191
# informative priors  189

# These are essentially the same

#=======================================

# B. Dispersion
# (Just for model with informative priors)

# Poisson GLMs assume that the mean and variance of the response variable 
# increase at the same rate. This assumption must be confirmed. 
# Overdispersion means that a Poisson distribution does not adequately 
# model the variance and is not appropriate for the analysis.

# We will implement the 7-step protocol of Zuur et al. (2017)
# to determine whether the Poisson GLM is over- or under-dispersed.

# 1. Apply Poisson GLM in INLA
# 2. Simulate a set of regression parameters from the posterior distribution
# 3. Calculate expected values
# 4. Simulate count data using the rpois function
# 5. Calculate Pearson residuals for simulated data
# 6. Repeat steps 2-5 x 1000 times
# 7. Compare sums of squared Pearson residuals with observed data


# 1. Refit the model with the 'config = TRUE' option,
# which permits us to simulate regression parameters.

I09 <- inla(f09, family = "poisson", data = ga,
            control.compute = list(config = TRUE),
            control.predictor = list(compute = TRUE), 
            control.fixed = list(mean.intercept = -3.4,
                                 prec.intercept = 0.6^(-2),
                                 mean = list(ph = 0.5,
                                             sl = 0.1),
                                 prec = list(ph = 0.2^(-2),
                                             sl = 0.05^(-2))))

# 2. Simulate regression parameters using the inla.posterior.sample 
# function to simulate from the model. The output is 
# stored in the Sim object.

set.seed(1966)
SimData <- inla.posterior.sample(n = 1, result = I09)

# This gives an object with 60 rows for this 
# data set and model. The first 57 
# rows are simulated values for eta = X * beta, 
# where X is the matrix with covariates, and 
# the last 3 rows are simulated regression parameters. 
# This is just one set of simulated values 
# (i.e. n = 1 in the function above). 

SimData[[1]]$latent

# 3: Calculate predicted values

BetaRows <- 58:60  #Last 3 rows are the betas
BetaRows <- (nrow(SimData[[1]]$latent) - 2):nrow(SimData[[1]]$latent)
Betas    <- SimData[[1]]$latent[BetaRows]
Betas  #These are the simulated betas (from the joint distribution)

# Now get the corresponding X matrix so that INLA can calculate X * Betas.
X <- model.matrix(~ ph + sl, data = ga)

# And now we can calculate: exp(X * Betas)
mu <- exp(X %*% Betas)

# 4. Simulate count data
Ysim <- rpois(n = nrow(ga), lambda = mu)

# We now have 57 simulated values that 
# correspond to the observed covariate values. 

# 5: Calculate summary statistic

mu1 <- I09$summary.fitted.values[,"mean"]  
Es <- (Ysim - mu1) /sqrt(mu1)
N    <- nrow(ga) 
Npar <- length(BetaRows)  
round(sum(Es^2) / (N - Npar),1)
# This should be close to 1.

# Step 6: Repeat steps 2 – 5 x 1000 times

SimData <- inla.posterior.sample(n = 1000, result = I09)

# We have 1000 simulated mean values from the model. 
# Extract the simulated betas and predict plate number.

N      <- nrow(ga)
Ysim   <- matrix(nrow = N, ncol = 1000)
mu.sim <- matrix(nrow = N, ncol = 1000)
for (i in 1:1000){
  Betas      <- SimData[[i]]$latent[BetaRows]
  mu.sim[,i] <- exp(X %*% Betas)
  Ysim[,i]   <- rpois(n = N, lambda = exp(X %*% Betas))
}

# 7: Compare simulation results and observed data

E1  <- (ga$plate - mu1) /sqrt(mu1)
Dispersion <- sum(E1^2) / (N - Npar)

# We calculate the same statistic for 
# each of the simulated data sets.
Dispersion.sim <- NULL
for(i in 1:1000){
  e2 <- (Ysim[,i] - mu.sim[,i]) /sqrt(mu.sim[,i])
  Dispersion.sim[i] <- sum(e2^2) / (N - Npar)
}

# Plot the simulated results
# Fig. 5.9
par(mfrow = c(1,1), mar = c(5,5,2,2), cex.lab = 1.3)
hist(Dispersion.sim,
     xlab = "Model dispersion",
     ylab = "Frequency",
     main = "",
     xlim = c(0, 2.5),
     breaks = 16)

# And identify the dispersion for the original Poisson GLM
points(x = Dispersion, 
       y = 0, 
       pch = 19, 
       cex = 2.5, 
       col = 1)

# The dispersion statistic for the 
# observed data is at the lower end of 
# the simulated data sets, indicating that 
# the variation in the observed data is lower 
# than is expected for a Poisson distribution, 
# which indicates mild under-dispersion.

#=======================================

# C. Posterior predictive check

I09.pred <- inla(f09, family = "poisson", data = ga,
              control.predictor = list(link = 1,
                                    compute = TRUE),
                 control.compute = list(dic = TRUE, 
                                        cpo = TRUE),
            control.fixed = list(mean.intercept = -3.4,
                                 prec.intercept = 0.6^(-2),
                                 mean = list(ph = 0.5,
                                             sl = 0.1),
                                 prec = list(ph = 0.2^(-2),
                                             sl = 0.05^(-2))))

ppp <- vector(mode = "numeric", length = nrow(ga))
for(i in (1:nrow(ga))) {
  ppp[i] <- inla.pmarginal(q = ga$plate[i],
                    marginal = I09.pred$marginals.fitted.values[[i]])
}

# Fig. 5.10 
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

# Problem.

#=======================================

# D. Cross-Validation Model Checking

# Use CPO and PIT
sum(I09.pred$cpo$failure)
# [1] 0
# 0 indicates CPO is reliable

#Extract pit values
PIT <- (I09.pred$cpo$pit)

#And plot
Pit1 <- ggplot() +
geom_histogram(aes(PIT), binwidth = 0.08, 
                        colour = "black", fill = "gray88") +
xlab("PIT") +
ylab("Frequency") +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))


Pit2 <- ggplot(mapping = aes(sample = I09.pred$cpo$pit)) +
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
# Fig. 5.11
ggarrange(Pit1, Pit2,
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

#PIT values are uniform

#=======================================

# E. Bayesian residual analysis

# Plot residuals versus fitted values & each covariate in the model

# Obtain fitted values
# (Ensure  compute = TRUE in `control.predictor`)

Fit <- I09.pred$summary.fitted.values[, "mean"]

# Calculate residuals
Res <- ga$plate - Fit
ResPlot <- cbind.data.frame(Fit,Res,ga$ph,ga$sl)

# Plot residuals against fitted
FigA <- ggplot(ResPlot, aes(x=Fit, y=Res)) + 
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
FigB <- ggplot(ResPlot, aes(x=ga$ph, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("pH") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

FigC <- ggplot(ResPlot, aes(x=ga$sl, y=Res)) + 
  geom_point(shape = 19, size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylab("") + xlab("SL (mm)") +
  theme(text = element_text(size=13)) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, 
                                    colour = "black", size = 1)) +
  theme(strip.background = element_rect
        (fill = "white", color = "white", size = 1))

# Combine plots
# Fig 5.12
ggarrange(FigA, FigB, FigC,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)

# Residuals analysis indicates no violation of model assumptions

# =======================================
  
# F. Prior sensitivity analysis

# Model with priors unchanged
I9.inform <- inla(f09, family = "poisson", data = ga,
                  control.predictor = list(link = 1,
                                        compute = TRUE),
                     control.compute = list(dic = TRUE, 
                                            cpo = TRUE),
                 control.fixed = list(mean.intercept = -3.4,
                                      prec.intercept = 0.6^(-2),
                                      mean = list(ph = 0.5,
                                                  sl = 0.1),
                                      prec = list(ph = 0.2^(-2),
                                                  sl = 0.05^(-2))))

Betas.Inform <- I9.inform$summary.fixed[,c("mean", "sd", 
                                           "0.025quant", 
                                           "0.975quant")] 
round(Betas.Inform, digits = 2)

#              mean  sd   0.025quant 0.975quant
# (Intercept) -3.74 0.43      -4.59      -2.89
# ph           0.44 0.10       0.25       0.62
# sl           0.05 0.02       0.02       0.09

#== 1. Increase priors by 20% ==
#   Intercept from -3.40 to -2.72 (var 0.36 to 0.43)
#   ph from 0.50 to 0.60 (var 0.040 to 0.0480)
#   sl from 0.10 to 0.12 (var 0.0025 to 0.003)

#Re-run model and obtain estimates of betas
I9.plus20 <- inla(f09, family = "poisson", data = ga,
                  control.predictor = list(link = 1,
                                        compute = TRUE),
                  control.compute = list(dic = TRUE, 
                                         cpo = TRUE),
                  control.fixed = list(mean.intercept = -2.72,
                                       prec.intercept = 0.66^(-2),
                                       mean = list(ph = 0.6,
                                                   sl = 0.12),
                                       prec = list(ph = 0.22^(-2),
                                                   sl = 0.055^(-2))))

Betas.plus20 <- I9.plus20$summary.fixed[,c("mean", "sd", 
                                           "0.025quant", 
                                           "0.975quant")] 
round(Betas.plus20, digits = 2)

#             mean    sd   0.025quant 0.975quant
# (Intercept) -3.49 0.45      -4.38      -2.60
# ph           0.42 0.10       0.23       0.62
# sl           0.05 0.02       0.01       0.09

#== 2. Decrease priors by 20% ==
#   Intercept from -3.40 to -4.08 (var 0.36 to 0.29)
#   ph from 0.50 to 0.40 (var 0.040 to 0.0320)
#   sl from 0.10 to 0.08 (var 0.0025 to 0.002)

#Re-run model and obtain estimates of betas
I9.minus20 <- inla(f09, family = "poisson", data = ga,
                   control.predictor = list(link = 1,
                                         compute = TRUE),
                  control.compute = list(dic = TRUE, 
                                         cpo = TRUE),
                  control.fixed = list(mean.intercept = -4.08,
                                       prec.intercept = 0.54^(-2),
                                       mean = list(ph = 0.4,
                                                   sl = 0.08),
                                       prec = list(ph = 0.179^(-2),
                                                   sl = 0.045^(-2))))

Betas.minus20 <- I9.minus20$summary.fixed[,c("mean", "sd", 
                                             "0.025quant", 
                                             "0.975quant")] 
round(Betas.minus20, digits = 2)

#             mean   sd   0.025quant 0.975quant
# (Intercept) -4.04 0.41      -4.83      -3.24
# ph           0.45 0.09       0.27       0.63
# sl           0.06 0.02       0.02       0.09

#=========================
# Plot posterior distributions of the alternative models

# Model intercept (Beta1)
PostBeta1.I9.inform  <- as.data.frame(I9.inform$marginals.fixed$`(Intercept)`)
PostBeta1.I9.plus20  <- as.data.frame(I9.plus20$marginals.fixed$`(Intercept)`)
PostBeta1.I9.minus20 <- as.data.frame(I9.minus20$marginals.fixed$`(Intercept)`)

beta1.sens <- ggplot() +
  geom_line(data = PostBeta1.I9.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid") +
geom_line(data = PostBeta1.I9.plus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44") +
geom_line(data = PostBeta1.I9.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44") +
xlab("Intercept") +
ylab("Density") +
xlim(-6,-1) +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))



# Male ph (Beta2)
PostBeta1.I9.inform  <- as.data.frame(I9.inform$marginals.fixed$`ph`)
PostBeta1.I9.plus20  <- as.data.frame(I9.plus20$marginals.fixed$`ph`)
PostBeta1.I9.minus20 <- as.data.frame(I9.minus20$marginals.fixed$`ph`)

beta2.sens <- ggplot() +
geom_line(data = PostBeta1.I9.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid") +
geom_line(data = PostBeta1.I9.plus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44") +
geom_line(data = PostBeta1.I9.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44") +
xlab("Slope for pH") +
ylab("Density") +
xlim(0,0.9) +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, 
                                           colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))

# Supplementary feeding (Beta3)
PostBeta1.I9.inform  <- as.data.frame(I9.inform$marginals.fixed$`sl`)
PostBeta1.I9.plus20  <- as.data.frame(I9.plus20$marginals.fixed$`sl`)
PostBeta1.I9.minus20 <- as.data.frame(I9.minus20$marginals.fixed$`sl`)

beta3.sens <- ggplot() +
geom_line(data = PostBeta1.I9.inform,
                   aes(y = y, x = x), lwd = 0.8, linetype = "solid") +
geom_line(data = PostBeta1.I9.plus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dashed", colour = "gray44") +
geom_line(data = PostBeta1.I9.minus20,
                   aes(y = y, x = x), lwd = 0.8, 
                   linetype = "dotted", colour = "gray44") +
xlab("Slope for SL") +
ylab("Density") +
xlim(-0.05,0.15) +
theme(text = element_text(size=13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
theme(strip.background = element_rect
               (fill = "white", color = "white", size = 1))

# Combine plots
# Fig. 5.13
 ggarrange(beta1.sens, beta2.sens, beta3.sens,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)

# Changes to posterior distributions are negligible

#=======================================

# 8. INTERPRET AND PRESENT MODEL OUTPUT

#=======================================

# The final model is: 
Final <- inla(f09, family = "poisson", data = ga,
                     control.predictor = list(link = 1,
                                           compute = TRUE),
                        control.compute = list(dic = TRUE, 
                                               cpo = TRUE),
                        control.fixed = list(mean.intercept = -3.4,
                                             prec.intercept = 0.6^(-2),
                                             mean = list(ph = 0.5,
                                                         sl = 0.1),
                                             prec = list(ph = 0.2^(-2),
                                                         sl = 0.05^(-2))))

# Posterior mean values and 95% CI for fixed effects
BetasFinal <- Final$summary.fixed[,c("mean", "sd", 
                                     "0.025quant", 
                                     "0.975quant")] 
round(BetasFinal, digits = 3)

#              mean  sd   0.025quant 0.975quant
# (Intercept) -3.74 0.43      -4.59      -2.89
# ph           0.44 0.10       0.25       0.62
# sl           0.05 0.02       0.02       0.09


## Model interpretation

# Positive effect of ph on plate number
# Positive effect of sl on plate number


#=======================================

# 9. VISUALISE THE RESULTS

#=======================================

#Plot figure for ph
MyData <- expand.grid(
  ph = seq(5.2, 8.7, length = 50),
  sl = 40)

# 2. Make a design matrix
Xmat <- model.matrix(~ ph + sl, data = MyData)

Xmat <- as.data.frame(Xmat)

lcb <- inla.make.lincombs(Xmat)

# Re-run the model in R-INLA using the combined data set, ensuring
# that `compute = TRUE` is selected in the `control.predictor` argument
Final.Pred <- inla(f09, family = "poisson", data = ga,
                       lincomb = lcb,
                  control.inla = list(lincomb.derived.only = TRUE),
             control.predictor = list(compute = TRUE),
                 control.fixed = list(mean.intercept = -3.4,
                                      prec.intercept = 0.6^(-2),
                                      mean = list(ph = 0.5,
                                                  sl = 0.1),
                                      prec = list(ph = 0.2^(-2),
                                                  sl = 0.05^(-2))))


# Run loop to get mu, selo and seup
Pred.marg <- Final.Pred$marginals.lincomb.derived
for (i in 1:nrow(MyData)){
  MyData$mu[i]  <- inla.emarginal(exp, Pred.marg[[i]])
  lo.up <- inla.qmarginal(c(0.025, 0.975), 
                          inla.tmarginal(exp, Pred.marg[[i]]))
  MyData$selo[i] <- lo.up[1]
  MyData$seup[i] <- lo.up[2]    	
}               


# Plot
phfig <- ggplot()  +
geom_jitter(data = ga, 
                     aes(y = plate, x = ph),
                     shape = 19, size = 2.2,
                     height = 0.25, width = 0.25, alpha = 0.6) +
xlab("Loch pH")  +
ylab("Posterior mean plate number")  +
ylim(0,13) +# + xlim(25,45) +
theme(text = element_text(size = 13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
theme(strip.background = element_rect (fill = "white", color = "white", size = 1)) +
geom_line(data = MyData, aes(x = ph, y = mu), size = 1) +
geom_ribbon(data = MyData, aes(x = ph, ymax = seup, ymin = selo), alpha = 0.5)


#=======================================

#Plot figure for sl
MyData <- expand.grid(
  ph = 7,
  sl = seq(25, 43, length = 50))

# 2. Make a design matrix
Xmat <- model.matrix(~ ph + sl, data = MyData)

Xmat <- as.data.frame(Xmat)

lcb <- inla.make.lincombs(Xmat)

# Re-run the model in R-INLA using the combined data set, ensuring
# that `compute = TRUE` is selected in the `control.predictor` argument
Final.Pred <- inla(f09, family = "poisson", data = ga,
                   lincomb = lcb,
                   control.inla = list(lincomb.derived.only = TRUE),
                   control.predictor = list(compute = TRUE),
                   control.fixed = list(mean.intercept = -3.4,
                                        prec.intercept = 0.6^(-2),
                                        mean = list(ph = 0.5,
                                                    sl = 0.1),
                                        prec = list(ph = 0.2^(-2),
                                                    sl = 0.05^(-2))))


# Run loop to get mu, selo and seup
for (i in 1:nrow(MyData)){
  MyData$mu[i]  <- inla.emarginal(exp, Pred.marg[[i]])
  lo.up <- inla.qmarginal(c(0.025, 0.975), 
                          inla.tmarginal(exp, Pred.marg[[i]]))
  MyData$selo[i] <- lo.up[1]
  MyData$seup[i] <- lo.up[2]    	
}               


# Plot
slfig <- ggplot() +
geom_jitter(data = ga, 
                     aes(y = plate, x = sl),
                     shape = 19, size = 2.2,
                     height = 0.25, width = 0.25, alpha = 0.6) +
xlab("Mean population SL (mm)")  +
ylab("")  +
ylim(0,13) + # + xlim(25,45)
theme(text = element_text(size = 13))  +
theme(panel.background = element_blank()) +
theme(panel.border = element_rect(fill = NA, colour = "black", size = 1)) +
theme(strip.background = element_rect (fill = "white", color = "white", size = 1)) +
geom_line(data = MyData, aes(x = sl, y = mu), size = 1) +
geom_ribbon(data = MyData, aes(x = sl, ymax = seup, ymin = selo), alpha = 0.5)


# Combine plots
# Fig. 5.14

ggarrange(phfig, slfig, 
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)


#===================END===================#
