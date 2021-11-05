#=======================================

# R code for Chapter 3 of "Bayesian GLMs in R for Ecology" 
# by Mark Warren & Carl Smith

#=======================================

#Load packages
library(arm)
library(car)
library(ggplot2)
library(GGally)
library(lattice)
library(lawstat)
library(outliers)
library(tidyverse)
library(rlang)
library(gridExtra)

#=======================================

#Import the data
cyan <- read.table(file = "cyanistes.txt", header = TRUE, 
                   dec = ".", stringsAsFactors = TRUE)

#Use 'str' to inspect the data frame and understand 
# what sort of data are in each column 
str(cyan)

# 'data.frame':	438 obs. of  7 variables:
# $ id    : Factor w/ 377 levels "B1","B10","B111",..: 1 2 3 4 5 6 7 8 9 10 ...
# $ zone  : Factor w/ 9 levels "B","C","CP","E",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ year  : int  2002 2001 2002 2002 2001 2001 2002 2002 2002 2002 ...
# $ multi : Factor w/ 2 levels "no","yes": 2 2 1 1 1 1 1 2 1 2 ...
# $ height: num  2.7 2 2 1.9 2.1 2.5 1.9 2.2 2.4 2.3 ...
# $ day   : int  8 25 10 6 28 23 10 10 11 7 ...
# $ depth : num  0.33 0.25 0.2 0.25 0.33 0.4 0.2 0.2 0.2 0.33 ...

# Are there missing values?
colSums(is.na(cyan))

# id   zone   year  multi  height  day  depth 
# 0      0      0     0      0      0     0 
# No missing data

#=======================================

# OUTLIERS

# Fig 3.1 a boxplot of nest depths for each woodland zone

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

# Then plot
cyan %>%
    ggplot(aes(y = depth, x = zone)) +
    labs(y = "Nest depth", x = "Woodland zone") +
    geom_boxplot(fill = "red1") +
    My_theme +
    theme(panel.border = element_rect(colour = "black", fill=NA, size = 1))

# Fig. 3.2 a multi-panel dotchart to view all 
# continuous variables simultaneously.

# A function for dotplots
multi_dotplot <- function(filename, Xvar, Yvar){
  filename %>%
    ggplot(aes(x = {{Xvar}})) +
    geom_point(aes(y = {{Yvar}})) +
    theme_bw() +
    My_theme +
    coord_flip() +
    labs(x = "Order of Data")}

#Add observations
cyan <- cyan %>%
  mutate(order = seq(1:nrow(cyan)))

#Choose continuous variables to plot
p1 <- multi_dotplot(cyan, order, height)
p2 <- multi_dotplot(cyan, order, day)
p3 <- multi_dotplot(cyan, order, depth)

#Plot as a grid
grid.arrange(p1, p2, p3, nrow = 1)
  
# Use Grubbs' test to assess whether a value that is 
# farthest (above or below) the mean is an outlier

#For nest depth:
grubbs.test(cyan$depth, type = 10) 
      # type 10 is used to detect only one outlier

# Grubbs test for one outlier
# 
# data:  cyan$depth
# G = 4.13755, U = 0.96074, p-value = 0.006479
# alternative hypothesis: highest value 0.75 is an outlier

#For nest box height there are two potential outliers:
grubbs.test(cyan$height, type = 11) 
        # type 11 is used to detect two outliers simultaneously

# Grubbs test for one outlier
# 
# data:  cyan$height
# G = 8.64561, U = 0.91371, p-value < 2.2e-16
# alternative hypothesis: 0.6 and 3.5 are outliers

# 2.2e-16 is 0.000000000000000022,
# In a paper you should write P <0.001

#For day:
grubbs.test(cyan$day, type = 10)

# data:  cyan$day
# G = 2.41610, U = 0.98661, p-value = 1
# alternative hypothesis: highest value 37 is an outlier

#=======================================

#NORMALITY AND HOMOGENEITY OF DEPENDENT VARIABLE

#Fig 3.3 - Frequency polygon plot
cyan %>% ggplot(aes(depth)) +
  geom_freqpoly( bins = 7) +
  labs(x = "Nest depth", y = "Frequency") +
    My_theme +
  theme(panel.border = element_rect(colour = "black", 
                                    fill=NA, size = 1))

#Shapiro-Wilk test for deviation from normality
shapiro.test(cyan$depth)

# data:  cyan$depth
# W = 0.89872, p-value < 2.2e-16

# p<0.05 indicates non-normality, but does not accommodate covariates

# Homogeneity of variance: even distribution of values around the mean

# Homogeneity of response variable depth by covariate height

# mean of response variable depth
mean(cyan$depth)
# 0.33

# mean of covariate height
mean(cyan$height)
# 2.19

# Fig 3.4 - height and depth by year
cyan %>% ggplot(aes(x = height, y = depth)) +
  geom_jitter(shape = 19, size = 2.5, height = 0.05, 
              width = 0.05, alpha = 0.5) +
  geom_hline(yintercept=0.33, linetype="dashed") +
  geom_vline(xintercept=2.19, linetype="dashed") +
  labs(x = "Nest box height (m)", y = "Nest depth") +
  xlim(0,4) + ylim(0,0.8) +
  My_theme +
  theme(panel.border = element_rect(colour = "black", 
                                    fill=NA, size = 1))

# Use modified Levene's test (based on median) - Brown & Forsythe test
leveneTest(cyan$depth,
           cyan$zone,
           location = c("median"), 
           trim.alpha = 0.25)

# Levene's Test for Homogeneity of Variance (center = median: c("median"))
#        Df F value Pr(>F)
# group   8  0.6051 0.7738
#       429

#=======================================

# CALCULATE NUMBER OF ZEROS

# What is the percentage of zeros for nest depth?

sum(cyan$depth == 0) * 100 / nrow(cyan)
#0% zeros

#=======================================

#COLLINEARITY
#
# Fig. 3.5: a summary using the ggpairs command from the GGally library
cyan %>% 
    ggpairs(columns = c("year", "height", "day"), aes(alpha = 0.8), 
            lower = list(combo = wrap("facethist", binwidth = 2))) + 
    My_theme

#Association between day and year

#=======================================
#PLOT RELATIONSHIPS

# Any evidence of an interaction between depth and zone?
# Fig. 3.6
ggplot(cyan, aes(x = height, y = depth)) +
  geom_point(colour = 'gray44') +
  geom_smooth(method = 'lm', colour = 'black', se = FALSE) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, size = 1)) +
  theme(strip.background = element_rect(fill = "white", 
                                        color = "white", size = 1)) +
  theme(text = element_text(size=16)) +
  facet_wrap(~zone)

# Possible interaction here? Slopes appear to vary with zone.

#=======================================
# 
# The data exploration showed:
#   
# 1.	An outlier in the response variable depth. 
#     There were additional outliers in the explanatory variable height.
# 2.	A non-normally distributed but homogenous response variable.
# 3.	No zeros in the response variable.
# 4.	Collinearity between the variables day and year.
# 5.	A potential interaction between depth and zone.
# 6.	Dependency in the data due to repeated sampling of the same nest boxes.
#
##################################### END





