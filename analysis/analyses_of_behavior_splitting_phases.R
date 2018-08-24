################################################################## Code info ###############################################################################################################################################################################################################


# Experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Gilles Pourtois, Ernst Koster, & Søren Andersen) (*: co-first authors)
# Code written by: Ivan Grahek & Antonio Schettino (2016-2018)
# Description: Code for the analysis of behavioral data for Experiment 1 of the SSVEP - reward project. 

################################################################## Importing data & first steps ###############################################################################################################################################################################################################


# Clear environemnt and import data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
# clear the console
cat("\014") 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2,brms,rstan, tidyverse, here)
# set seed
set.seed(42) 
# import data
data.raw = read.csv(file = here("data","Data_behavior_exp1_48pps.csv"),header=TRUE,na.strings="NaN")

# Prepare the dataset------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Adding and renaming variables 
# rename EventType variable
names(data.raw)[names(data.raw) == "EventType"] = "MovedDots"
# add a variable with the name of the attended color instead of a numbers
data.raw$AttendedColor = ifelse(data.raw$AttendedColor==1,"red","blue")
# add a variable saying which color was linked with High_Rew (even numbers - blue was High_Rew)
data.raw$RewardedColor = ifelse(data.raw$ParticipantNo%%2==0,"blue","red") 
# add a variable with the name of the moved color instead of a numbers
data.raw$MovedDots = ifelse(data.raw$MovedDots==1,"red","blue") 
# split experimental phases into 6 isntead of 3 phases (trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext)
data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,100,200,300,400,500,600),labels=c("Bsln1","Bsln2","Acq1","Acq2","Ext1","Ext2")) 
# split experimental phases into 3 phases (trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext)
#data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,200,400,600),labels=c("Bsln","Acq","Ext")) # trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext

### Convert variables to be used in analyses into factors
data.raw[c("ParticipantNo", "AttendedColor","RewardedColor", "MovedDots", "ExpPhase" )] = 
  lapply(data.raw[c("ParticipantNo", "AttendedColor","RewardedColor", "MovedDots", "ExpPhase" )], factor)

### Create variables needed for the accuracy analyses
# count hits, false alarms, misses, correct rejections, and RT separately for each participant (their calculation is done in Matlab: see DataProcessing.m)
data.final = ddply(data.raw,.(ParticipantNo,ExpPhase,AttendedColor,RewardedColor,MovedDots),summarize,
                  numtrials=length(which(Response!=99)), # number of trials per condition (anything that is not 99 or any other number that we're not using)
                  Hits=length(which(Response==1)), # hits: attended color moved, correct response
                  FAs=length(which(Response==2)), # false alarms: attended color did not move, (wrong) response
                  Misses=length(which(Response==0)), # misses: attended color moved, no response
                  CRs=length(which(Response==3)), # correct rejections: attended color did not move, no response
                  mean.RT=mean(RT,na.rm=TRUE)) # mean RT per condition


################################################################## Calculate accuracy and RTs per condition ###############################################################################################################################################################################################################

# Prepare the data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Calculate Hits and False alarms
# Hits are calculated for each participant in each condition on trials when they are attending the color that moved. 
# False alarms are  calculated for each participant in each condition on trials when they are attending the color that didn't move (the unattended color moved, but they responded)  
# Here we create the same number of hits & fas for each of the two conditions (moved attended or not)
data.final = ddply(data.final, .(ParticipantNo,ExpPhase,AttendedColor), transform, 
                 Hits = Hits[MovedDots==AttendedColor],
                 FAs = FAs[MovedDots!=AttendedColor])

# Keep only trials on which the attended color moved (we can do behavioral analysis only on those)
data.final = subset(data.final,MovedDots==AttendedColor)

### Calculate d'
# use loglinear transformation: add 0.5 to Hits, FAs, Misses, and CRs (Hautus, 1995, Behavior Research Methods, Instruments, & Computers),
# which is preferred over the 1/2N rule (Macmillan & Kaplan, 1985, Psychological Bulletin) because it results in less biased estimates of d'.
data.final =  ddply(data.final,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,numtrials),summarize,
                      tot.Hits=Hits+.5, # hits
                      tot.FAs=FAs+.5, # false alarms
                      tot.Misses=(numtrials-tot.Hits)+.5, # misses
                      tot.CRs=(numtrials-tot.FAs)+.5, # correct rejections
                      Hit.Rate=tot.Hits/(tot.Hits+tot.Misses), # hit rate
                      FA.Rate=tot.FAs/(tot.FAs+tot.CRs), # false alarm rate
                      dprime=qnorm(Hit.Rate)-qnorm(FA.Rate),
                      Hits.RTs=mean(mean.RT,na.rm=TRUE)) # mean RTs
                      

# Create the final dataframe for accuracy and RTs ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or Low_Rewed color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"High_Rew","Low_Rew")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)

################################################################## Plotting ###############################################################################################################################################################################################################

# # Plot Hit rates------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Prepare the dataset
data.plot = data.final

# rename variables
colnames(data.plot)[colnames(data.plot)=="ExpPhase"] <- "Reward phase"
colnames(data.plot)[colnames(data.plot)=="Condition"] <- "Reward probability"

# rename conditions
data.plot$`Reward phase` = recode(data.plot$`Reward phase`,
                                "Acq1" = "Acquisition1",
                                "Acq2" = "Acquisition2",
                                "Bsln1" = "Baseline1",
                                "Bsln2" = "Baseline2",
                                "Ext1" = "Extinction1",
                                "Ext2" = "Extinction2")

data.plot$`Reward probability` = recode(data.plot$`Reward probability`,
                                      "High_Rew" = "High",
                                      "Low_Rew" = "Low")




  # Pirate plot
  pirateplot(formula=Hit.Rate ~ `Reward phase` + `Reward probability`, # dependent~independent variables
             data=data.plot, # data frame
             main = "Hit rates",
             ylim=c(0.1,0.9), # y-axis: limits
             ylab="Hit Rate", # y-axis: label
             theme=0, # preset theme (0: use your own)
             point.col="black", # points: color
             point.o=.3, # points: opacity (0-1)
             avg.line.col="black", # average line: color
             avg.line.lwd=2, # average line: line width
             avg.line.o=1, # average line: opacity (0-1)
             bean.b.col="black", # bean border, color
             bean.lwd=0.6, # bean border, line width
             bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
             bean.b.o=0.3, # bean border, opacity (0-1)
             bean.f.col="gray", # bean filling, color
             bean.f.o=.1, # bean filling, opacity (0-1)
             cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
             gl.col="gray", # gridlines: color
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=1, # axis labels: size
             cex.axis=1, # axis numbers: size
             cex.names = 1,
             bty="l", # plot box type
             back.col="white") # background, color
             
  
  
#  Plot Reaction times------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

  # Pirate plot
  pirateplot(formula = Hits.RTs ~ `Reward phase` + `Reward probability`, # dependent~independent variables
             data=data.plot, # data frame
             ylim=c(400,700), # y-axis: limits
             ylab="Reaction time", # y-axis: label
             main = "Reaction times",
             theme=0, # preset theme (0: use your own)
             point.col="black", # points: color
             point.o=.3, # points: opacity (0-1)
             avg.line.col="black", # average line: color
             avg.line.lwd=2, # average line: line width
             avg.line.o=1, # average line: opacity (0-1)
             bean.b.col="black", # bean border, color
             bean.lwd=0.6, # bean border, line width
             bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
             bean.b.o=0.3, # bean border, opacity (0-1)
             bean.f.col="gray", # bean filling, color
             bean.f.o=.1, # bean filling, opacity (0-1)
             cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
             gl.col="gray", # gridlines: color
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=1, # axis labels: size
             cex.axis=1, # axis numbers: size
             cex.names = 1,
             bty="l", # plot box type
             back.col="white") # background, color
  
  
################################################################## Stats ###############################################################################################################################################################################################################

# brms reaction times------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the working directory where to save the models
setwd(here("brms_models"))

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln1")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")


#Interaction model
# model.full.RT = brm(Hits.RTs ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
#                  data=data.final,
#                  family=gaussian(),
#                  warmup = 2000,
#                  iter = 10000,
#                  save_all_pars = TRUE,
#                  control = list(adapt_delta = 0.99),
#                  cores = 4,
#                  sample_prior = TRUE)
# saveRDS(model.full.RT,file="model.full.RT.6phases.rds")

# # read in the models and comparisons

model.full.RT = readRDS("model.full.RT.6phases.rds")
# compare.waic = readRDS("compare.RT.waic")


# Analyzing the posterior and differences between conditions
post = posterior_samples(model.full.RT, "^b")


################################################ Baseline 1 ####

######### High reward
Baseline1_High = post[["b_Intercept"]]
######### Low reward
Baseline1_Low = post[["b_Intercept"]] + 
  post[["b_ConditionLow_Rew"]] 

################################################ Baseline 2 ####

######### High reward
Baseline2_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseBsln2"]] 
######### Low reward
Baseline2_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseBsln2"]] +
  post[["b_ConditionLow_Rew"]] +
  post[["b_ExpPhaseBsln2:ConditionLow_Rew"]]

################################################ Acquistion 1

######### High reward
Acquisition1_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq1"]] 
######### Low reward
Acquisition1_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq1"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseAcq1:ConditionLow_Rew"]]

################################################ Acquistion 2

######### High reward
Acquisition2_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq2"]] 
######### Low reward
Acquisition2_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq2"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseAcq2:ConditionLow_Rew"]]

################################################ Extinction 1

######### High reward
Extinction1_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt1"]] 
######### Low reward
Extinction1_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt1"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseExt1:ConditionLow_Rew"]]

################################################ Extinction 2

######### High reward
Extinction2_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt2"]] 
######### Low reward
Extinction2_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt2"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseExt2:ConditionLow_Rew"]]


# Difference between high and low reward in baseline
Diff_Rew_Bsln = Baseline_High - Baseline_Low
plotPost(Diff_Rew_Bsln, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)

# Difference between high and low reward in acquisition 
Diff_Rew_Acq = Acquisition_High - Acquisition_Low
plotPost(Diff_Rew_Acq, xlab = "", col = "#b3cde0", cex = 1, showCurve = FALSE, compVal = 0)

# Difference between high and low reward in extinction
Diff_Rew_Ext = Extinction_High - Extinction_Low
plotPost(Diff_Rew_Ext, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)


########### plotting the posterior

# make a data frame
posterior_conditions = melt(data.frame(Baseline1_High, Baseline2_High, Baseline1_Low, Baseline2_Low, Acquisition1_High, Acquisition1_Low, Acquisition2_High, Acquisition2_Low, Extinction1_High, Extinction1_Low, Extinction2_High, Extinction2_Low))

posterior_conditions =  posterior_conditions %>% separate(variable, c("Reward Phase", "Reward Probability"), "_", extra = "merge")

names(posterior_conditions)[3] = "Reaction time"

# Pirate plot
pirateplot(formula = `Reaction time` ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
           data=posterior_conditions, # data frame
           main="Reaction times", # main title
           ylim=c(485,600), # y-axis: limits
           ylab="Reaction time", # y-axis: label
           theme=0, # preset theme (0: use your own)
           avg.line.col="black", # average line: color
           avg.line.lwd=2, # average line: line width
           avg.line.o=1, # average line: opacity (0-1)
           bean.b.col="black", # bean border, color
           bean.lwd=0.6, # bean border, line width
           bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
           bean.b.o=0.3, # bean border, opacity (0-1)
           bean.f.col="gray", # bean filling, color
           bean.f.o=.1, # bean filling, opacity (0-1)
           cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
           gl.col="gray", # gridlines: color
           gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
           cex.lab=1, # axis labels: size
           cex.axis=1, # axis numbers: size
           cex.names = 1,
           sortx = "sequential",
           bty="l", # plot box type
           back.col="white") # background, color
           

# brms accuracy------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln1")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")


#Interaction model
# model.full.Acc = brm(Hit.Rate ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
#                  data=data.final,
#                  family=gaussian(),
#                  warmup = 2000,
#                  iter = 10000,
#                  save_all_pars = TRUE,
#                  control = list(adapt_delta = 0.99),
#                  cores = 4,
#                  sample_prior = TRUE)
# saveRDS(model.full.Acc,file="model.full.Acc.6phases.rds")

#read in the models and comparisons
model.full.Acc = readRDS("model.full.Acc.6phases.rds")

# Analyzing the posterior and differences between conditions

post = posterior_samples(model.full.Acc, "^b")


################################################ Baseline 1 ####

######### High reward
Baseline1_High = post[["b_Intercept"]]
######### Low reward
Baseline1_Low = post[["b_Intercept"]] + 
  post[["b_ConditionLow_Rew"]] 

################################################ Baseline 2 ####

######### High reward
Baseline2_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseBsln2"]] 
######### Low reward
Baseline2_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseBsln2"]] +
  post[["b_ConditionLow_Rew"]] +
  post[["b_ExpPhaseBsln2:ConditionLow_Rew"]]

################################################ Acquistion 1

######### High reward
Acquisition1_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq1"]] 
######### Low reward
Acquisition1_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq1"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseAcq1:ConditionLow_Rew"]]

################################################ Acquistion 2

######### High reward
Acquisition2_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq2"]] 
######### Low reward
Acquisition2_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq2"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseAcq2:ConditionLow_Rew"]]

################################################ Extinction 1

######### High reward
Extinction1_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt1"]] 
######### Low reward
Extinction1_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt1"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseExt1:ConditionLow_Rew"]]

################################################ Extinction 2

######### High reward
Extinction2_High = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt2"]] 
######### Low reward
Extinction2_Low = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt2"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseExt2:ConditionLow_Rew"]]


# Difference between high and low reward in baseline
Diff_Rew_Bsln = Baseline_High - Baseline_Low
plotPost(Diff_Rew_Bsln, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)

# Difference between high and low reward in acquisition 
Diff_Rew_Acq = Acquisition_High - Acquisition_Low
plotPost(Diff_Rew_Acq, xlab = "", col = "#b3cde0", cex = 1, showCurve = FALSE, compVal = 0)

# Difference between high and low reward in extinction
Diff_Rew_Ext = Extinction_High - Extinction_Low
plotPost(Diff_Rew_Ext, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)


########### plotting the posterior

# make a data frame
posterior_conditions = melt(data.frame(Baseline1_High, Baseline2_High, Baseline1_Low, Baseline2_Low, Acquisition1_High, Acquisition1_Low, Acquisition2_High, Acquisition2_Low, Extinction1_High, Extinction1_Low, Extinction2_High, Extinction2_Low))

posterior_conditions =  posterior_conditions %>% separate(variable, c("Reward Phase", "Reward Probability"), "_", extra = "merge")

names(posterior_conditions)[3] = "Hit rate"

# Pirate plot
pirateplot(formula = `Hit rate` ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
           data=posterior_conditions, # data frame
           main="Hit rates", # main title
           ylim=c(0.5, 0.7), # y-axis: limits
           ylab="Hit rate", # y-axis: label
           theme=0, # preset theme (0: use your own)
           avg.line.col="black", # average line: color
           avg.line.lwd=2, # average line: line width
           avg.line.o=1, # average line: opacity (0-1)
           bean.b.col="black", # bean border, color
           bean.lwd=0.6, # bean border, line width
           bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
           bean.b.o=0.3, # bean border, opacity (0-1)
           bean.f.col="gray", # bean filling, color
           bean.f.o=.1, # bean filling, opacity (0-1)
           cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
           gl.col="gray", # gridlines: color
           gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
           cex.lab=1, # axis labels: size
           cex.axis=1, # axis numbers: size
           cex.names = 1,
           sortx = "sequential",
           bty="l", # plot box type
           back.col="white") # background, color


