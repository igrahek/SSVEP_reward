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
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2,brms,rstan)
# set seed
set.seed(42) 
# set directory
setwd("C:/Users/igrahek/Documents/Studies/SSVEP Reward - Soren & Antonio/Experiment 1/SSVEP and reward/")
# import data
data.raw = read.csv(file="./data/Data_behavior_exp1_48pps.csv",header=TRUE,na.strings="NaN") 

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
#data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,100,200,300,400,500,600),labels=c("Bsln1","Bsln2","Acq1","Acq2","Ext1","Ext2")) 
# split experimental phases into 3 phases (trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext)
data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,200,400,600),labels=c("Bsln","Acq","Ext")) # trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext

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
                      

# Handle outliers------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Determine outliers and cut them
# outliers based on hit rate at any condition
#crit = .6 # minimum 60% hit rate in any condition .6
# select participants below the criterion
#criterion = subset(ddply(data.final,.(ParticipantNo),summarize,mean.Hit.Rate=mean(Hit.Rate)),mean.Hit.Rate<crit)$Participant # minimum 60% hit rate across all conditions

#criterion = subset(data.final,data.final$Hit.Rate<.3)$Participant # minimum 60% hit rate across all conditions

# eliminate ouotliers from data frame
#data.final = data.final[!data.final$ParticipantNo %in% unique(criterion),] 

# Create the final dataframe for accuracy and RTs ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or Low_Rewed color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"High_Rew","Low_Rew")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)

################################################################## Plotting ###############################################################################################################################################################################################################

# # Plot Hit rates------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 
  # Pirate plot
  pirateplot(formula=Hit.Rate~ExpPhase+Condition, # dependent~independent variables
             data=data.final, # data frame
             main="Hit rates", # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(0,1), # y-axis: limits
             ylab="Hit Rate", # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=ParticipantNo, # ID variable
             theme=0, # preset theme (0: use your own)
             # theme settings
             # pal="xman", # color palette [see piratepal(palette="all")]
             point.col="black", # points: color
             point.o=.3, # points: opacity (0-1)
             avg.line.col="black", # average line: color
             avg.line.lwd=2, # average line: line width
             avg.line.o=1, # average line: opacity (0-1)
             bar.b.col=NULL, # bars, border: color
             bar.lwd=0, # bars, border: line width
             bar.b.o=0, # bars, border: opacity (0-1)
             bar.f.col=NULL, # bars, filling: color
             bar.f.o=0, # bars, filling: opacity (0-1)
             inf.b.col="black", # inference band, border: color
             inf.lwd=0.1, # inference band, border: line width
             inf.b.o=1, # inference band, border: opacity (0-1)
             inf.f.col="black", # inference band, filling: color
             inf.f.o=0, # inference band, filling: opacity (0-1)
             bean.b.col="black", # bean border, color
             bean.lwd=0.6, # bean border, line width
             bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
             bean.b.o=0.3, # bean border, opacity (0-1)
             bean.f.col="gray", # bean filling, color
             bean.f.o=.1, # bean filling, opacity (0-1)
             cap.beans=TRUE, # max and min values of bean densities are capped at the limits found in the data
             # quant=c(.1,.9), # quantiles (e.g., 10th and 90th)
             # quant.col="black", # quantiles, line: color
             # quant.length=.7, # quantiles, horizontal line length
             # quant.lwd=2, # quantiles, line width
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color

#  Plot Reaction times------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# 
  # Pirate plot
  pirateplot(formula=Hits.RTs~ExpPhase+Condition, # dependent~independent variables
             data=data.final, # data frame
             main="Reaction times", # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(400,700), # y-axis: limits
             ylab="Reaction time", # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=ParticipantNo, # ID variable
             theme=0, # preset theme (0: use your own)
             # theme settings
             # pal="xman", # color palette [see piratepal(palette="all")]
             point.col="black", # points: color
             point.o=.3, # points: opacity (0-1)
             avg.line.col="black", # average line: color
             avg.line.lwd=2, # average line: line width
             avg.line.o=1, # average line: opacity (0-1)
             bar.b.col=NULL, # bars, border: color
             bar.lwd=0, # bars, border: line width
             bar.b.o=0, # bars, border: opacity (0-1)
             bar.f.col=NULL, # bars, filling: color
             bar.f.o=0, # bars, filling: opacity (0-1)
             inf.b.col="black", # inference band, border: color
             inf.lwd=0.1, # inference band, border: line width
             inf.b.o=1, # inference band, border: opacity (0-1)
             inf.f.col="black", # inference band, filling: color
             inf.f.o=0, # inference band, filling: opacity (0-1)
             bean.b.col="black", # bean border, color
             bean.lwd=0.6, # bean border, line width
             bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
             bean.b.o=0.3, # bean border, opacity (0-1)
             bean.f.col="gray", # bean filling, color
             bean.f.o=.1, # bean filling, opacity (0-1)
             cap.beans=TRUE, # max and min values of bean densities are capped at the limits found in the data
             # quant=c(.1,.9), # quantiles (e.g., 10th and 90th)
             # quant.col="black", # quantiles, line: color
             # quant.length=.7, # quantiles, horizontal line length
             # quant.lwd=2, # quantiles, line width
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color
#   
################################################################## Stats ###############################################################################################################################################################################################################

# BF ANOVA------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# number of MonteCarlo iterations (default: 10000)
# num.iter=10000 
# 
# # DV: dprime
# # 3 x 2 rANOVA
# # Assuming a medium prior d~Cauchy(0,.707):
# rANOVA_dPrime.bf.med <- anovaBF(dprime~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
# rANOVA_dPrime.bf.med
# 
# # DV: Hit.Rate
# # 3 x 2 rANOVA
# # Assuming a medium prior d~Cauchy(0,.707):
# rANOVA_HitRate.bf.med <- anovaBF(Hit.Rate~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
# rANOVA_HitRate.bf.med
# 
# # DV: FA rate
# # 3 x 2 rANOVA
# # Assuming a medium prior d~Cauchy(0,.707):
# rANOVA_FArate.bf.med <- anovaBF(FA.Rate~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
# rANOVA_FArate.bf.med
# 
# # DV: RTs for Hits 
# # 3 x 2 rANOVA
# 
# # Assuming a medium prior d~Cauchy(0,.707):
# rANOVA_RT.bf.med <- anovaBF(Hits.RTs~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
# rANOVA_RT.bf.med
# 

# brms reaction times------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the working directory where to save the models
setwd("./brms_models")

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")


# Null model
model.null.RT = brm(Hits.RTs ~ 1 + (1|ParticipantNo),
                 data=data.final,
                 family=gaussian(),
                 warmup = 2000,
                 iter = 10000,
                 save_all_pars = TRUE,
                 control = list(adapt_delta = 0.99),
                 cores = 4)
saveRDS(model.null.RT,file="nullmodel.RT.rds")

# ExpPhase model
model.expphase.RT = brm(Hits.RTs ~ ExpPhase + (ExpPhase|ParticipantNo),
                 data=data.final,
                 family=gaussian(),
                 warmup = 2000,
                 iter = 10000,
                 save_all_pars = TRUE,
                 control = list(adapt_delta = 0.99),
                 cores = 4)
saveRDS(model.expphase.RT,file="expphasemodel.RT.rds")

# Condition model
model.condition.RT = brm(Hits.RTs ~ Condition + (Condition|ParticipantNo),
                     data=data.final,
                     family=gaussian(),
                     warmup = 2000,
                     iter = 10000,
                     save_all_pars = TRUE,
                     control = list(adapt_delta = 0.99),
                     cores = 4)
saveRDS(model.condition.RT,file="model.condition.RT.rds")

# # Two main effects model
model.twomaineffects.RT = brm(Hits.RTs ~ ExpPhase + Condition + (ExpPhase + Condition|ParticipantNo),
                 data=data.final,
                 family=gaussian(),
                 warmup = 2000,
                 iter = 10000,
                 save_all_pars = TRUE,
                 control = list(adapt_delta = 0.99),
                 cores = 4)
saveRDS(model.twomaineffects.RT,file="model.twomaineffects.RT.rds")

#Interaction model
model.full.RT = brm(Hits.RTs ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
                 data=data.final,
                 family=gaussian(),
                 warmup = 2000,
                 iter = 10000,
                 save_all_pars = TRUE,
                 control = list(adapt_delta = 0.99),
                 cores = 4)
saveRDS(model.full.RT,file="model.full.RT.rds")

# # read in the models and comparisons
# model.null.RT = readRDS("nullmodel.RT.rds")
# model.condition.RT = readRDS("model.condition.RT.rds")
# model.expphase.RT = readRDS("expphasemodel.RT.rds")
# model.twomaineffects.RT = readRDS("model.twomaineffects.RT.rds")
# model.full.RT = readRDS("model.full.RT.rds")
#compare.loo = readRDS("compare.RT.loo")
#compare.waic = readRDS("compare.RT.waic")

#WAIC
compare.RT.waic = WAIC(model.null.RT,model.condition.RT,model.expphase.RT,model.twomaineffects.RT,model.full.RT, compare = FALSE)
saveRDS(compare.RT.waic,file="compare.RT.waic")

# #LOO crossvalidation
# compare.RT.loo = LOO(model.null,model.condition,model.expphase,model.twomaineffects,model.full,reloo=TRUE) #,reloo=TRUE
# saveRDS(compare.RT.loo,file="compare.RT.loo")


# brms accuracy------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")


# Null model
model.null.Acc = brm(Hit.Rate ~ 1 + (1|ParticipantNo),
                 data=data.final,
                 family=gaussian(),
                 warmup = 2000,
                 iter = 10000,
                 save_all_pars = TRUE,
                 control = list(adapt_delta = 0.99),
                 cores = 4)
saveRDS(model.null.Acc,file="nullmodel.Acc.rds")

# ExpPhase model
model.expphase.Acc = brm(Hit.Rate ~ ExpPhase + (ExpPhase|ParticipantNo),
                     data=data.final,
                     family=gaussian(),
                     warmup = 2000,
                     iter = 10000,
                     save_all_pars = TRUE,
                     control = list(adapt_delta = 0.99),
                     cores = 4)
saveRDS(model.expphase.Acc,file="expphasemodel.Acc.rds")

# Condition model
model.condition.Acc = brm(Hit.Rate ~ Condition + (Condition|ParticipantNo),
                      data=data.final,
                      family=gaussian(),
                      warmup = 2000,
                      iter = 10000,
                      save_all_pars = TRUE,
                      control = list(adapt_delta = 0.99),
                      cores = 4)
saveRDS(model.expphase.Acc,file="model.condition.Acc.rds")

# # Two main effects model
model.twomaineffects.Acc = brm(Hit.Rate ~ ExpPhase + Condition + (ExpPhase + Condition|ParticipantNo),
                           data=data.final,
                           family=gaussian(),
                           warmup = 2000,
                           iter = 10000,
                           save_all_pars = TRUE,
                           control = list(adapt_delta = 0.99),
                           cores = 4)
saveRDS(model.twomaineffects.Acc,file="model.twomaineffects.Acc.rds")

#Interaction model
model.full.Acc = brm(Hit.Rate ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
                 data=data.final,
                 family=gaussian(),
                 warmup = 2000,
                 iter = 10000,
                 save_all_pars = TRUE,
                 control = list(adapt_delta = 0.99),
                 cores = 4)
saveRDS(model.full.Acc,file="model.full.Acc.rds")

#read in the models and comparisons
# model.null.Acc = readRDS("nullmodel.Acc.rds")
# model.condition.Acc = readRDS("model.condition.Acc.rds")
# model.expphase.Acc = readRDS("expphasemodel.Acc.rds")
# model.twomaineffects.Acc = readRDS("model.twomaineffects.Acc.rds")
#model.full.Acc = readRDS("model.full.Acc_nocorrelationvarying.rds")
#compare.loo.Acc = readRDS("compare.Acc.loo")
# compare.waic.Acc = readRDS("compare.Acc.waic")

#WAIC
compare.Acc.waic = WAIC(model.null.Acc,model.condition.Acc,model.expphase.Acc,model.twomaineffects.Acc,model.full.Acc, compare = FALSE)
saveRDS(compare.Acc.waic,file="compare.Acc.waic")

# #LOO crossvalidation
# compare.RT.loo = LOO(model.null,model.condition,model.expphase,model.twomaineffects,model.full,reloo=TRUE) #,reloo=TRUE
# saveRDS(compare.RT.loo,file="compare.Acc.loo")
