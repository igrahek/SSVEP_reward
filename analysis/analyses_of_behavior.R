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
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2)
# set seed
set.seed(42) 
# import data
data.raw = read.csv(file="./data/Data_behavior_exp1_48pps.csv",header=TRUE,na.strings="NaN") 

# Prepare the dataset------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Adding and renaming variables 
# rename EventType variable
data.raw = rename(data.raw,c(EventType="MovedDots")) 
# add a variable with the name of the attended color instead of a numbers
data.raw$AttendedColor = ifelse(data.raw$AttendedColor==1,"red","blue")
# add a variable saying which color was linked with high reward (even numbers - blue was high reward)
data.raw$RewardedColor = ifelse(data.raw$ParticipantNo%%2==0,"blue","red") 
# add a variable with the name of the moved color instead of a numbers
data.raw$MovedDots = ifelse(data.raw$MovedDots==1,"red","blue") 
# split experimental phases into 6 isntead of 3 phases (trial 0-200: baseline; trial 201-400: acquisition; trial 401-600: extinction)
#data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,100,200,300,400,500,600),labels=c("baseline1","baseline2","acquisition1","acquisition2","extinction1","extinction2")) 
# split experimental phases into 3 phases (trial 0-200: baseline; trial 201-400: acquisition; trial 401-600: extinction)
data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,200,400,600),labels=c("baseline","acquisition","extinction")) # trial 0-200: baseline; trial 201-400: acquisition; trial 401-600: extinction

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
# eliminate ouotliers from data frame
#data.final = data.final[!data.final$ParticipantNo %in% unique(criterion),] 

# Create the final dataframe for accuracy and RTs ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or low rewarded color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"RewAtt","NotRewAtt")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)


################################################################## STATS ###############################################################################################################################################################################################################

# number of MonteCarlo iterations (default: 10000)
num.iter=10000 

# DV: dprime
# 3 x 2 rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_dPrime.bf.med <- anovaBF(dprime~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_dPrime.bf.med

# DV: Hit.Rate
# 3 x 2 rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_HitRate.bf.med <- anovaBF(Hit.Rate~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_HitRate.bf.med

# DV: FA rate
# 3 x 2 rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_FArate.bf.med <- anovaBF(FA.Rate~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_FArate.bf.med

# DV: RTs for Hits 
# 3 x 2 rANOVA

# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_RT.bf.med <- anovaBF(Hits.RTs~ExpPhase*Condition+ParticipantNo,data=data.final,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_RT.bf.med

#### brms analysis  ####
library(brms)
# referencing for easier interpretation
data.summaryRT.wide.ssj$ExpPhase=relevel(data.summaryRT.wide.ssj$ExpPhase,ref="acquisition")
data.summaryRT.wide.ssj$Condition=relevel(data.summaryRT.wide.ssj$Condition,ref="RewAtt")

#data.summaryRT.wide.ssj = within(data.summaryRT.wide.ssj,Condition=relevel(Condition,ref="RewAtt")) # reference Phase to high reward
# Null model
model.null <-brm(Hits.RTs~1+(1|ParticipantNo),
                 data=data.summaryRT.wide.ssj,
                 family=exgaussian(),
                 warmup = 2000, #200 & 10000
                 iter = 10000)
saveRDS(model.null,file="nullmodel.RT.rds")

# ExpPhase model
model.expphase <-brm(Hits.RTs~ExpPhase+(1|ParticipantNo),
                 data=data.summaryRT.wide.ssj,
                 family=exgaussian(),
                 warmup = 2000, #200 & 10000
                 iter = 10000)
saveRDS(model.expphase,file="expphasemodel.RT.rds")

# Condition model
model.condition <-brm(Hits.RTs~Condition+(1|ParticipantNo),
                     data=data.summaryRT.wide.ssj,
                     family=exgaussian(),
                     warmup = 2000, #200 & 10000
                     iter = 10000)
saveRDS(model.expphase,file="model.condition.RT.rds")

# # Two main effects model
model.twomaineffects <-brm(Hits.RTs~ExpPhase+Condition+(1|ParticipantNo),
                 data=data.summaryRT.wide.ssj,
                 family=exgaussian(),
                 warmup = 2000, #200 & 10000
                 iter = 10000)
saveRDS(model.twomaineffects,file="model.twomaineffects.RT.rds")

#Interaction model
model.full <-brm(Hits.RTs~ExpPhase+Condition+ExpPhase*Condition+(1|ParticipantNo),
                 data=data.summaryRT.wide.ssj,
                 family=exgaussian(),
                 warmup = 2000, #200 & 10000
                 iter = 10000)
saveRDS(model.full,file="model.full.RT.rds")

#Interaction model with random effect of condition 
model.full.random.cond <-brm(Hits.RTs~ExpPhase+Condition+ExpPhase*Condition+(Condition|ParticipantNo),
                        data=data.summaryRT.wide.ssj,
                        family=exgaussian(),
                        warmup = 2000, #200 & 10000
                        iter = 10000)
saveRDS(model.full.random.cond,file="model.full.random.cond.RT.rds")

#Interaction model with random effects of condition and experiment phase
model.full.random <-brm(Hits.RTs~ExpPhase+Condition+ExpPhase*Condition+(Condition*ExpPhase|ParticipantNo),
                 data=data.summaryRT.wide.ssj,
                 family=exgaussian(),
                 warmup = 2000, #200 & 10000
                 iter = 10000)
saveRDS(model.full.random,file="model.full.random.RT.rds")

#Interaction model with random effect of exp phase 
model.full.random.expphase <-brm(Hits.RTs~ExpPhase+Condition+ExpPhase*Condition+(ExpPhase|ParticipantNo),
                             data=data.summaryRT.wide.ssj,
                             family=exgaussian(),
                             warmup = 2000, #200 & 10000
                             iter = 10000)
saveRDS(model.full.random.expphase,file="model.full.random.expphase.RT.rds")

#Interaction model with random effects of condition and experiment phase
model.full.random1 <-brm(Hits.RTs~ExpPhase+Condition+ExpPhase*Condition+(ExpPhase*Condition|ParticipantNo),
                        data=data.summaryRT.wide.ssj,
                        family=exgaussian(),
                        warmup = 2000, #200 & 10000
                        iter = 10000)
saveRDS(model.full.random1,file="model.full.random1.RT.rds")

#Interaction model with random effects of condition and experiment phase
model.full.random2 <-brm(Hits.RTs~ExpPhase+Condition+ExpPhase*Condition+(ExpPhase*Condition|ParticipantNo),
                         data=data.summaryRT.wide.ssj,
                         family=exgaussian(),
                         warmup = 2000, #200 & 10000
                         iter = 10000)
saveRDS(model.full.random2,file="model.full.random2.RT.rds")

#WAIC
compare.RT.waic <- WAIC(model.null,model.condition,model.expphase,model.twomaineffects,model.full.random.expphase,model.full.random.cond,model.full.random1,model.full.random2)
saveRDS(compare.RT.waic,file="compare.RT.waic")

# #LOO crossvalidation
compare.RT.loo <- LOO(model.null,model.condition,model.expphase,model.twomaineffects,model.full,reloo=TRUE)
saveRDS(compare.RT.loo,file="compare.RT.loo")



# read in the models and comparisons
model.null <- readRDS("nullmodel.RT.rds")
model.condition <- readRDS("model.condition.RT.rds")
model.expphase <- readRDS("expphasemodel.RT.rds")
model.twomaineffects <- readRDS("model.twomaineffects.RT.rds")
model.full <- readRDS("model.RT.full.rds")
compare.loo = readRDS("compare.RT.loo")
compare.waic = readRDS("compare.RT.waic")

launch_shinystan(model.full.random1)

# Plot marginal effects for each predictor
plot(marginal_effects(model.full.random2),ask=FALSE)


