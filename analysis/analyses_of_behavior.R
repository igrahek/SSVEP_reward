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
crit = .6 # minimum 60% hit rate in any condition .6
# select participants below the criterion
criterion = subset(ddply(data.final,.(ParticipantNo),summarize,mean.Hit.Rate=mean(Hit.Rate)),mean.Hit.Rate<crit)$Participant # minimum 60% hit rate across all conditions
# eliminate ouotliers from data frame
data.final = data.final[!data.final$ParticipantNo %in% unique(criterion),] 

# Create the final dataframe for accuracy and RTs ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or low rewarded color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"RewAtt","NotRewAtt")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)


################################################################## STATS ###############################################################################################################################################################################################################

dataSDT.ssj$Condition <- as.factor(dataSDT.ssj$Condition) # convert to factor
dataSDT.ssj$RewardedColor <- as.factor(dataSDT.ssj$RewardedColor) # convert to factor
data.summaryRT.wide.ssj$Condition <- as.factor(data.summaryRT.wide.ssj$Condition) # convert to factor
data.summaryRT.wide.ssj$RewardedColor <- as.factor(data.summaryRT.wide.ssj$RewardedColor) # convert to factor
num.iter=10000 # number of MonteCarlo iterations (default: 10000)

# DV: dprime

#Deleting one participant with a missing RT average for one condition
# not necessary anymore: VP14 is already considered outlier
# dataSDT.ssj <- subset(dataSDT.ssj,ParticipantNo!=14)


# 3 x 2 rANOVA
# Assuming a medium prior d~Cauchy(0,.707):


rANOVA_dPrime.bf.med <- anovaBF(dprime~ExpPhase*Condition+ParticipantNo,data=dataSDT.ssj,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_dPrime.bf.med

# DV: Hit.Rate
# 3 x 2 rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_HitRate.bf.med <- anovaBF(Hit.Rate~ExpPhase*Condition+ParticipantNo,data=dataSDT.ssj,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_HitRate.bf.med

HitRate=subset(dataSDT.ssj,select=c(ParticipantNo,ExpPhase,Condition,Hit.Rate))
HitRate=dcast(HitRate,ParticipantNo~ExpPhase+Condition, value.var="Hit.Rate")
write.csv2(HitRate,file="Hit Rates.csv")

# DV: FA rate
# 3 x 2 rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_FArate.bf.med <- anovaBF(FA.Rate~ExpPhase*Condition+ParticipantNo,data=dataSDT.ssj,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_FArate.bf.med

# DV: RTs for Hits 
# 3 x 2 rANOVA
#Deleting one participant with a missing RT average for one condition
# not necessary anymore: VP14 is already considered outlier

#data.summaryRT.wide.ssj=subset(data.summaryRT.wide.ssj, data.summaryRT.wide.ssj$ParticipantNo!=14 & data.summaryRT.wide.ssj$ParticipantNo!=34)

# Assuming a medium prior d~Cauchy(0,.707):
rANOVA_RT.bf.med <- anovaBF(Hits.RTs~ExpPhase*Condition+ParticipantNo,data=data.summaryRT.wide.ssj,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_RT.bf.med

#### lmBF analysis ####
#### lmBF only random effect of subject ####

## Null model
bfNull = lmBF(Hits.RTs ~ ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj)
## Reward magnitude
bfCondition = lmBF(Hits.RTs ~ Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj)
## Block
bfExpPhase = lmBF(Hits.RTs ~ ExpPhase + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj)
## Main effects
bfMainEffects = lmBF(Hits.RTs ~ ExpPhase + Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj)
## Interaction
bfInteraction = lmBF(Hits.RTs ~ ExpPhase + Condition + ExpPhase:Condition + ParticipantNo, whichRandom = "ParticipantNo", data = data.summaryRT.wide.ssj)

bfCondition/bfNull
bfExpPhase/bfNull
bfMainEffects/bfNull
bfInteraction/bfNull

(bfExpPhase/bfNull)/(bfCondition/bfNull)
(bfExpPhase/bfNull)/(bfMainEffects/bfNull)
(bfExpPhase/bfNull)/(bfInteraction/bfNull)

chains = posterior(bfMainEffects, iterations = 10000)
summary(chains) 
plot(chains)

#### lmBF all random effects ####
## Null model
bfNull = lmBF(Hits.RTs ~ 1 + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj)
## Reward magnitude
bfCondition = lmBF(Hits.RTs ~ Condition + ParticipantNo + Condition:ParticipantNo, whichRandom = c("ParticipantNo","Condition:ParticipantNo"),data = data.summaryRT.wide.ssj)
## Block
bfExpPhase = lmBF(Hits.RTs ~ ExpPhase + ParticipantNo + ExpPhase:ParticipantNo, whichRandom = c("ParticipantNo","ExpPhase:ParticipantNo"),data = data.summaryRT.wide.ssj)
## Main effects
bfMainEffects = lmBF(Hits.RTs ~ ExpPhase + Condition + ParticipantNo + Condition:ParticipantNo + ExpPhase:ParticipantNo, whichRandom = c("ParticipantNo","Condition:ParticipantNo","ExpPhase:ParticipantNo"),data = data.summaryRT.wide.ssj)
## Interaction
bfInteraction = lmBF(Hits.RTs ~ ExpPhase + Condition + ExpPhase:Condition + ParticipantNo + Condition:ExpPhase:ParticipantNo, whichRandom = c("ParticipantNo","Condition:ExpPhase:ParticipantNo"), data = data.summaryRT.wide.ssj)

bfCondition/bfNull
bfExpPhase/bfNull
bfMainEffects/bfNull
bfInteraction/bfNull

(bfExpPhase/bfNull)/(bfCondition/bfNull)
(bfExpPhase/bfNull)/(bfMainEffects/bfNull)
(bfExpPhase/bfNull)/(bfInteraction/bfNull)

chains = posterior(bfMainEffects, iterations = 10000)
summary(chains) 
plot(chains)

#### lmBF sensitivity analysis ####

# priors for fixed effects
#create a list of priors
prior=seq(0.1,1,0.1)
#initialize a data frame to store the results
sensitivity = as.data.frame(prior)


# calculate BFs for RTs

sensitivity.1 = ddply(sensitivity,.(prior),summarise,
                      #Exp Phase versus the Null model
                      ExpPhase.vs.Null.model=extractBF(
                        lmBF(Hits.RTs ~ ExpPhase + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior)/
                        lmBF(Hits.RTs ~ ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior))$bf,
                      #Condition versus the Null model
                      Condition.vs.Null.model=extractBF(
                        lmBF(Hits.RTs ~ Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior)/
                        lmBF(Hits.RTs ~ ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior))$bf,
                      #Main effects versus the Null model
                      MainEffects.vs.Null.model=extractBF(
                        lmBF(Hits.RTs ~ ExpPhase + Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior)/
                        lmBF(Hits.RTs ~ ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior))$bf,
                      #Full model versus the Null model
                      Full.vs.Null.model=extractBF(
                        lmBF(Hits.RTs ~ ExpPhase + Condition + ExpPhase:Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior)/
                        lmBF(Hits.RTs ~ ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior))$bf,
                      #Main effects versus the Exp Phase model
                      MainEffects.vs.Full.model=extractBF(
                        lmBF(Hits.RTs ~ ExpPhase + Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior)/
                        lmBF(Hits.RTs ~ ExpPhase + Condition + ExpPhase:Condition + ParticipantNo, whichRandom = "ParticipantNo",data = data.summaryRT.wide.ssj,rscaleFixed = prior))$bf)

plot(prior,sensitivity.1$Condition.vs.Null.model)

# calculate BFs for Hit rates

sensitivity.1 = ddply(sensitivity,.(prior),summarise,
                      #Exp Phase versus the Null model
                      ExpPhase.vs.Null.model=extractBF(
                        lmBF(Hit.Rate ~ ExpPhase + ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior)/
                        lmBF(Hit.Rate ~ ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior))$bf,
                      #Condition versus the Null model
                      Condition.vs.Null.model=extractBF(
                        lmBF(Hit.Rate ~ Condition + ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior)/
                        lmBF(Hit.Rate ~ ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior))$bf,
                      #Main effects versus the Null model
                      MainEffects.vs.Null.model=extractBF(
                        lmBF(Hit.Rate ~ ExpPhase + Condition + ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior)/
                        lmBF(Hit.Rate ~ ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior))$bf,
                      #Full model versus the Null model
                      Full.vs.Null.model=extractBF(
                        lmBF(Hit.Rate ~ ExpPhase + Condition + ExpPhase:Condition + ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior)/
                        lmBF(Hit.Rate ~ ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior))$bf,
                      #Main effects versus the Exp Phase model
                      MainEffects.vs.Full.model=extractBF(
                        lmBF(Hit.Rate ~ ExpPhase + Condition + ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior)/
                        lmBF(Hit.Rate ~ ExpPhase + Condition + ExpPhase:Condition + ParticipantNo, whichRandom = "ParticipantNo",data = dataSDT.ssj,rscaleFixed = prior))$bf)

plot(prior,sensitivity.1$ExpPhase.vs.Null.model)

##### lme analysis ####
library(lme4)
m.null = lmer(Hits.RTs~1+(1|ParticipantNo),data.summaryRT.wide.ssj)
m.expphase = lmer(Hits.RTs~ExpPhase+(1|ParticipantNo),data.summaryRT.wide.ssj)
m.condition = lmer(Hits.RTs~Condition+(1|ParticipantNo),data.summaryRT.wide.ssj)
m.twomaineffects = lmer(Hits.RTs~Condition+ExpPhase+(1|ParticipantNo),data.summaryRT.wide.ssj)
m.full = lmer(Hits.RTs~ExpPhase*Condition+(1|ParticipantNo),data.summaryRT.wide.ssj)
m.full.random.condition = lmer(Hits.RTs~ExpPhase*Condition+(Condition|ParticipantNo),data.summaryRT.wide.ssj)
m.full.random.expphase = lmer(Hits.RTs~ExpPhase*Condition+(ExpPhase|ParticipantNo),data.summaryRT.wide.ssj)
m.full.random.both = lmer(Hits.RTs~ExpPhase*Condition+(Condition*ExpPhase|ParticipantNo),data.summaryRT.wide.ssj)

anova(m.null, m.expphase,m.condition, m.twomaineffects,m.full,m.full.random.condition,m.full.random.expphase) 



##plot
library(ggplot2)
tmp <- as.data.frame(confint(glht(m.full.random.both))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()
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


#### export data ####

#export only exp phase
RTs = ddply(data.summaryRT.wide.ssj,.(ParticipantNo,ExpPhase),summarize,
                              Hits.RTs1 = mean(Hits.RTs,na.rm=TRUE))
                             
RTs$Hits.RTs=RTs$Hits.RTs1

RTs=subset(RTs,select=c(ParticipantNo,ExpPhase,Hits.RTs))
RTs=dcast(RTs,ParticipantNo~ExpPhase, value.var="Hits.RTs")
write.csv2(RTs,file="RTsOnlyExpPhaseP.csv")

# export only the acquisition phase
RTs=subset(subset(data.summaryRT.wide.ssj,ExpPhase=="acquisition"),select=c(ParticipantNo,Condition,Hits.RTs))
RTs=dcast(RTs,ParticipantNo~Condition, value.var="Hits.RTs")
write.csv2(RTs,file="RTsOnlyAcquisition.csv")

# DV: RTs for Hits - rewarded color
# 3-way rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
dataRTsBayesRew = subset(data.summaryRT.wide.ssj,Condition=="RewAtt")
rANOVA_RT.bf.med.Rew <- anovaBF(Hits.RTs~ExpPhase+ParticipantNo,data=dataRTsBayesRew,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_RT.bf.med.Rew

# DV: RTs for Hits - not rewarded color
# 3-way rANOVA
# Assuming a medium prior d~Cauchy(0,.707):
dataRTsBayesNotRew = subset(data.summaryRT.wide.ssj,Condition=="NotRewAtt")
rANOVA_RT.bf.med.NotRew <- anovaBF(Hits.RTs~ExpPhase+ParticipantNo,data=dataRTsBayesNotRew,iterations=num.iter,whichRandom="ParticipantNo",rscaleRandom="nuisance",rscaleFixed=sqrt(2)/2)
rANOVA_RT.bf.med.NotRew


# #Export data
# 
Wide_hitrate=subset(dataSDT.ssj,select=c(ParticipantNo,Hit.Rate,Condition,ExpPhase))
Wide_hitrate=dcast(Wide_hitrate,ParticipantNo~Condition+ExpPhase, value.var="Hit.Rate")
write.csv2(Wide_hitrate,file="Wide_hitrate48pps_splitBlocks.csv")

Wide_dprime=subset(dataSDT.ssj,select=c(ParticipantNo,dprime,Condition,ExpPhase))
Wide_dprime=dcast(Wide_dprime,ParticipantNo~Condition+ExpPhase, value.var="dprime")
write.csv2(Wide_dprime,file="Wide_dprime48pps_splitBlocks.csv")

Wide_FARate=subset(dataSDT.ssj,select=c(ParticipantNo,FA.Rate,Condition,ExpPhase))
Wide_FARate=dcast(Wide_FARate,ParticipantNo~Condition+ExpPhase, value.var="FA.Rate")
write.csv2(Wide_FARate,file="Wide_FAs48pps_splitBlocks.csv")
# 
Wide_RTs=subset(data.summaryRT.wide.ssj,select=c(ParticipantNo,Hits.RTs,Condition,ExpPhase))
Wide_RTs=dcast(Wide_RTs,ParticipantNo~Condition+ExpPhase, value.var="Hits.RTs")
write.csv2(Wide_RTs,file="Wide_RTs_48pps.csv")
