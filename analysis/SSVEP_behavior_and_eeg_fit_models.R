
#### Code info ####

# Experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Ernst Koster, & Søren Andersen) (*: co-first authors)
# Code written by: Ivan Grahek & Antonio Schettino (2016-2019)
# Description: Code for data preprocessing and the analysis of EEG and behavioral data for Experiment 1 of the SSVEP - reward project. 


#### EEG  ####

# Clear the environement and import the data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
# clear the console
cat("\014") 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(reshape2,yarrr,BayesFactor,plyr,ez,schoRsch,brms,lme4, BEST, brmstools, here, tidyverse)
# set seed
set.seed(42) 
# Set working directory
setwd(here())
# import data
data.raw = read.csv(file = here("EEG_preprocessing/movement","grandAverage_amplitudes.csv"),header=TRUE,na.strings="NaN") 
# Prepare the dataset------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Rename columns
colnames(data.raw) = c("Subject","Frequency","BslnRedAttendedNomov","BslnBlueAttendedNomov","AcqRedAttendedNomov","AcqBlueAttendedNomov","ExtRedAttendedNomov","ExtBlueAttendedNomov","BslnRedAttendedMov","BslnBlueAttendedMov","AcqRedAttendedMov","AcqBlueAttendedMov","ExtRedAttendedMov","ExtBlueAttendedMov")

# # Reshape to long format - take only no movement trials
data = melt(data.raw,id.vars=c("Subject","Frequency"),
            measure.vars=c("BslnRedAttendedNomov","BslnBlueAttendedNomov","AcqRedAttendedNomov","AcqBlueAttendedNomov","ExtRedAttendedNomov","ExtBlueAttendedNomov","BslnRedAttendedMov","BslnBlueAttendedMov","AcqRedAttendedMov","AcqBlueAttendedMov","ExtRedAttendedMov","ExtBlueAttendedMov"),
            variable.name="Condition",value.name="Amplitude")

# Sort the new dataframe by participant name
data = data[order(data$Subject),]

# Split the variable Condition based on capital letters
data$Condition = gsub("(?!^)(?=[[:upper:]])", " ", data$Condition, perl=T)

# Split the variable condition into multiple variables
Conditions = colsplit(data$Condition, pattern="\\s+",names = c('ExpPhase', 'ColorMoved',"attended","no","moved"))

# Add the variable defining which color is rewarded based on the participant number
data$RewardedColor = ifelse(data$Subject%%2==0,"Blue","Red") # if participant number is even, blue was rewarded

# Add the Conditions needed to the dataset
data$ExpPhase = Conditions[,1]
data$AttendedColor = Conditions[,2]
data$Movement = Conditions[,4]

# Switch the Frequency to the color
data$RecordedFrequency = ifelse(data$Frequency==10,"Blue","Red") # if the recorded frequency is 10Hz assign Blue (color flickering at 10Hz), otherwise assign Red (color flickering at 12Hz)

# Make a new condition based on the attended color and the rewarded color
data$Condition = ifelse(data$AttendedColor==data$RewardedColor, "High_Rew","Low_Rew")

# Make a new condition based on the attended color and the recorded frequency
data$Attention = ifelse(data$AttendedColor==data$RecordedFrequency, "Att","NotAtt")

# Make a new condition based the Condition and the Attention
data$RecordingAndCondition = with(data, paste0(Condition,"_",Attention))

# Select variables which we want to keep
data = subset(data, select=c("Subject","RewardedColor","ExpPhase","AttendedColor","Condition","RecordedFrequency","Attention","RecordingAndCondition","Movement","Amplitude"))

# Sort the data 
data = data[with(data, order(Subject)), ]

# Make a new variable with mean amplitude across all conditions for each participant and each frequency
data = ddply(data,.(Subject,RecordedFrequency),transform,
             MeanAmplitude = mean(Amplitude[ExpPhase=="Bsln"],na.rm=TRUE),
             SDAmplitude =   sd(Amplitude,na.rm=TRUE))

# Divide amplitudes in each Subject, Frequency, and Condition by the Mean Amplitude (normalization)
data$Amplitude = data$Amplitude/data$MeanAmplitude

# Calculate the reward index - High reward minus Low reward
data.reward = ddply(data, .(Subject,ExpPhase,Attention), transform, Reward = Amplitude[Condition=="High_Rew"]-Amplitude[Condition=="Low_Rew"])
# Delete the Attention column and rows which are not necessary (indexes repeated twice)
data.reward = subset(data.reward,Condition=="High_Rew") #keep only Att as it is equal to NotAtt
data.reward$Condition = NULL  #drop the Condition column

# Sort the data 
data.reward$ExpPhase = factor(data.reward$ExpPhase, levels = c("Bsln","Acq","Ext"))
data.reward = data.reward[order(data.reward$Subject,data.reward$Attention,data.reward$ExpPhase),]

# Convert variables to be used in analyses into factors
data[c("Subject", "Condition","ExpPhase", "RewardedColor", "Attention", "RecordingAndCondition")] = 
  lapply(data[c("Subject", "Condition","ExpPhase", "RewardedColor", "Attention", "RecordingAndCondition")], factor)

data.diff[c("Subject", "Condition","ExpPhase", "RecordingAndCondition")] = 
  lapply(data.diff[c("Subject", "Condition","ExpPhase",  "RecordingAndCondition")], factor)

# Save the final data into a new variable and delete the subject with the noisy EEG data
data_all = subset(data, Subject!=14)

# Fitting the models------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Select the data
data = data_all

# Set the working directory where to save the models
setwd(here("brms_models"))

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Set the intercept model
data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High_Rew")
data$Attention=relevel(data$Attention,ref="Att")

# Prior for the intercept only model
prior = c(
  prior(normal(1, 2), class = Intercept)) #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2


# Null model
null = brm(Amplitude ~ 1 + (1|Subject),
                              data=data,
                              family=gaussian(),
                              prior = prior,
                              warmup = 2000,
                              iter = 6000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(null,file="null.EEG.alltrials.rds")

# Priors for the models with slopes
prior = c(
  prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(normal(0, 3), class = b)) #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2


# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                                  data=data,
                                  family=gaussian(),
                                  prior = prior,
                                  warmup = 2000,
                                  iter = 6000,
                                  save_all_pars = TRUE,
                                  control = list(adapt_delta = 0.99),
                                  cores = 4,
                                  sample_prior = TRUE)
saveRDS(expphase,file="expphase.EEG.alltrials.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Attention|Subject),
                                   data=data,
                                   family=gaussian(),
                                   prior = prior,  
                                   warmup = 2000,
                                   iter = 6000,
                                   save_all_pars = TRUE,
                                   control = list(adapt_delta = 0.99),
                                   cores = 4,
                                   sample_prior = TRUE)
saveRDS(attention,file="attention.EEG.alltrials.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + (ExpPhase + Attention|Subject),
                                        data=data,
                                        family=gaussian(),
                                        prior = prior,
                                        warmup = 2000,
                                        iter = 6000,
                                        save_all_pars = TRUE,
                                        control = list(adapt_delta = 0.99),
                                        cores = 4,
                                        sample_prior = TRUE)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.alltrials.rds")

# Interaction between phase and attention
phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + (ExpPhase + Attention|Subject),
                        data=data,
                        family=gaussian(),
                        prior = prior,
                        warmup = 2000,
                        iter = 6000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99),
                        cores = 4,
                        sample_prior = TRUE)
saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.alltrials.rds")

# Interaction between expphase and reward magnitude plus attention
rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (ExpPhase + Attention + Condition|Subject),
                         data=data,
                         family=gaussian(),
                         prior = prior,
                         warmup = 2000,
                         iter = 6000,
                         save_all_pars = TRUE,
                         control = list(adapt_delta = 0.99),
                         cores = 4,
                         sample_prior = TRUE)
saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.alltrials.rds")

# Full model
full = brm(Amplitude ~ Condition * ExpPhase * Attention + (ExpPhase + Attention + Condition|Subject),
                              data=data,
                              family=gaussian(),
                              prior = prior,
                              warmup = 2000,
                              iter = 6000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(full,file="full.EEG.alltrials.rds")


# WAIC
compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.alltrials.rds")

# Bayesian R2
#Null
bR2.null.EEG = bayes_R2(null)
saveRDS(bR2.null.EEG,file="bR2.null.EEG.alltrials")
#ExpPhase
bR2.expphase.EEG = bayes_R2(expphase)
saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG.alltrials")
#Attention
bR2.attention.EEG = bayes_R2(attention)
saveRDS(bR2.attention.EEG,file="bR2.attention.EEG.alltrials")
#Phase and attention
bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG.alltrials")
#Phase and attention interaction
bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG.alltrials")
#Reward times phase plus attention
bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG.alltrials")
#Full
bR2.full.EEG = bayes_R2(full)
saveRDS(bR2.full.EEG,file="bR2.full.EEG.alltrials")









#### Behavior ####

# Clear the environement and import the data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
# clear the console
cat("\014") 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2,brms, broom, tidyverse, brmstools, BEST, knitr, here,psych)
# set seed
set.seed(42) 
# set directory
setwd(here())
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
#data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,100,200,300,400,500,600),labels=c("Bsln1","Bsln2","Acq1","Acq2","Ext1","Ext2")) 
# split experimental phases into 3 phases (trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext)
data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,200,400,600),labels=c("Bsln","Acq","Ext")) # trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext

### Convert variables to be used in analyses into factors
data.raw[c("ParticipantNo", "AttendedColor","RewardedColor", "MovedDots", "ExpPhase" )] = 
  lapply(data.raw[c("ParticipantNo", "AttendedColor","RewardedColor", "MovedDots", "ExpPhase" )], factor)

### Create variables needed for the accuracy analyses
# count hits, false alarms, misses, correct rejections, and RT separately for each participant (their calculation is done in Matlab: see DataProcessing.m)
data.final = ddply(data.raw,.(ParticipantNo,ExpPhase,AttendedColor,RewardedColor,MovedDots),summarise,
                   numtrials=length(which(Response!=99)), # number of trials per condition (anything that is not 99 or any other number that we're not using)
                   Hits=length(which(Response==1)), # hits: attended color moved, correct response
                   FAs=length(which(Response==2)), # false alarms: attended color did not move, (wrong) response
                   Misses=length(which(Response==0)), # misses: attended color moved, no response
                   CRs=length(which(Response==3)), # correct rejections: attended color did not move, no response
                   mean.RT=mean(RT,na.rm=TRUE)) # mean RT per condition



# Exclude subjects without EEG
missing_eeg = c(5,10,13,27)
data.final = data.final[!data.final$ParticipantNo %in% missing_eeg,]
# Calculate the hits and false alarms ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
data.final =  ddply(data.final,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,numtrials),summarise,
                    tot.Hits=Hits+.5, # hits
                    tot.FAs=FAs+.5, # false alarms
                    tot.Misses=(numtrials-Hits)+.5, # misses
                    tot.CRs=(numtrials-FAs)+.5, # correct rejections
                    Hit.Rate=tot.Hits/(tot.Hits+tot.Misses), # hit rate
                    FA.Rate=tot.FAs/(tot.FAs+tot.CRs), # false alarm rate
                    dprime_by_hand=qnorm(Hit.Rate)-qnorm(FA.Rate),
                    Hits.RTs=mean(mean.RT,na.rm=TRUE)) # mean RTs

# Calculate SDT indices with psycho
indices = psycho::dprime(data.final$tot.Hits, data.final$tot.FAs, data.final$tot.Misses, data.final$tot.CRs) 

data.final = cbind(data.final, indices) 


### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or Low_Rewed color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"High_Rew","Low_Rew")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)


# Exclude subjects without EEG
missing_eeg = c(5,10,13,27,14)
data.final = data.final[!data.final$ParticipantNo %in% missing_eeg,]

# Fitting the models------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Reaction times------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the working directory where to save the models
setwd(here("brms_models"))

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")

# Set the prior for the intercept only model
prior = c(
  prior(normal(500, 200), class = Intercept)) # A wide prior sensible for this type of task


# Null model
model.null.RT = brm(Hits.RTs ~ 1 + (1|ParticipantNo),
                    data=data.final,
                    family=gaussian(),
                    prior = prior,
                    iter = 6000,
                    save_all_pars = TRUE,
                    control = list(adapt_delta = 0.99,max_treedepth = 15),
                    cores = 4,
                    sample_prior = TRUE,
                    inits = 0)
saveRDS(model.null.RT,file="nullmodel.RT.rds")

# Set the priors for the models with slope
prior = c(
  prior(normal(500, 200), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 200), class = b)) # a wide prior

# ExpPhase model
model.expphase.RT = brm(Hits.RTs ~ ExpPhase + (ExpPhase|ParticipantNo),
                        data=data.final,
                        family=gaussian(),
                        prior = prior,
                        iter = 6000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99,max_treedepth = 15),
                        cores = 4,
                        sample_prior = TRUE,
                        inits = 0)
saveRDS(model.expphase.RT,file="expphasemodel.RT.rds")

#Interaction model
model.full.RT = brm(Hits.RTs ~ ExpPhase * Condition + (ExpPhase + Condition|ParticipantNo),
                    data=data.final,
                    family=gaussian(),
                    prior = prior,
                    iter = 6000,
                    save_all_pars = TRUE,
                    control = list(adapt_delta = 0.99,max_treedepth = 15),
                    cores = 4,
                    sample_prior = TRUE,
                    inits = 0)
saveRDS(model.full.RT,file="model.full.RT.rds")



#WAIC
compare.RT.waic = WAIC(model.null.RT,model.expphase.RT,model.full.RT, comapre = TRUE)
saveRDS(compare.RT.waic,file="compare.RT.waic")


# Bayesian R2
#Null
bR2.null.RT = bayes_R2(model.null.RT)
saveRDS(bR2.null.RT,file="bR2.null.RT")
#ExpPhase
bR2.expphase.RT = bayes_R2(model.expphase.RT)
saveRDS(bR2.expphase.RT,file="bR2.expphase.RT")
#Full
bR2.full.RT = bayes_R2(model.full.RT)
saveRDS(bR2.full.RT,file="bR2.full.RT")


# Accuracy (d prime)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set the working directory where to save the models
setwd(here("brms_models"))
# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")

# Set the priors for the model with only intercept
prior = c(
  prior(normal(1.8, 1), class = Intercept)) # based on Andersen & Mueller, 2011


# Null model
model.null.Acc.dprime = brm(dprime ~ 1 + (1|ParticipantNo),
                            data=data.final,
                            family=gaussian(),
                            prior = prior,
                            iter = 6000,
                            save_all_pars = TRUE,
                            control = list(adapt_delta = 0.99,max_treedepth = 15),
                            cores = 4,
                            sample_prior = TRUE,
                            inits = 0)
saveRDS(model.null.Acc.dprime,file="nullmodel.Acc.dprime.rds")

# Set the priors for the models with slopes
prior = c(
  prior(normal(1.8, 1), class = Intercept), # based on Andersen & Mueller, 2011
  prior(normal(0, 2), class = b)) # a wide prior

# ExpPhase model
model.expphase.Acc.dprime = brm(dprime ~ ExpPhase + (ExpPhase|ParticipantNo),
                                data=data.final,
                                family=gaussian(),
                                prior = prior,
                                iter = 6000,
                                save_all_pars = TRUE,
                                control = list(adapt_delta = 0.99,max_treedepth = 15),
                                cores = 4,
                                sample_prior = TRUE,
                                inits = 0)
saveRDS(model.expphase.Acc.dprime,file="expphasemodel.Acc.dprime.rds")

#Interaction model
model.full.Acc.dprime = brm(dprime ~ ExpPhase * Condition + (ExpPhase + Condition|ParticipantNo),
                            data=data.final,
                            family=gaussian(),
                            prior = prior,
                            iter = 6000,
                            save_all_pars = TRUE,
                            control = list(adapt_delta = 0.99,max_treedepth = 15),
                            cores = 4,
                            sample_prior = TRUE,
                            inits = 0)
saveRDS(model.full.Acc.dprime,file="model.full.Acc.dprime.rds")


#WAIC
compare.Acc.dprime.waic = WAIC(model.null.Acc.dprime,model.expphase.Acc.dprime,model.full.Acc.dprime, compare = TRUE)
saveRDS(compare.Acc.dprime.waic,file="compare.Acc.dprime.waic")

# Weighted waic
compare.Acc.dprime.waic.weights = model_weights(model.null.Acc.dprime,model.expphase.Acc.dprime,model.full.Acc.dprime, weights = "waic")
saveRDS(compare.Acc.dprime.waic.weights,file="compare.Acc.dprime.waic.weights")

# Bayesian R2
#Null
bR2.null.Acc.dprime = bayes_R2(model.null.Acc.dprime)
saveRDS(bR2.null.Acc.dprime,file="bR2.null.Acc.dprime")
#ExpPhase
bR2.expphase.Acc.dprime = bayes_R2(model.expphase.Acc.dprime)
saveRDS(bR2.expphase.Acc.dprime,file="bR2.expphase.Acc.dprime")
#Full
bR2.full.Acc.dprime = bayes_R2(model.full.Acc.dprime)
saveRDS(bR2.full.Acc.dprime,file="bR2.full.Acc.dprime")





#### Behavior - splitting phases (supplementary analyses) #########

# Clear the environement and import the data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
# clear the console
cat("\014") 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2,brms, broom, tidyverse, brmstools, BEST, knitr, here,psych)
# set seed
set.seed(42) 
# set directory
setwd(here())
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
# data.raw$ExpPhase = cut(data.raw$Trial,breaks=c(0,200,400,600),labels=c("Bsln","Acq","Ext")) # trial 0-200: Bsln; trial 201-400: Acq; trial 401-600: Ext

### Convert variables to be used in analyses into factors
data.raw[c("ParticipantNo", "AttendedColor","RewardedColor", "MovedDots", "ExpPhase" )] = 
  lapply(data.raw[c("ParticipantNo", "AttendedColor","RewardedColor", "MovedDots", "ExpPhase" )], factor)

### Create variables needed for the accuracy analyses
# count hits, false alarms, misses, correct rejections, and RT separately for each participant (their calculation is done in Matlab: see DataProcessing.m)
data.final = ddply(data.raw,.(ParticipantNo,ExpPhase,AttendedColor,RewardedColor,MovedDots),summarise,
                   numtrials=length(which(Response!=99)), # number of trials per condition (anything that is not 99 or any other number that we're not using)
                   Hits=length(which(Response==1)), # hits: attended color moved, correct response
                   FAs=length(which(Response==2)), # false alarms: attended color did not move, (wrong) response
                   Misses=length(which(Response==0)), # misses: attended color moved, no response
                   CRs=length(which(Response==3)), # correct rejections: attended color did not move, no response
                   mean.RT=mean(RT,na.rm=TRUE)) # mean RT per condition



# Exclude subjects without EEG
missing_eeg = c(5,10,13,27)
data.final = data.final[!data.final$ParticipantNo %in% missing_eeg,]
# Calculate the hits and false alarms------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
data.final =  ddply(data.final,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,numtrials),summarise,
                    tot.Hits=Hits+.5, # hits
                    tot.FAs=FAs+.5, # false alarms
                    tot.Misses=(numtrials-Hits)+.5, # misses
                    tot.CRs=(numtrials-FAs)+.5, # correct rejections
                    Hit.Rate=tot.Hits/(tot.Hits+tot.Misses), # hit rate
                    FA.Rate=tot.FAs/(tot.FAs+tot.CRs), # false alarm rate
                    dprime_by_hand=qnorm(Hit.Rate)-qnorm(FA.Rate),
                    Hits.RTs=mean(mean.RT,na.rm=TRUE)) # mean RTs

# Calculate SDT indices with psycho
indices = psycho::dprime(data.final$tot.Hits, data.final$tot.FAs, data.final$tot.Misses, data.final$tot.CRs) 

data.final = cbind(data.final, indices) 


### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or Low_Rewed color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"High_Rew","Low_Rew")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)


# Exclude subjects without EEG
missing_eeg = c(5,10,13,27,14)
data.final = data.final[!data.final$ParticipantNo %in% missing_eeg,]



# Reaction times------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the working directory where to save the models
setwd(here("brms_models"))

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln1")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")

# Set the prior for the intercept only model
prior = c(
  prior(normal(500, 200), class = Intercept)) # A wide prior sensible for this type of task


# Null model
model.null.RT = brm(Hits.RTs ~ 1 + (1|ParticipantNo),
                    data=data.final,
                    family=gaussian(),
                    prior = prior,
                    iter = 6000,
                    save_all_pars = TRUE,
                    control = list(adapt_delta = 0.99,max_treedepth = 15),
                    cores = 4,
                    sample_prior = TRUE,
                    inits = 0)
saveRDS(model.null.RT,file="nullmodel.RT.split.rds")

# Set the priors for the models with slope
prior = c(
  prior(normal(500, 200), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 200), class = b)) # a wide prior

# ExpPhase model
model.expphase.RT = brm(Hits.RTs ~ ExpPhase + (ExpPhase|ParticipantNo),
                        data=data.final,
                        family=gaussian(),
                        prior = prior,
                        iter = 6000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99,max_treedepth = 15),
                        cores = 4,
                        sample_prior = TRUE,
                        inits = 0)
saveRDS(model.expphase.RT,file="expphasemodel.RT.split.rds")

#Interaction model
model.full.RT = brm(Hits.RTs ~ ExpPhase * Condition + (ExpPhase + Condition|ParticipantNo),
                    data=data.final,
                    family=gaussian(),
                    prior = prior,
                    iter = 6000,
                    save_all_pars = TRUE,
                    control = list(adapt_delta = 0.99,max_treedepth = 15),
                    cores = 4,
                    sample_prior = TRUE,
                    inits = 0)
saveRDS(model.full.RT,file="model.full.RT.split.rds")


#WAIC
compare.RT.waic = WAIC(model.null.RT,model.expphase.RT,model.full.RT, comapre = TRUE)
saveRDS(compare.RT.waic,file="compare.RT.split.waic")


# Weighted waic
compare.RT.waic.weights = model_weights(model.null.RT,model.expphase.RT,model.full.RT, weights = "waic")
saveRDS(compare.RT.waic.weights,file="compare.RT.split.waic.weights")


# Bayesian R2
#Null
bR2.null.RT = bayes_R2(model.null.RT)
saveRDS(bR2.null.RT,file="bR2.null.RT.split")
#ExpPhase
bR2.expphase.RT = bayes_R2(model.expphase.RT)
saveRDS(bR2.expphase.RT,file="bR2.expphase.RT.split")
#Full
bR2.full.RT = bayes_R2(model.full.RT)
saveRDS(bR2.full.RT,file="bR2.full.RT.split")

# Accuracy (d prime)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set the working directory where to save the models
setwd(here("brms_models"))
# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln1")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")

# Set the priors for the model with only intercept
prior = c(
  prior(normal(1.8, 1), class = Intercept)) # based on Andersen & Mueller, 2011


# Null model
model.null.Acc.dprime = brm(dprime ~ 1 + (1|ParticipantNo),
                            data=data.final,
                            family=gaussian(),
                            prior = prior,
                            iter = 6000,
                            save_all_pars = TRUE,
                            control = list(adapt_delta = 0.99,max_treedepth = 15),
                            cores = 4,
                            sample_prior = TRUE,
                            inits = 0)
saveRDS(model.null.Acc.dprime,file="nullmodel.Acc.dprime.split.rds")

# Set the priors for the models with slopes
prior = c(
  prior(normal(1.8, 1), class = Intercept), # based on Andersen & Mueller, 2011
  prior(normal(0, 2), class = b)) # a wide prior

# ExpPhase model
model.expphase.Acc.dprime = brm(dprime ~ ExpPhase + (ExpPhase|ParticipantNo),
                                data=data.final,
                                family=gaussian(),
                                prior = prior,
                                iter = 6000,
                                save_all_pars = TRUE,
                                control = list(adapt_delta = 0.99,max_treedepth = 15),
                                cores = 4,
                                sample_prior = TRUE,
                                inits = 0)
saveRDS(model.expphase.Acc.dprime,file="expphasemodel.Acc.dprime.split.rds")

#Interaction model
model.full.Acc.dprime = brm(dprime ~ ExpPhase * Condition + (ExpPhase + Condition|ParticipantNo),
                            data=data.final,
                            family=gaussian(),
                            prior = prior,
                            iter = 6000,
                            save_all_pars = TRUE,
                            control = list(adapt_delta = 0.99,max_treedepth = 15),
                            cores = 4,
                            sample_prior = TRUE,
                            inits = 0)
saveRDS(model.full.Acc.dprime,file="model.full.Acc.dprime.split.rds")

#WAIC
compare.Acc.dprime.waic = WAIC(model.null.Acc.dprime,model.expphase.Acc.dprime,model.full.Acc.dprime, compare = TRUE)
saveRDS(compare.Acc.dprime.waic,file="compare.Acc.dprime.split.waic")

# Weighted waic
compare.Acc.dprime.waic.weights = model_weights(model.null.Acc.dprime,model.expphase.Acc.dprime,model.full.Acc.dprime, weights = "waic")
saveRDS(compare.Acc.dprime.waic.weights,file="compare.Acc.dprime.split.waic.weights")

# Bayesian R2
#Null
bR2.null.Acc.dprime = bayes_R2(model.null.Acc.dprime)
saveRDS(bR2.null.Acc.dprime,file="bR2.null.Acc.dprime.split")
#ExpPhase
bR2.expphase.Acc.dprime = bayes_R2(model.expphase.Acc.dprime)
saveRDS(bR2.expphase.Acc.dprime,file="bR2.expphase.Acc.dprime.split")
#Full
bR2.full.Acc.dprime = bayes_R2(model.full.Acc.dprime)
saveRDS(bR2.full.Acc.dprime,file="bR2.full.Acc.dprime.split")
