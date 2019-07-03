
# Code info ###############################################################################################################################################################################################################


# Experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Gilles Pourtois, Ernst Koster, & Søren Andersen) (*: co-first authors)
# Code written by: Ivan Grahek & Antonio Schettino (2016-2019)
# Description: Code for the analysis of EEG data for Experiment 1 of the SSVEP - reward project. 


# Importing data & first steps ###############################################################################################################################################################################################################


# Clear the environement and import data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# Calculate the attention indexes - Selectivity (attended-unattended) & total enhancement (attended+unattended) (Andersen & Muller, 2010, PNAS)
data.diff = ddply(data, .(Subject,ExpPhase,Condition), transform, Selectivity = Amplitude[Attention=="Att"]-Amplitude[Attention=="NotAtt"],TotalEnhancement=Amplitude[Attention=="Att"]+Amplitude[Attention=="NotAtt"])
# Delete the Attention column and rows which are not necessary (indexes repeated twice)
data.diff = subset(data.diff,Attention=="Att") #keep only Att as it is equal to NotAtt
data.diff$Attention = NULL  #drop the Attention column

# Sort the data 
data.diff$ExpPhase = factor(data.diff$ExpPhase, levels = c("Bsln","Acq","Ext"))
data.diff = data.diff[order(data.diff$Subject,data.diff$Condition,data.diff$ExpPhase),]

# Calculate the reward index - High reward minus Low reward
data.reward = ddply(data, .(Subject,ExpPhase,Attention), transform, Reward = Amplitude[Condition=="High_Rew"]-Amplitude[Condition=="Low_Rew"])
# Delete the Attention column and rows which are not necessary (indexes repeated twice)
data.reward = subset(data.reward,Condition=="High_Rew") #keep only Att as it is equal to NotAtt
data.reward$Condition = NULL  #drop the Condition column

# Sort the data 
data.reward$ExpPhase = factor(data.reward$ExpPhase, levels = c("Bsln","Acq","Ext"))
data.reward = data.reward[order(data.reward$Subject,data.reward$Attention,data.reward$ExpPhase),]

hist(subset(data.reward,data.reward$Attention=="Att" & data.reward$ExpPhase=="Acq" )$Reward)
hist(subset(data.reward,data.reward$Attention=="NotAtt")$Reward)

# Convert variables to be used in analyses into factors
data[c("Subject", "Condition","ExpPhase", "RewardedColor", "Attention", "RecordingAndCondition")] = 
  lapply(data[c("Subject", "Condition","ExpPhase", "RewardedColor", "Attention", "RecordingAndCondition")], factor)

data.diff[c("Subject", "Condition","ExpPhase", "RecordingAndCondition")] = 
  lapply(data.diff[c("Subject", "Condition","ExpPhase",  "RecordingAndCondition")], factor)

# Save the final data into a new variable
data_all = data

# All trials------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
  prior(normal(1, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)

# Null model
null = brm(Amplitude ~ 1 + (1|Subject),
                              data=data,
                              family=gaussian(),
                              prior = prior,
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(null,file="null.EEG.alltrials.rds")

# Priors for the models with slopes
prior = c(
  prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(normal(0, 3), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)

# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                                  data=data,
                                  family=gaussian(),
                                  prior = prior,
                                  warmup = 2000,
                                  iter = 10000,
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
                                   iter = 10000,
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
                                        iter = 10000,
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
                        iter = 10000,
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
                         iter = 10000,
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
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(full,file="full.EEG.alltrials.rds")


# WAIC
compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.alltrials.rds")

# Weighted waic
# compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")

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



# Only movement trials------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Take only the movement trials
data = subset(data_all, Movement=="Mov")

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
  prior(normal(1, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)

# Null model
null = brm(Amplitude ~ 1 + (1|Subject),
           data=data,
           family=gaussian(),
           prior = prior,
           warmup = 2000,
           iter = 10000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4,
           sample_prior = TRUE)
saveRDS(null,file="null.EEG.movementtrials.rds")

# Priors for the models with slopes
prior = c(
  prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(normal(0, 3), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)

# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
               data=data,
               family=gaussian(),
               prior = prior,
               warmup = 2000,
               iter = 10000,
               save_all_pars = TRUE,
               control = list(adapt_delta = 0.99),
               cores = 4,
               sample_prior = TRUE)
saveRDS(expphase,file="expphase.EEG.movementtrials.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Attention|Subject),
                data=data,
                family=gaussian(),
                prior = prior,  
                warmup = 2000,
                iter = 10000,
                save_all_pars = TRUE,
                control = list(adapt_delta = 0.99),
                cores = 4,
                sample_prior = TRUE)
saveRDS(attention,file="attention.EEG.movementtrials.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + (ExpPhase + Attention|Subject),
                        data=data,
                        family=gaussian(),
                        prior = prior,
                        warmup = 2000,
                        iter = 10000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99),
                        cores = 4,
                        sample_prior = TRUE)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.movementtrials.rds")

# Interaction between phase and attention
phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + (ExpPhase + Attention|Subject),
                                    data=data,
                                    family=gaussian(),
                                    prior = prior,
                                    warmup = 2000,
                                    iter = 10000,
                                    save_all_pars = TRUE,
                                    control = list(adapt_delta = 0.99),
                                    cores = 4,
                                    sample_prior = TRUE)
saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.movementtrials.rds")

# Interaction between expphase and reward magnitude plus attention
rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (ExpPhase + Attention + Condition|Subject),
                              data=data,
                              family=gaussian(),
                              prior = prior,
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.movementtrials.rds")

# Full model
full = brm(Amplitude ~ Condition * ExpPhase * Attention + (ExpPhase + Attention + Condition|Subject),
           data=data,
           family=gaussian(),
           prior = prior,
           warmup = 2000,
           iter = 10000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4,
           sample_prior = TRUE)
saveRDS(full,file="full.EEG.movementtrials.rds")


# WAIC
compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.movementtrials.rds")

# Weighted waic
# compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")

# Bayesian R2
#Null
bR2.null.EEG = bayes_R2(null)
saveRDS(bR2.null.EEG,file="bR2.null.EEG.movementtrials")
#ExpPhase
bR2.expphase.EEG = bayes_R2(expphase)
saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG.movementtrials")
#Attention
bR2.attention.EEG = bayes_R2(attention)
saveRDS(bR2.attention.EEG,file="bR2.attention.EEG.movementtrials")
#Phase and attention
bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG.movementtrials")
#Phase and attention interaction
bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG.movementtrials")
#Reward times phase plus attention
bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG.movementtrials")
#Full
bR2.full.EEG = bayes_R2(full)
saveRDS(bR2.full.EEG,file="bR2.full.EEG.movementtrials")



# Only no movement trials------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Take only the movement trials
data = subset(data_all, Movement=="Nomov")

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
  prior(normal(1, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)

# Null model
null = brm(Amplitude ~ 1 + (1|Subject),
           data=data,
           family=gaussian(),
           prior = prior,
           warmup = 2000,
           iter = 10000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4,
           sample_prior = TRUE)
saveRDS(null,file="null.EEG.nomovementtrials.rds")

# Priors for the models with slopes
prior = c(
  prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(normal(0, 3), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)

# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
               data=data,
               family=gaussian(),
               prior = prior,
               warmup = 2000,
               iter = 10000,
               save_all_pars = TRUE,
               control = list(adapt_delta = 0.99),
               cores = 4,
               sample_prior = TRUE)
saveRDS(expphase,file="expphase.EEG.nomovementtrials.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Attention|Subject),
                data=data,
                family=gaussian(),
                prior = prior,  
                warmup = 2000,
                iter = 10000,
                save_all_pars = TRUE,
                control = list(adapt_delta = 0.99),
                cores = 4,
                sample_prior = TRUE)
saveRDS(attention,file="attention.EEG.nomovementtrials.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + (ExpPhase + Attention|Subject),
                        data=data,
                        family=gaussian(),
                        prior = prior,
                        warmup = 2000,
                        iter = 10000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99),
                        cores = 4,
                        sample_prior = TRUE)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.nomovementtrials.rds")

# Interaction between phase and attention
phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + (ExpPhase + Attention|Subject),
                                    data=data,
                                    family=gaussian(),
                                    prior = prior,
                                    warmup = 2000,
                                    iter = 10000,
                                    save_all_pars = TRUE,
                                    control = list(adapt_delta = 0.99),
                                    cores = 4,
                                    sample_prior = TRUE)
saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.nomovementtrials.rds")

# Interaction between expphase and reward magnitude plus attention
rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (ExpPhase + Attention + Condition|Subject),
                              data=data,
                              family=gaussian(),
                              prior = prior,
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.nomovementtrials.rds")

# Full model
full = brm(Amplitude ~ Condition * ExpPhase * Attention + (ExpPhase + Attention + Condition|Subject),
           data=data,
           family=gaussian(),
           prior = prior,
           warmup = 2000,
           iter = 10000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4,
           sample_prior = TRUE)
saveRDS(full,file="full.EEG.nomovementtrials.rds")


# WAIC
compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.nomovementtrials.rds")

# Weighted waic
# compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")

# Bayesian R2
#Null
bR2.null.EEG = bayes_R2(null)
saveRDS(bR2.null.EEG,file="bR2.null.EEG.nomovementtrials")
#ExpPhase
bR2.expphase.EEG = bayes_R2(expphase)
saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG.nomovementtrials")
#Attention
bR2.attention.EEG = bayes_R2(attention)
saveRDS(bR2.attention.EEG,file="bR2.attention.EEG.nomovementtrials")
#Phase and attention
bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG.nomovementtrials")
#Phase and attention interaction
bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG.nomovementtrials")
#Reward times phase plus attention
bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG.nomovementtrials")
#Full
bR2.full.EEG = bayes_R2(full)
saveRDS(bR2.full.EEG,file="bR2.full.EEG.nomovementtrials")