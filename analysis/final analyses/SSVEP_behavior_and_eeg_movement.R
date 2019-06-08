
# Code info ###############################################################################################################################################################################################################


# Experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Gilles Pourtois, Ernst Koster, & Søren Andersen) (*: co-first authors)
# Code written by: Ivan Grahek & Antonio Schettino (2016-2018)
# Description: Code for the analysis of EEG data for Experiment 1 of the SSVEP - reward project. 


# Importing data & first steps ###############################################################################################################################################################################################################


# Clear environemnt and import data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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


# Average the movement and no movement trials
# data.raw = ddply(data.raw,.(Subject,Frequency),transform,
#                  BslnRedAttended = (BslnRedAttended + BslnRedAttendedMov)/2,
#                  BslnBlueAttended = (BslnBlueAttended + BslnBlueAttendedMov)/2,
#                  AcqRedAttended = (AcqRedAttended + AcqRedAttendedMov)/2,
#                  AcqBlueAttended = (AcqBlueAttended + AcqBlueAttendedMov)/2,
#                  ExtRedAttended = (ExtRedAttended + ExtRedAttendedMov)/2,
#                  ExtBlueAttended = (ExtBlueAttended + ExtBlueAttendedMov)/2)

# Select only no movement trials
# data.raw = dplyr::select(data.raw,"Subject","Frequency","BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended")

# Select only  movement trials
data.raw = dplyr::select(data.raw,"Subject","Frequency","BslnRedAttendedMov","BslnBlueAttendedMov","AcqRedAttendedMov","AcqBlueAttendedMov","ExtRedAttendedMov","ExtBlueAttendedMov")

colnames(data.raw) = c("Subject","Frequency","BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended")

# # Reshape to long format - take only no movement trials
data = melt(data.raw,id.vars=c("Subject","Frequency"),
            measure.vars=c("BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended"),
            variable.name="Condition",value.name="Amplitude")

# Reshape to long format - take only movement trials
# data = melt(data.raw,id.vars=c("Subject","Frequency"),
#             measure.vars=c("BslnRedAttendedMov","BslnBlueAttendedMov","AcqRedAttendedMov","AcqBlueAttendedMov","ExtRedAttendedMov","ExtBlueAttendedMov"),
#             variable.name="Condition",value.name="Amplitude")

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

# Switch the Frequency to the color
data$RecordedFrequency = ifelse(data$Frequency==10,"Blue","Red") # if the recorded frequency is 10Hz assign Blue (color flickering at 10Hz), otherwise assign Red (color flickering at 12Hz)

# Make a new condition based on the attended color and the rewarded color
data$Condition = ifelse(data$AttendedColor==data$RewardedColor, "High_Rew","Low_Rew")

# Make a new condition based on the attended color and the recorded frequency
data$Attention = ifelse(data$AttendedColor==data$RecordedFrequency, "Att","NotAtt")

# Make a new condition based the Condition and the Attention
data$RecordingAndCondition = with(data, paste0(Condition,"_",Attention))

# Select variables which we want to keep
data = subset(data, select=c("Subject","RewardedColor","ExpPhase","AttendedColor","Condition","RecordedFrequency","Attention","RecordingAndCondition","Amplitude"))

# Sort the data 
data = data[with(data, order(Subject)), ]

# Normalize the two frequencies
# Make a new variable with mean amplitude across all conditions for each participant and each frequency  !!! originally did not include ExpPhase below !!!
# data = ddply(data,.(Subject,RecordedFrequency),transform,
#                     MeanAmplitude = mean(Amplitude,na.rm=TRUE),
#                     SDAmplitude =   sd(Amplitude,na.rm=TRUE))

# Make a new variable with mean amplitude across all conditions for each participant and each frequency
data = ddply(data,.(Subject,RecordedFrequency),transform,
             MeanAmplitude = mean(Amplitude[ExpPhase=="Bsln"],na.rm=TRUE),
             SDAmplitude =   sd(Amplitude,na.rm=TRUE))

#MeanAmplitude = mean(Amplitude[ExpPhase=="Bsln"],na.rm=TRUE),   [ExpPhase=="Bsln"]

# Divide amplitudes in each Subject, Frequency, and Condition by the Mean Amplitude
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

# Plotting ###############################################################################################################################################################################################################

# prepare data for plotting
dataPlot = data

# rename variables
colnames(dataPlot)[colnames(dataPlot)=="ExpPhase"] <- "Reward phase"
colnames(dataPlot)[colnames(dataPlot)=="Condition"] <- "Reward probability"

# rename conditions
dataPlot$`Reward phase` = recode(dataPlot$`Reward phase`,
                                  "Acq" = "Acquisition",
                                  "Bsln" = "Baseline",
                                  "Ext" = "Extinction")

dataPlot$`Reward probability` = recode(dataPlot$`Reward probability`,
                                        "High_Rew" = "High",
                                        "Low_Rew" = "Low")

#order
dataPlot$`Reward phase` = factor(dataPlot$`Reward phase`, levels = c("Baseline","Acquisition","Extinction"))
dataPlot = dataPlot[order(dataPlot$Attention,dataPlot$`Reward phase`,dataPlot$`Reward probability`),]

plottingConditions = c("Attended","Unattended" )
for (i in 1:length(plottingConditions)){
  
  if(plottingConditions[i]=="Attended"){dataAmplitudePlot=subset(dataPlot,Attention=="Att")}
  
  if(plottingConditions[i]=="Unattended"){dataAmplitudePlot=subset(dataPlot,Attention=="NotAtt")}  

# Pirate plot

    pirateplot(formula = Amplitude ~ `Reward phase` + `Reward probability`, # dependent~independent variables
             data=dataAmplitudePlot, # data frame
             main=plottingConditions[i], # main title
             ylim=c(0.2,2.2), # y-axis: limits
             ylab=expression(paste("Amplitude (",mu,"V)")), # y-axis: label
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
}
  
pirateplot(formula = Selectivity ~ ExpPhase + Condition, # dependent~independent variables
           data=data.diff, # data frame
           main='Selectivity', # main title
           ylim=c(-0.5,1), # y-axis: limits
           ylab=expression(paste("Amplitude (",mu,"V)")), # y-axis: label
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


#### lme4 ####

library(lme4)
library(sjstats)
library(sjPlot)
library(jtools)
library(MASS)
# Set the intercept model
# data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
# data$Condition=relevel(data$Condition,ref="High_Rew")
# data$Attention=relevel(data$Attention,ref="Att")


data$ExpPhase = ordered(data$ExpPhase, levels = c("Bsln", "Acq", "Ext"))  # first contrast is facilitaion and the second one is interference
contrasts(data$ExpPhase) <- contr.sdif(3)
colnames(attr(data$ExpPhase, "contrasts")) =  c("Acq_min_Bsln", "Ext_min_Acq") # Change the name of the con

data$Attention = ordered(data$Attention, levels = c("NotAtt", "Att"))  # first contrast is facilitaion and the second one is interference
contrasts(data$Attention) <- contr.sdif(2)
colnames(attr(data$Attention, "contrasts")) =  c("Att_min_NotAtt") # Change the name of the con

data$Condition = ordered(data$Condition, levels = c("Low_Rew", "High_Rew"))  # first contrast is facilitaion and the second one is interference
contrasts(data$Condition) <- contr.sdif(2)
colnames(attr(data$Condition, "contrasts")) =  c("High_Rew_min_Low_Rew") # Change the name of the con


# Contrast coding
# data$Condition = ifelse(data$Condition == "High_Rew", 0.5, -0.5)
# data$Attention = ifelse(data$Attention == "Att", 0.5, -0.5)



# Null model
null = lmer(Amplitude ~ 1 + (1|Subject),
           data=data,REML = FALSE)

tab_model(null, show.re.var = T, show.icc = FALSE, show.ci = 0.95)


# Exp phase model
expphase = lmer(Amplitude ~ ExpPhase + (1|Subject),
            data=data,REML = FALSE)

tab_model(expphase, show.re.var = T, show.icc = FALSE, show.ci = 0.95)

# Attention model
attention = lmer(Amplitude ~ Attention + (1|Subject),
                data=data,REML = FALSE)

tab_model(attention, show.re.var = T, show.icc = FALSE, show.ci = 0.95)



# Two main effects - phase and attention
phaseANDattention = lmer(Amplitude ~ ExpPhase + Attention + (1|Subject),
                 data=data,REML = FALSE)

tab_model(phaseANDattention, show.re.var = T, show.icc = FALSE, show.ci = 0.95)

# Interaction between phase and attention
phaseANDattention_interaction = lmer(Amplitude ~ ExpPhase * Attention + (1|Subject),
                         data=data,REML = FALSE)

tab_model(phaseANDattention_interaction, show.re.var = T, show.icc = FALSE, show.ci = 0.95)

# Interaction between expphase and reward magnitude plus attention
rewardTimesPhasePlusAtt = lmer(Amplitude ~ Condition * ExpPhase + Attention + (1|Subject),
                                     data=data,REML = FALSE)

tab_model(rewardTimesPhasePlusAtt, show.re.var = T, show.icc = FALSE, show.ci = 0.95)

# Full model
full = lmer(Amplitude ~ Condition * ExpPhase * Attention + (1|Subject),
                               data=data,REML = FALSE)

tab_model(full, show.re.var = T, show.icc = FALSE, show.ci = 0.95)

cat_plot(full, pred = "ExpPhase", modx = "Condition",mod2 = "Attention",interval = T)


# 
# # Set the intercept model
# data.diff$ExpPhase=relevel(data.diff$ExpPhase,ref="Bsln")
# data.diff$Condition=relevel(data.diff$Condition,ref="High_Rew")
# 
# 


# # Selectivity
# phaseANDattention = brm(Selectivity ~ ExpPhase * Condition + (ExpPhase * Condition|Subject),
#                         data=data.diff,
#                         family=gaussian(),
#                         warmup = 2000,
#                         iter = 10000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99),
#                         cores = 4,
#                         sample_prior = TRUE)
# saveRDS(phaseANDattention,file="full.selectivity.EEG.allsubs.rds")


# brms three factors------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the working directory where to save the models
setwd(here("brms_models"))

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Modelling the effects of phase, attention, and reward magnitude - All subjects

# Set the intercept model
data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High_Rew")
data$Attention=relevel(data$Attention,ref="Att")

# Contrast coding
# data$Condition = ifelse(data$Condition == "High_Rew", 0.5, -0.5)
# data$Attention = ifelse(data$Attention == "Att", 0.5, -0.5)

# Null model
null = brm(Amplitude ~ 1 + (1|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(null,file="null.EEG.allsub.rds")

# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                                  data=data,
                                  family=gaussian(),
                                  warmup = 2000,
                                  iter = 10000,
                                  save_all_pars = TRUE,
                                  control = list(adapt_delta = 0.99),
                                  cores = 4,
                                  sample_prior = TRUE)
saveRDS(expphase,file="expphase.EEG.allsubs.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Attention|Subject),
                                   data=data,
                                   family=gaussian(),
                                   warmup = 2000,
                                   iter = 10000,
                                   save_all_pars = TRUE,
                                   control = list(adapt_delta = 0.99),
                                   cores = 4,
                                   sample_prior = TRUE)
saveRDS(attention,file="attention.EEG.allsubs.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + (ExpPhase + Attention|Subject),
                                        data=data,
                                        family=gaussian(),
                                        warmup = 2000,
                                        iter = 10000,
                                        save_all_pars = TRUE,
                                        control = list(adapt_delta = 0.99),
                                        cores = 4,
                                        sample_prior = TRUE)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs.rds")

# Interaction between phase and attention
phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + (ExpPhase * Attention|Subject),
                        data=data,
                        family=gaussian(),
                        warmup = 2000,
                        iter = 10000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99),
                        cores = 4,
                        sample_prior = TRUE)
saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.allsubs.rds")

# Interaction between expphase and reward magnitude plus attention
rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (Condition * ExpPhase + Attention|Subject),
                         data=data,
                         family=gaussian(),
                         warmup = 2000,
                         iter = 10000,
                         save_all_pars = TRUE,
                         control = list(adapt_delta = 0.99),
                         cores = 4,
                         sample_prior = TRUE)
saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.allsubs.rds")

# Full model
full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4,
                              sample_prior = TRUE)
saveRDS(full,file="full.EEG.allsubs.rds")



# Set the intercept model
data.diff$ExpPhase=relevel(data.diff$ExpPhase,ref="Bsln")
data.diff$Condition=relevel(data.diff$Condition,ref="High_Rew")




# # Selectivity
# phaseANDattention = brm(Selectivity ~ ExpPhase * Condition + (ExpPhase * Condition|Subject),
#                         data=data.diff,
#                         family=gaussian(),
#                         warmup = 2000,
#                         iter = 10000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99),
#                         cores = 4,
#                         sample_prior = TRUE)
# saveRDS(phaseANDattention,file="full.selectivity.EEG.allsubs.rds")

#Trial
#full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject/Trial),

# read in the models and comparisons
# null = readRDS("null.EEG.allsub.rds")
# attention = readRDS("attention.EEG.allsubs.rds")
# expphase = readRDS("expphase.EEG.allsubs.rds")
# phaseANDattention = readRDS("phaseANDattention.EEG.allsubs.rds")
# phaseANDattention_interaction = readRDS("phaseANDattention_interaction.EEG.allsubs.rds")
# rewardTimesPhasePlusAtt = readRDS("rewardTimesPhasePlusAtt.EEG.allsubs.rds")
# full = readRDS("full.EEG.allsubs.rds")
# waic = readRDS("compare.EEG.waic.allsubs.rds")

# color_scheme_set("viridis")
# pp_check(full, type = "stat_grouped", nsamples = 100, group = "ExpPhase")

# WAIC
compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs.rds")

# Weighted waic
# compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")

# Bayesian R2
#Null
bR2.null.EEG = bayes_R2(null)
saveRDS(bR2.null.EEG,file="bR2.null.EEG")
#ExpPhase
bR2.expphase.EEG = bayes_R2(expphase)
saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG")
#Attention
bR2.attention.EEG = bayes_R2(attention)
saveRDS(bR2.attention.EEG,file="bR2.attention.EEG")
#Phase and attention
bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG")
#Phase and attention interaction
bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG")
#Reward times phase plus attention
bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG")
#Full
bR2.full.EEG = bayes_R2(full)
saveRDS(bR2.full.EEG,file="bR2.full.EEG")

# null = readRDS("null.EEG.allsub.rds")
# attention = readRDS("attention.EEG.allsubs.rds")
# expphase = readRDS("expphase.EEG.allsubs.rds")
# phaseANDattention = readRDS("phaseANDattention.EEG.allsubs.rds")
# rewardTimesPhasePlusAtt = readRDS("rewardTimesPhasePlusAtt.EEG.allsubs.rds")
# full = readRDS("full.EEG.allsubs.rds")


# Analyzing the posterior and differences between conditions

post = posterior_samples(full, "^b")


################################################ Baseline ####

##################### Attended

######### High reward
Baseline_High_Attended = post[["b_Intercept"]]
######### Low reward
Baseline_Low_Attended = post[["b_Intercept"]] + 
  post[["b_ConditionLow_Rew"]] 

##################### Not Attended

######### High reward
Baseline_High_NotAttended = post[["b_Intercept"]] + 
  post[["b_AttentionNotAtt"]]
######### Low reward
Baseline_Low_NotAttended = post[["b_Intercept"]] + 
  post[["b_AttentionNotAtt"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ConditionLow_Rew:AttentionNotAtt"]]

################################################ Acquistion

##################### Attended

######### High reward
Acquisition_High_Attended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq"]] 
######### Low reward
Acquisition_Low_Attended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ConditionLow_Rew:ExpPhaseAcq"]]

##################### Not Attended

######### High reward
Acquisition_High_NotAttended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq"]] + 
  post[["b_AttentionNotAtt"]] + 
  post[["b_ExpPhaseAcq:AttentionNotAtt"]]
######### Low reward
Acquisition_Low_NotAttended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseAcq"]] + 
  post[["b_AttentionNotAtt"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseAcq:AttentionNotAtt"]] +
  post[["b_ConditionLow_Rew:ExpPhaseAcq"]] + 
  post[["b_ConditionLow_Rew:ExpPhaseAcq:AttentionNotAtt"]]

################################################ Extinction

##################### Attended

######### High reward
Extinction_High_Attended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt"]] 
######### Low reward
Extinction_Low_Attended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ConditionLow_Rew:ExpPhaseExt"]]

##################### Not Attended

######### High reward
Extinction_High_NotAttended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt"]] + 
  post[["b_AttentionNotAtt"]] + 
  post[["b_ExpPhaseExt:AttentionNotAtt"]]
######### Low reward
Extinction_Low_NotAttended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt"]] + 
  post[["b_AttentionNotAtt"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseExt:AttentionNotAtt"]] +
  post[["b_ConditionLow_Rew:ExpPhaseExt"]] + 
  post[["b_ConditionLow_Rew:ExpPhaseExt:AttentionNotAtt"]]



### Plotting the posterior ###

# make a data frame

posterior_conditions = melt(data.frame(Baseline_High_Attended, Baseline_High_NotAttended, Baseline_Low_Attended, Baseline_Low_NotAttended, Acquisition_High_Attended, Acquisition_High_NotAttended, Acquisition_Low_Attended, Acquisition_Low_NotAttended, Extinction_High_Attended, Extinction_High_NotAttended, Extinction_Low_Attended, Extinction_Low_NotAttended))

posterior_conditions =  posterior_conditions %>% separate(variable, c("Reward Phase", "Reward Probability", "Attention"), "_", extra = "merge")

posterior_conditions$Attention = recode(posterior_conditions$Attention,
                                        "Attended" = "Attended",
                                        "NotAttended" = "Unattended")

names(posterior_conditions)[4] = "Amplitude"


#order
#dataPlot$`Reward phase` = factor(dataPlot$`Reward phase`, levels = c("Baseline","Acquisition","Extinction"))
#dataPlot = dataPlot[order(dataPlot$Attention,dataPlot$`Reward phase`,dataPlot$`Reward probability`),]


plottingConditions = c("Attended","Unattended" )
for (i in 1:length(plottingConditions)){
  
  if(plottingConditions[i]=="Attended"){dataAmplitudePlot=subset(posterior_conditions,Attention=="Attended")}
  
  if(plottingConditions[i]=="Unattended"){dataAmplitudePlot=subset(posterior_conditions,Attention=="Unattended")}  
  
  # Pirate plot
  
  pirateplot(formula = Amplitude ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
             data=dataAmplitudePlot, # data frame
             main=plottingConditions[i], # main title
             ylim=c(0.7,1.2), # y-axis: limits
             ylab=expression(paste("Amplitude (",mu,"V)")), # y-axis: label
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
}

#Check the difference between high and low reward in acquisition attended

# Difference between high and low reward in acquisition attended
Diff_Rew_Acq_Att = Acquisition_High_Attended - Acquisition_Low_Attended
plotPost(Diff_Rew_Acq_Att, xlab = "", col = "#b3cde0", cex = 1, showCurve = FALSE, compVal = 0)

mean(Acquisition_High_Attended>Acquisition_Low_Attended)


#Check the difference between high and low reward in baseline attended

# Difference between high and low reward in acquisition attended
Diff_Rew_Bsln_Att = Baseline_High_Attended - Baseline_Low_Attended
plotPost(Diff_Rew_Bsln_Att, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)





