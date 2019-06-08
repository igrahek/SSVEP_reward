################################################################## Code info ###############################################################################################################################################################################################################


# Experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Gilles Pourtois, Ernst Koster, & Søren Andersen) (*: co-first authors)
# Code written by: Ivan Grahek & Antonio Schettino (2016-2018)
# Description: Code for the analysis of EEG data for Experiment 1 of the SSVEP - reward project. 


################################################################## Importing data & first steps ###############################################################################################################################################################################################################


# Clear environemnt and import data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# clear the environment
rm(list=ls()) 
# clear the console
cat("\014") 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(reshape2,yarrr,BayesFactor,plyr,ez,schoRsch,brms,rstan, lme4, BEST, brmstools, here, tidyverse, Hmisc, car,sjstats,jtools)
# set seed
set.seed(42) 
#set.seed(32) 

# Set working directory
setwd(here())
# import data
data.raw = read.csv(file = here("data","singleTrial_amplitudes_movement_and_nomovement.csv"),header=TRUE,na.strings="NaN")
# data.raw = read.csv(file = here("data","singleTrial_amplitudes.csv"),header=TRUE,na.strings="NaN") 

# Prepare the dataset------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Reshape to long format
# data = melt(data.raw,id.vars=c("Subject","Frequency"),
#              measure.vars=c("BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended"),
#              variable.name="Condition",value.name="Amplitude") 

# Sort the new dataframe by participant name
# data = data[order(data$Subject),]

# Split the variable Condition based on capital letters
# data$Condition = gsub("(?!^)(?=[[:upper:]])", " ", data$Condition, perl=T)

# Split the variable condition into multiple variables
# Conditions = colsplit(data$Condition, pattern="\\s+",names = c('ExpPhase', 'ColorMoved',"attended","no","moved"))

data = data.raw

# Clean the subject name variable
data$participant = gsub('VP', '', data$participant)
data$participant = as.numeric(data$participant)

# Change the names of the variables
names(data)[names(data) == "participant"] = "Subject"
names(data)[names(data) == "amplitude"] = "Amplitude"
names(data)[names(data) == "frequency"] = "Frequency"


# Add new variables based on the condition
data$ExpPhase[data$condition == 1 | data$condition == 2 | data$condition == 11 | data$condition == 12]="Bsln"
data$ExpPhase[data$condition == 3 | data$condition == 4 | data$condition == 13 | data$condition == 14]="Acq"
data$ExpPhase[data$condition == 5 | data$condition == 6 | data$condition == 15 | data$condition == 16]="Ext"

data$AttendedColor[data$condition == 1 | data$condition == 3 | data$condition == 5 | data$condition == 11 |data$condition == 13 | data$condition == 15]="Red"
data$AttendedColor[data$condition == 2 | data$condition == 4 | data$condition == 6 | data$condition == 12 |data$condition == 14 | data$condition == 16]="Blue"

data$Movement[data$condition == 1 | data$condition == 2 | data$condition == 3 | data$condition == 4 |data$condition == 5 | data$condition == 6]="NoMovement"
data$Movement[data$condition == 11 | data$condition == 12 | data$condition == 13 | data$condition == 14 |data$condition == 15 | data$condition == 16]="Movement"

data$ExpPhase[data$condition == 3 | data$condition == 4 | data$condition == 13 | data$condition == 14]="Acq"
data$ExpPhase[data$condition == 5 | data$condition == 6 | data$condition == 15 | data$condition == 16]="Ext"

# Add the variable defining which color is rewarded based on the participant number
data$RewardedColor = ifelse(data$Subject%%2==0,"Blue","Red") # if participant number is even, blue was rewarded

# Switch the Frequency to the color
data$RecordedFrequency = ifelse(data$Frequency==10,"Blue","Red") # if the recorded frequency is 10Hz assign Blue (color flickering at 10Hz), otherwise assign Red (color flickering at 12Hz)

# Make a new condition based on the attended color and the rewarded color
data$Condition = ifelse(data$AttendedColor==data$RewardedColor, "High_Rew","Low_Rew")

# Make a new condition based on the attended color and the recorded frequency
data$Attention = ifelse(data$AttendedColor==data$RecordedFrequency, "Att","NotAtt")

# Make a new condition based the Condition and the Attention
data$RecordingAndCondition = with(data, paste0(Condition,"_",Attention))

# Select variables which we want to keep
data = subset(data, select=c("Subject","RewardedColor","ExpPhase","AttendedColor","Condition","RecordedFrequency","Attention","RecordingAndCondition","Amplitude","Movement"))

# Sort the data 
data = data[with(data, order(Subject)), ]

# Normalize the two frequencies
#Make a new variable with mean amplitude across all conditions for each participant and each frequency  !!! originally did not include ExpPhase below !!!
# data = ddply(data,.(Subject,RecordedFrequency),transform,
#                     MeanAmplitude = mean(Amplitude[ExpPhase=="Bsln"],na.rm=TRUE),
#                     SDAmplitude =   sd(Amplitude,na.rm=TRUE))

data = ddply(data,.(Subject,RecordedFrequency),transform,
             MeanAmplitude = mean(Amplitude,na.rm=TRUE),
             SDAmplitude =   sd(Amplitude,na.rm=TRUE))

#MeanAmplitude = mean(Amplitude[ExpPhase=="Baseline"],na.rm=TRUE),   [ExpPhase=="Bsln"]

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


################################################################## Add Questionnaires ###############################################################################################################################################################################################################

# import data
questionnaires = read.csv2(file = here("data","Questionnaires.csv"),header=TRUE,na.strings="NaN") 

#### Exclude participants with no EEG data ####

missing_eeg = c(35,39,44)
questionnaires = questionnaires[!questionnaires$Subject %in% missing_eeg,] 

#### Calculate gender ####
gender = ddply(questionnaires,.(Subject),summarise,
               Gender = unique(Sex))
table(unlist(gender$Gender))


#### Calculate age ####
age = ddply(questionnaires,.(Subject),summarise,
            Age = unique(Age))
age$Age=as.numeric(age$Age)
median(age$Age)

#### Calculate BDI ####
bdi = ddply(questionnaires,.(Subject),summarise,
            BDI = sum(ItemdisplayBDI3.RESP,na.rm=TRUE))

#### Calculate BIS BAS ####

bisbas = subset(questionnaires, Procedure.Block. == "QBISBAS", select = c(Subject,ItemQ,ItemdisplayBISBAS.RESP,ItemsBISBAS))

# Recode reversed items (2 and 22)
bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 2] = recode(bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 2], "0=3; 1=2; 2=1; 3=0")
bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 22] = recode(bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 22], "0=3; 1=2; 2=1; 3=0")

# Put on the scale from 1-4 instead of 0-3
bisbas$ItemdisplayBISBAS.RESP = bisbas$ItemdisplayBISBAS.RESP + 1

bisbas = ddply(bisbas,.(Subject),summarise,
               BIS = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 2],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 8],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 13],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 16],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 19],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 22],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 24]),
               
               DRIVE = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 3],
                           ItemdisplayBISBAS.RESP[ItemsBISBAS == 9],
                           ItemdisplayBISBAS.RESP[ItemsBISBAS == 12],
                           ItemdisplayBISBAS.RESP[ItemsBISBAS == 21]),
               
               FUN = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 5],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 10],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 15],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 20]),
               
               REWARD = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 4],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 7],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 14],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 18],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 23]),
               
               BAS = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 3],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 9],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 12],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 21],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 5],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 10],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 15],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 20],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 4],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 7],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 14],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 18],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 23]))

questionnaires = merge(age,gender)
questionnaires = merge(questionnaires,bdi)
questionnaires = merge(questionnaires,bisbas)

# Add information that we know about the missing participant's data
newrow = rep(NA, ncol(questionnaires))
questionnaires = rbind(questionnaires,newrow )
questionnaires$Gender[40] = "male"
questionnaires$Subject[40] = 1

# Impute the missing data
questionnaires$BDI = impute(questionnaires$BDI)
questionnaires$BAS = impute(questionnaires$BAS)
questionnaires$REWARD = impute(questionnaires$REWARD)
questionnaires$Age = impute(questionnaires$Age)
questionnaires$DRIVE = impute(questionnaires$DRIVE)
questionnaires$FUN = impute(questionnaires$FUN)
questionnaires$BIS = impute(questionnaires$BIS)


questionnaires$BDI[40] = median(questionnaires$BDI, na.rm = TRUE)
questionnaires$BAS[40] = median(questionnaires$BAS, na.rm = TRUE)
questionnaires$REWARD[40] = median(questionnaires$REWARD, na.rm = TRUE)



#names(questionnaires)[names(questionnaires) == "Subject"] = "ParticipantNo"

# Merge with the behavioral data
data=merge(data, questionnaires, all.x=TRUE, sort=FALSE)



# Add the participant for which we only know the name
#data.final$Gender[data.final$ParticipantNo==1] = "male"

# Center the questionnaire data
data$BDI = scale(data$BDI, scale= FALSE, center = TRUE)
data$BAS = scale(data$BAS, scale= FALSE, center = TRUE)
data$REWARD = scale(data$REWARD, scale= FALSE, center = TRUE)


################################################################## Plotting ###############################################################################################################################################################################################################

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


# brms movement and no movement trials------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# # Set the working directory where to save the models
# setwd(here("brms_models"))
# 
# #help stan run faster
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# # Modelling the effects of phase, attention, and reward magnitude - All subjects
# 
# # Set the intercept model
# data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
# data$Condition=relevel(data$Condition,ref="High_Rew")
# data$Attention=relevel(data$Attention,ref="Att")
# 
# # Contrast coding
# # data$Condition = ifelse(data$Condition == "High_Rew", 0.5, -0.5)
# # data$Attention = ifelse(data$Attention == "Att", 0.5, -0.5)
# 
# # # Prior for the intercept only model
# # prior = c(
# #   prior(normal(1, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
# #   prior(student_t(3, 0,3), class = sd),
# #   prior(student_t(3, 0,3), class = sigma)
# # )
# # 
# # # Null model
# # null = brm(Amplitude ~ 1 + (1|Subject),
# #            data=data,
# #            family=gaussian(),
# #            prior = prior,
# #            iter = 6000,
# #            save_all_pars = TRUE,
# #            control = list(adapt_delta = 0.99,max_treedepth = 15),
# #            cores = 4,
# #            sample_prior = TRUE,
# #            inits = 0)
# # saveRDS(null,file="null.EEG.allsub.rds")
# 
# # Priors for the models with slopes
# prior = c(
#   prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
#   prior(normal(0, 3), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
#   prior(student_t(3, 0,3), class = sd),
#   prior(student_t(3, 0,3), class = sigma)
#   )
# 
# # # Exp phase model
# # expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
# #                data=data,
# #                family=gaussian(),
# #                prior = prior,
# #                iter = 6000,
# #                save_all_pars = TRUE,
# #                control = list(adapt_delta = 0.99,max_treedepth = 15),
# #                cores = 4,
# #                sample_prior = TRUE,
# #                inits = 0)
# # saveRDS(expphase,file="expphase.EEG.allsubs.rds")
# # 
# # # Attention model
# # attention = brm(Amplitude ~ Attention + (Attention|Subject),
# #                 data=data,
# #                 family=gaussian(),
# #                 prior = prior,
# #                 iter = 6000,
# #                 save_all_pars = TRUE,
# #                 control = list(adapt_delta = 0.99,max_treedepth = 15),
# #                 cores = 4,
# #                 sample_prior = TRUE,
# #                 inits = 0)
# # saveRDS(attention,file="attention.EEG.allsubs.rds")
# # 
# # # Two main effects - phase and attention
# # phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + ( ExpPhase + Attention|Subject),
# #                         data=data,
# #                         family=gaussian(),
# #                         prior = prior,
# #                         iter = 6000,
# #                         save_all_pars = TRUE,
# #                         control = list(adapt_delta = 0.99,max_treedepth = 15),
# #                         cores = 4,
# #                         sample_prior = TRUE,
# #                         inits = 0)
# # saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs.rds")
# # 
# # # Interaction between phase and attention
# # phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + ( ExpPhase * Attention|Subject),
# #                         data=data,
# #                         family=gaussian(),
# #                         prior = prior,
# #                         iter = 6000,
# #                         save_all_pars = TRUE,
# #                         control = list(adapt_delta = 0.99,max_treedepth = 15),
# #                         cores = 4,
# #                         sample_prior = TRUE,
# #                         inits = 0)
# # saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.allsubs.rds")
# # 
# # # Interaction between expphase and reward magnitude plus attention
# # rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (Condition * ExpPhase + Attention|Subject),
# #                          data=data,
# #                          family=gaussian(),
# #                          prior = prior,
# #                          iter = 6000,
# #                          save_all_pars = TRUE,
# #                          control = list(adapt_delta = 0.99,max_treedepth = 15),
# #                          cores = 4,
# #                          sample_prior = TRUE,
# #                          inits = 0)
# # saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.allsubs.rds")
# 
# # Full model
# 
# 
# full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),   
#                              data=data,
#                              family=gaussian(),
#                              prior = prior,
#                              iter = 6000,
#                              save_all_pars = TRUE,
#                              control = list(adapt_delta = 0.99),
#                              cores = 4,
#                              sample_prior = TRUE,
#                              inits = 0)
# saveRDS(full,file="full.EEG.allsubs_old_normalization.rds")

# library(lme4)
# library(sjPlot)
# m1=lmer(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),data = data)
# tab_model(m1)
# cat_plot(m1, pred = "ExpPhase",modx = "Condition",mod2 = "Attention",interval = T)
# 
# m2=lmer(Amplitude ~ Condition * ExpPhase * Attention + (Condition + ExpPhase + Attention|Subject),data = subset(data,data$Movement=="Movement"))
# tab_model(m2)
# cat_plot(m2, pred = "ExpPhase",modx = "Condition",mod2 = "Attention",interval = T)
# 
# m3=lmer(Amplitude ~ Condition * ExpPhase * Attention + (Condition + ExpPhase + Attention|Subject),data = subset(data,data$Movement=="NoMovement"))
# tab_model(m3)
# cat_plot(m3, pred = "ExpPhase",modx = "Condition",mod2 = "Attention",interval = T)

# full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),   
#            data=data,
#            prior = prior,
#            family=gaussian(),
#            warmup = 2000,
#            iter = 10000,
#            save_all_pars = TRUE,
#            control = list(adapt_delta = 0.99),
#            cores = 4,
#            sample_prior = TRUE)
# saveRDS(full,file="full.EEG.allsubs.rds")
# 
# # Depression model
# depression = brm(Amplitude ~ Condition * ExpPhase * Attention * BDI + (Condition * ExpPhase * Attention * BDI|Subject),   
#            data=data,
#            family=gaussian(),
#            warmup = 2000,
#            iter = 10000,
#            save_all_pars = TRUE,
#            control = list(adapt_delta = 0.99),
#            cores = 4,
#            sample_prior = TRUE)
# saveRDS(depression,file="full.EEG.depression.allsubs.rds")
# 
# # BAS model
# bas = brm(Amplitude ~ Condition * ExpPhase * Attention * BAS + (Condition * ExpPhase * Attention * BAS|Subject),   
#            data=data,
#            family=gaussian(),
#            warmup = 2000,
#            iter = 10000,
#            save_all_pars = TRUE,
#            control = list(adapt_delta = 0.99),
#            cores = 4,
#            sample_prior = TRUE)
# saveRDS(bas,file="full.EEG.bas.allsubs.rds")
# 
# # REWARD model
# reward = brm(Amplitude ~ Condition * ExpPhase * Attention * REWARD + (Condition * ExpPhase * Attention * REWARD|Subject),   
#            data=data,
#            family=gaussian(),
#            warmup = 2000,
#            iter = 10000,
#            save_all_pars = TRUE,
#            control = list(adapt_delta = 0.99),
#            cores = 4,
#            sample_prior = TRUE)
# saveRDS(reward,file="full.EEG.reward.allsubs.rds")



# # Set the intercept model
# data.diff$ExpPhase=relevel(data.diff$ExpPhase,ref="Bsln")
# data.diff$Condition=relevel(data.diff$Condition,ref="High_Rew")
# 
# 
# 
# 
# #Selectivity
# 
# # Priors for the models with slopes
# prior = c(
#   prior(normal(0, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
#   prior(normal(0, 2), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
#   prior(student_t(3, 0,3), class = sd),
#   prior(student_t(3, 0,3), class = sigma)
# )
# 
# # # apply contrasts
# contrasts(data.diff$Condition) <- contr.sdif(2)
# contrasts(data.diff$ExpPhase) <- contr.sdif(3)
# 
# 
# # fit the model
# phaseANDattention = brm(Selectivity ~ ExpPhase * Condition + (1|Subject),
#                         data=data.diff,
#                         family=gaussian(),
#                         prior = prior,
#                         iter = 6000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99,max_treedepth = 15),
#                         cores = 4,
#                         sample_prior = TRUE,
#                         inits = 0)
# saveRDS(phaseANDattention,file="full.selectivity.EEG.allsubs.rds")
# 
# 
# marginal_effects(phaseANDattention)
# summary(phaseANDattention)

# hyp = hypothesis(phaseANDattention, "ExpPhase2M1:Condition2M1 > 0")
# 
# 
# cat_plot(phaseANDattention, pred = "ExpPhase", modx = "Condition",interval = T)
# 



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
# compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
# saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs.rds")

# # Questionnaires WAIC
# compare.EEG.waic.questionnaires = WAIC(full, depression, bas, reward, compare = TRUE)
# saveRDS(compare.EEG.waic.questionnaires,file="compare.EEG.waic.allsubs.questionnaires.rds")

# Weighted waic
# compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")

# Bayesian R2
#Null
# bR2.null.EEG = bayes_R2(null)
# saveRDS(bR2.null.EEG,file="bR2.null.EEG")
# #ExpPhase
# bR2.expphase.EEG = bayes_R2(expphase)
# saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG")
# #Attention
# bR2.attention.EEG = bayes_R2(attention)
# saveRDS(bR2.attention.EEG,file="bR2.attention.EEG")
# #Phase and attention
# bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
# saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG")
# #Phase and attention interaction
# bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
# saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG")
# #Reward times phase plus attention
# bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
# saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG")
# #Full
# bR2.full.EEG = bayes_R2(full)
# saveRDS(bR2.full.EEG,file="bR2.full.EEG")
# #Depression
# bR2.full.EEG.depression = bayes_R2(depression)
# saveRDS(bR2.full.EEG.depression,file="bR2.full.EEG.depression")
# #BAS
# bR2.full.EEG.BAS = bayes_R2(bas)
# saveRDS(bR2.full.EEG.BAS,file="bR2.full.EEG.BAS")
# #REWARD
# bR2.full.EEG.REWARD = bayes_R2(reward)
# saveRDS(bR2.full.EEG.REWARD,file="bR2.full.EEG.REWARD")

# null = readRDS("null.EEG.allsub.rds")
# attention = readRDS("attention.EEG.allsubs.rds")
# expphase = readRDS("expphase.EEG.allsubs.rds")
# phaseANDattention = readRDS("phaseANDattention.EEG.allsubs.rds")
# rewardTimesPhasePlusAtt = readRDS("rewardTimesPhasePlusAtt.EEG.allsubs.rds")
# full = readRDS("full.EEG.allsubs.rds")









# brms movement trials------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# # Prior for the intercept only model
# prior = c(
#   prior(normal(1, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
#   prior(student_t(3, 0,3), class = sd),
#   prior(student_t(3, 0,3), class = sigma)
# )
# 
# # Null model
# null = brm(Amplitude ~ 1 + (1|Subject),
#            data=subset(data,Movement=="Movement"),
#            family=gaussian(),
#            prior = prior,
#            iter = 6000,
#            save_all_pars = TRUE,
#            control = list(adapt_delta = 0.99,max_treedepth = 15),
#            cores = 4,
#            sample_prior = TRUE,
#            inits = 0)
# saveRDS(null,file="null.EEG.allsub.movement.rds")
# 
# Priors for the models with slopes
prior = c(
  prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(normal(0, 3), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)
# 
# # Exp phase model
# expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
#                data=subset(data,Movement=="Movement"),
#                family=gaussian(),
#                prior = prior,
#                iter = 6000,
#                save_all_pars = TRUE,
#                control = list(adapt_delta = 0.99,max_treedepth = 15),
#                cores = 4,
#                sample_prior = TRUE,
#                inits = 0)
# saveRDS(expphase,file="expphase.EEG.allsubs.movement.rds")
# 
# # Attention model
# attention = brm(Amplitude ~ Attention + (Attention|Subject),
#                 data=subset(data,Movement=="Movement"),
#                 family=gaussian(),
#                 prior = prior,
#                 iter = 6000,
#                 save_all_pars = TRUE,
#                 control = list(adapt_delta = 0.99,max_treedepth = 15),
#                 cores = 4,
#                 sample_prior = TRUE,
#                 inits = 0)
# saveRDS(attention,file="attention.EEG.allsubs.movement.rds")
# 
# # Two main effects - phase and attention
# phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + ( ExpPhase + Attention|Subject),
#                         data=subset(data,Movement=="Movement"),
#                         family=gaussian(),
#                         prior = prior,
#                         iter = 6000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99,max_treedepth = 15),
#                         cores = 4,
#                         sample_prior = TRUE,
#                         inits = 0)
# saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs.movement.rds")
# 
# # Interaction between phase and attention
# phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + ( ExpPhase * Attention|Subject),
#                                     data=subset(data,Movement=="Movement"),
#                                     family=gaussian(),
#                                     prior = prior,
#                                     iter = 6000,
#                                     save_all_pars = TRUE,
#                                     control = list(adapt_delta = 0.99,max_treedepth = 15),
#                                     cores = 4,
#                                     sample_prior = TRUE,
#                                     inits = 0)
# saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.allsubs.movement.rds")
# 
# # Interaction between expphase and reward magnitude plus attention
# rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (Condition * ExpPhase + Attention|Subject),
#                               data=subset(data,Movement=="Movement"),
#                               family=gaussian(),
#                               prior = prior,
#                               iter = 6000,
#                               save_all_pars = TRUE,
#                               control = list(adapt_delta = 0.99,max_treedepth = 15),
#                               cores = 4,
#                               sample_prior = TRUE,
#                               inits = 0)
# saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.allsubs.movement.rds")

# Full model


full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),   
           data=subset(data,Movement=="Movement"),
           family=gaussian(),
           prior = prior,
           iter = 6000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4,
           sample_prior = TRUE,
           inits = 0)
saveRDS(full,file="full.EEG.allsubs.movement_old_normalization.rds")

# WAIC
# compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
# saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs.movement.rds")

# # Questionnaires WAIC
# compare.EEG.waic.questionnaires = WAIC(full, depression, bas, reward, compare = TRUE)
# saveRDS(compare.EEG.waic.questionnaires,file="compare.EEG.waic.allsubs.questionnaires.rds")

# Weighted waic
# compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")

# Bayesian R2
# #Null
# bR2.null.EEG = bayes_R2(null)
# saveRDS(bR2.null.EEG,file="bR2.null.EEG.movement")
# #ExpPhase
# bR2.expphase.EEG = bayes_R2(expphase)
# saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG.movement")
# #Attention
# bR2.attention.EEG = bayes_R2(attention)
# saveRDS(bR2.attention.EEG,file="bR2.attention.EEG.movement")
# #Phase and attention
# bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
# saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG.movement")
# #Phase and attention interaction
# bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
# saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG.movement")
# #Reward times phase plus attention
# bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
# saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG.movement")
# #Full
# bR2.full.EEG = bayes_R2(full)
# saveRDS(bR2.full.EEG,file="bR2.full.EEG.movement")
# brms No movement trials------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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

# # Prior for the intercept only model
# prior = c(
#   prior(normal(1, 2), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
#   prior(student_t(3, 0,3), class = sd),
#   prior(student_t(3, 0,3), class = sigma)
# )
# 
# # Null model
# null = brm(Amplitude ~ 1 + (1|Subject),
#            data=subset(data,Movement=="NoMovement"),
#            family=gaussian(),
#            prior = prior,
#            iter = 6000,
#            save_all_pars = TRUE,
#            control = list(adapt_delta = 0.99,max_treedepth = 15),
#            cores = 4,
#            sample_prior = TRUE,
#            inits = 0)
# saveRDS(null,file="null.EEG.allsub.nomovement.rds")
# 
# Priors for the models with slopes
prior = c(
  prior(normal(1, 3), class = Intercept), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(normal(0, 3), class = b), #based on Anderson & Mueller, 2010 and on logic - max attention effect (removal of the unattended stimulus) == 2
  prior(student_t(3, 0,3), class = sd),
  prior(student_t(3, 0,3), class = sigma)
)
# 
# # Exp phase model
# expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
#                data=subset(data,Movement=="NoMovement"),
#                family=gaussian(),
#                prior = prior,
#                iter = 6000,
#                save_all_pars = TRUE,
#                control = list(adapt_delta = 0.99,max_treedepth = 15),
#                cores = 4,
#                sample_prior = TRUE,
#                inits = 0)
# saveRDS(expphase,file="expphase.EEG.allsubs.nomovement.rds")
# 
# # Attention model
# attention = brm(Amplitude ~ Attention + (Attention|Subject),
#                 data=subset(data,Movement=="NoMovement"),
#                 family=gaussian(),
#                 prior = prior,
#                 iter = 6000,
#                 save_all_pars = TRUE,
#                 control = list(adapt_delta = 0.99,max_treedepth = 15),
#                 cores = 4,
#                 sample_prior = TRUE,
#                 inits = 0)
# saveRDS(attention,file="attention.EEG.allsubs.nomovement.rds")
# 
# # Two main effects - phase and attention
# phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + ( ExpPhase + Attention|Subject),
#                         data=subset(data,Movement=="NoMovement"),
#                         family=gaussian(),
#                         prior = prior,
#                         iter = 6000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99,max_treedepth = 15),
#                         cores = 4,
#                         sample_prior = TRUE,
#                         inits = 0)
# saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs.nomovement.rds")
# 
# # Interaction between phase and attention
# phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + ( ExpPhase * Attention|Subject),
#                                     data=subset(data,Movement=="NoMovement"),
#                                     family=gaussian(),
#                                     prior = prior,
#                                     iter = 6000,
#                                     save_all_pars = TRUE,
#                                     control = list(adapt_delta = 0.99,max_treedepth = 15),
#                                     cores = 4,
#                                     sample_prior = TRUE,
#                                     inits = 0)
# saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.allsubs.nomovement.rds")
# 
# # Interaction between expphase and reward magnitude plus attention
# rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (Condition * ExpPhase + Attention|Subject),
#                               data=subset(data,Movement=="NoMovement"),
#                               family=gaussian(),
#                               prior = prior,
#                               iter = 6000,
#                               save_all_pars = TRUE,
#                               control = list(adapt_delta = 0.99,max_treedepth = 15),
#                               cores = 4,
#                               sample_prior = TRUE,
#                               inits = 0)
# saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.allsubs.nomovement.rds")

# Full model


full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),   
           data=subset(data,Movement=="NoMovement"),
           family=gaussian(),
           prior = prior,
           iter = 6000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4,
           sample_prior = TRUE,
           inits = 0)
saveRDS(full,file="full.EEG.allsubs.nomovement_old_normalization.rds")

# # WAIC
# compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
# saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs.nomovement.rds")
# 
# # # Questionnaires WAIC
# # compare.EEG.waic.questionnaires = WAIC(full, depression, bas, reward, compare = TRUE)
# # saveRDS(compare.EEG.waic.questionnaires,file="compare.EEG.waic.allsubs.questionnaires.rds")
# 
# # Weighted waic
# # compare.EEG.waic.weights = model_weights(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, weights = "waic")
# # saveRDS(compare.EEG.waic.weights,file="compare.EEG.waic.weights")
# 
# # Bayesian R2
# #Null
# bR2.null.EEG = bayes_R2(null)
# saveRDS(bR2.null.EEG,file="bR2.null.EEG.nomovement")
# #ExpPhase
# bR2.expphase.EEG = bayes_R2(expphase)
# saveRDS(bR2.expphase.EEG,file="bR2.expphase.EEG.nomovement")
# #Attention
# bR2.attention.EEG = bayes_R2(attention)
# saveRDS(bR2.attention.EEG,file="bR2.attention.EEG.nomovement")
# #Phase and attention
# bR2.phaseANDattention.EEG = bayes_R2(phaseANDattention)
# saveRDS(bR2.phaseANDattention.EEG,file="bR2.phaseANDattention.EEG.nomovement")
# #Phase and attention interaction
# bR2.phaseANDattention_interaction.EEG = bayes_R2(phaseANDattention_interaction)
# saveRDS(bR2.phaseANDattention_interaction.EEG,file="bR2.phaseANDattention_interaction.EEG.nomovement")
# #Reward times phase plus attention
# bR2.rewardTimesPhasePlusAtt.EEG = bayes_R2(rewardTimesPhasePlusAtt)
# saveRDS(bR2.rewardTimesPhasePlusAtt.EEG,file="bR2.rewardTimesPhasePlusAtt.EEG.nomovement")
# #Full
# bR2.full.EEG = bayes_R2(full)
# saveRDS(bR2.full.EEG,file="bR2.full.EEG.nomovement")




















### Behavior ####

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
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2,brms,rstan, tidyverse, here,car, Hmisc,psych)
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
# data.final = ddply(data.raw,.(ParticipantNo,ExpPhase,AttendedColor,RewardedColor,MovedDots),summarise,
#                   numtrials=length(which(Response!=99)), # number of trials per condition (anything that is not 99 or any other number that we're not using)
#                   Hits=length(which(Response==1)), # hits: attended color moved, correct response
#                   FAs=length(which(Response==2)), # false alarms: attended color did not move, (wrong) response
#                   Misses=length(which(Response==0)), # misses: attended color moved, no response
#                   CRs=length(which(Response==3)), # correct rejections: attended color did not move, no response
#                   mean.RT=mean(RT,na.rm=TRUE)) # mean RT per condition


### Create variables needed for the trial-by-trial accuracy analyses
# count hits, false alarms, misses, correct rejections, and RT separately for each participant (their calculation is done in Matlab: see DataProcessing.m)

# First we take out the responses that are correct rejections
data.final = subset(data.raw, data.raw$Response!=3 & data.raw$Response!=0)

data.final$Response = ifelse(data.final$Response==1,"Hit","FA")
data.final$Response = ifelse(data.final$Response=="Hit",1,0)
data.final$Response = as.factor(data.final$Response)


# For the correct vs. incorrect analysis
data.final.corr = subset(data.raw, data.raw$Response!=3)
data.final.corr = subset(data.final.corr,MovedDots==AttendedColor)

data.final.corr$Response = ifelse(data.final.corr$Response==1,"Correct","Incorrect")
data.final.corr$Response = ifelse(data.final.corr$Response=="Correct",1,0)
data.final.corr$Response = as.factor(data.final.corr$Response)

# Then we turn the responses coded as hits (1) and false alarms ()

# data.final = ddply(data.raw,.(ParticipantNo,ExpPhase,AttendedColor,RewardedColor,MovedDots,Trial),summarise,
#                    Response = if()
#                    numtrials=length(which(Response!=99)), # number of trials per condition (anything that is not 99 or any other number that we're not using)
#                    Hits=length(which(Response==1)), # hits: attended color moved, correct response
#                    FAs=length(which(Response==2)), # false alarms: attended color did not move, (wrong) response
#                    Misses=length(which(Response==0)), # misses: attended color moved, no response
#                    CRs=length(which(Response==3)), # correct rejections: attended color did not move, no response
#                    mean.RT=mean(RT,na.rm=TRUE)) # mean RT per condition


################################################################## Calculate accuracy and RTs per condition ###############################################################################################################################################################################################################

# Prepare the data------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Calculate Hits and False alarms
# Hits are calculated for each participant in each condition on trials when they are attending the color that moved. 
# False alarms are  calculated for each participant in each condition on trials when they are attending the color that didn't move (the unattended color moved, but they responded)  
# Here we create the same number of hits & fas for each of the two conditions (moved attended or not)
# data.final = ddply(data.final, .(ParticipantNo,ExpPhase,AttendedColor), transform, 
#                  Hits = Hits[MovedDots==AttendedColor],
#                  FAs = FAs[MovedDots!=AttendedColor])

# Keep only trials on which the attended color moved (we can do behavioral analysis only on those)
# data.final = subset(data.final,MovedDots==AttendedColor)

### Calculate d'
# use loglinear transformation: add 0.5 to Hits, FAs, Misses, and CRs (Hautus, 1995, Behavior Research Methods, Instruments, & Computers),
# which is preferred over the 1/2N rule (Macmillan & Kaplan, 1985, Psychological Bulletin) because it results in less biased estimates of d'.
# data.final =  ddply(data.final,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,numtrials),summarise,
#                       tot.Hits=Hits+.5, # hits
#                       tot.FAs=FAs+.5, # false alarms
#                       tot.Misses=(numtrials-tot.Hits)+.5, # misses
#                       tot.CRs=(numtrials-tot.FAs)+.5, # correct rejections
#                       Hit.Rate=tot.Hits/(tot.Hits+tot.Misses), # hit rate
#                       FA.Rate=tot.FAs/(tot.FAs+tot.CRs), # false alarm rate
#                       dprime=qnorm(Hit.Rate)-qnorm(FA.Rate),
#                       Hits.RTs=mean(mean.RT,na.rm=TRUE)) # mean RTs
# 
# data.final =  ddply(data.final,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,numtrials),summarise,
#                     tot.Hits=Hits+.5, # hits
#                     tot.FAs=FAs+.5, # false alarms
#                     tot.Misses=(numtrials-Hits)+.5, # misses
#                     tot.CRs=(numtrials-FAs)+.5, # correct rejections
#                     #Hit.Rate=tot.Hits/(tot.Hits+tot.Misses), # hit rate
#                     #FA.Rate=tot.FAs/(tot.FAs+tot.CRs), # false alarm rate
#                     #dprime=qnorm(Hit.Rate)-qnorm(FA.Rate),
#                     Hits.RTs=mean(mean.RT,na.rm=TRUE)) # mean RTs

# Calculate SDT indices with psycho
# indices = psycho::dprime(data.final$tot.Hits, data.final$tot.FAs, data.final$tot.Misses, data.final$tot.CRs) 
# 
# data.final = cbind(data.final, indices)                      

# Handle outliers------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Determine outliers and cut them
# outliers based on hit rate at any condition
#crit = .6 # minimum 60% hit rate in any condition .6
# select participants below the criterion
#criterion = subset(ddply(data.final,.(ParticipantNo),summarise,mean.Hit.Rate=mean(Hit.Rate)),mean.Hit.Rate<crit)$Participant # minimum 60% hit rate across all conditions

#criterion = subset(data.final,data.final$dprime<0)$Participant # minimum 60% hit rate across all conditions

# eliminate ouotliers from data frame
#data.final = data.final[!data.final$ParticipantNo %in% unique(criterion),] 

# Create the final dataframe for accuracy and RTs ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### Create a final dataframe for accuracy and RTs analyses
# add a new variable specifying whether the participant is attending the high or Low_Rewed color
data.final$Condition = ifelse(data.final$RewardedColor==data.final$AttendedColor,"High_Rew","Low_Rew")
# make this variable a factor for further analyses
data.final$Condition = factor(data.final$Condition)

# add a new variable specifying whether the participant is attending the high or Low_Rewed color
data.final.corr$Condition = ifelse(data.final.corr$RewardedColor==data.final.corr$AttendedColor,"High_Rew","Low_Rew")
# make this variable a factor for further analyses
data.final.corr$Condition = factor(data.final.corr$Condition)


################################################################## Add Questionnaires ###############################################################################################################################################################################################################

# import data
questionnaires = read.csv2(file = here("data","Questionnaires.csv"),header=TRUE,na.strings="NaN") 

#### Exclude participants with no EEG data ####

missing_eeg = c(35,39,44)
questionnaires = questionnaires[!questionnaires$Subject %in% missing_eeg,] 

#### Calculate gender ####
gender = ddply(questionnaires,.(Subject),summarise,
               Gender = unique(Sex))
table(unlist(gender$Gender))


#### Calculate age ####
age = ddply(questionnaires,.(Subject),summarise,
            Age = unique(Age))
age$Age=as.numeric(age$Age)
median(age$Age)

#### Calculate BDI ####
bdi = ddply(questionnaires,.(Subject),summarise,
            BDI = sum(ItemdisplayBDI3.RESP,na.rm=TRUE))

#### Calculate BIS BAS ####

bisbas = subset(questionnaires, Procedure.Block. == "QBISBAS", select = c(Subject,ItemQ,ItemdisplayBISBAS.RESP,ItemsBISBAS))

# Recode reversed items (2 and 22)
bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 2] = recode(bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 2], "0=3; 1=2; 2=1; 3=0")
bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 22] = recode(bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 22], "0=3; 1=2; 2=1; 3=0")

# Put on the scale from 1-4 instead of 0-3
bisbas$ItemdisplayBISBAS.RESP = bisbas$ItemdisplayBISBAS.RESP + 1

bisbas = ddply(bisbas,.(Subject),summarise,
               BIS = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 2],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 8],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 13],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 16],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 19],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 22],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 24]),
               
               DRIVE = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 3],
                           ItemdisplayBISBAS.RESP[ItemsBISBAS == 9],
                           ItemdisplayBISBAS.RESP[ItemsBISBAS == 12],
                           ItemdisplayBISBAS.RESP[ItemsBISBAS == 21]),
               
               FUN = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 5],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 10],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 15],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 20]),
               
               REWARD = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 4],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 7],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 14],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 18],
                            ItemdisplayBISBAS.RESP[ItemsBISBAS == 23]),
               
               BAS = sum(ItemdisplayBISBAS.RESP[ItemsBISBAS == 3],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 9],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 12],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 21],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 5],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 10],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 15],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 20],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 4],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 7],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 14],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 18],
                         ItemdisplayBISBAS.RESP[ItemsBISBAS == 23]))

questionnaires = merge(age,gender)
questionnaires = merge(questionnaires,bdi)
questionnaires = merge(questionnaires,bisbas)

# Add information that we know about the missing participant's data
newrow = rep(NA, ncol(questionnaires))
questionnaires = rbind(questionnaires,newrow )
questionnaires$Gender[40] = "male"
questionnaires$Subject[40] = 1

# Impute the missing data
questionnaires$BDI = impute(questionnaires$BDI)
questionnaires$BAS = impute(questionnaires$BAS)
questionnaires$REWARD = impute(questionnaires$REWARD)
questionnaires$Age = impute(questionnaires$Age)
questionnaires$DRIVE = impute(questionnaires$DRIVE)
questionnaires$FUN = impute(questionnaires$FUN)
questionnaires$BIS = impute(questionnaires$BIS)


questionnaires$BDI[40] = median(questionnaires$BDI, na.rm = TRUE)
questionnaires$BAS[40] = median(questionnaires$BAS, na.rm = TRUE)
questionnaires$REWARD[40] = median(questionnaires$REWARD, na.rm = TRUE)



names(questionnaires)[names(questionnaires) == "Subject"] = "ParticipantNo"

# Merge with the behavioral data
data.final=merge(data.final, questionnaires, all.x=TRUE, sort=FALSE)



# Add the participant for which we only know the name
#data.final$Gender[data.final$ParticipantNo==1] = "male"

# Center the questionnaire data
data.final$BDI = scale(data.final$BDI, scale= FALSE, center = TRUE)
data.final$BAS = scale(data.final$BAS, scale= FALSE, center = TRUE)
data.final$REWARD = scale(data.final$REWARD, scale= FALSE, center = TRUE)

################################################################## Plotting ###############################################################################################################################################################################################################

# # Plot Hit rates------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# # Prepare the dataset
# data.plot = data.final
# 
# # rename variables
# colnames(data.plot)[colnames(data.plot)=="ExpPhase"] <- "Reward phase"
# colnames(data.plot)[colnames(data.plot)=="Condition"] <- "Reward probability"
# 
# # rename conditions
# data.plot$`Reward phase` = recode(data.plot$`Reward phase`,
#                                 "Acq" = "Acquisition",
#                                 "Bsln" = "Baseline",
#                                 "Ext" = "Extinction")
# 
# data.plot$`Reward probability` = recode(data.plot$`Reward probability`,
#                                       "High_Rew" = "High",
#                                       "Low_Rew" = "Low")
# 
# 
# 
# 
#   # Pirate plot
#   pirateplot(formula=Hit.Rate ~ `Reward phase` + `Reward probability`, # dependent~independent variables
#              data=data.plot, # data frame
#              main = "Hit rates",
#              ylim=c(0.1,0.9), # y-axis: limits
#              ylab="Hit Rate", # y-axis: label
#              theme=0, # preset theme (0: use your own)
#              point.col="black", # points: color
#              point.o=.3, # points: opacity (0-1)
#              avg.line.col="black", # average line: color
#              avg.line.lwd=2, # average line: line width
#              avg.line.o=1, # average line: opacity (0-1)
#              bean.b.col="black", # bean border, color
#              bean.lwd=0.6, # bean border, line width
#              bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
#              bean.b.o=0.3, # bean border, opacity (0-1)
#              bean.f.col="gray", # bean filling, color
#              bean.f.o=.1, # bean filling, opacity (0-1)
#              cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
#              gl.col="gray", # gridlines: color
#              gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
#              cex.lab=1, # axis labels: size
#              cex.axis=1, # axis numbers: size
#              cex.names = 1,
#              bty="l", # plot box type
#              back.col="white") # background, color


# # Plot False alarms------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Prepare the dataset
# data.plot = data.final
# 
# # rename variables
# colnames(data.plot)[colnames(data.plot)=="ExpPhase"] <- "Reward phase"
# colnames(data.plot)[colnames(data.plot)=="Condition"] <- "Reward probability"
# 
# # rename conditions
# data.plot$`Reward phase` = recode(data.plot$`Reward phase`,
#                                   "Acq" = "Acquisition",
#                                   "Bsln" = "Baseline",
#                                   "Ext" = "Extinction")
# 
# data.plot$`Reward probability` = recode(data.plot$`Reward probability`,
#                                         "High_Rew" = "High",
#                                         "Low_Rew" = "Low")




# Pirate plot
pirateplot(formula=FA.Rate ~ `Reward phase` + `Reward probability`, # dependent~independent variables
           data=data.plot, # data frame
           main = "False alarm rates",
           ylim=c(-0.1,0.7), # y-axis: limits
           ylab="False alarms", # y-axis: label
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



# # Plot D prime------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Prepare the dataset
# data.plot = data.final
# 
# # rename variables
# colnames(data.plot)[colnames(data.plot)=="ExpPhase"] <- "Reward phase"
# colnames(data.plot)[colnames(data.plot)=="Condition"] <- "Reward probability"
# 
# # rename conditions
# data.plot$`Reward phase` = recode(data.plot$`Reward phase`,
#                                   "Acq" = "Acquisition",
#                                   "Bsln" = "Baseline",
#                                   "Ext" = "Extinction")
# 
# data.plot$`Reward probability` = recode(data.plot$`Reward probability`,
#                                         "High_Rew" = "High",
#                                         "Low_Rew" = "Low")




# Pirate plot
# pirateplot(formula=dprime ~ `Reward phase` + `Reward probability`, # dependent~independent variables
#            data=data.plot, # data frame
#            main = "Dprime",
#            #ylim=c(-0.5,5), # y-axis: limits
#            ylab="Dprime", # y-axis: label
#            theme=0, # preset theme (0: use your own)
#            point.col="black", # points: color
#            point.o=.3, # points: opacity (0-1)
#            avg.line.fun = median,
#            avg.line.col="black", # average line: color
#            avg.line.lwd=2, # average line: line width
#            avg.line.o=1, # average line: opacity (0-1)
#            bean.b.col="black", # bean border, color
#            bean.lwd=0.6, # bean border, line width
#            bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
#            bean.b.o=0.3, # bean border, opacity (0-1)
#            bean.f.col="gray", # bean filling, color
#            bean.f.o=.1, # bean filling, opacity (0-1)
#            cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
#            gl.col="gray", # gridlines: color
#            gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
#            cex.lab=1, # axis labels: size
#            cex.axis=1, # axis numbers: size
#            cex.names = 1,
#            bty="l", # plot box type
#            back.col="white") # background, color



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

# # Set the working directory where to save the models
# setwd(here("brms_models"))
# 
# #help stan run faster
# rstan_options(auto_write = TRUE)
# options(mc.cores = parallel::detectCores())
# 
# # referencing for easier interpretation
# data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
# data.final$Condition=relevel(data.final$Condition,ref="High_Rew")
# 
# # # Set the prior for the intercept only model
# # prior = c(
# #   prior(normal(500, 200), class = Intercept)) # A wide prior sensible for this type of task
# # 
# # 
# # # Null model
# # model.null.RT = brm(RT ~ 1 + (1|ParticipantNo),
# #                  data=subset(data.final,data.final$Response=="Hit"),
# #                  family=exgaussian(),
# #                  prior = prior,
# #                  iter = 6000,
# #                  save_all_pars = TRUE,
# #                  control = list(adapt_delta = 0.99),
# #                  cores = 4,
# #                  sample_prior = TRUE,
# #                  inits = 0)
# # saveRDS(model.null.RT,file="nullmodel.RT.rds")
# # 
# # 
# # 
# # Set the priors for the models with slope
# prior = c(
#   prior(normal(500, 200), class = Intercept), # A wide prior sensible for this type of task
#   prior(normal(0, 200), class = b)) # a wide prior
# # 
# # # ExpPhase model
# # model.expphase.RT = brm(RT ~ ExpPhase + (ExpPhase|ParticipantNo),
# #                         data=subset(data.final,data.final$Response=="Hit"),
# #                         family=exgaussian(),
# #                         prior = prior,
# #                         iter = 6000,
# #                         save_all_pars = TRUE,
# #                         control = list(adapt_delta = 0.99),
# #                         cores = 4,
# #                         sample_prior = TRUE,
# #                         inits = 0)
# # saveRDS(model.expphase.RT,file="expphasemodel.RT.rds")
# 
# #Interaction model
# model.full.RT = brm(RT ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
#                     data=subset(data.final,data.final$Response=="Hit"),
#                     family=exgaussian(),
#                     prior = prior,
#                     iter = 6000,
#                     save_all_pars = TRUE,
#                     control = list(adapt_delta = 0.99),
#                     cores = 4,
#                     sample_prior = TRUE,
#                     inits = 0)
# saveRDS(model.full.RT,file="model.full.RT.rds")
# 
# # #Depression model
# # model.full.RT.depression = brm(Hits.RTs ~ ExpPhase * Condition * BDI + (ExpPhase * Condition * BDI|ParticipantNo),
# #                     data=data.final,
# #                     family=gaussian(),
# #                     prior = prior,
# #                     iter = 6000,
# #                     save_all_pars = TRUE,
# #                     control = list(adapt_delta = 0.99),
# #                     cores = 4,
# #                     sample_prior = TRUE)
# # saveRDS(model.full.RT.depression,file="model.full.RT.depression.rds")
# # 
# # #BAS model
# # model.full.RT.BAS = brm(Hits.RTs ~ ExpPhase * Condition * BAS + (ExpPhase * Condition * BAS|ParticipantNo),
# #                     data=data.final,
# #                     family=gaussian(),
# #                     prior = prior,
# #                     iter = 6000,
# #                     save_all_pars = TRUE,
# #                     control = list(adapt_delta = 0.99),
# #                     cores = 4,
# #                     sample_prior = TRUE)
# # saveRDS(model.full.RT.BAS,file="model.full.RT.BAS.rds")
# # 
# # #REWARD model
# # model.full.RT.REWARD = brm(Hits.RTs ~ ExpPhase * Condition * REWARD + (ExpPhase * Condition * REWARD|ParticipantNo),
# #                     data=data.final,
# #                     family=gaussian(),
# #                     prior = prior,
# #                     iter = 6000,
# #                     save_all_pars = TRUE,
# #                     control = list(adapt_delta = 0.99),
# #                     cores = 4,
# #                     sample_prior = TRUE)
# # saveRDS(model.full.RT.REWARD,file="model.full.RT.REWARD.rds")
# 
# # # read in the models and comparisons
# model.null.RT = readRDS("nullmodel.RT.rds")
# model.expphase.RT = readRDS("expphasemodel.RT.rds")
# #model.full.RT = readRDS("model.full.RT.rds")
# # compare.waic = readRDS("compare.RT.waic")
# 
# #WAIC
# compare.RT.waic = WAIC(model.null.RT,model.expphase.RT,model.full.RT, comapre = TRUE)
# saveRDS(compare.RT.waic,file="compare.RT.waic")
# 
# # # WAIC Questionnaires
# # compare.RT.waic.questionnaires = WAIC(model.full.RT,model.full.RT.depression, model.full.RT.REWARD, model.full.RT.BAS, comapre = TRUE)
# # saveRDS(compare.RT.waic.questionnaires,file="compare.RT.questionnaires.waic")
# 
# # Weighted waic
# compare.RT.waic.weights = model_weights(model.null.RT,model.expphase.RT,model.full.RT, weights = "waic")
# saveRDS(compare.RT.waic.weights,file="compare.RT.waic.weights")
# 
# # # Questionnaires Weighted waic
# # compare.RT.waic.weights.questionnaires = model_weights(model.full.RT,model.full.RT.depression, model.full.RT.REWARD, model.full.RT.BAS, comapre = TRUE)
# # saveRDS(compare.RT.waic.weights.questionnaires,file="compare.RT.waic.questionnaires.weights")
# 
# # Bayesian R2
# #Null
# bR2.null.RT = bayes_R2(model.null.RT)
# saveRDS(bR2.null.RT,file="bR2.null.RT")
# #ExpPhase
# bR2.expphase.RT = bayes_R2(model.expphase.RT)
# saveRDS(bR2.expphase.RT,file="bR2.expphase.RT")
# #Full
# bR2.full.RT = bayes_R2(model.full.RT)
# saveRDS(bR2.full.RT,file="bR2.full.RT")
# #Depression
# bR2.full.RT.depression = bayes_R2(model.full.RT.depression)
# saveRDS(bR2.full.RT.depression,file="bR2.full.RT.depression")
# #BAS
# bR2.full.RT.BAS = bayes_R2(model.full.RT.BAS)
# saveRDS(bR2.full.RT.BAS,file="bR2.full.RT.BAS")
# #REWARD
# bR2.full.RT.REWARD = bayes_R2(model.full.RT.REWARD)
# saveRDS(bR2.full.RT.REWARD,file="bR2.full.RT.REWARD")

# # Analyzing the posterior and differences between conditions
# 
# post = posterior_samples(model.full.RT, "^b")
# 
# 
# ################################################ Baseline ####
# 
# ######### High reward
# Baseline_High = post[["b_Intercept"]]
# ######### Low reward
# Baseline_Low = post[["b_Intercept"]] + 
#   post[["b_ConditionLow_Rew"]] 
# 
# ################################################ Acquistion
# 
# ######### High reward
# Acquisition_High = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseAcq"]] 
# ######### Low reward
# Acquisition_Low = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseAcq"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ExpPhaseAcq:ConditionLow_Rew"]]
# 
# ################################################ Extinction
# 
# ######### High reward
# Extinction_High = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseExt"]] 
# ######### Low reward
# Extinction_Low = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseExt"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ExpPhaseExt:ConditionLow_Rew"]]
# 
# 
# # Difference between high and low reward in baseline
# Diff_Rew_Bsln = Baseline_High - Baseline_Low
# plotPost(Diff_Rew_Bsln, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)
# 
# # Difference between high and low reward in acquisition 
# Diff_Rew_Acq = Acquisition_High - Acquisition_Low
# plotPost(Diff_Rew_Acq, xlab = "", col = "#b3cde0", cex = 1, showCurve = FALSE, compVal = 0)
# 
# # Difference between high and low reward in extinction
# Diff_Rew_Ext = Extinction_High - Extinction_Low
# plotPost(Diff_Rew_Ext, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)
# 
# 
# ########### plotting the posterior
# 
# # make a data frame
# posterior_conditions = melt(data.frame(Baseline_High, Baseline_Low, Acquisition_High, Acquisition_Low, Extinction_High, Extinction_Low))
# 
# posterior_conditions =  posterior_conditions %>% separate(variable, c("Reward Phase", "Reward Probability"), "_", extra = "merge")
# 
# names(posterior_conditions)[3] = "Reaction time"
# 
# # Pirate plot
# pirateplot(formula = `Reaction time` ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
#            data=posterior_conditions, # data frame
#            main="Reaction times", # main title
#            ylim=c(485,600), # y-axis: limits
#            ylab="Reaction time", # y-axis: label
#            theme=0, # preset theme (0: use your own)
#            avg.line.col="black", # average line: color
#            avg.line.lwd=2, # average line: line width
#            avg.line.o=1, # average line: opacity (0-1)
#            bean.b.col="black", # bean border, color
#            bean.lwd=0.6, # bean border, line width
#            bean.lty=1, # bean border, line type (1: solid; 2:dashed; 3: dotted; ...)
#            bean.b.o=0.3, # bean border, opacity (0-1)
#            bean.f.col="gray", # bean filling, color
#            bean.f.o=.1, # bean filling, opacity (0-1)
#            cap.beans=FALSE, # max and min values of bean densities are capped at the limits found in the data
#            gl.col="gray", # gridlines: color
#            gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
#            cex.lab=1, # axis labels: size
#            cex.axis=1, # axis numbers: size
#            cex.names = 1,
#            sortx = "sequential",
#            bty="l", # plot box type
#            back.col="white") # background, color


# # brms accuracy (hit rates)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# # Set the working directory where to save the models
setwd(here("brms_models"))
# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")
# 
# Set the prior for the intercept only model
prior = c(
  prior(normal(0.5, 1), class = Intercept)) # A wide prior sensible for this type of task

# Null model
model.null.Acc = brm(Response ~ 1 + (1|ParticipantNo),
                     data=data.final,
                     family=bernoulli(),
                     prior = prior,
                     iter = 6000,
                     save_all_pars = TRUE,
                     control = list(adapt_delta = 0.99),
                     cores = 4,
                     sample_prior = TRUE,
                     inits = 0)
saveRDS(model.null.Acc,file="nullmodel.Acc.rds")

# Set the priors for the models with slope
prior = c(
  prior(normal(0.5, 1), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 1), class = b)) # a wide prior

# ExpPhase model
model.expphase.Acc = brm(Response ~ ExpPhase + (ExpPhase|ParticipantNo),
                         data=data.final,
                         family=bernoulli(),
                         prior = prior,
                         iter = 6000,
                         save_all_pars = TRUE,
                         control = list(adapt_delta = 0.99),
                         cores = 4,
                         sample_prior = TRUE,
                         inits = 0)
saveRDS(model.expphase.Acc,file="expphasemodel.Acc.rds")

#Interaction model
model.full.Acc = brm(Response ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
                     data=data.final,
                     family=bernoulli(),
                     prior = prior,
                     iter = 6000,
                     save_all_pars = TRUE,
                     control = list(adapt_delta = 0.99),
                     cores = 4,
                     sample_prior = TRUE,
                     inits = 0)
saveRDS(model.full.Acc,file="model.full.Acc.rds")

# 
# #Depression model
# model.full.Acc.depression = brm(Hit.Rate ~ ExpPhase * Condition * BDI + (ExpPhase * Condition * BDI|ParticipantNo),
#                      data=data.final,
#                      family=gaussian(),
#                      prior = prior,
#                      iter = 6000,
#                      save_all_pars = TRUE,
#                      control = list(adapt_delta = 0.99),
#                      control = list(max_treedepth = 15),
#                      cores = 4,
#                      sample_prior = TRUE)
# saveRDS(model.full.Acc,file="model.full.Acc.depression.rds")
# 
# #BAS model
# model.full.Acc.BAS = brm(Hit.Rate ~ ExpPhase * Condition * BAS + (ExpPhase * Condition * BAS|ParticipantNo),
#                                 data=data.final,
#                                 family=gaussian(),
#                          prior = prior,
#                          iter = 6000,
#                          save_all_pars = TRUE,
#                          control = list(adapt_delta = 0.99),
#                          control = list(max_treedepth = 15),
#                          cores = 4,
#                          sample_prior = TRUE)
# saveRDS(model.full.Acc,file="model.full.Acc.BAS.rds")
# 
# #REWARD model
# model.full.Acc.REWARD = brm(Hit.Rate ~ ExpPhase * Condition * REWARD + (ExpPhase * Condition * REWARD|ParticipantNo),
#                          data=data.final,
#                          family=gaussian(),
#                          prior = prior,
#                          iter = 6000,
#                          save_all_pars = TRUE,
#                          control = list(adapt_delta = 0.99),
#                          control = list(max_treedepth = 15),
#                          cores = 4,
#                          sample_prior = TRUE)
# saveRDS(model.full.Acc,file="model.full.Acc.REWARD.rds")
# 
# #read in the models and comparisons
#  # model.null.Acc = readRDS("nullmodel.Acc.rds")
#  # model.expphase.Acc = readRDS("expphasemodel.Acc.rds")
#  # model.full.Acc = readRDS("model.full.Acc.rds")
# # compare.waic.Acc = readRDS("compare.Acc.waic")
# 
# #WAIC
compare.Acc.waic = WAIC(model.null.Acc,model.expphase.Acc,model.full.Acc, compare = TRUE)
saveRDS(compare.Acc.waic,file="compare.Acc.waic")
# 
# #WAIC questionnaires
# compare.Acc.waic.questionnaires = WAIC(model.full.Acc,model.full.Acc.depression,model.full.Acc.BAS,model.full.Acc.REWARD, compare = TRUE)
# saveRDS(compare.Acc.waic.questionnaires,file="compare.Acc.waic.questionnaires")
# 
# # Weighted waic
compare.Acc.waic.weights = model_weights(model.null.Acc,model.expphase.Acc,model.full.Acc, weights = "waic")
saveRDS(compare.Acc.waic.weights,file="compare.Acc.waic.weights")
# 
# # Weighted waic questionnaires
# compare.Acc.waic.questionnaires.weights = WAIC(model.full.Acc,model.full.Acc.depression,model.full.Acc.BAS,model.full.Acc.REWARD, compare = TRUE)
# saveRDS(compare.Acc.waic.questionnaires.weights,file="compare.Acc.waic.questionnaires.weights")
# 
# # Bayesian R2
#Null
bR2.null.Acc = bayes_R2(model.null.Acc)
saveRDS(bR2.null.Acc,file="bR2.null.Acc")
#ExpPhase
bR2.expphase.Acc = bayes_R2(model.expphase.Acc)
saveRDS(bR2.expphase.Acc,file="bR2.expphase.Acc")
#Full
bR2.full.Acc = bayes_R2(model.full.Acc)
saveRDS(bR2.full.Acc,file="bR2.full.Acc")
# #Depression
# bR2.full.Acc.depression = bayes_R2(model.full.Acc.depression)
# saveRDS(bR2.full.Acc.depression,file="bR2.full.Acc.depression")
# #BAS
# bR2.full.Acc.BAS = bayes_R2(model.full.Acc.BAS)
# saveRDS(bR2.full.Acc.BAS,file="bR2.full.Acc.BAS")
# #REWARD
# bR2.full.Acc.REWARD = bayes_R2(model.full.Acc.REWARD)
# saveRDS(bR2.full.Acc.REWARD,file="bR2.full.Acc.REWARD")

# # Analyzing the posterior and differences between conditions
# 
# post = posterior_samples(model.full.Acc, "^b")
# 
# 
