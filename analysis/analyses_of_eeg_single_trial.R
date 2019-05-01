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

# Add new variable names based on the condition
data$ExpPhase[data$condition == 1 | data$condition == 2 | data$condition == 11 | data$condition == 12]="Bsln"
data$ExpPhase[data$condition == 3 | data$condition == 4 | data$condition == 13 | data$condition == 14]="Acq"
data$ExpPhase[data$condition == 5 | data$condition == 6 | data$condition == 15 | data$condition == 16]="Ext"

  data.raw$Rew_eff_fb[data.raw$feedbackOrder == 2 & data.raw$IsRewarded == "Reward"]="Reward"

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
#Make a new variable with mean amplitude across all conditions for each participant and each frequency  !!! originally did not include ExpPhase below !!!
data = ddply(data,.(Subject,RecordedFrequency),transform,
                    MeanAmplitude = mean(Amplitude[ExpPhase=="Bsln"],na.rm=TRUE),
                    SDAmplitude =   sd(Amplitude,na.rm=TRUE))

# data = ddply(data,.(Subject,RecordedFrequency),transform,
#              MeanAmplitude = mean(Amplitude,na.rm=TRUE),
#              SDAmplitude =   sd(Amplitude,na.rm=TRUE))

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
           iter = 6000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99,max_treedepth = 15),
           cores = 4,
           sample_prior = TRUE,
           inits = 0)
saveRDS(null,file="null.EEG.allsub.rds")

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
               iter = 6000,
               save_all_pars = TRUE,
               control = list(adapt_delta = 0.99,max_treedepth = 15),
               cores = 4,
               sample_prior = TRUE,
               inits = 0)
saveRDS(expphase,file="expphase.EEG.allsubs.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Attention|Subject),
                data=data,
                family=gaussian(),
                prior = prior,
                iter = 6000,
                save_all_pars = TRUE,
                control = list(adapt_delta = 0.99,max_treedepth = 15),
                cores = 4,
                sample_prior = TRUE,
                inits = 0)
saveRDS(attention,file="attention.EEG.allsubs.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + ( ExpPhase + Attention|Subject),
                        data=data,
                        family=gaussian(),
                        prior = prior,
                        iter = 6000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99,max_treedepth = 15),
                        cores = 4,
                        sample_prior = TRUE,
                        inits = 0)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs.rds")

# Interaction between phase and attention
phaseANDattention_interaction = brm(Amplitude ~ ExpPhase * Attention + ( ExpPhase + Attention|Subject),
                        data=data,
                        family=gaussian(),
                        prior = prior,
                        iter = 6000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99,max_treedepth = 15),
                        cores = 4,
                        sample_prior = TRUE,
                        inits = 0)
saveRDS(phaseANDattention_interaction,file="phaseANDattention_interaction.EEG.allsubs.rds")

# Interaction between expphase and reward magnitude plus attention
rewardTimesPhasePlusAtt = brm(Amplitude ~ Condition * ExpPhase + Attention + (Condition + ExpPhase + Attention|Subject),
                         data=data,
                         family=gaussian(),
                         prior = prior,
                         iter = 6000,
                         save_all_pars = TRUE,
                         control = list(adapt_delta = 0.99,max_treedepth = 15),
                         cores = 4,
                         sample_prior = TRUE,
                         inits = 0)
saveRDS(rewardTimesPhasePlusAtt,file="rewardTimesPhasePlusAtt.EEG.allsubs.rds")

# Full model


full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition + ExpPhase + Attention|Subject),   
                             data=data,
                             family=gaussian(),
                             prior = prior,
                             iter = 6000,
                             save_all_pars = TRUE,
                             control = list(adapt_delta = 0.99,max_treedepth = 15),
                             cores = 4,
                             sample_prior = TRUE,
                             inits = 0)
saveRDS(full,file="full.EEG.allsubs.rds")



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
compare.EEG.waic = WAIC(null, expphase, attention, phaseANDattention, phaseANDattention_interaction, rewardTimesPhasePlusAtt, full, compare = TRUE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs.rds")

# # Questionnaires WAIC
# compare.EEG.waic.questionnaires = WAIC(full, depression, bas, reward, compare = TRUE)
# saveRDS(compare.EEG.waic.questionnaires,file="compare.EEG.waic.allsubs.questionnaires.rds")

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


# Analyzing the posterior and differences between conditions

# post = posterior_samples(full, "^b")

# 
# ################################################ Baseline ####
# 
# ##################### Attended
# 
# ######### High reward
# Baseline_High_Attended = post[["b_Intercept"]]
# ######### Low reward
# Baseline_Low_Attended = post[["b_Intercept"]] + 
#   post[["b_ConditionLow_Rew"]] 
# 
# ##################### Not Attended
# 
# ######### High reward
# Baseline_High_NotAttended = post[["b_Intercept"]] + 
#   post[["b_AttentionNotAtt"]]
# ######### Low reward
# Baseline_Low_NotAttended = post[["b_Intercept"]] + 
#   post[["b_AttentionNotAtt"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ConditionLow_Rew:AttentionNotAtt"]]
# 
# ################################################ Acquistion
# 
# ##################### Attended
# 
# ######### High reward
# Acquisition_High_Attended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseAcq"]] 
# ######### Low reward
# Acquisition_Low_Attended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseAcq"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ConditionLow_Rew:ExpPhaseAcq"]]
# 
# ##################### Not Attended
# 
# ######### High reward
# Acquisition_High_NotAttended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseAcq"]] + 
#   post[["b_AttentionNotAtt"]] + 
#   post[["b_ExpPhaseAcq:AttentionNotAtt"]]
# ######### Low reward
# Acquisition_Low_NotAttended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseAcq"]] + 
#   post[["b_AttentionNotAtt"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ExpPhaseAcq:AttentionNotAtt"]] +
#   post[["b_ConditionLow_Rew:ExpPhaseAcq"]] + 
#   post[["b_ConditionLow_Rew:ExpPhaseAcq:AttentionNotAtt"]]
# 
# ################################################ Extinction
# 
# ##################### Attended
# 
# ######### High reward
# Extinction_High_Attended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseExt"]] 
# ######### Low reward
# Extinction_Low_Attended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseExt"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ConditionLow_Rew:ExpPhaseExt"]]
# 
# ##################### Not Attended
# 
# ######### High reward
# Extinction_High_NotAttended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseExt"]] + 
#   post[["b_AttentionNotAtt"]] + 
#   post[["b_ExpPhaseExt:AttentionNotAtt"]]
# ######### Low reward
# Extinction_Low_NotAttended = post[["b_Intercept"]] + 
#   post[["b_ExpPhaseExt"]] + 
#   post[["b_AttentionNotAtt"]] + 
#   post[["b_ConditionLow_Rew"]] + 
#   post[["b_ExpPhaseExt:AttentionNotAtt"]] +
#   post[["b_ConditionLow_Rew:ExpPhaseExt"]] + 
#   post[["b_ConditionLow_Rew:ExpPhaseExt:AttentionNotAtt"]]
# 
# 
# 
# ### Plotting the posterior ###
# 
# # make a data frame
# 
# posterior_conditions = melt(data.frame(Baseline_High_Attended, Baseline_High_NotAttended, Baseline_Low_Attended, Baseline_Low_NotAttended, Acquisition_High_Attended, Acquisition_High_NotAttended, Acquisition_Low_Attended, Acquisition_Low_NotAttended, Extinction_High_Attended, Extinction_High_NotAttended, Extinction_Low_Attended, Extinction_Low_NotAttended))
# 
# posterior_conditions =  posterior_conditions %>% separate(variable, c("Reward Phase", "Reward Probability", "Attention"), "_", extra = "merge")
# 
# posterior_conditions$Attention = recode(posterior_conditions$Attention,
#                                         "Attended" = "Attended",
#                                         "NotAttended" = "Unattended")
# 
# names(posterior_conditions)[4] = "Amplitude"
# 
# 
# #order
# #dataPlot$`Reward phase` = factor(dataPlot$`Reward phase`, levels = c("Baseline","Acquisition","Extinction"))
# #dataPlot = dataPlot[order(dataPlot$Attention,dataPlot$`Reward phase`,dataPlot$`Reward probability`),]
# 
# 
# plottingConditions = c("Attended","Unattended" )
# for (i in 1:length(plottingConditions)){
#   
#   if(plottingConditions[i]=="Attended"){dataAmplitudePlot=subset(posterior_conditions,Attention=="Attended")}
#   
#   if(plottingConditions[i]=="Unattended"){dataAmplitudePlot=subset(posterior_conditions,Attention=="Unattended")}  
#   
#   # Pirate plot
#   
#   pirateplot(formula = Amplitude ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
#              data=dataAmplitudePlot, # data frame
#              main=plottingConditions[i], # main title
#              ylim=c(0.7,1.2), # y-axis: limits
#              ylab=expression(paste("Amplitude (",mu,"V)")), # y-axis: label
#              theme=0, # preset theme (0: use your own)
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
#              sortx = "sequential",
#              bty="l", # plot box type
#              back.col="white") # background, color
# }
# 
# #Check the difference between high and low reward in acquisition attended
# 
# # Difference between high and low reward in acquisition attended
# Diff_Rew_Acq_Att = Acquisition_High_Attended - Acquisition_Low_Attended
# plotPost(Diff_Rew_Acq_Att, xlab = "", col = "#b3cde0", cex = 1, showCurve = FALSE, compVal = 0)
# 
# mean(Acquisition_High_Attended>Acquisition_Low_Attended)
# 
# 
# #Check the difference between high and low reward in baseline attended
# 
# # Difference between high and low reward in acquisition attended
# Diff_Rew_Bsln_Att = Baseline_High_Attended - Baseline_Low_Attended
# plotPost(Diff_Rew_Bsln_Att, xlab = "", col = "#b3cde0", showCurve = FALSE, cex = 1, compVal = 0)





