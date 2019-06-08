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

# Set the working directory where to save the models
setwd(here("brms_models"))

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# referencing for easier interpretation
data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
data.final$Condition=relevel(data.final$Condition,ref="High_Rew")

# # Set the prior for the intercept only model
# prior = c(
#   prior(normal(500, 200), class = Intercept)) # A wide prior sensible for this type of task
# 
# 
# # Null model
# model.null.RT = brm(RT ~ 1 + (1|ParticipantNo),
#                  data=subset(data.final,data.final$Response=="Hit"),
#                  family=exgaussian(),
#                  prior = prior,
#                  iter = 6000,
#                  save_all_pars = TRUE,
#                  control = list(adapt_delta = 0.99),
#                  cores = 4,
#                  sample_prior = TRUE,
#                  inits = 0)
# saveRDS(model.null.RT,file="nullmodel.RT.rds")
# 
# 
# 
# Set the priors for the models with slope
prior = c(
  prior(normal(500, 200), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 200), class = b)) # a wide prior
# 
# # ExpPhase model
# model.expphase.RT = brm(RT ~ ExpPhase + (ExpPhase|ParticipantNo),
#                         data=subset(data.final,data.final$Response=="Hit"),
#                         family=exgaussian(),
#                         prior = prior,
#                         iter = 6000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99),
#                         cores = 4,
#                         sample_prior = TRUE,
#                         inits = 0)
# saveRDS(model.expphase.RT,file="expphasemodel.RT.rds")

#Interaction model
model.full.RT = brm(RT ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
                    data=subset(data.final,data.final$Response=="Hit"),
                    family=exgaussian(),
                    prior = prior,
                    iter = 6000,
                    save_all_pars = TRUE,
                    control = list(adapt_delta = 0.99),
                    cores = 4,
                    sample_prior = TRUE,
                    inits = 0)
saveRDS(model.full.RT,file="model.full.RT.rds")

# #Depression model
# model.full.RT.depression = brm(Hits.RTs ~ ExpPhase * Condition * BDI + (ExpPhase * Condition * BDI|ParticipantNo),
#                     data=data.final,
#                     family=gaussian(),
#                     prior = prior,
#                     iter = 6000,
#                     save_all_pars = TRUE,
#                     control = list(adapt_delta = 0.99),
#                     cores = 4,
#                     sample_prior = TRUE)
# saveRDS(model.full.RT.depression,file="model.full.RT.depression.rds")
# 
# #BAS model
# model.full.RT.BAS = brm(Hits.RTs ~ ExpPhase * Condition * BAS + (ExpPhase * Condition * BAS|ParticipantNo),
#                     data=data.final,
#                     family=gaussian(),
#                     prior = prior,
#                     iter = 6000,
#                     save_all_pars = TRUE,
#                     control = list(adapt_delta = 0.99),
#                     cores = 4,
#                     sample_prior = TRUE)
# saveRDS(model.full.RT.BAS,file="model.full.RT.BAS.rds")
# 
# #REWARD model
# model.full.RT.REWARD = brm(Hits.RTs ~ ExpPhase * Condition * REWARD + (ExpPhase * Condition * REWARD|ParticipantNo),
#                     data=data.final,
#                     family=gaussian(),
#                     prior = prior,
#                     iter = 6000,
#                     save_all_pars = TRUE,
#                     control = list(adapt_delta = 0.99),
#                     cores = 4,
#                     sample_prior = TRUE)
# saveRDS(model.full.RT.REWARD,file="model.full.RT.REWARD.rds")

# # read in the models and comparisons
model.null.RT = readRDS("nullmodel.RT.rds")
model.expphase.RT = readRDS("expphasemodel.RT.rds")
  #model.full.RT = readRDS("model.full.RT.rds")
# compare.waic = readRDS("compare.RT.waic")

#WAIC
compare.RT.waic = WAIC(model.null.RT,model.expphase.RT,model.full.RT, comapre = TRUE)
saveRDS(compare.RT.waic,file="compare.RT.waic")

# # WAIC Questionnaires
# compare.RT.waic.questionnaires = WAIC(model.full.RT,model.full.RT.depression, model.full.RT.REWARD, model.full.RT.BAS, comapre = TRUE)
# saveRDS(compare.RT.waic.questionnaires,file="compare.RT.questionnaires.waic")

# Weighted waic
compare.RT.waic.weights = model_weights(model.null.RT,model.expphase.RT,model.full.RT, weights = "waic")
saveRDS(compare.RT.waic.weights,file="compare.RT.waic.weights")

# # Questionnaires Weighted waic
# compare.RT.waic.weights.questionnaires = model_weights(model.full.RT,model.full.RT.depression, model.full.RT.REWARD, model.full.RT.BAS, comapre = TRUE)
# saveRDS(compare.RT.waic.weights.questionnaires,file="compare.RT.waic.questionnaires.weights")

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
# # brms accuracy (hit rates vs. everything else)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# # Set the working directory where to save the models
setwd(here("brms_models"))
# referencing for easier interpretation
data.final.corr$ExpPhase=relevel(data.final.corr$ExpPhase,ref="Bsln")
data.final.corr$Condition=relevel(data.final.corr$Condition,ref="High_Rew")
# 
# Set the prior for the intercept only model
prior = c(
  prior(normal(0.5, 1), class = Intercept)) # A wide prior sensible for this type of task

# Null model
model.null.Cor = brm(Response ~ 1 + (1|ParticipantNo),
                     data=data.final.corr,
                     family=bernoulli(),
                     prior = prior,
                     iter = 6000,
                     save_all_pars = TRUE,
                     control = list(adapt_delta = 0.99),
                     cores = 4,
                     sample_prior = TRUE,
                     inits = 0)
saveRDS(model.null.Cor,file="nullmodel.Cor.rds")

# Set the priors for the models with slope
prior = c(
  prior(normal(0.5, 1), class = Intercept), # A wide prior sensible for this type of task
  prior(normal(0, 1), class = b)) # a wide prior

# ExpPhase model
model.expphase.Cor = brm(Response ~ ExpPhase + (ExpPhase|ParticipantNo),
                         data=data.final.corr,
                         family=bernoulli(),
                         prior = prior,
                         iter = 6000,
                         save_all_pars = TRUE,
                         control = list(adapt_delta = 0.99),
                         cores = 4,
                         sample_prior = TRUE,
                         inits = 0)
saveRDS(model.expphase.Cor,file="expphasemodel.Cor.rds")

#Interaction model
model.full.Cor = brm(Response ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
                     data=data.final.corr,
                     family=bernoulli(),
                     prior = prior,
                     iter = 6000,
                     save_all_pars = TRUE,
                     control = list(adapt_delta = 0.99),
                     cores = 4,
                     sample_prior = TRUE,
                     inits = 0)
saveRDS(model.full.Cor,file="model.full.Cor.rds")


# data.x = data.final.corr
# data.x$Response = ifelse(data.x$Response==1,0,1)
# data.x$Response = as.numeric(data.x$Response)
# # 
# library(lme4)
# x = glmer(Response~ExpPhase*Condition + (1|ParticipantNo),data=data.x, family = "binomial")
# summary(x)
# 
# 
# cat_plot(x,pred="ExpPhase", modx = "Condition")



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
compare.Cor.waic = WAIC(model.null.Cor,model.expphase.Cor,model.full.Cor, compare = TRUE)
saveRDS(compare.Cor.waic,file="compare.Cor.waic")
# 
# #WAIC questionnaires
# compare.Acc.waic.questionnaires = WAIC(model.full.Acc,model.full.Acc.depression,model.full.Acc.BAS,model.full.Acc.REWARD, compare = TRUE)
# saveRDS(compare.Acc.waic.questionnaires,file="compare.Acc.waic.questionnaires")
# 
# # Weighted waic
compare.Cor.waic.weights = model_weights(model.null.Cor,model.expphase.Cor,model.full.Cor, weights = "waic")
saveRDS(compare.Cor.waic.weights,file="compare.Cor.waic.weights")
# 
# # Weighted waic questionnaires
# compare.Acc.waic.questionnaires.weights = WAIC(model.full.Acc,model.full.Acc.depression,model.full.Acc.BAS,model.full.Acc.REWARD, compare = TRUE)
# saveRDS(compare.Acc.waic.questionnaires.weights,file="compare.Acc.waic.questionnaires.weights")
# 
# # Bayesian R2
#Null
bR2.null.Cor = bayes_R2(model.null.Cor)
saveRDS(bR2.null.Cor,file="bR2.null.Cor")
#ExpPhase
bR2.expphase.Cor = bayes_R2(model.expphase.Cor)
saveRDS(bR2.expphase.Cor,file="bR2.expphase.Cor")
#Full
bR2.full.Corr = bayes_R2(model.full.Cor)
saveRDS(bR2.full.Cor,file="bR2.full.Cor")

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
# names(posterior_conditions)[3] = "Hit rate"
# 
# # Pirate plot
# pirateplot(formula = `Hit rate` ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
#            data=posterior_conditions, # data frame
#            main="Hit rates", # main title
#            ylim=c(0.5, 0.7), # y-axis: limits
#            ylab="Hit rate", # y-axis: label
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
# 
# 
# 
# # brms accuracy (false alarms)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# # Set the working directory where to save the models
# setwd(here("brms_models"))
# # referencing for easier interpretation
# data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
# data.final$Condition=relevel(data.final$Condition,ref="High_Rew")
# 
# 
# # Null model
# model.null.Acc.FA = brm(FA.Rate ~ 1 + (1|ParticipantNo),
#                      data=data.final,
#                      family=gaussian(),
#                      warmup = 2000,
#                      iter = 10000,
#                      save_all_pars = TRUE,
#                      control = list(adapt_delta = 0.99),
#                      cores = 4,
#                      sample_prior = TRUE)
# saveRDS(model.null.Acc.FA,file="nullmodel.Acc.FA.rds")
# 
# # ExpPhase model
# model.expphase.Acc.FA = brm(FA.Rate ~ ExpPhase + (ExpPhase|ParticipantNo),
#                          data=data.final,
#                          family=gaussian(),
#                          warmup = 2000,
#                          iter = 10000,
#                          save_all_pars = TRUE,
#                          control = list(adapt_delta = 0.99),
#                          cores = 4,
#                          sample_prior = TRUE)
# saveRDS(model.expphase.Acc.FA,file="expphasemodel.Acc.FA.rds")
# 
# #Interaction model
# model.full.Acc.FA = brm(FA.Rate ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
#                      data=data.final,
#                      family=gaussian(),
#                      warmup = 2000,
#                      iter = 10000,
#                      save_all_pars = TRUE,
#                      control = list(adapt_delta = 0.99),
#                      cores = 4,
#                      sample_prior = TRUE)
# saveRDS(model.full.Acc.FA,file="model.full.Acc.FA.rds")
# 
# #read in the models and comparisons
# # model.null.Acc.FA = readRDS("nullmodel.Acc.FA.rds")
# # model.expphase.Acc.FA = readRDS("expphasemodel.Acc.FA.rds")
# # model.full.Acc.FA = readRDS("model.full.Acc.FA.rds")
# # compare.waic.Acc = readRDS("compare.Acc.waic")
# 
# #WAIC
# compare.Acc.FA.waic = WAIC(model.null.Acc.FA,model.expphase.Acc.FA,model.full.Acc.FA, compare = TRUE)
# saveRDS(compare.Acc.FA.waic,file="compare.Acc.FA.waic")
# 
# # Weighted waic
# compare.Acc.FA.waic.weights = model_weights(model.null.Acc.FA,model.expphase.Acc.FA,model.full.Acc.FA, weights = "waic")
# saveRDS(compare.Acc.FA.waic.weights,file="compare.Acc.FA.waic.weights")
# 
# # Bayesian R2
# #Null
# bR2.null.Acc.FA = bayes_R2(model.null.Acc.FA)
# saveRDS(bR2.null.Acc.FA,file="bR2.null.Acc.FA")
# #ExpPhase
# bR2.expphase.Acc.FA = bayes_R2(model.expphase.Acc.FA)
# saveRDS(bR2.expphase.Acc.FA,file="bR2.expphase.Acc.FA")
# #Full
# bR2.full.Acc.FA = bayes_R2(model.full.Acc.FA)
# saveRDS(bR2.full.Acc.FA,file="bR2.full.Acc.FA")
# 
# # Analyzing the posterior and differences between conditions
# 
# post = posterior_samples(model.full.Acc.FA, "^b")
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
# names(posterior_conditions)[3] = "FA rate"
# 
# # Pirate plot
# pirateplot(formula = `FA rate` ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
#            data=posterior_conditions, # data frame
#            main="FA rates", # main title
#            ylim=c(0.5, 0.7), # y-axis: limits
#            ylab="FA rate", # y-axis: label
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
# 
# 
# 
# brms accuracy (d prime)------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Set the working directory where to save the models
# setwd(here("brms_models"))
# # referencing for easier interpretation
# data.final$ExpPhase=relevel(data.final$ExpPhase,ref="Bsln")
# data.final$Condition=relevel(data.final$Condition,ref="High_Rew")
# 
# # Set the priors for the model with only intercept
# prior = c(
#   prior(normal(1.8, 1), class = Intercept)) # based on Andersen & Mueller, 2011
# 
# 
# # Null model
# model.null.Acc.dprime = brm(dprime ~ 1 + (1|ParticipantNo),
#                         data=data.final,
#                         family=gaussian(),
#                         prior = prior,
#                         iter = 6000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99),
#                         cores = 4,
#                         sample_prior = TRUE,
#                         inits = 0)
# saveRDS(model.null.Acc.dprime,file="nullmodel.Acc.dprime.rds")
# 
# # Set the priors for the models with slopes
# prior = c(
#   prior(normal(1.8, 1), class = Intercept), # based on Andersen & Mueller, 2011
#   prior(normal(0, 2), class = b)) # a wide prior
# 
# # ExpPhase model
# model.expphase.Acc.dprime = brm(dprime ~ ExpPhase + (ExpPhase|ParticipantNo),
#                             data=data.final,
#                             family=gaussian(),
#                             prior = prior,
#                             iter = 6000,
#                             save_all_pars = TRUE,
#                             control = list(adapt_delta = 0.99),
#                             cores = 4,
#                             sample_prior = TRUE,
#                             inits = 0)
# saveRDS(model.expphase.Acc.dprime,file="expphasemodel.Acc.dprime.rds")
# 
# #Interaction model
# model.full.Acc.dprime = brm(dprime ~ ExpPhase * Condition + (ExpPhase * Condition|ParticipantNo),
#                         data=data.final,
#                         family=gaussian(),
#                         prior = prior,
#                         iter = 6000,
#                         save_all_pars = TRUE,
#                         control = list(adapt_delta = 0.99),
#                         cores = 4,
#                         sample_prior = TRUE,
#                         inits = 0)
# saveRDS(model.full.Acc.dprime,file="model.full.Acc.dprime.rds")
# 


#read in the models and comparisons
# model.null.Acc.FA = readRDS("nullmodel.Acc.FA.rds")
# model.expphase.Acc.FA = readRDS("expphasemodel.Acc.FA.rds")
# model.full.Acc.FA = readRDS("model.full.Acc.FA.rds")
# compare.waic.Acc = readRDS("compare.Acc.waic")

# #WAIC
# compare.Acc.dprime.waic = WAIC(model.null.Acc.dprime,model.expphase.Acc.dprime,model.full.Acc.dprime, compare = TRUE)
# saveRDS(compare.Acc.dprime.waic,file="compare.Acc.dprime.waic")
# 
# # Weighted waic
# compare.Acc.dprime.waic.weights = model_weights(model.null.Acc.dprime,model.expphase.Acc.dprime,model.full.Acc.dprime, weights = "waic")
# saveRDS(compare.Acc.dprime.waic.weights,file="compare.Acc.dprime.waic.weights")
# 
# # Bayesian R2
# #Null
# bR2.null.Acc.dprime = bayes_R2(model.null.Acc.dprime)
# saveRDS(bR2.null.Acc.dprime,file="bR2.null.Acc.dprime")
# #ExpPhase
# bR2.expphase.Acc.dprime = bayes_R2(model.expphase.Acc.dprime)
# saveRDS(bR2.expphase.Acc.dprime,file="bR2.expphase.Acc.dprime")
# #Full
# bR2.full.Acc.dprime = bayes_R2(model.full.Acc.dprime)
# saveRDS(bR2.full.Acc.dprime,file="bR2.full.Acc.dprime")
# 
# # Analyzing the posterior and differences between conditions
# 
# post = posterior_samples(model.full.Acc.dprime, "^b")
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
# names(posterior_conditions)[3] = "Dprime"
# 
# # Pirate plot
# pirateplot(formula = `Dprime` ~ `Reward Phase` + `Reward Probability`, # dependent~independent variables
#            data=posterior_conditions, # data frame
#            main="D prime", # main title
#            ylim=c(1, 2.5), # y-axis: limits
#            ylab="FA rate", # y-axis: label
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
# 
# 
