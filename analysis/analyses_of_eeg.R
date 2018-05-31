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
pacman::p_load(reshape2,yarrr,BayesFactor,plyr,ez,schoRsch,brms,lme4)
# set seed
set.seed(42) 
# import data
data.raw = read.csv2(file="C:/Users/igrahek/Documents/Studies/SSVEP Reward - Soren & Antonio/Experiment 1/SSVEP and reward/data/amplitudes_rewardBoth_wholeSample.csv",header=TRUE,na.strings="NaN") #only good behavior: amplitudes_rewardBoth.csv  # full sample: amplitudes_rewardBoth_wholeSample.csv POORPERF_amplitude_rewardBoth.csv

# Prepare the dataset------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Reshape to long format
data = melt(data.raw,id.vars=c("Subject","Frequency"),
             measure.vars=c("BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended"),
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
# Make a new variable with mean amplitude across all conditions for each participant and each frequency 
data = ddply(data,.(Subject,RecordedFrequency),transform,
                    MeanAmplitude = mean(Amplitude,na.rm=TRUE),
                    SDAmplitude =   sd(Amplitude,na.rm=TRUE))

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

################################################################## Plotting ###############################################################################################################################################################################################################

# Plot amplitude across experiment phases------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

plottingConditions = c("Amplitude both reward conditions","Amplitude High_Rew","Amplitude Low_Rew" )
for (i in 1:length(plottingConditions)){
  
  if(plottingConditions[i]=="Amplitude both reward conditions"){
    
    #Average over the reward condition and order data for plotting
    dataAmplitudePlot=ddply(data,.(Subject,Attention,ExpPhase),plyr::summarize,Amplitude=mean(Amplitude,na.rm=TRUE)) 
    
    #Order the data again in order to be able to plot Condition in the same order as in the other plots
    dataAmplitudePlot$Attention = factor(dataAmplitudePlot$Attention, levels = c("NotAtt","Att"))
    dataAmplitudePlot$ExpPhase = factor(dataAmplitudePlot$ExpPhase, levels = c("Bsln","Acq","Ext"))
    dataAmplitudePlot = dataAmplitudePlot[order(dataAmplitudePlot$Subject,dataAmplitudePlot$ExpPhase,dataAmplitudePlot$Attention),]}
  
  if(plottingConditions[i]=="Amplitude High_Rew"){dataAmplitudePlot=subset(data,Condition=="High_Rew")}
  
  if(plottingConditions[i]=="Amplitude Low_Rew"){
    
    #Order the data again in order to be able to plot Condition in the same order as in the other plots
    dataAmplitudePlot=subset(data,Condition=="Low_Rew")
    dataAmplitudePlot$Attention = factor(dataAmplitudePlot$Attention, levels = c("NotAtt","Att"))
    dataAmplitudePlot$ExpPhase = factor(dataAmplitudePlot$ExpPhase, levels = c("Bsln","Acq","Ext"))
    dataAmplitudePlot = dataAmplitudePlot[order(dataAmplitudePlot$Subject,dataAmplitudePlot$ExpPhase,dataAmplitudePlot$Attention),] #order the data again in order to be able to plot Condition in the same order as in the other plots
    }  

# Pirate plot

  pirateplot(formula=Amplitude~Attention+ExpPhase, # dependent~independent variables
             data=dataAmplitudePlot, # data frame
             main=plottingConditions[i], # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(0.2,2.2), # y-axis: limits
             ylab=expression(paste("Amplitude (",mu,"V)")), # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=Subject, # ID variable
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
             sortx="sequential",
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color
}
  

# Plot the selectivity index across experiment phases and reward magnitude------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Pirate plot

pirateplot(formula=Selectivity ~ Condition + ExpPhase, # dependent~independent variables
           data=data.diff, # data frame
           main="Selectivity index", # main title
           xlim=NULL, # x-axis: limits
           xlab="", # x-axis: label
           ylim=c(-1,1.5), # y-axis: limits
           ylab=expression(paste("Amplitude diff (",mu,"V)")), # y-axis: label
           inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
           hdi.iter=5000, # number of iterations for estimation of HDI
           inf.within=Subject, # ID variable
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
           sortx="sequential",
           gl.col="gray", # gridlines: color
           gl.lwd=c(.75,0), # gridlines: line width
           gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
           cex.lab=0.8, # axis labels: size
           cex.axis=1, # axis numbers: size
           bty="l", # plot box type
           back.col="white") # background, color

# Plot the total enhancement index across experiment phases and reward magnitude------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

    # Pirate plot
    
    pirateplot(formula=TotalEnhancement ~ ExpPhase + Condition, # dependent~independent variables
               data=data.diff, # data frame
               main="Total enhancement index", # main title
               xlim=NULL, # x-axis: limits
               xlab="", # x-axis: label
               ylim=c(1,3.5), # y-axis: limits
               ylab=expression(paste("Amplitude sum (",mu,"V)")), # y-axis: label
               inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
               hdi.iter=5000, # number of iterations for estimation of HDI
               inf.within=Subject, # ID variable
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
               sortx="sequential",
               gl.col="gray", # gridlines: color
               gl.lwd=c(.75,0), # gridlines: line width
               gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
               cex.lab=0.8, # axis labels: size
               cex.axis=1, # axis numbers: size
               bty="l", # plot box type
               back.col="white") # background, color
    
##### STATS ####

# # Set the parameters 
# num.iter=10000 # number of MonteCarlo iterations (default: 10000)
# # select scaling factor r of Cauchy(0,r) prior on standardized effect sizes
# medprior=sqrt(2)/2 # medium prior
# 
# 
# 
# # Amplitude - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) X 2(high vs. low reward) 
# ssVEPamp.bf.amplitude <- anovaBF(Amplitude ~ ExpPhase * Condition  * Attention + Subject,data=data,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
# sort(ssVEPamp.bf.amplitude)
# 
# # Amplitude - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) 
# ssVEPamp.bf.amplitude <- anovaBF(Amplitude ~ ExpPhase * Attention + Subject,data=data,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
# sort(ssVEPamp.bf.amplitude)
# 
# 
# # Selectivity index - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) X 2(high vs. low reward) 
# ssVEPamp.bf.selectivity <- anovaBF(Selectivity ~ ExpPhase * Condition + Subject,data=data.diff,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
# sort(ssVEPamp.bf.selectivity)
# 
# # Total enhancement index - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) X 2(high vs. low reward) 
# ssVEPamp.bf.totalenhancement <- anovaBF(TotalEnhancement ~ ExpPhase * Condition + Subject,data=data.diff,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
# sort(ssVEPamp.bf.totalenhancement)
# 

# brms three factors------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Set the working directory where to save the models
setwd("./brms_models")

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Modelling the effects of phase, attention, and reward magnitude - All subjects

# Set the intercept model
data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High_Rew")
data$Attention=relevel(data$Attention,ref="Att")

# Null model
null = brm(Amplitude ~ 1 + (1|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4)
saveRDS(null,file="null.EEG.allsub.rds")

# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                                  data=data,
                                  family=gaussian(),
                                  warmup = 2000,
                                  iter = 10000,
                                  save_all_pars = TRUE,
                                  control = list(adapt_delta = 0.99),
                                  cores = 4)
saveRDS(expphase,file="expphase.EEG.allsubs.rds")

# Condition model
condition = brm(Amplitude ~ Condition + (Condition|Subject),
                                   data=data,
                                   family=gaussian(),
                                   warmup = 2000,
                                   iter = 10000,
                                   save_all_pars = TRUE,
                                   control = list(adapt_delta = 0.99),
                                   cores = 4)
saveRDS(condition,file="condition.EEG.allsubs.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Attention|Subject),
                                   data=data,
                                   family=gaussian(),
                                   warmup = 2000,
                                   iter = 10000,
                                   save_all_pars = TRUE,
                                   control = list(adapt_delta = 0.99),
                                   cores = 4)
saveRDS(attention,file="attention.EEG.allsubs.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + (ExpPhase + Attention|Subject),
                                        data=data,
                                        family=gaussian(),
                                        warmup = 2000,
                                        iter = 10000,
                                        save_all_pars = TRUE,
                                        control = list(adapt_delta = 0.99),
                                        cores = 4)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs.rds")

# Two main effects - reward magnitude and attention
rewardANDattention = brm(Amplitude ~ Condition + Attention + (Condition + Attention|Subject),
                                           data=data,
                                           family=gaussian(),
                                           warmup = 2000,
                                           iter = 10000,
                                           save_all_pars = TRUE,
                                           control = list(adapt_delta = 0.99),
                                           cores = 4)
saveRDS(rewardANDattention,file="rewardANDattention.EEG.allsubs.rds")

# Three main effects
threemain = brm(Amplitude ~ Condition + ExpPhase + Attention + (Condition + ExpPhase + Attention|Subject),
                                          data=data,
                                          family=gaussian(),
                                          warmup = 2000,
                                          iter = 10000,
                                          save_all_pars = TRUE,
                                          control = list(adapt_delta = 0.99),
                                          cores = 4)
saveRDS(threemain,file="threemain.EEG.allsubs.rds")

# Full model
full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE,
                              control = list(adapt_delta = 0.99),
                              cores = 4)
saveRDS(full,file="full.EEG.allsubs.rds")

# read in the models and comparisons
null = readRDS("null.EEG.allsub.rds")
condition = readRDS("condition.EEG.allsubs.rds")
attention = readRDS("attention.EEG.allsubs.rds")
expphase = readRDS("expphase.EEG.allsubs.rds")
rewardANDattention = readRDS("rewardANDattention.EEG.allsubs.rds")
phaseANDattention = readRDS("phaseANDattention.EEG.allsubs.rds")
threemain = readRDS("threemain.EEG.allsubs.rds")
full = readRDS("full.EEG.allsubs.rds")

#WAIC
compare.EEG.waic = WAIC(null, condition, expphase, attention, phaseANDattention, rewardANDattention, threemain, full, compare = FALSE)
saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs.rds")
#LOO crossvalidation
compare.EEG.loo = LOO(null, condition, expphase, attention, threemain, reloo = TRUE, compare = FALSE)
saveRDS(compare.EEG.loo,file="compare.EEG.loo.allsubs.rds")

# Sample from the posterior
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
  post[["b_AttentionNotAtt"]]
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
  post[["b_AttentionNotAtt"]]
######### Low reward
Extinction_Low_NotAttended = post[["b_Intercept"]] + 
  post[["b_ExpPhaseExt"]] + 
  post[["b_AttentionNotAtt"]] + 
  post[["b_ConditionLow_Rew"]] + 
  post[["b_ExpPhaseExt:AttentionNotAtt"]] +
  post[["b_ConditionLow_Rew:ExpPhaseExt"]] + 
  post[["b_ConditionLow_Rew:ExpPhaseExt:AttentionNotAtt"]]


#Check the difference between high and low reward in acquisition attended

# Difference between high and low reward in acquisition attended
Diff_Rew_Acq_Att = Acquisition_High_Attended - Acquisition_Low_Attended
plotPost(Diff_Rew_Acq_Att, xlab = "", col = "#b3cde0", cex = 1, showCurve = TRUE, ROPE = c(0,1))

# Evidence ratio for the hypothesis that the high rewarded condition is higher than the low rewarded condition for attended
h1 = hypothesis(post, "0 > b_ConditionLow_Rew  + b_ConditionLow_Rew:ExpPhaseAcq")
print(h1)
plot(h1)



########### All varying effects in all models ###########
# Set the working directory where to save the models
setwd("./brms_models")

#help stan run faster
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Modelling the effects of phase, attention, and reward magnitude - All subjects

# Set the intercept model
data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High_Rew")
data$Attention=relevel(data$Attention,ref="Att")

# Null model
null = brm(Amplitude ~ 1 + (Condition * ExpPhase * Attention|Subject),
           data=data,
           family=gaussian(),
           warmup = 2000,
           iter = 10000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4)
saveRDS(null,file="null.EEG.allsub_allvaryinginallmodels.rds")

# Exp phase model
expphase = brm(Amplitude ~ ExpPhase + (Condition * ExpPhase * Attention|Subject),
               data=data,
               family=gaussian(),
               warmup = 2000,
               iter = 10000,
               save_all_pars = TRUE,
               control = list(adapt_delta = 0.99),
               cores = 4)
saveRDS(expphase,file="expphase.EEG.allsubs_allvaryinginallmodels.rds")

# Condition model
condition = brm(Amplitude ~ Condition + (Condition * ExpPhase * Attention|Subject),
                data=data,
                family=gaussian(),
                warmup = 2000,
                iter = 10000,
                save_all_pars = TRUE,
                control = list(adapt_delta = 0.99),
                cores = 4)
saveRDS(condition,file="condition.EEG.allsubs_allvaryinginallmodels.rds")

# Attention model
attention = brm(Amplitude ~ Attention + (Condition * ExpPhase * Attention|Subject),
                data=data,
                family=gaussian(),
                warmup = 2000,
                iter = 10000,
                save_all_pars = TRUE,
                control = list(adapt_delta = 0.99),
                cores = 4)
saveRDS(attention,file="attention.EEG.allsubs_allvaryinginallmodels.rds")

# Two main effects - phase and attention
phaseANDattention = brm(Amplitude ~ ExpPhase + Attention + (Condition * ExpPhase * Attention|Subject),
                        data=data,
                        family=gaussian(),
                        warmup = 2000,
                        iter = 10000,
                        save_all_pars = TRUE,
                        control = list(adapt_delta = 0.99),
                        cores = 4)
saveRDS(phaseANDattention,file="phaseANDattention.EEG.allsubs_allvaryinginallmodels.rds")

# Two main effects - reward magnitude and attention
rewardANDattention = brm(Amplitude ~ Condition + Attention + (Condition * ExpPhase * Attention|Subject),
                         data=data,
                         family=gaussian(),
                         warmup = 2000,
                         iter = 10000,
                         save_all_pars = TRUE,
                         control = list(adapt_delta = 0.99),
                         cores = 4)
saveRDS(rewardANDattention,file="rewardANDattention.EEG.allsubs_allvaryinginallmodels.rds")

# Three main effects
threemain = brm(Amplitude ~ Condition + ExpPhase + Attention + (Condition * ExpPhase * Attention|Subject),
                data=data,
                family=gaussian(),
                warmup = 2000,
                iter = 10000,
                save_all_pars = TRUE,
                control = list(adapt_delta = 0.99),
                cores = 4)
saveRDS(threemain,file="threemain.EEG.allsubs_allvaryinginallmodels.rds")

# Full model
full = brm(Amplitude ~ Condition * ExpPhase * Attention + (Condition * ExpPhase * Attention|Subject),
           data=data,
           family=gaussian(),
           warmup = 2000,
           iter = 10000,
           save_all_pars = TRUE,
           control = list(adapt_delta = 0.99),
           cores = 4)
saveRDS(full,file="full.EEG.allsubs_allvaryinginallmodels.rds")

model.full.threefactors = readRDS("full.EEG.allsubs_allvaryinginallmodels.rds")
model.threemain = readRDS("threemain.EEG.allsubs_allvaryinginallmodels.rds")
model.null = readRDS("null.EEG.allsub_allvaryinginallmodels.rds")
model.attention = readRDS("attention.EEG.allsubs_allvaryinginallmodels.rds")

#WAIC
compare.EEG.waic = WAIC(model.null, model.attention, model.threemain, model.full.threefactors, compare = FALSE) #phaseANDattention, rewardANDattention, null, condition, expphase, attention, 
saveRDS(compare.EEG.waic,file="compare.EEG.waic.allsubs_allvaryinginallmodels.rds")
