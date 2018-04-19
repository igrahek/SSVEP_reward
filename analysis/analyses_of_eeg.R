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
data.raw = read.csv2(file="./data/amplitudes_rewardBoth_wholeSample.csv",header=TRUE,na.strings="NaN") #only good behavior: amplitudes_rewardBoth.csv

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
data$Condition = ifelse(data$AttendedColor==data$RewardedColor, "High Reward Attended","Low Reward Attended")

# Make a new condition based on the attended color and the recorded frequency
data$ConditionRecording = ifelse(data$AttendedColor==data$RecordedFrequency, "AttRec","NotAttRec")

# Make a new condition based the Condition and the ConditionRecording
data$RecordingAndCondition = with(data, paste0(Condition,"_",ConditionRecording))

# Select variables which we want to keep
data = subset(data, select=c("Subject","RewardedColor","ExpPhase","AttendedColor","Condition","RecordedFrequency","ConditionRecording","RecordingAndCondition","Amplitude"))

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
data.diff = ddply(data, .(Subject,ExpPhase,Condition), transform, Selectivity = Amplitude[ConditionRecording=="AttRec"]-Amplitude[ConditionRecording=="NotAttRec"],TotalEnhancement=Amplitude[ConditionRecording=="AttRec"]+Amplitude[ConditionRecording=="NotAttRec"])
# Delete the ConditionRecording column and rows which are not necessary (indexes repeated twice)
data.diff = subset(data.diff,ConditionRecording=="AttRec") #keep only AttRec as it is equal to NotAttRec
data.diff$ConditionRecording = NULL  #drop the ConditionRecording column

# Sort the data 
data.diff$ExpPhase = factor(data.diff$ExpPhase, levels = c("Bsln","Acq","Ext"))
data.diff = data.diff[order(data.diff$Subject,data.diff$Condition,data.diff$ExpPhase),]

# Convert variables to be used in analyses into factors
data[c("Subject", "Condition","ExpPhase", "RewardedColor", "ConditionRecording", "RecordingAndCondition")] = 
  lapply(data.raw[c("Subject", "Condition","ExpPhase", "RewardedColor", "ConditionRecording", "RecordingAndCondition")], factor)

data.dif[c("Subject", "Condition","ExpPhase", "RecordingAndCondition")] = 
  lapply(data.raw[c("Subject", "Condition","ExpPhase",  "RecordingAndCondition")], factor)

################################################################## Plotting ###############################################################################################################################################################################################################

# Plot amplitude across experiment phases------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

plottingConditions = c("Amplitude both reward conditions","Amplitude high reward attended","Amplitude low reward attended" )
for (i in 1:length(plottingConditions)){
  
  if(plottingConditions[i]=="Amplitude both reward conditions"){
    
    #Average over the reward condition and order data for plotting
    dataAmplitudePlot=ddply(data,.(Subject,ConditionRecording,ExpPhase),plyr::summarize,Amplitude=mean(Amplitude,na.rm=TRUE)) 
    
    #Order the data again in order to be able to plot Condition in the same order as in the other plots
    dataAmplitudePlot$ConditionRecording = factor(dataAmplitudePlot$ConditionRecording, levels = c("NotAttRec","AttRec"))
    dataAmplitudePlot$ExpPhase = factor(dataAmplitudePlot$ExpPhase, levels = c("Bsln","Acq","Ext"))
    dataAmplitudePlot = dataAmplitudePlot[order(dataAmplitudePlot$Subject,dataAmplitudePlot$ExpPhase,dataAmplitudePlot$ConditionRecording),]}
  
  if(plottingConditions[i]=="Amplitude high reward attended"){dataAmplitudePlot=subset(data,Condition=="High Reward Attended")}
  
  if(plottingConditions[i]=="Amplitude low reward attended"){
    
    #Order the data again in order to be able to plot Condition in the same order as in the other plots
    dataAmplitudePlot=subset(data,Condition=="Low Reward Attended")
    dataAmplitudePlot$ConditionRecording = factor(dataAmplitudePlot$ConditionRecording, levels = c("NotAttRec","AttRec"))
    dataAmplitudePlot$ExpPhase = factor(dataAmplitudePlot$ExpPhase, levels = c("Bsln","Acq","Ext"))
    dataAmplitudePlot = dataAmplitudePlot[order(dataAmplitudePlot$Subject,dataAmplitudePlot$ExpPhase,dataAmplitudePlot$ConditionRecording),] #order the data again in order to be able to plot Condition in the same order as in the other plots
    }  

# Pirate plot

  pirateplot(formula=Amplitude~ConditionRecording+ExpPhase, # dependent~independent variables
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

# Set the parameters 
num.iter=10000 # number of MonteCarlo iterations (default: 10000)
# select scaling factor r of Cauchy(0,r) prior on standardized effect sizes
medprior=sqrt(2)/2 # medium prior



# Amplitude - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) X 2(high vs. low reward) 
ssVEPamp.bf.amplitude <- anovaBF(Amplitude ~ ExpPhase * Condition  *ConditionRecording + Subject,data=data,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
sort(ssVEPamp.bf.amplitude)


# Selectivity index - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) X 2(high vs. low reward) 
ssVEPamp.bf.selectivity <- anovaBF(Selectivity ~ ExpPhase * Condition + Subject,data=data.diff,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
sort(ssVEPamp.bf.selectivity)

# Total enhancement index - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) X 2(high vs. low reward) 
ssVEPamp.bf.totalenhancement <- anovaBF(TotalEnhancement ~ ExpPhase * Condition + Subject,data=data.dif,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
sort(ssVEPamp.bf.totalenhancement)


# brms------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Modelling the effects of phase, attention, and reward magnitude - All subjects

data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High Reward Attended")
data$ConditionRecording=relevel(data$ConditionRecording,ref="AttRec")

# Null model
model.null.threefactors <-brm(Amplitude ~ 1 + (1|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE)
saveRDS(model.null.threefactors,file="model.null.threefactors.EEG.allsub.rds")

# Exp phase model
model.expphase.threefactors <-brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                                  data=data,
                                  family=gaussian(),
                                  warmup = 2000,
                                  iter = 10000,
                                  save_all_pars = TRUE)
saveRDS(model.expphase.threefactors,file="model.expphase.threefactors.EEG.allsubs.rds")

# Condition model
model.condition.threefactors <-brm(Amplitude ~ Condition + (Condition|Subject),
                                   data=data,
                                   family=gaussian(),
                                   warmup = 2000,
                                   iter = 10000,
                                   save_all_pars = TRUE)
saveRDS(model.condition.threefactors,file="model.condition.threefactors.EEG.allsubs.rds")

# Attention model
model.attention.threefactors <-brm(Amplitude ~ ConditionRecording + (ConditionRecording|Subject),
                                   data=data,
                                   family=gaussian(),
                                   warmup = 2000,
                                   iter = 10000,
                                   save_all_pars = TRUE)
saveRDS(model.attention.threefactors,file="model.attention.threefactors.EEG.allsubs.rds")

# Two main effects - phase and attention
model.phaseANDattention.threefactors <-brm(Amplitude ~ ExpPhase + ConditionRecording + (ExpPhase + ConditionRecording|Subject),
                                        data=data,
                                        family=gaussian(),
                                        warmup = 2000,
                                        iter = 10000,
                                        save_all_pars = TRUE)
saveRDS(model.phaseANDattention.threefactors,file="model.phaseANDattention.threefactors.EEG.allsubs.rds")

# Two main effects - reward magnitude and attention
model.rewardmagnitudeANDattention.threefactors <-brm(Amplitude ~ Condition + ConditionRecording + (Condition + ConditionRecording|Subject),
                                           data=data,
                                           family=gaussian(),
                                           warmup = 2000,
                                           iter = 10000,
                                           save_all_pars = TRUE)
saveRDS(model.rewardmagnitudeANDattention.threefactors,file="model.rewardmagnitudeANDattention.threefactors.EEG.allsubs.rds")

# Three main effects
model.threemaineffects.threefactors <-brm(Amplitude ~ Condition + ExpPhase + ConditionRecording + (Condition + ExpPhase + ConditionRecording|Subject),
                                          data=data,
                                          family=gaussian(),
                                          warmup = 2000,
                                          iter = 10000,
                                          save_all_pars = TRUE)
saveRDS(model.threemaineffects.threefactors,file="model.threemaineffects.threefactors.EEG.allsubs.rds")

# Full model
model.full.threefactors <-brm(Amplitude ~ Condition * ExpPhase * ConditionRecording + (Condition * ExpPhase * ConditionRecording|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE)
saveRDS(model.full.threefactors,file="model.full.threefactors.EEG.allsubs.rds")


#LOO crossvalidation
compare.threefactors.EEG.loo <- LOO(model.null.threefactors,model.condition.threefactors,model.attention.threefactors,model.expphase.threefactors,model.phaseANDattention.threefactors,model.rewardmagnitudeANDattention.threefactors,model.threemaineffects.threefactors,model.full.threefactors)
saveRDS(compare.threefactors.EEG.loo,file="compare.threefactors.EEG.loo.allsubs.rds")
#WAIC
compare.threefactors.EEG.waic <- WAIC(model.null.threefactors,model.condition.threefactors,model.attention.threefactors,model.expphase.threefactors,model.phaseANDattention.threefactors,model.rewardmagnitudeANDattention.threefactors,model.threemaineffects.threefactors,model.full.threefactors)
saveRDS(compare.threefactors.EEG.waic,file="compare.threefactors.EEG.waic.allsubs.rds")

# Plot chains
plot(model.threemaineffects.threefactors, pars = parnames(model.threemaineffects.threefactors)[1:5])

# Plot parameter estimates
pairs(fit, pars = parnames(model.threemaineffects.threefactors)[1:5], exact_match = TRUE)

# Plot marginal effects for each predictor
plot(marginal_effects(model.full.threefactors),ask=FALSE)

pp_check(model.full.threefactors)

bayes_R2(model.attention.threefactors)
bayes_R2(model.condition.threefactors)
bayes_R2(model.expphase.threefactors)
bayes_R2(model.twomaineffects.threefactors)
bayes_R2(model.threemaineffects.threefactors)
bayes_R2(model.full.threefactors)

bayes_factor(model.threemaineffects.threefactors,model.attention.threefactors)
bayes_factor(model.threemaineffects.threefactors,model.null.threefactors)
bayes_factor(model.attention.threefactors,model.null.threefactors)
bayes_factor(model.twomaineffects.threefactors,model.null.threefactors)
bayes_factor(model.full.threefactors,model.null.threefactors)
bayes_factor(model.full.threefactors,model.attention.threefactors)

# read in the models and comparisons
model.null.threefactors = readRDS("model.null.threefactors.EEG.allsubs.rds")
model.condition.threefactors = readRDS("model.condition.threefactors.EEG.allsubs.rds")
model.attention.threefactors = readRDS("model.attention.threefactors.EEG.allsubs.rds")
model.expphase.threefactors = readRDS("model.expphase.threefactors.EEG.rds")
model.twomaineffects.threefactors = readRDS("model.twomaineffects.threefactors.EEG.allsubs.rds")
model.threemaineffects.threefactors = readRDS("model.threemaineffects.threefactors.EEG.allsubs.rds")
model.full.threefactors = readRDS("model.full.threefactors.EEG.allsubs.rds")
compare.threefactors.EEG.loo = readRDS("compare.threefactors.EEG.loo.allsubs.rds")
compare.threefactors.EEG.waic = readRDS("compare.threefactors.EEG.waic.allsubs.rds")
