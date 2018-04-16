# experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Gilles Pourtois, Ernst Koster, & S?ren Andersen)
# (*: co-first authors)
# code contributions: Antonio Schettino, Ivan Grahek
# analysis of behavioral data
#

##### Importing data & first steps ####
# dev.off() # clear plots
rm(list=ls()) # clear environment
cat("\014") # clear console

# load packages and install them if they're not installed. the pacman package will
# automatically check for the requested packages and download & install them if they are not on the computer.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(reshape2,yarrr,BayesFactor,plyr,ez,schoRsch,brms,lme4)

set.seed(42) # the Answer to the Ultimate Question of Life, the Universe and Everything
wd <- ("C:/Users/igrahek/Desktop/All participants/FSAReward-master/Exp1_fullsample/") # Antonio's work directory
# wd <- ("C:/Users/igrahek/Google Drive/Work computer/FSAReward/Analysis/EEG/Exp1") # Ivan's work directory
setwd(wd) # set work directory
#data.raw <- read.csv2("amplitudes_rewardBoth.csv",header=TRUE,na.strings="NaN") # read data #only good behavior: amplitudes_rewardBoth.csv  amplitudes_rewardBoth_wholeSample.csv

data.raw <- read.csv2("amplitudes_rewardBoth_wholeSample.csv",header=TRUE,na.strings="NaN") # read data #only good behavior: amplitudes_rewardBoth.csv  amplitudes_rewardBoth_wholeSample.csv


# Reshape to long format
data <- melt(data.raw,id.vars=c("Subject","Frequency"),
             measure.vars=c("BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended"),
             variable.name="Condition",value.name="Amplitude") 

# data1 <- melt(data.raw1,id.vars=c("Subject","Frequency"),
#              measure.vars=c("BslnRedAttended","BslnBlueAttended","AcqRedAttended","AcqBlueAttended","ExtRedAttended","ExtBlueAttended"),
#              variable.name="Condition",value.name="Amplitude") 

# #Exclude participants without full behavioral RT data
# behavioral.data = c(1,3,7,8,9,11,12,16,18,19,21,22,23,28,29,30,32,33,36,37,40,42,43,46,47,48) # 26 participants, criterion is 0.6
# # behavioral.data = c(1,  3,  4,  6,  7,  8,  9,  11, 12, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 28, 29, 30, 32, 33, 36, 37, 38, 40, 41, 42, 43, 45, 46, 47, 48) # 37 participants, criterion is 0.5
#behavioral.data = c(4,8,12,11,14,15,17,18,20,31,34,41) # eclude people with at least one condition below .5 ##4,8,12 are excluded to have an equal number of blue and red
#data = data[!data$Subject %in% behavioral.data, ]

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

# Z score the amplitudes in each Subject, Frequency, and Condition 
#data$Amplitude = (data$Amplitude-data$MeanAmplitude)/data$SDAmplitude

# Calculate the attention indexes - Selectivity (attended-unattended) & total enhancement (attended+unattended) (Andersen & Muller, 2010, PNAS)
data.diff = ddply(data, .(Subject,ExpPhase,Condition), transform, Selectivity = Amplitude[ConditionRecording=="AttRec"]-Amplitude[ConditionRecording=="NotAttRec"],TotalEnhancement=Amplitude[ConditionRecording=="AttRec"]+Amplitude[ConditionRecording=="NotAttRec"])
# Delete the ConditionRecording column and rows which are not necessary (indexes repeated twice)
data.diff = subset(data.diff,ConditionRecording=="AttRec") #keep only AttRec as it is equal to NotAttRec
data.diff$ConditionRecording = NULL  #drop the ConditionRecording column

# Sort the data 
data.diff$ExpPhase = factor(data.diff$ExpPhase, levels = c("Bsln","Acq","Ext"))
data.diff = data.diff[order(data.diff$Subject,data.diff$Condition,data.diff$ExpPhase),]


##### PLOT ATTENDED vs. NOT ATTENDED #####


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
  tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image
  
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
  
  dev.off()
}

##### PLOT EXP PHASE vs. REWARD SELECTIVITY INDEX#####
plottingConditions = c("Selectivity both colors","Selectivity blue rewarded","Selectivity red rewarded" )
for (i in 1:length(plottingConditions)){
  if(plottingConditions[i]=="Selectivity both colors"){dataAmplitudePlot=data.diff}
  if(plottingConditions[i]=="Selectivity blue rewarded"){dataAmplitudePlot=subset(data.diff,RewardedColor=="Blue")}
  if(plottingConditions[i]=="Selectivity red rewarded"){dataAmplitudePlot=subset(data.diff,RewardedColor=="Red")}
  

# Pirate plot
tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image

pirateplot(formula=Selectivity~Condition+ExpPhase, # dependent~independent variables
           data=dataAmplitudePlot, # data frame
           main=plottingConditions[i], # main title
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

dev.off()
}

##### PLOT EXP PHASE vs. REWARD TOTAL ENHANCEMENT INDEX#####
  plottingConditions = c("Total Enhancement both colors","Total Enhancement blue rewarded","Total Enhancement red rewarded" )
  for (i in 1:length(plottingConditions)){
    if(plottingConditions[i]=="Total Enhancement both colors"){dataAmplitudePlot=data.diff}
    if(plottingConditions[i]=="Total Enhancement blue rewarded"){dataAmplitudePlot=subset(data.diff,RewardedColor=="Blue")}
    if(plottingConditions[i]=="Total Enhancement red rewarded"){dataAmplitudePlot=subset(data.diff,RewardedColor=="Red")}
    
    
    # Pirate plot
    tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image
    
    pirateplot(formula=TotalEnhancement~ExpPhase+Condition, # dependent~independent variables
               data=dataAmplitudePlot, # data frame
               main=plottingConditions[i], # main title
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
    
    dev.off()
  }
##### STATS ####

# Set the parameters 
data$Subject <- as.factor(data$Subject) # convert to factor
data$Condition <- as.factor(data$Condition) # convert to factor
data$ExpPhase <- as.factor(data$ExpPhase) # convert to factor
data$RewardedColor <- as.factor(data$RewardedColor) # convert to factor
data$ConditionRecording <- as.factor(data$ConditionRecording) # convert to factor
data$RecordingAndCondition <- as.factor(data$RecordingAndCondition) # convert to factor
data.diff$Subject <- as.factor(data.diff$Subject) # convert to factor
data.diff$Condition <- as.factor(data.diff$Condition) # convert to factor
data.diff$ExpPhase <- as.factor(data.diff$ExpPhase) # convert to factor
data.diff$RecordingAndCondition <- as.factor(data.diff$RecordingAndCondition) # convert to factor
num.iter=10000 # number of MonteCarlo iterations (default: 10000)
# select scaling factor r of Cauchy(0,r) prior on standardized effect sizes
medprior=sqrt(2)/2 # medium prior



# Test the effect of attention - Bayesian ANOVA 3(Exp Phase) X 2(Attended recorded vs. NotAttended recorded) 

# medium prior
# models: 
# numerator
#         (1) main effect of ExpPhase
#         (2) main effect of Condition Recording
#         (3) both main effects
#         (4) full (main effects & interaction)
# denominator
#         (0) null model (only random effect of Subject, considered as nuisance)
dataEffectOfAttention = ddply(data,.(Subject,ExpPhase,ConditionRecording),summarize,
             Amplitude1 = mean(Amplitude,na.rm=TRUE),
             SD =   sd(Amplitude,na.rm=TRUE))
dataEffectOfAttention$Amplitude=dataEffectOfAttention$Amplitude1
dataEffectOfAttention$ExpPhase <- as.factor(dataEffectOfAttention$ExpPhase) # convert to factor
dataEffectOfAttention$ConditionRecording <- as.factor(dataEffectOfAttention$ConditionRecording) # convert to factor
dataEffectOfAttention$Subject <- as.factor(dataEffectOfAttention$Subject) # convert to factor

ssVEPamp.bf.attention <- anovaBF(Amplitude~ExpPhase*ConditionRecording+Subject,data=dataEffectOfAttention,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
ssVEPamp.bf.attention

ssVEPamp.bf.attention.highreward <- anovaBF(Amplitude~ExpPhase*ConditionRecording+Subject,data=subset(data,Condition=="High Reward Attended"),iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
ssVEPamp.bf.attention.highreward

ssVEPamp.bf.attention.lowreward <- anovaBF(Amplitude~ExpPhase*ConditionRecording+Subject,data=subset(data,Condition=="Low Reward Attended"),iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
ssVEPamp.bf.attention.lowreward

ssVEPamp.bf.attention <- anovaBF(Amplitude~ExpPhase*Condition*ConditionRecording+Subject,data=data,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
sort(ssVEPamp.bf.attention)


ssVEPamp.bf.attention <- anovaBF(Amplitude~ExpPhase*Condition*ConditionRecording+Subject,data=data,iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
sort(ssVEPamp.bf.attention)

# length(unique(subset(data,RewardedColor=="Red")$Subject))
# length(unique(subset(data,RewardedColor=="Blue")$Subject))
# average over rewarded color
# analysis
ssVEPamp.bf.attention.acquisition <- anovaBF(Amplitude~Condition*ConditionRecording+Subject,data=subset(data,ExpPhase=="Acq"),iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
ssVEPamp.bf.attention.acquisition

ssVEPamp.bf.attention.acquisition <- anovaBF(Amplitude~Condition*ConditionRecording+Subject,data=subset(data,ExpPhase=="Acq"&RewardedColor=="Blue"),iterations=num.iter,whichRandom="Subject",rscaleRandom="nuisance",rscaleFixed=medprior)
ssVEPamp.bf.attention.acquisition

#Export the effect of attention averaged over reward conditions
Amplitude=subset(dataEffectOfAttention,select=c(Subject,ExpPhase,ConditionRecording,Amplitude))
Amplitude=dcast(Amplitude,Subject~ExpPhase+ConditionRecording, value.var="Amplitude")
write.csv2(Amplitude,file="Amplitude Averaged over reward.csv")

#Export the effect of attention for high reward attended trials
Amplitude=subset(subset(data,Condition=="High Reward Attended"),select=c(Subject,ExpPhase,ConditionRecording,Amplitude))
Amplitude=dcast(Amplitude,Subject~ExpPhase+ConditionRecording, value.var="Amplitude")
write.csv2(Amplitude,file="AmplitudeHighRewardAtt.csv")

#Export the effect of attention for low reward attended trials
Amplitude=subset(subset(data,Condition=="Low Reward Attended"),select=c(Subject,ExpPhase,ConditionRecording,Amplitude))
Amplitude=dcast(Amplitude,Subject~ExpPhase+ConditionRecording, value.var="Amplitude")
write.csv2(Amplitude,file="AmplitudeLowRewardAtt.csv")

#Export the data for aquisition
Amplitude=subset(subset(data,ExpPhase=="Acq"),select=c(Subject,Condition,ConditionRecording,Amplitude))
Amplitude=dcast(Amplitude,Subject~Condition+ConditionRecording, value.var="Amplitude")
write.csv2(Amplitude,file="OnlyAcquisition.csv")

#Export the data for aquisition
Amplitude=subset(subset(data.diff,ExpPhase=="Acq"),select=c(Subject,Condition,Selectivity))
Amplitude=dcast(Amplitude,Subject~Condition, value.var="Selectivity")
write.csv2(Amplitude,file="OnlyAcquisitionSelectivity.csv")


# Export data
Amplitude=dcast(subset(data,RewardedColor=="Red"),Subject~ExpPhase+Condition+ConditionRecording, value.var="Amplitude")
write.csv2(Amplitude,file="AmplitudeAllConditionsRedRewarded.csv")

# Export data
Amplitude=dcast(data,Subject~ExpPhase+Condition+ConditionRecording, value.var="Amplitude")
write.csv2(Amplitude,file="AmplitudeAllConditionsPoorPerformenrs.csv")

#### lmBF ####

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Only random effect of subject - Attention & reward block ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Null
bfNull = lmBF(Amplitude ~ 1 + Subject, data = dataEffectOfAttention, whichRandom = "Subject")
## Main effect of ExpPhase 
bfExpPhase = lmBF(Amplitude ~ ExpPhase + Subject, data = dataEffectOfAttention, whichRandom = "Subject")
## Main effect of ConditionRecording 
bfConditionRecording = lmBF(Amplitude ~ ConditionRecording + Subject, data = dataEffectOfAttention, whichRandom = "Subject")
## Two main effects
bfTwoMain = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject, data = dataEffectOfAttention, whichRandom = "Subject")
## Interaction
bfInteraction = lmBF(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase:ConditionRecording + Subject, data = dataEffectOfAttention, whichRandom = "Subject")

bfConditionRecording/bfInteraction
bfExpPhase/bfNull
bfTwoMain/bfNull
bfInteraction/bfNull

(bfConditionRecording/bfNull)/(bfExpPhase/bfNull)
(bfConditionRecording/bfNull)/(bfExpPhase/bfNull)
(bfConditionRecording/bfNull)/(bfTwoMain/bfNull)
(bfConditionRecording/bfNull)/(bfInteraction/bfNull)

chains = posterior(bfInteraction, iterations = 100000)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### All Random effects - Attention & reward block ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Null 
bfNull = lmBF(Amplitude ~ 1 + Subject, data = dataEffectOfAttention, whichRandom = "Subject")
## Main effect of ConditionRecording 
bfConditionRecording.rand = lmBF(Amplitude ~ ConditionRecording + Subject + ConditionRecording:Subject, data = dataEffectOfAttention, whichRandom = c("Subject","ConditionRecording:Subject"))
## Main effect of ExpPhase 
bfExpPhase.rand = lmBF(Amplitude ~ ExpPhase + Subject + ExpPhase:Subject, data = dataEffectOfAttention, whichRandom = c("Subject","ExpPhase:Subject"))
## Main effects 
bfTwoMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject + ExpPhase:Subject + ConditionRecording:Subject, data = dataEffectOfAttention, whichRandom = c("Subject","ExpPhase:Subject", "ConditionRecording:Subject"))
## Interaction 
bfInteraction.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase:ConditionRecording + Subject + ExpPhase:ConditionRecording:Subject, data = dataEffectOfAttention, whichRandom = c("Subject","ExpPhase:ConditionRecording:Subject"))

bfConditionRecording.rand/bfNull
bfExpPhase.rand/bfNull
bfTwoMain.rand/bfNull
bfInteraction.rand/bfNull

(bfConditionRecording.rand/bfNull)/(bfExpPhase.rand/bfNull)
(bfConditionRecording.rand/bfNull)/(bfExpPhase.rand/bfNull)
(bfConditionRecording.rand/bfNull)/(bfTwoMain.rand/bfNull)
(bfConditionRecording.rand/bfNull)/(bfInteraction.rand/bfNull)

chains = posterior(bfInteraction.rand, iterations = 10000)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Only random effects of subject - Attention, reward block, and reward magnitude ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## Null
bfNull = lmBF(Amplitude ~ 1 + Subject, data = data, whichRandom = "Subject")
## Attention
bfAttention.rand = lmBF(Amplitude ~ ConditionRecording + Subject, data = data, whichRandom = "Subject")
## Two Main effects 
bfTwoMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject, data = data, whichRandom = "Subject")
## Three Main effects 
bfThreeMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Condition + Subject, data = data, whichRandom = "Subject")
## Interaction 
bfInteractionThree.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Condition + ExpPhase:ConditionRecording:Condition + Subject, data = data, whichRandom = "Subject")

bfAttention.rand/bfThreeMain.rand
bfTwoMain.rand/bfNull
bfThreeMain.rand/bfNull
bfInteractionThree.rand/bfNull

(bfAttention.rand/bfNull)/(bfTwoMain.rand/bfNull)
(bfAttention.rand/bfNull)/(bfThreeMain.rand/bfNull)
(bfAttention.rand/bfNull)/(bfInteractionThree.rand/bfNull)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#### All random effects - Attention, reward block, and reward magnitude ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Null
bfNull = lmBF(Amplitude ~ 1 + Subject, data = data, whichRandom = "Subject")
## Attention
bfAttention.rand = lmBF(Amplitude ~ ConditionRecording + Subject + ConditionRecording:Subject, data = data, whichRandom = c("Subject","ConditionRecording:Subject"))
## Two Main effects 
bfTwoMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject + ExpPhase:Subject + ConditionRecording:Subject + ExpPhase*Subject, data = data, whichRandom = c("Subject","ExpPhase:Subject","ConditionRecording:Subject"))
## Three Main effects 
bfThreeMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Condition + Subject + ExpPhase:Subject + Condition:Subject + ConditionRecording:Subject, data = data, whichRandom = c("Subject","ExpPhase:Subject","ConditionRecording:Subject","Condition:Subject"))
## Interaction 
bfInteractionThree.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Condition + ExpPhase:ConditionRecording:Condition + Subject + ExpPhase:ConditionRecording:Condition:Subject, data = data, whichRandom = c("Subject","ExpPhase:ConditionRecording:Condition:Subject"))



bfAttention.rand/bfThreeMain.rand
bfTwoMain.rand/bfNull
bfThreeMain.rand/bfNull
bfInteractionThree.rand/bfNull

(bfThreeMain.rand/bfNull)/(bfAttention.rand/bfNull)
(bfThreeMain.rand/bfNull)/(bfTwoMain.rand/bfNull)
(bfThreeMain.rand/bfNull)/(bfInteractionThree.rand/bfNull)


bfThreeMain.rand/bfAttention.rand
chains = posterior(bfThreeMain.rand, iterations = 10000)
plot(chains[,2:3])
summary(chains)


#### Random effect of Condition - Attention, reward block, and reward magnitude ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
## Null
bfNull = lmBF(Amplitude ~ 1 + Subject, data = data, whichRandom = "Subject")
## Attention
bfAttention.rand = lmBF(Amplitude ~ ConditionRecording + Subject + Condition:Subject, data = data, whichRandom = c("Subject","Condition:Subject"))
## Two Main effects 
bfTwoMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject + Condition:Subject + ExpPhase*Subject, data = data, whichRandom = c("Subject","Condition:Subject"))
## Three Main effects 
bfThreeMain.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Condition + Subject + Condition:Subject, data = data, whichRandom = c("Subject","Condition:Subject"))
## Interaction 
bfInteractionThree.rand = lmBF(Amplitude ~ ExpPhase + ConditionRecording + Condition + Subject + Condition:Subject, data = data, whichRandom = c("Subject","Condition:Subject"))



bfAttention.rand/bfThreeMain.rand
bfTwoMain.rand/bfNull
bfThreeMain.rand/bfNull
bfInteractionThree.rand/bfNull

(bfThreeMain.rand/bfNull)/(bfAttention.rand/bfNull)
(bfThreeMain.rand/bfNull)/(bfTwoMain.rand/bfNull)
(bfThreeMain.rand/bfNull)/(bfInteractionThree.rand/bfNull)

#### Sensitivity analysis ####
# priors for fixed effects
#create a list of priors
prior=seq(0.1,1,0.1)
#initialize a data frame to store the results
sensitivity = as.data.frame(prior)

# Analysis of the main effects of attention and phase

dataEffectOfAttention = ddply(data,.(Subject,ExpPhase,ConditionRecording),summarize,
                              Amplitude1 = mean(Amplitude,na.rm=TRUE),
                              SD =   sd(Amplitude,na.rm=TRUE))
dataEffectOfAttention$Amplitude=dataEffectOfAttention$Amplitude1
dataEffectOfAttention$ExpPhase <- as.factor(dataEffectOfAttention$ExpPhase) # convert to factor
dataEffectOfAttention$ConditionRecording <- as.factor(dataEffectOfAttention$ConditionRecording) # convert to factor
dataEffectOfAttention$Subject <- as.factor(dataEffectOfAttention$Subject) # convert to factor

# Sensitivity

data.selectivity=dataEffectOfAttention
sensitivity.1 = ddply(sensitivity,.(prior),summarise,
                      #Exp Phase versus the Null model
                      ExpPhase.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #ConditionRecording versus the Null model
                      ConditionRecording.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #Main effects versus the Null model
                      MainEffects.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #Full model versus the Null model
                      Full.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase:ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #Main effects versus the Exp Phase model
                      MainEffects.vs.Full.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase:ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf)


# Analysis separated by reward magnitude
data.selectivity=subset(data,Condition=="Low Reward Attended")
# calculate BFs for RTs

sensitivity.1 = ddply(sensitivity,.(prior),summarise,
                      #Exp Phase versus the Null model
                      ExpPhase.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #ConditionRecording versus the Null model
                      ConditionRecording.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #Main effects versus the Null model
                      MainEffects.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #Full model versus the Null model
                      Full.vs.Null.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase:ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf,
                      #Main effects versus the Exp Phase model
                      MainEffects.vs.Full.model=extractBF(
                        lmBF(Amplitude ~ ExpPhase + ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior)/
                          lmBF(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase:ConditionRecording + Subject, whichRandom = "Subject",data = data.selectivity,rscaleFixed = prior))$bf)

plot(prior,sensitivity.1$Condition.vs.Null.model)

#### brms ####
dataEffectOfAttention = ddply(data,.(Subject,ExpPhase,ConditionRecording),summarize,
                              Amplitude1 = mean(Amplitude,na.rm=TRUE),
                              SD =   sd(Amplitude,na.rm=TRUE))
dataEffectOfAttention$Amplitude=dataEffectOfAttention$Amplitude1
dataEffectOfAttention$ExpPhase <- as.factor(dataEffectOfAttention$ExpPhase) # convert to factor
dataEffectOfAttention$ConditionRecording <- as.factor(dataEffectOfAttention$ConditionRecording) # convert to factor
dataEffectOfAttention$Subject <- as.factor(dataEffectOfAttention$Subject) # convert to factor

#lme4

m.null = lmer(Amplitude~1+(1|Subject),dataEffectOfAttention)
m.expphase = lmer(Amplitude~ExpPhase+(1|Subject),dataEffectOfAttention)
m.attention = lmer(Amplitude~ConditionRecording+(1|Subject),dataEffectOfAttention)
m.twomaineffects = lmer(Amplitude~ConditionRecording+ExpPhase+(1|Subject),dataEffectOfAttention)
m.interaction = lmer(Amplitude~ConditionRecording*ExpPhase+(1|Subject),dataEffectOfAttention)
anova(m.null, m.expphase, m.attention, m.twomaineffects,m.interaction)



m.full = lmer(Amplitude~ExpPhase*ConditionRecording+(1|Subject),dataEffectOfAttention)
m.full.rand.attention=lmer(Amplitude~ExpPhase*ConditionRecording+(ConditionRecording|Subject),dataEffectOfAttention)
m.full.rand.expphase=lmer(Amplitude~ExpPhase*ConditionRecording+(ExpPhase|Subject),dataEffectOfAttention)
m.full.rand.both=lmer(Amplitude~ExpPhase*ConditionRecording+(ConditionRecording*ExpPhase|Subject),dataEffectOfAttention)
anova(m.null, m.expphase, m.attention, m.twomaineffects) 

#lme4 all factors mixed

m.null = lmer(Amplitude~1+(1|Subject),data)
m.expphase = lmer(Amplitude ~ ExpPhase +(1|Subject), data)
m.attention = lmer(Amplitude ~ ConditionRecording + (1|Subject), data)
m.twomaineffects = lmer(Amplitude ~ ConditionRecording + Condition +(1|Subject),data)
m.twomaineffects1 = lmer(Amplitude ~ Condition + ExpPhase +(1|Subject),data)
m.twomaineffects2 = lmer(Amplitude ~ ConditionRecording + ExpPhase +(1|Subject),data)
m.threemaineffects = lmer(Amplitude ~ ConditionRecording + ExpPhase + Condition + (1|Subject),data)
m.full = lmer(Amplitude ~ ExpPhase * ConditionRecording * Condition + (1|Subject),data)

anova(m.null,m.expphase,m.attention,m.twomaineffects,m.twomaineffects1,m.twomaineffects2,m.threemaineffects,m.full)

m.full2 = lmer(Amplitude~ExpPhase*ConditionRecording*Condition+(RewardedColor*Condition|Subject),data)
m.full3 = lmer(Amplitude~ExpPhase*ConditionRecording*Condition+(Condition|Subject),data)
anova(m.null, m.expphase,m.attention,m.twomaineffects,m.threemaineffects,m.full,m.full2,m.full3) 

m.twomaineffects.rand = lmer(Amplitude~ConditionRecording+ExpPhase+(RewardedColor|Subject),data)
m.twomaineffects = lmer(Amplitude~ConditionRecording+ExpPhase+(1|Subject),data)
anova(m.twomaineffects.rand,m.twomaineffects)

tmp <- as.data.frame(confint(glht(m.threemaineffects))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()
#### Modelling the effects of phase and attention #### 
# referencing for easier interpretation
dataEffectOfAttention$ExpPhase=relevel(data.summaryRT.wide.ssj$ExpPhase,ref="Acq")
dataEffectOfAttention$Condition=relevel(data.summaryRT.wide.ssj$Condition,ref="AttRec")

# Null model
model.null <-brm(Amplitude~1+(1|Subject),
                 data=dataEffectOfAttention,
                 family=gaussian(),
                 warmup = 500, #200 & 10000
                 iter = 2000)
saveRDS(model.null,file="model.null.EEG.rds")
# Main effect of Exp Phase
model.expphase <-brm(Amplitude ~ ExpPhase + (1|Subject),
                 data=dataEffectOfAttention,
                 family=gaussian(),
                 warmup = 500,
                 iter = 2000)
saveRDS(model.expphase,file="model.expphase.EEG.rds")
# Main effect of attention
model.attention <-brm(Amplitude ~ ConditionRecording + (1|Subject),
                     data=dataEffectOfAttention,
                     family=gaussian(),
                     warmup = 500,
                     iter = 2000)
saveRDS(model.attention,file="model.attention.EEG.rds")
# Two main effects
model.twomaineffects <-brm(Amplitude ~ ExpPhase + ConditionRecording + (1|Subject),
                           data=dataEffectOfAttention,
                           family=gaussian(),
                           warmup = 500,
                           iter = 2000)
saveRDS(model.twomaineffects,file="model.twomaineffects.EEG.rds")
# Full model
model.full <-brm(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase * ConditionRecording + (1|Subject),
                           data=dataEffectOfAttention,
                           family=gaussian(),
                           warmup = 500,
                           iter = 2000)
saveRDS(model.full,file="model.full.EEG.rds")

# Full model with random effects
model.full.rand <-brm(Amplitude ~ ExpPhase + ConditionRecording + ExpPhase * ConditionRecording + (ExpPhase * ConditionRecording|Subject),
                 data=dataEffectOfAttention,
                 family=gaussian(),
                 warmup = 500,
                 iter = 2000)
saveRDS(model.full.rand,file="model.full.rand.EEG.rds")



# read in the models and comparisons
model.null <- readRDS("model.null.EEG.rds")
model.expphase <- readRDS("model.expphase.EEG.rds")
model.attention <- readRDS("model.attention.EEG.rds")
model.twomaineffects <- readRDS("model.twomaineffects.EEG.rds")
model.full <- readRDS("model.full.EEG.rds")
model.full.rand <- readRDS("model.full.rand.EEG.rds")
compare.loo = readRDS("compare.EEG.loo.rds")
compare.waic = readRDS("compare.EEG.waic.rds")

launch_shinystan(model.full.rand)

# Plot marginal effects for each predictor
plot(marginal_effects(model.full.rand),ask=FALSE)

# dataEffectOfAttention$ConditionRecording = ifelse(ConditionRecording=="AttRec",-.5,.5)
# dataEffectOfAttention$ExpPhase = ifelse(ExpPhase=="Bsln",-.5,.5)


#LOO crossvalidation
compare.EEG.loo <- LOO(model.null,model.expphase,model.attention,model.twomaineffects,model.full,reloo=TRUE)
saveRDS(compare.EEG.loo,file="compare.EEG.loo")
#WAIC
compare.EEG.waic <- WAIC(model.null,model.expphase,model.attention,model.twomaineffects,model.full,model.full.rand)
saveRDS(compare.EEG.waic,file="compare.EEG.waic")

# Plot marginal effects for each predictor
plot(marginal_effects(model.full),ask=FALSE)

#Save data
data.brms=subset(subset(dataEffectOfAttention,select=c(Subject,ExpPhase,ConditionRecording,Amplitude)))
write.csv2(data.brms,file="Amplitudes.brms.csv")


# library(tidyverse)
# model.null %>%
#   plot(
#     combo = c("hist", "trace"), widths = c(1, 1.5),
#     theme = theme_bw(base_size = 10)
#   )


#### Modelling the effects of phase, attention, and reward magnitude ####
library(brms)

data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High Reward Attended")
data$ConditionRecording=relevel(data$ConditionRecording,ref="AttRec")

# Null model
model.null.threefactors <-brm(Amplitude ~ 1 + (1|Subject),
                 data=data,
                 family=gaussian(),
                 save_all_pars = TRUE)
saveRDS(model.null.threefactors,file="model.null.threefactors.EEG.rds")

# Exp phase model
model.expphase.threefactors <-brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                              data=data,
                              family=gaussian(),
                              save_all_pars = TRUE)
saveRDS(model.expphase.threefactors,file="model.expphase.threefactors.EEG.rds")

# Condition model
model.condition.threefactors <-brm(Amplitude ~ Condition + (Condition|Subject),
                                  data=data,
                                  family=gaussian(),
                                  save_all_pars = TRUE)
saveRDS(model.condition.threefactors,file="model.condition.threefactors.EEG.rds")

# Attention model
model.attention.threefactors <-brm(Amplitude ~ ConditionRecording + (ConditionRecording|Subject),
                                  data=data,
                                  family=gaussian(),
                                  save_all_pars = TRUE)
saveRDS(model.attention.threefactors,file="model.attention.threefactors.EEG.rds")

# Two main effects
model.twomaineffects.threefactors <-brm(Amplitude ~ ExpPhase + ConditionRecording + (ExpPhase + ConditionRecording|Subject),
                           data=data,
                           family=gaussian(),
                           save_all_pars = TRUE)
saveRDS(model.twomaineffects.threefactors,file="model.twomaineffects.threefactors.EEG.rds")

# Three main effects
model.threemaineffects.threefactors <-brm(Amplitude ~ Condition + ExpPhase + ConditionRecording + (Condition + ExpPhase + ConditionRecording|Subject),
                           data=data,
                           family=gaussian(),
                           save_all_pars = TRUE)
saveRDS(model.threemaineffects.threefactors,file="model.threemaineffects.threefactors.EEG.rds")

# Full model
model.full.threefactors <-brm(Amplitude ~ Condition + ExpPhase + ConditionRecording + Condition * ExpPhase * ConditionRecording + (Condition * ExpPhase * ConditionRecording|Subject),
                 data=data,
                 family=gaussian(),
                 save_all_pars = TRUE)
saveRDS(model.full.threefactors,file="model.full.threefactors.EEG.rds")


#LOO crossvalidation
compare.threefactors.EEG.loo <- LOO(model.null.threefactors,model.condition.threefactors,model.attention.threefactors,model.expphase.threefactors,model.twomaineffects.threefactors,model.threemaineffects.threefactors,model.full.threefactors)
saveRDS(compare.threefactors.EEG.loo,file="compare.threefactors.EEG.loo.rds")
#WAIC
compare.threefactors.EEG.waic <- WAIC(model.null.threefactors,model.condition.threefactors,model.attention.threefactors,model.expphase.threefactors,model.twomaineffects.threefactors,model.threemaineffects.threefactors,model.full.threefactors)
saveRDS(compare.threefactors.EEG.waic,file="compare.threefactors.EEG.waic.rds")

# Plot chains
plot(model.threemaineffects.threefactors, pars = parnames(model.threemaineffects.threefactors)[1:5])

# Plot parameter estimates
pairs(fit, pars = parnames(model.threemaineffects.threefactors)[1:5], exact_match = TRUE)

# Plot marginal effects for each predictor
plot(marginal_effects(model.threemaineffects.threefactors),ask=FALSE)

pp_check(model.threemaineffects.threefactors)

pp_check(model.threemaineffects.threefactors,nsamples=NULL,type="stat_2d",group=c("ConditionRecording"))

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
model.null.threefactors <- readRDS("model.null.threefactors.EEG.rds")
model.condition.threefactors <- readRDS("model.condition.threefactors.EEG.rds")
model.attention.threefactors <- readRDS("model.attention.threefactors.EEG.rds")
model.expphase.threefactors <- readRDS("model.expphase.threefactors.EEG.rds")
model.twomaineffects.threefactors <- readRDS("model.twomaineffects.threefactors.EEG.rds")
model.threemaineffects.threefactors = readRDS("model.threemaineffects.threefactors.EEG.rds")
model.full.threefactors <- readRDS("model.full.threefactors.EEG.rds")
compare.threefactors.EEG.loo = readRDS("compare.threefactors.EEG.loo.rds")
compare.threefactors.EEG.waic = readRDS("compare.threefactors.EEG.waic.rds")

#### Modelling the effects of phase, attention, and reward magnitude - All subjects ####
library(brms)

data$ExpPhase=relevel(data$ExpPhase,ref="Bsln")
data$Condition=relevel(data$Condition,ref="High Reward Attended")
data$ConditionRecording=relevel(data$ConditionRecording,ref="AttRec")

# Null model
model.null.threefactors <-brm(Amplitude ~ 1 + (1|Subject),
                              data=data,
                              family=gaussian(),
                              save_all_pars = TRUE)
saveRDS(model.null.threefactors,file="model.null.threefactors.EEG.allsub.rds")

# Exp phase model
model.expphase.threefactors <-brm(Amplitude ~ ExpPhase + (ExpPhase|Subject),
                                  data=data,
                                  family=gaussian(),
                                  save_all_pars = TRUE)
saveRDS(model.expphase.threefactors,file="model.expphase.threefactors.EEG.allsubs.rds")

# Condition model
model.condition.threefactors <-brm(Amplitude ~ Condition + (Condition|Subject),
                                   data=data,
                                   family=gaussian(),
                                   save_all_pars = TRUE)
saveRDS(model.condition.threefactors,file="model.condition.threefactors.EEG.allsubs.rds")

# Attention model
model.attention.threefactors <-brm(Amplitude ~ ConditionRecording + (ConditionRecording|Subject),
                                   data=data,
                                   family=gaussian(),
                                   save_all_pars = TRUE)
saveRDS(model.attention.threefactors,file="model.attention.threefactors.EEG.allsubs.rds")

# Two main effects
model.twomaineffects.threefactors <-brm(Amplitude ~ ExpPhase + ConditionRecording + (ExpPhase + ConditionRecording|Subject),
                                        data=data,
                                        family=gaussian(),
                                        save_all_pars = TRUE)
saveRDS(model.twomaineffects.threefactors,file="model.twomaineffects.threefactors.EEG.allsubs.rds")

# Three main effects
model.threemaineffects.threefactors <-brm(Amplitude ~ Condition + ExpPhase + ConditionRecording + (Condition + ExpPhase + ConditionRecording|Subject),
                                          data=data,
                                          family=gaussian(),
                                          save_all_pars = TRUE)
saveRDS(model.threemaineffects.threefactors,file="model.threemaineffects.threefactors.EEG.allsubs.rds")

# Full model
model.full.threefactors <-brm(Amplitude ~ Condition + ExpPhase + ConditionRecording + Condition * ExpPhase * ConditionRecording + (Condition * ExpPhase * ConditionRecording|Subject),
                              data=data,
                              family=gaussian(),
                              warmup = 2000,
                              iter = 10000,
                              save_all_pars = TRUE)
saveRDS(model.full.threefactors,file="model.full.threefactors.EEG.allsubs.rds")


#LOO crossvalidation
compare.threefactors.EEG.loo <- LOO(model.null.threefactors,model.condition.threefactors,model.attention.threefactors,model.expphase.threefactors,model.twomaineffects.threefactors,model.threemaineffects.threefactors,model.full.threefactors)
saveRDS(compare.threefactors.EEG.loo,file="compare.threefactors.EEG.loo.allsubs.rds")
#WAIC
compare.threefactors.EEG.waic <- WAIC(model.null.threefactors,model.condition.threefactors,model.attention.threefactors,model.expphase.threefactors,model.twomaineffects.threefactors,model.threemaineffects.threefactors,model.full.threefactors)
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
model.null.threefactors <- readRDS("model.null.threefactors.EEG.allsubs.rds")
model.condition.threefactors <- readRDS("model.condition.threefactors.EEG.allsubs.rds")
model.attention.threefactors <- readRDS("model.attention.threefactors.EEG.allsubs.rds")
model.expphase.threefactors <- readRDS("model.expphase.threefactors.EEG.rds")
model.twomaineffects.threefactors <- readRDS("model.twomaineffects.threefactors.EEG.allsubs.rds")
model.threemaineffects.threefactors = readRDS("model.threemaineffects.threefactors.EEG.allsubs.rds")
model.full.threefactors <- readRDS("model.full.threefactors.EEG.allsubs.rds")
compare.threefactors.EEG.loo = readRDS("compare.threefactors.EEG.loo.allsubs.rds")
compare.threefactors.EEG.waic = readRDS("compare.threefactors.EEG.waic.allsubs.rds")
