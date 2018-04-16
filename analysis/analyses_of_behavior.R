# experiment: FSAReward (Ivan Grahek*, Antonio Schettino*, Gilles Pourtois, Ernst Koster, & Søren Andersen)
# (*: co-first authors)
# code contributions: Antonio Schettino, Ivan Grahek
# analysis of behavioral data
#

##### Importing data & first steps ####
# dev.off() # clear plots
rm(list=ls()) # clear environment
cat("\014") # clear console

#load packages and install them if they're not installed. the pacman package will 
#automatically check for the requested packages and download & install them if they are not on the computer.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,Rmisc,yarrr,BayesFactor,reshape2)

set.seed(42) # the Answer to the Ultimate Question of Life, the Universe and Everything
# wd <- ("D:/OneDrive - UGent/FSAReward/Analysis/Behavior/Exp1/") # Antonio's work directory
wd <- ("C:/Users/igrahek/Desktop/All participants/FSAReward-master/Exp1_fullsample/") # Ivan's work directory
setwd(wd) # set work directory
data.raw <- read.csv("Data_behavior_exp1_48pps.csv",header=TRUE,na.strings="NaN") # read data
#data.raw <- read.csv("allssj.csv",header=TRUE,na.strings="NaN") # read data

data.raw <- rename(data.raw,c(EventType="MovedDots")) # rename EventType variable

data.raw$AttendedColor <- ifelse(data.raw$AttendedColor==1,"red","blue") # attended color: 1 -> red; 2 -> blue
data.raw$RewardedColor <- ifelse(data.raw$ParticipantNo%%2==0,"blue","red") # if participant number is even, blue was rewarded
data.raw$MovedDots <- ifelse(data.raw$MovedDots==1,"red","blue") # color that moved: 1 -> red; 2 -> blue
#data.raw$ExpPhase <- cut(data.raw$Trial,breaks=c(0,100,200,300,400,500,600),labels=c("baseline1","baseline2","acquisition1","acquisition2","extinction1","extinction2")) # trial 0-200: baseline; trial 201-400: acquisition; trial 401-600: extinction
data.raw$ExpPhase <- cut(data.raw$Trial,breaks=c(0,200,400,600),labels=c("baseline","acquisition","extinction")) # trial 0-200: baseline; trial 201-400: acquisition; trial 401-600: extinction
data.raw$ParticipantNo <- factor(data.raw$ParticipantNo) # convert to factor

# count hits, false alarms, misses, correct rejections, and RT separately for each participant
# (their calculation is done in Matlab: see DataProcessing.m)
data.ssj <- ddply(data.raw,.(ParticipantNo,ExpPhase,AttendedColor,RewardedColor,MovedDots),summarize,
                  numtrials=length(which(Response!=99)), # number of trials per condition (anything that is not 99 or any other number that we're not using)
                  Hits=length(which(Response==1)), # hits: attended color moved, correct response
                  FAs=length(which(Response==2)), # false alarms: attended color did not move, (wrong) response
                  Misses=length(which(Response==0)), # misses: attended color moved, no response
                  CRs=length(which(Response==3)), # correct rejections: attended color did not move, no response
                  mean.RT=mean(RT,na.rm=TRUE)) # mean RT per condition

#### Calculate d' prime' #########
data.ssj = ddply(data.ssj, .(ParticipantNo,ExpPhase,AttendedColor), transform, 
                 Hits = Hits[MovedDots==AttendedColor],
                 FAs = FAs[MovedDots!=AttendedColor])

data.ssj = subset(data.ssj,MovedDots==AttendedColor)

# calculate d'
# use loglinear transformation: add 0.5 to Hits, FAs, Misses, and CRs (Hautus, 1995, Behavior Research Methods, Instruments, & Computers),
# which is preferred over the 1/2N rule (Macmillan & Kaplan, 1985, Psychological Bulletin) because it results in less biased estimates of d'.
dataSDT.ssj <-  ddply(data.ssj,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,numtrials),summarize,
                      tot.Hits=Hits+.5, # hits
                      tot.FAs=FAs+.5, # false alarms
                      tot.Misses=(numtrials-tot.Hits)+.5, # misses
                      tot.CRs=(numtrials-tot.FAs)+.5, # correct rejections
                      Hit.Rate=tot.Hits/(tot.Hits+tot.Misses), # hit rate
                      FA.Rate=tot.FAs/(tot.FAs+tot.CRs), # false alarm rate
                      dprime=qnorm(Hit.Rate)-qnorm(FA.Rate)) # d' (see Pallier, 2002)

# Outliers
# Outliers based on hit rate at any condition
crit <- .6 # minimum 60% hit rate in any condition .6
# # Delete participants below the criterion
criterion <- subset(ddply(dataSDT.ssj,.(ParticipantNo),summarize,mean.Hit.Rate=mean(Hit.Rate)),mean.Hit.Rate<crit)$Participant # minimum 60% hit rate across all conditions

# The most principled, we would eliminate 11 participants: 4, 6, 14, 15, 17, 20, 24, 25, 26, 31, 34.
# criterion <- dataSDT.ssj[dataSDT.ssj$Hit.Rate<crit,]$Participant # which participants are marked as outliers
dataSDT.ssj <- dataSDT.ssj[!dataSDT.ssj$ParticipantNo %in% unique(criterion),] # eliminate ouotliers from data frame

# summary d'
summary.dataSDT.dprime <- summarySEwithin(dataSDT.ssj,measurevar="dprime",withinvars=c("ExpPhase","RewardedColor","AttendedColor"),idvar="ParticipantNo",na.rm=TRUE)

# summary d', add condition
dataSDT.ssj$Condition <- ifelse(dataSDT.ssj$RewardedColor==dataSDT.ssj$AttendedColor,"RewAtt","NotRewAtt")

##### PLOTS #####

### summary of d' for merged colors and separately
plottingConditions = c("d' both colors","d' blue rewarded","d' red rewarded" )
for (i in 1:length(plottingConditions)){
  if(plottingConditions[i]=="d' both colors"){datadPrimePlot=dataSDT.ssj}
  if(plottingConditions[i]=="d' blue rewarded"){datadPrimePlot=subset(dataSDT.ssj,RewardedColor=="blue")}
  if(plottingConditions[i]=="d' red rewarded"){datadPrimePlot=subset(dataSDT.ssj,RewardedColor=="red")}
  
  tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image
  
  # Pirate plot
  pirateplot(formula=dprime~ExpPhase+Condition, # dependent~independent variables
             data=datadPrimePlot, # data frame
             main=plottingConditions[i], # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(0,4), # y-axis: limits
             ylab="d'", # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=ParticipantNo, # ID variable
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
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color
  
  dev.off()
}

# summary of Hit rates for merged colors and separately
plottingConditions = c("Hit Rate both colors","Hit Rate blue rewarded","Hit Rate red rewarded" )
for (i in 1:length(plottingConditions)){
  if(plottingConditions[i]=="Hit Rate both colors"){dataHitRatePlot=dataSDT.ssj}
  if(plottingConditions[i]=="Hit Rate blue rewarded"){dataHitRatePlot=subset(dataSDT.ssj,RewardedColor=="blue")}
  if(plottingConditions[i]=="Hit Rate red rewarded"){dataHitRatePlot=subset(dataSDT.ssj,RewardedColor=="red")}
  
  tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image
  
  # Pirate plot
  pirateplot(formula=Hit.Rate~ExpPhase+Condition, # dependent~independent variables
             data=dataHitRatePlot, # data frame
             main=plottingConditions[i], # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(0,1), # y-axis: limits
             ylab="Hit Rate", # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=ParticipantNo, # ID variable
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
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color
  
  dev.off()
}

# summary of False alarms for merged colors and separately
plottingConditions = c("FAs both colors","FAs blue rewarded","FAs red rewarded" )
for (i in 1:length(plottingConditions)){
  if(plottingConditions[i]=="FAs both colors"){dataFAsPlot=dataSDT.ssj}
  if(plottingConditions[i]=="FAs blue rewarded"){dataFAsPlot=subset(dataSDT.ssj,RewardedColor=="blue")}
  if(plottingConditions[i]=="FAs red rewarded"){dataFAsPlot=subset(dataSDT.ssj,RewardedColor=="red")}
  
  tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image
  
  # Pirate plot
  pirateplot(formula=FA.Rate~ExpPhase+Condition, # dependent~independent variables
             data=dataFAsPlot, # data frame
             main=plottingConditions[i], # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(0,.8), # y-axis: limits
             ylab="FA Rate", # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=ParticipantNo, # ID variable
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
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color
  
  dev.off()
}

##### Calculate RTs #####

# Clean outliers (mean +/- 2 SD)
# LowerBoundary <- mean(data.ssj[data.ssj$AttendedColor==data.ssj$MovedDots,]$mean.RT)-2.5*sd(data.ssj[data.ssj$AttendedColor==data.ssj$MovedDots,]$mean.RT)
# UpperBoundary <- mean(data.ssj[data.ssj$AttendedColor==data.ssj$MovedDots,]$mean.RT)+2.5*sd(data.ssj[data.ssj$AttendedColor==data.ssj$MovedDots,]$mean.RT)
# data.ssj$mean.RT <- ifelse(data.ssj$mean.RT<LowerBoundary,NA,data.ssj$mean.RT)
# data.ssj$mean.RT <- ifelse(data.ssj$mean.RT>UpperBoundary,NA,data.ssj$mean.RT)


# Calculate RTs
data.summaryRT.wide.ssj <- ddply(data.ssj,.(ParticipantNo,ExpPhase,RewardedColor,AttendedColor,MovedDots),summarize,
                                 Hits.RTs=mean(mean.RT,na.rm=TRUE)) # mean RTs of hits

# Delete RTs for False alarms
data.summaryRT.wide.ssj <- subset(data.summaryRT.wide.ssj,AttendedColor==MovedDots)

# Make a new variable "Condition" with two levels - RewardedAttended and NotRewardedAttended                                 
data.summaryRT.wide.ssj$Condition <- ifelse(data.summaryRT.wide.ssj$AttendedColor==data.summaryRT.wide.ssj$RewardedColor,"RewAtt","NotRewAtt")

# Outliers - based on previously calculated criterion (see d' section)
#data.summaryRT.wide.ssj <- data.summaryRT.wide.ssj[!data.summaryRT.wide.ssj$ParticipantNo %in% unique(criterion),]




###### CALCULATING DIFFERENCE SCORES ######
# #Calculate difference scores between rewarded and non rewarded trials for each participant
# #Go into wide format
# data.summaryRT.wide.ssj=subset(data.summaryRT.wide.ssj,select=c(ParticipantNo,Hits.RTs,Condition,ExpPhase))
# data.summaryRT.wide.ssj.diff=dcast(data.summaryRT.wide.ssj,ParticipantNo~Condition+ExpPhase, value.var="Hits.RTs")
# 
# data.summaryRT.wide.ssj.diff$BaselineToAcquisition = data.summaryRT.wide.ssj.diff$RewardedAttended_baseline - data.summaryRT.wide.ssj.diff$RewardedAttended_acquisition
# data.summaryRT.wide.ssj.diff$BaselineToExtinction = data.summaryRT.wide.ssj.diff$RewardedAttended_baseline - data.summaryRT.wide.ssj.diff$RewardedAttended_extinction
# write.csv2(data.summaryRT.wide.ssj.diff,file = "data.summaryRT.wide.ssj.diff.csv")
######ADDING QUESTIONNAIRES#########
# 
# #Add BDI (depression) and BISBAS scores to each participant
# BDIdata = read.csv2("BDI.csv",header=TRUE,na.strings="NaN") # read data
# BDIdata = arrange(BDIdata,Participant)#sort BDI data by participant name
# #Delete participants which we don't consider based on the accuracy cutoff
# # Sets accuracy cutoff based on what you selected before
# if (condition == "All_Subjects"){
#   BDIdata = BDIdata
# }else if (condition == "60_Acc"){
#   BDIdata = BDIdata[!BDIdata$Participant %in% criterion60,]
# }else if (condition == "75_Acc"){
#   BDIdata = BDIdata[!BDIdata$Participant %in% criterion75,]
# }
# 
# #Add individual differences to the difference scores
# diff.final.data_sub = cbind(diff.final.data_sub,BDIdata) #merge the two dataframes


##### PLOTS #####

### summary of RTs of hits (blue + red colors)
plottingConditions = c("RTs Hits both colors","RTs Hits blue rewarded","RTs Hits red rewarded" )
for (i in 1:length(plottingConditions)){
  if(plottingConditions[i]=="RTs Hits both colors"){dataRTPlot=data.summaryRT.wide.ssj}
  if(plottingConditions[i]=="RTs Hits blue rewarded"){dataRTPlot=subset(data.summaryRT.wide.ssj,RewardedColor=="blue")}
  if(plottingConditions[i]=="RTs Hits red rewarded"){dataRTPlot=subset(data.summaryRT.wide.ssj,RewardedColor=="red")}
  
  tiff(file = paste(wd,"graphs/",plottingConditions[i],".tiff",sep="")) # save the plot as .tiff image
  
  # Pirate plot
  pirateplot(formula=Hits.RTs~ExpPhase+Condition, # dependent~independent variables
             data=dataRTPlot, # data frame
             main=plottingConditions[i], # main title
             xlim=NULL, # x-axis: limits
             xlab="", # x-axis: label
             ylim=c(400,700), # y-axis: limits
             ylab="Reaction time", # y-axis: label
             inf.method="hdi", # type of inference: 95% Bayesian Highest Density Intervals
             hdi.iter=5000, # number of iterations for estimation of HDI
             inf.within=ParticipantNo, # ID variable
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
             gl.col="gray", # gridlines: color
             gl.lwd=c(.75,0), # gridlines: line width
             gl.lty=2, # gridlines: line type (1: solid; 2:dashed; 3: dotted; ...)
             cex.lab=0.8, # axis labels: size
             cex.axis=1, # axis numbers: size
             bty="l", # plot box type
             back.col="white") # background, color
  
  dev.off()
}


##### STATS #####

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
