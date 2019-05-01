#### Clear environemnt and import data ####

# clear the environment
rm(list=ls()) 
# clear the console
cat("\014") 
#load packages and install them if they're not installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,reshape2,tidyverse,here,car)
# set seed
set.seed(42) 
# set directory
setwd(here())
# import data
questionnaires = read.csv2(file = here("data","Questionnaires.csv"),header=TRUE,na.strings="NaN") 

#### Exclude participants with no EEG data ####

missing_eeg = c(35,39,44)
questionnaires = questionnaires[!questionnaires$Subject %in% missing_eeg,] 

#### Calculate gender ####
gender = ddply(questionnaires,.(Subject),summarize,
               Gender = unique(Sex))
table(unlist(gender$Gender))


#### Calculate age ####
age = ddply(questionnaires,.(Subject),summarize,
               Age = unique(Age))
age$Age=as.numeric(age$Age)
median(age$Age)

#### Calculate BDI ####
bdi = ddply(questionnaires,.(Subject),summarize,
            BDI = sum(ItemdisplayBDI3.RESP,na.rm=TRUE))

#### Calculate BIS BAS ####

bisbas = subset(questionnaires, Procedure.Block. == "QBISBAS", select = c(Subject,ItemQ,ItemdisplayBISBAS.RESP,ItemsBISBAS))

# Recode reversed items (2 and 22)
bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 2] = recode(bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 2], "0=3; 1=2; 2=1; 3=0")
bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 22] = recode(bisbas$ItemdisplayBISBAS.RESP[bisbas$ItemsBISBAS == 22], "0=3; 1=2; 2=1; 3=0")

# Put on the scale from 1-4 instead of 0-3
bisbas$ItemdisplayBISBAS.RESP = bisbas$ItemdisplayBISBAS.RESP + 1

bisbas = ddply(bisbas,.(Subject),summarize,
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

