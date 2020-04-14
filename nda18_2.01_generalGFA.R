#########################################
#
# NDA 18 2.01 General GFA Script 
#
# modified based on NDA 18.0 now only brain GFA
#
# Last modified 10.31.2019
#
#########################################

library(GFA)
library(DMwR)

# support function for data chunks from:
# https://stackoverflow.com/questions/5577221/how-can-i-load-an-object-into-a-variable-name-that-i-specify-from-an-r-data-file

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

centerData <- function(mydata){
  df <- mydata
  df_means=t(apply(df,1,mean))
  df_sds=t(apply(df,1,sd))
  df_z=sweep(sweep(df,1,df_means,"-"),1,df_sds,"/")
  return(list(df_z,df_means,df_sds))
}

# Directories with the data set:
mydir <- paste0("/home/librad.laureateinstitute.org/mpaulus/ABCD_data/NDA18_2.01/")
# mydir <- paste0("/Users/mpaulus/Dropbox (Personal)/Private/RDataAnalysis/ABCD_Data/NDA2.01/")
setwd(mydir)

# File components:
myfile <- c("nda18_2.01_")
GFAtext <- c("preGFA_")
dateext <- c("_10.31.2019")

# Steps:
# Read in Data Chunks

# Free Surfer Data
# The baseline data are with eventname = baseline_year_1_arm_1

chunkfile <- paste0(mydir,myfile,"FS_cort_subcort",".RDATA")
MPPbrain <- loadRData(chunkfile)
fsnames <- names(MPPbrain)

# Cortical Thickness
thickvars <- fsnames[grep("smri_thick_cort.desikan",fsnames)]
# Remove the mean
thickvars <- thickvars[-grep("mean",thickvars)]

# Sulcal Depth
sulcvars <- fsnames[grep("smri_sulc_cort.desikan",fsnames)]
sulcvars <- sulcvars[-grep("mean",sulcvars)]

# Gray Matter Volume
volvars <- fsnames[grep("smri_vol_cort.desikan",fsnames)]
volvars <- volvars[-grep("total",volvars)]

# Cortical Surface 
areavars <- fsnames[grep("smri_area_cort.desikan",fsnames)]
areavars <- areavars[-grep("total",areavars)]
  
# Subcortical Volumes
subcortvars <- fsnames[grep("smri_vol_subcort",fsnames)]
subcortvars <- subcortvars[-grep("volume",subcortvars)]
subcortvars <- subcortvars[-grep("wholebrain",subcortvars)]
subcortvars <- subcortvars[-grep("ventricle",subcortvars)]
subcortvars <- subcortvars[-grep("_csf",subcortvars)]
subcortvars <- subcortvars[-grep("lesion",subcortvars)]

# Symptom Data - Parent
chunkfile <- paste0(mydir,myfile,"CBCL",".RDATA")
MPPcbcl <- loadRData(chunkfile)
fsnames <- names(MPPcbcl)
cbclvars <- fsnames[grep("cbcl_scr",fsnames)]

# Symptom Data - Youth
chunkfile <- paste0(mydir,myfile,"bis",".RDATA")
MPPbis <- loadRData(chunkfile)
fsnames <- names(MPPbis)
bisvars <- fsnames[grep("bis_y_ss",fsnames)]

chunkfile <- paste0(mydir,myfile,"upps",".RDATA")
MPPupps <- loadRData(chunkfile)
fsnames <- names(MPPupps)
uppsvars <- fsnames[grep("upps_y_ss",fsnames)]

# Cognition Data
chunkfile <- paste0(mydir,myfile,"cog",".RDATA")
MPPcog <- loadRData(chunkfile)
fsnames <- names(MPPcog)
cogvars <- fsnames[-grep("src_",fsnames)]
cogvars <- cogvars[-grep("eventname",cogvars)]

# Screen Media Data
chunkfile <- paste0(mydir,myfile,"SMAn",".RDATA")
MPPsmaN <- loadRData(chunkfile)
fsnames <- names(MPPsma)
smavars <- fsnames[-grep("eventname",fsnames)]
smavars <- smavars[-grep("src_",smavars)]

# This is the SMA variable set with the X-rated and mature-rated items removed
smavarsR <- smavars[-grep("13",smavars)]
smavarsR <- smavarsR[-grep("14",smavarsR)]

# Environmental Data
chunkfile <- paste0(mydir,myfile,"socs",".RDATA")
MPPsocs <- loadRData(chunkfile)
fsnames <- names(MPPsocs)
socvars <- fsnames[-grep("src_",fsnames)]
socvars <- socvars[-grep("eventname",socvars)]

# Get demographic variables
chunkfile <- paste0(mydir,myfile,"demovars",".RDATA")
MPPdemo <- loadRData(chunkfile)

# Labels for different variable sets:
cbclabels <- c("Anxious/Depressed","Withdrawn","Somatic Sx","Social Problems","Thought Problems","Attention Problems","Rule Breaking","Aggressive Behavior","Internalizing","Externalizing","Total Problems")
bislabels <- c("BIS Total","BAS Reward Reactivity","BAS Drive","BAS Fun Seeking")
uppslabels <- c("Negative Urgency","Lack of Planning","Sensation Seeking","Positive Urgency","Lack of Perseverance")
coglabels <- c("Picture Vocabulary","Flanker Tes","List Sorting","Card Sorting","Pattern Comparison","Picture Sequence","Oral Reading Recog","Fluid Composite","Crystallized Composite","Cognition Total","RAVLT Short Delay","RAVLT Long Delay","WISC-V Matrix Reasoning")
soclabels <- c("P: Family Conflict","Y: Family Conflict","P: CPBRI Acceptance","SC: CPBRI Acceptance",
               "P: Prosocial Behavior","Y: Prosocial","Y: Parental Monitoring")
smalabels <- c("WD TV/Movie","WD Videos","WD Games","WD Texting","WD Social Network","WD Chat",
               "WE TV/Movie","WE Videos","WE Games","WE Texting","WE Social Network","WE Chat",
               "Mature Videogame","R-rated Movie")
smalabelsR <- c("WD TV/Movie","WD Videos","WD Games","WD Texting","WD Social Network","WD Chat",
               "WE TV/Movie","WE Videos","WE Games","WE Texting","WE Social Network","WE Chat")

# Get the labels from the desikan atlas csv file

# Desikan labels:
myall <- paste0(mydir,"ABCD_FS_labels",".csv")
desikan <- read.csv(myall,header=TRUE)

# Reorder variables to fit location and side and label:
thickvars_r <-thickvars[c(desikan$Order_orig)]
thicklabels <- paste(desikan[,2],desikan[,4],"TH",sep=" ")
sulcvars_r <- sulcvars[c(desikan$Order_orig)]
sulclabels <- paste(desikan[,2],desikan[,4],"SD",sep=" ")
volvars_r <- volvars[c(desikan$Order_orig)]
vollabels <- paste(desikan[,2],desikan[,4],"GV",sep=" ")
areavars_r <- areavars[c(desikan$Order_orig)]
arealabels <- paste(desikan[,2],desikan[,4],"AR",sep=" ")

# I am here and will need to modify the code from here on down:

# Combine Data File: Cognition, Psychopathology, Environment, SMA
MPPall <- merge(MPPcog[c("src_subject_id","eventname",cogvars)],
                MPPcbcl[c("src_subject_id","eventname",cbclvars)],by=c("src_subject_id","eventname"),all=FALSE)
MPPall <- merge(MPPall,
                MPPsocs[c("src_subject_id","eventname",socvars)],by=c("src_subject_id","eventname"),all=FALSE)
MPPall <- merge(MPPall,
                MPPsmaN[c("src_subject_id","eventname",smavarsR)],by=c("src_subject_id","eventname"),all=FALSE)

# select only baseline and only quality controlled data
MPPall_b <- subset(MPPall,(MPPall$eventname=="baseline_year_1_arm_1"))

# Prep Data for GFA
indepvars <- c(cogvars,cbclvars,socvars,smavarsR)

#Original Variables
COG <- as.matrix(MPPall_b[,c(cogvars)])
CBCL <- as.matrix(MPPall_b[,c(cbclvars)])
SOC <- as.matrix(MPPall_b[,c(socvars)])
SMA <- as.matrix(MPPall_b[,c(smavarsR)])

tmp <- cbind(COG,CBCL,SOC,SMA)

# Imputation goes here, most of the imputation software would not work on the server R version
# Instead I use DMwR
# https://stats.stackexchange.com/questions/61110/knn-imputation-r-packages

impute_tmp <- knnImputation(tmp)

#iCOG <- knnImputation(COG)
#iCBCL <- knnImputation(CBCL)
#iSOC <- knnImputation(SOC)
#iSMA <- knnImputation(SMA)

# generate a list of variables:

MY <- list(impute_tmp[,cogvars],
           impute_tmp[,cbclvars],
           impute_tmp[,socvars],
           impute_tmp[,smavarsR])

# MY <- list(iCOG,iCBCL,iSOC,iSMA)

# File details:
GFAtext <- c("COG_CBCL_SOC_SMA")
dateext <- c("_10.31.2019")

# Save the data frame
MYdf <- data.frame(MY)

myall <- paste0(mydir,myfile,GFAtext,"data",".RData",sep="")
save(MYdf, file=myall)

# Normalize all features
mynorm <- normalizeData(MY, type="scaleFeatures")

# set up the GFA defaults

# Get the default options
opts <- getDefaultOpts()
opts$vrbose <- 0

# number of data for posterior vector:
opts$iter.saved = 100

startK = length(c(indepvars))

# Run GFA

set.seed(TEMP);
res <- list()
for(i in TEMP:TEMP){
  print(i)
  myall <- paste0(mydir,myfile,GFAtext,i,dateext,".RData",sep="")
  print(myall)
  
  res[[i]] <- gfa(mynorm$train, K=startK, opts=opts)
  
  # Save as an interim variable
  myres <- res[[i]]
  
  # Save the GFA results
  # Write the result to a file
  myall <- paste0(mydir,myfile,GFAtext,i,dateext,".RData",sep="")
  save(myres, file=myall)
}

# Save robust components
# rob <- robustComponents(res)
# myall <- paste0(mydir,myfile,GFAtext,"Robust",dateext,".RData",sep="")
# save(rob, file=myall)
