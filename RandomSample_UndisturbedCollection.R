#Subset Routine of 60 random sounds from each collection

setwd("/Volumes/G-RAID with Thunderbolt/SoundEcology_R/METRICSsoundBank/")
dir()

#Sample from Arizona
#set a seed
set.seed(1234) 

#Get complete list of all files in directory to be sampled
filesA <- list.files(path = "Undisturbed_Collections/Arizona")

#Run random sample and get the sub file list
subfilesA <- files[sample(1:length(filesA), 60, replace=FALSE)]

#copy those files to a new path
sapply(subfilesA, function(x) file.copy(from=paste("Undisturbed_Collections/Arizona/",x,sep=""), to="RandomSample_Undisturbed/AZ_RS"))

#set new seed
set.seed(1234) 

# Sample from Maine
#Get complete list of all files in directory to be sampled
filesM <- list.files(path = "Undisturbed_Collections/Maine")

#Run random sample and get the sub file list
subfiles <- files[sample(1:length(filesM), 60, replace=FALSE)]

#copy those files to a new path
sapply(subfilesM, function(x) file.copy(from=paste("Undisturbed_Collections/Maine/",x,sep=""), to="RandomSample_Undisturbed/ME_RS"))

#set new seed
set.seed(1234)

#Sample from Costa Rica 2011
#Get complete list of all files in directory to be sampled
filesL <- list.files(path = "Undisturbed_Collections/LaSelva2011")

#Run random sample and get the sub file list
subfiles <- files[sample(1:length(files), 60, replace=FALSE)]

#copy those files to a new path
sapply(subfiles, function(x) file.copy(from=paste("Undisturbed_Collections/LaSelva2011/",x,sep=""), to="RandomSample_Undisturbed/CR_RS"))
