###Script to process the metadata to a proper format

#Load file paths for inputs
DATADIR <- '../data/IBS_Mayo_secondrun/'
metadata <- paste(DATADIR, "Metadata_IBS.csv", sep='')
taxfp <- paste(DATADIR,'taxatable.txt',sep='')
medsfp <- paste(DATADIR, "Baseline_meds.txt", sep='')
condense <- paste(DATADIR, "condensed_metadata_file_v3.csv", sep='')

#Load the metadata and meds files
meta <- read.delim(metadata, 
                   sep=',', 
                   head=T, 
                   comment="",
                   quote="")

meds <- read.delim(medsfp,
                   sep="\t",
                   head=T,
                   comment="",
                   quote="")

cons <- read.delim(condense,
                   sep=",",
                   head=T,
                   comment="",
                   quote="",
                   check.names = F)

#Reformat the sample IDs
for(r in 1:nrow(meta)){
  sample <- paste(meta[r,"ID_on_tube"], meta[r,"Timepoint"], sep=".T.")
  meta$SampleID[r] <- sample
}
m <- meta
for(r in 1:nrow(meta)){
  sample <- paste(meta[r,"study_id"], meta[r,"ID_on_tube"], sep=".")
  meta$study_vial[r] <- sample
}


#Read in the OTU table
x <- t(read.delim(taxfp, row=1)) #532 samples by 2783 OTUs
taxa_names <- colnames(x)

#clean duplicates from the run that was messed up
rownames(x) <- gsub(".S[0-9]+.R1.001","",rownames(x));      # Clean old plate IDs
x <- data.frame(x[order(rownames(x)),]);              # Sort nicely by sample ID, OLDs will be second
x$samples <- rownames(x)
x$samples <- gsub(".OLD","",x$samples); # remove OLD and NPM
x$samples <- gsub(".NPM", "", x$samples); #532
x2 <- x[match(unique(x$samples), x$samples),] #keep first instance of dups
rownames(x2) <- x2$samples    # Keep new runs only
x <- x2[,!colnames(x2)=="samples"]
rownames(x) <- gsub("Study.ID.", "", rownames(x))

#Any duplicated rownames?
row.names(m)[which(duplicated(row.names(m)))] 

#Format the metadata
m <- as.matrix(m) #This is because there were duplicates in the map
rownames(m) <- m[,"SampleID"]
m <- as.data.frame(m, stringsAsFactors = F)

#Find missing samples
keeps <- intersect(rownames(x), rownames(m)) #456 samples overlap between the OTU table and metadata (original)
not_in_map <- rownames(x)[which(!rownames(x) %in% rownames(m))]
cat("These samples aren't in the map:")
not_in_map
#[1] "Negative.Control.H12.100" "Pos.Ctrl"                 "Positive.Control.B01.101"
# [4] "Positive.Control.G12.100" "10.T.6"                   "11.T.4"                  
# [7] "14.T.1"                   "14.T.2"                   "14.T.3"                  
# [10] "14.T.4"                   "14.T.5"                   "14.T.6"                  
# [13] "15.T.4"                   "16.T.2"                   "18.T.3"                  
# [16] "18.T.4"                   "52.T.5"                   "56.T.2"                  
# [19] "57.T.4"                   "6.T.5"                    "7.T.3"                   
# [22] "7.T.6"                    "73.T.4"                   "76.T.4"                  
# [25] "80.T.0"                   "81.T.0"                   "82.T.0"                  
# [28] "83.T.0"                   "84.T.0"                  
not_in_otu <- rownames(m)[which(!rownames(m) %in% rownames(x))]
cat("These samples aren't in the otu table:")
not_in_otu
#"33.T.1" "50.T.1" "51.T.0" "59.T.1" "63.T.4" "65.T.3" "65.T.4" "65.T.6"

#Fill in missing samples:
m2 <- as.data.frame(t(m), stringsAsFactors = F)

for(i in 5:length(not_in_map)){
  id <- not_in_map[i]
  id2 <- paste(strsplit(id, ".",fixed=T)[[1]][1], strsplit(id, ".",fixed=T)[[1]][2], "0", sep=".")
  if(id2 %in% colnames(m2)){
    m2$new <- m2[,id2]
  } else {
    m2 <- cbind(m2, c(rep(NA, nrow(m2))))
  }
  colnames(m2)[ncol(m2)] <- id
  m2[14:nrow(m2),id] <- NA
  m2["Timepoint",id] <- strsplit(id, ".", fixed=T)[[1]][3]
  m2["SampleID",id] <- id
  m2["study_vial",id] <- paste(m2["study_id", id], m2["Timepoint", id], sep=".")
}

m2 <- as.data.frame(t(m2))
keeps <- intersect(rownames(x), rownames(m2)) #Now 481 samples overlap between the OTU table and metadata
not_in_map <- rownames(x)[which(!rownames(x) %in% rownames(m2))]
cat("These samples aren't in the map:")
not_in_map #Now just pos and neg controls

#merge metadata with meds
merged_m2 <- merge(m2, meds, by= "SampleID", all=T)
merged_m2$q10_subj.x <- as.character(merged_m2$q10_subj.x)
merged_m2$yogurt_pro <- merged_m2$q10_subj.x
for(r in 1:nrow(merged_m2)){
  if(merged_m2$yogurt_pro[r] == "" | is.na(merged_m2$yogurt_pro[r])){
    merged_m2[r, "yogurt_pro"] <- "no"
  } else {
    merged_m2[r, "yogurt_pro"] <- "yes"
  }
}

#write to .csv
#NOTE: Go in and replace the "," with ";" in excel!!
write.csv(m2, "../data/IBS_Mayo_secondrun/Cleaned_Metadata.csv", quote=F, row.names=F)
write.csv(merged_m2,  "../data/IBS_Mayo_secondrun/Cleaned_Metadata_wMeds.csv", quote=F, row.names=F )

######################################################################

#Order the samples the same
rownames(cons) <- cons$SampleID
rownames(merged_m2) <- merged_m2$SampleID
rownames(cons) %in% rownames(merged_m2)
merged_m2 <- merged_m2[rownames(cons),]
rownames(cons) == rownames(merged_m2)
cols_keep <- c("SampleID", "Flare_timepoint", "timeline", "laxative", 
               "birthcontrol", "estrogen", "anti_viral", "pain_med",
               "proton_pump", "statin","allergy_med","anti_depress", 
               "cal_channel", "face_abx","probiotic", "baseline.medications", 
               "Antibiotics", "which_abx", "yogurt_pro")

m2_keeps <- merged_m2[, cols_keep]
cons2 <- merge(cons, m2_keeps, by="SampleID")


col_list <- c("laxative", 
              "birthcontrol", "estrogen", "anti_viral", "pain_med",
              "proton_pump", "statin","allergy_med","anti_depress", 
              "cal_channel", "face_abx","probiotic", "Antibiotics")

for(l in 1:length(col_list)){
  this <- col_list[l]
  for(r in 1:nrow(cons2)){
    cons2[,this] <- as.character(cons2[,this])
    if(cons2[r,this] == "" | is.na(cons2[r, this])){
      cons2[r,this] <- "no"
    } else {
      cons2[r,this] <- "yes"
    }
  }
}

cons2$DietID <- paste0(cons2$study_id, ".", cons2$Timepoint)

write.csv(cons2,  "../data/IBS_Mayo_secondrun/Final_Cleaned_Condensed_Metadata.csv", quote=F, row.names=F)
