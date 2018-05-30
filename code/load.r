#Set your working directory to be a new directory within this main one
library('vegan')
library('RColorBrewer')
library('ggplot2')
library('ggpubr')
library('reshape2')
library('plyr')
library('ape')
library('grid')
library('gridExtra')
library('cowplot')
library('ROCR')
library('flux')
library('randomForest')
library('robCompositions')
library('polycor')

#Create output dirs
source("../code/output.dirs.r")
source('../code/stats.r')
source('../code/util.r')
source('../code/risk.index.r')

#Load file paths for inputs
DATADIR <- '../data/IBS_Mayo_secondrun/'
metadata <- paste(DATADIR, "Final_Cleaned_Condensed_Metadata.csv", sep='')
metadatabiopsy <- paste(DATADIR, "mucosal_bacteria_condensed_metadata_file_v2.csv", sep="")
taxfp <- paste(DATADIR,'/shogun_output/id_update/taxatable.burst.capitalist.txt',sep='')
taxbiopyfp <- paste(DATADIR, '/shogun_output_biopsy/taxatable.burst.capitalist.txt', sep="")
taxALLfp <- paste(DATADIR, 'stool_biopsy_otus.txt', sep="")
modulefp <- paste(DATADIR, '/shogun_output/id_update/taxatable.strain.kegg.modules.txt', sep='')
module_mapfp <- paste(DATADIR, 'Kegg_ID_Map.txt', sep='')
keggfp <- paste(DATADIR, '/shogun_output/id_update/taxatable.strain.kegg.txt', sep='')
keggmapfp <- paste(DATADIR, "ko-enzyme-annotations.txt", sep='')
pathwayfp <- paste(DATADIR, '/shogun_output/id_update/taxatable.strain.kegg.pathways.txt', sep='')
trans1_fp <- paste(DATADIR, "NormalizedExpression_DEGs_T1.txt", sep ="")
trans2_fp <- paste(DATADIR, "NormalizedExpression_DEGs_T2.txt", sep ="")

#Load metadata (this was generated with metadata_prep.r)
m <- read.delim(metadata, 
                   sep=',', 
                   head=T, 
                   comment="",
                   row=1)
m$SampleID <- rownames(m)

m_bio <- read.delim(metadatabiopsy, 
                sep=',', 
                head=T, 
                comment="",
                quote="",
                row=1)
m_bio$SampleID <- rownames(m_bio)
m_bio$Timepoint <- gsub("nd", "", m_bio$Timepoint)
m_bio$Timepoint <- gsub("st", "", m_bio$Timepoint)

#Read in the OTU table, reformat the names of the samples, remove duplicates
x <- t(read.delim(taxfp, row=1, check.names=F, as.is =T)) #532 samples by 2783 OTUs
taxa_names <- colnames(x) #store this because it replaces ";" with "." later

x_bio <- t(read.delim(taxbiopyfp, row=1, check.names = F, as.is=T))
taxa_bio <- colnames(x_bio)

#Find samples are still not accounted for:
keeps <- intersect(rownames(x), rownames(m)) #481 samples overlap between the OTU table and metadata
not_in_map <- rownames(x)[which(!rownames(x) %in% rownames(m))]
cat("These samples aren't in the map:")
not_in_map
#[1] "Negative.Control.H12.100" "Pos.Ctrl" "Positive.Control.B01.101" "Positive.Control.G12.100" 

not_in_otu <- rownames(m)[which(!rownames(m) %in% rownames(x))]
cat("These samples aren't in the otu table:")
not_in_otu
#"33.T.1" "50.T.1" "51.T.0" "59.T.1" "63.T.4" "65.T.3" "65.T.4" "65.T.6" "58.T.0"

not_in_biospy <- rownames(m_bio)[which(!rownames(m_bio) %in% rownames(x_bio))]
patients_no_biopsy <- m_bio$ID_on_tube[which(!rownames(m_bio) %in% rownames(x_bio))]
cat("these many patients don't have biopsy sequences:")
unique(patients_no_biopsy)
cat("These samples aren't in the biopsy otu table:")
not_in_biospy
keeps <- intersect(rownames(x_bio), rownames(m_bio))
x_bio <- x_bio[keeps,]
m_bio <- m_bio[keeps,]

#remove family samples
m <- m[! is.na(m$study_id),]

#Keep only samples with 10000 counts
#original = 485 samples, 2343 OTUs 
x <- x[rowSums(x) > 10000,] #477 samples pass 
keeps <- intersect(rownames(x), rownames(m))
m <- m[keeps,]
x <- x[keeps,]

#Add aligned sequence counts to map
m$counts <- rowSums(x)
colnames(x) <- taxa_names #replace proper taxa names

m_bio$counts <- rowSums(x_bio)

#Remove archaea to sep table, archaea = 21
archaea <- x[,grep("Archaea", colnames(x))]
m$archaea <- rowSums(archaea)

archea_bio <- x_bio[,grep("Archaea", colnames(x_bio))]
m_bio$archea <- rowSums(archea_bio)

#Remove viruses/phage to sep table = 274
viruses <- x[,grep("Viruses", colnames(x))]
m$viruses <- rowSums(viruses)

viruses_bio <- x_bio[,grep("Viruses", colnames(x_bio))]
m_bio$viruses <- rowSums(viruses_bio)

#Keep only bacteria
xb <- x[,grep("Bacteria", colnames(x))] #now 2487 OTUs
m$bacteria <- rowSums(xb) 

# #Keep species name, or the lowest name listed for the taxon
x <- t(x)
x_bio <- t(x_bio)

##keep only species
#keeps <- grepl("s__", rownames(x)) & ! grepl("s__$", rownames(x))
#sum(keeps)
#x <- x[keeps,]

#remove these because some have them and some don't :(
# rownames(x) <- gsub("k__", "", rownames(x))
# rownames(x) <- gsub("p__", "", rownames(x))
# rownames(x) <- gsub("c__", "", rownames(x))
# rownames(x) <- gsub("o__", "", rownames(x))
# rownames(x) <- gsub("f__", "", rownames(x))
# rownames(x) <- gsub("g__", "", rownames(x))
# rownames(x) <- gsub("s__", "", rownames(x))
# rownames(x) <- gsub("t__", "", rownames(x))
# 

#Now split on delim ";"
split <- strsplit(rownames(x),";")
for(i in 1:length(split)){
  if(length(split[[i]]) > 6){
    rownames(x)[i] <- split[[i]][7] #If it has species, keep that
  } else {
    rownames(x)[i] <- paste(split[[i]][1:(length(split[[i]]))],collapse="_") #Or collapse down to most specific
  }
}
x <- aggregate(x, by=list(rownames(x)),sum) #aggregate them together, now 2209 taxa


split <- strsplit(rownames(x_bio),";")
for(i in 1:length(split)){
  if(length(split[[i]]) > 6){
    rownames(x_bio)[i] <- split[[i]][7] #If it has species, keep that
  } else {
    rownames(x_bio)[i] <- paste(split[[i]][1:(length(split[[i]]))],collapse="_") #Or collapse down to most specific
  }
}
x_bio <- aggregate(x_bio, by=list(rownames(x_bio)),sum) #aggregate them together, now 873 taxa *before = 1174

# ##remove square brackets
# #x$Group.1 <- gsub("\\[|\\]", "", x$Group.1)
# #x$Group.1 <- gsub("[.]", "", x$Group.1)
rownames(x) <- x$Group.1 #Aggregating store the rownames as a column
x <- x[,!colnames(x) == "Group.1"]
x <- t(x)

x.raw <- x #store raw
x <- sweep(x, 1, rowSums(x), '/') #store RA

rownames(x_bio) <- x_bio$Group.1
x_bio <- x_bio[,!colnames(x_bio) == "Group.1"]
x_bio <- t(x_bio)

x_bio_raw <- x_bio #store raw
x_bio <- sweep(x_bio_raw, 1, rowSums(x_bio_raw), '/') #store RA

rownames(x_bio) <- gsub(" ", "_", rownames(x_bio))
rownames(x_bio_raw) <- gsub(" ", "_", rownames(x_bio_raw))
#Write these to output so you can merge in QIIME!
file_name <- paste(DATADIR, "modified_stool_otus.txt", sep="")
cat("#OTUID\t", file=file_name)
write.table(t(x.raw), file=file_name, append=T, row.names = T, quote=F, sep="\t")

file_name <- paste(DATADIR, "modified_biopsy_otus.txt", sep="")
cat("#OTUID\t", file=file_name)
write.table(t(x_bio_raw), file=file_name, append=T, row.names = T, quote=F, sep="\t")

#Collapse by subject
m$ID_on_tube <- as.integer(m$ID_on_tube)
x2 <- x[!m$Timepoint == "Flare",] #Leave out Flares!
m2 <- m[!m$Timepoint == "Flare",]

#xc will be collapsed with no flares
xc <- apply(x2,2,function(xx) sapply(split(xx,m2$ID_on_tube),mean))
rownames(xc) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))
#xc will be collapsed with flares
xc_flares <- apply(x,2,function(xx) sapply(split(xx,m$ID_on_tube),mean))
rownames(xc_flares) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

#mc = no flares
mc <- m2[sapply(split(1:nrow(m2),m2$ID_on_tube),'[',1),,drop=TRUE]
rownames(mc) <- sprintf('Subject_%03d',sapply(split(m2$ID_on_tube,m2$ID_on_tube),'[',1))

#mc_flares = has flares
#mc = no flares
mc_flares <- m[sapply(split(1:nrow(m),m$ID_on_tube),'[',1),,drop=TRUE]
rownames(mc_flares) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

#This contains flares
xc.raw <- apply(x.raw,2,function(xx) sapply(split(xx,m$ID_on_tube),sum))
rownames(xc.raw) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

# drop rare bugs (show up in less than 10% of subjects)
# Goes to 1095 OTUs
rare.ix <- colMeans(xc > 0) < .10
x.raw <- x.raw[,!rare.ix]
xc.raw <- xc.raw[,!rare.ix]
xc_flares <- xc_flares[,!rare.ix]
x <- x[,!rare.ix]
xc <- xc[,!rare.ix]

rare.ix <- colMeans(x_bio > 0) < .05
x_bio_raw <- x_bio_raw[,!rare.ix] #450 taxa
x_bio <- x_bio[,!rare.ix]

#Load the full OTU table (combed in qiime and mod names, etc)
x_all <- t(read.delim(taxALLfp, row=1, skip=1, check.names=F, as.is =T)) #546 samples by 1545 OTUs
taxa_names <- colnames(x_all) #store this because it replaces ";" with "." later
m$sType <- m$Flare
m$sType <- as.character(m$sType)
m["sType"][is.na(m["sType"])] <- "Stool"
m_bio$sType <- rep("Biopsy", nrow(m_bio))

map_cols <- intersect(colnames(m), colnames(m_bio))
m_all <- data.frame(rbind(m[,map_cols], m_bio[,map_cols]))

keeps <- intersect(rownames(x_all), rownames(m_all))
x_all <- x_all[keeps,] #546 samples
m_all <- m_all[keeps,]

x_all.raw <- x_all #store raw
x_all <- sweep(x_all, 1, rowSums(x_all), '/') #store RA

# drop rare bugs (show up in less than 5 subjects)
# Goes to 1095 OTUs
rare.ix <- colSums(x_all > 0) < 5
x_all.raw <- x_all.raw[,!rare.ix]
x_all <- x_all[,!rare.ix] #1501 OTUs

# Load Modules
modules <- t(read.table(modulefp, header=T, sep='\t', row=1, comment='', check.names = F, as.is=T))
modules <- modules[rownames(m),]
modules <- modules[,colSums(modules)>0]
rare.ix <- colMeans(modules > 0) < .25
modules <- modules[,!rare.ix]

coef.variation <- function(xx) {
  sqrt(var(xx))/mean(xx)
}

modules <- data.frame(modules)
colnames(modules) <- gsub("\\.", "",colnames(modules))
co_varies <- sapply(modules[,colnames(modules)],coef.variation)
modules <- modules[,colnames(modules)[order(co_varies, decreasing = T)]]

#Collapse by subject
modules2 <- modules[!m$Timepoint == "Flare",] #Leave out Flares!

#Load module map
module_names <- read.table(module_mapfp, sep="\t", header=T, as.is=T)
module_names$module <- sapply(strsplit(as.character(module_names$KEGG_Label),'\ '), "[", 1)
module_names$BugBase_ID <- gsub(".*?_(.+)", "\\1", module_names$BugBase_ID)
rownames(module_names) <- module_names$module
module_names <- module_names[colnames(modules),]
colnames(modules) <- module_names$BugBase_ID
colnames(modules2) <- module_names$BugBase_ID

#xc will be collapsed with no flares
modulesc <- apply(modules2,2,function(xx) sapply(split(xx,m2$ID_on_tube),mean))
modulesc <- modulesc[,colSums(modulesc) > 0]
rownames(modulesc) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))
#xc will be collapsed with flares
modulesc_flares <- apply(modules,2,function(xx) sapply(split(xx,m$ID_on_tube),mean))
rownames(modulesc_flares) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

# Load KEGG
kegg <- t(read.table(keggfp, header=T, sep='\t', row=1, comment='', check.names=F, as.is=T)) #532 samples by 6669 kegg
kegg <- kegg[rownames(m),]
kegg <- kegg[,colSums(kegg)>1]
rare.ix <- colMeans(kegg > 0) < .25
kegg <- kegg[,!rare.ix] #now 532 samples by 4076 kegg

k_map <- read.table(keggmapfp, sep='\t', comment='') #6451 keggs
keeps <- intersect(colnames(kegg), k_map$V1)
k_map <- k_map[k_map$V1 %in% keeps,] #now 1876
#co_varies <- sapply(kegg[,colnames(kegg)],coef.variation) 
#kegg <- kegg[,colnames(kegg)[order(co_varies, decreasing = T)]] #order by covariance


#Collapse by subject
kegg2 <- kegg[!m$Timepoint == "Flare",] #Leave out Flares!
#keggc will be collapsed with no flares
keggc <- apply(kegg2,2,function(xx) sapply(split(xx,m2$ID_on_tube),mean))
keggc <- keggc[,colSums(keggc) > 0]
rownames(keggc) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))
#keggf_flares will be collapsed with flares
keggc_flares <- apply(kegg,2,function(xx) sapply(split(xx,m$ID_on_tube),mean))
rownames(keggc_flares) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

# gather differentiation indicies
ixc.hc <- mc$Cohort == "H"
ix.hc <- m$Cohort == "H" & is.na(m$Flare)
ixc.ibsc <- mc$Cohort == "C"
ix.ibsc <- m$Cohort == "C" & is.na(m$Flare)
ixc.ibsd <- mc$Cohort == "D"
ix.ibsd <- m$Cohort == "D" & is.na(m$Flare)
ixc.ibs <- mc$Cohort != "H"
ix.ibs <- m$Cohort != "H" & is.na(m$Flare)

ixcb.hc <- m_bio$Cohort == "H"
ixcb.ibs <- m_bio$Cohort != "H"
ixcb.ibsd <- m_bio$Cohort == "D"
ixcb.ibsc <- m_bio$Cohort == "C"

cols <- c("#cb1b4a", "#42aeb8", "#FDB316", "#c3c823", "#00797e", "#053058", "#aaada6", "#ae1848", "#368b90", "#2823c8",  "#ca9012", "#124cca", "#9ba11e", "#f1f2ec", "#d9d3df", "#348fbe", "#ff8340", "#ffAf40", "#bf503f", "#503fbf", "#951b72", "#b75f6d")
cols2 <- colorRampPalette(cols)
cols_ibs <- c("#FDB316", "#9ba11e")
cols_dh <- c("#42aeb8","#FDB316")
cols_ch <- c("#cb1b4a", "#FDB316")
cols_yn <- c("#951B72", "#AAADA6")

#There are extra samples in here that don't match
mc$IBS <- as.character(mc$Cohort)
mc[mc$IBS == "H", "IBS"] <- "Healthy"
mc[mc$IBS == "C" | mc$IBS == "D", "IBS"] <- "IBS"
m$IBS <- as.character(m$Cohort)
m[m$IBS == "H" & !is.na(m$IBS), "IBS"] <- "Healthy"
m[!m$IBS == "Healthy" & !is.na(m$IBS), "IBS"] <- "IBS"

mc$IBS <- factor(mc$IBS)
m$IBS <- factor(m$IBS)

m$Cohort <- as.character(m$Cohort)

#how many timelines are complete?
complete_timeline <- c(0, 0, 0)
names(complete_timeline) <- c("H", "C", "D")
for(i in 1:length(unique(m$study_id))){
  working <- unique(m$study_id)[i]
  if(! is.na(working)){
    if(length(which(m$study_id == working)) == 7){
      cohort <- as.character(m[which(m$study_id == working),"Cohort"])[1]
      complete_timeline[cohort] <- c(complete_timeline[cohort] + 1)
    }
  }
}
print("Complete Timelines:")
print(complete_timeline)

#How many have paired biopsies?
complete_timeline <- c(0, 0, 0)
names(complete_timeline) <- c("H", "C", "D")
for(i in 1:length(unique(m_bio$study_id))){
  working <- unique(m_bio$study_id)[i]
  if(! is.na(working)){
    if(length(which(m_bio$study_id == working)) == 2){
      cohort <- as.character(m_bio[which(m_bio$study_id == working),"Cohort"])[1]
      complete_timeline[cohort] <- c(complete_timeline[cohort] + 1)
    }
  }
}
print("2 biopsies:")
print(complete_timeline)

#what are the average number of samples each person has?
average_timeline <- list(NA, NA, NA)
names(average_timeline) <- c("H", "C", "D")
for(i in 1:length(unique(m$study_id))){
  working <- unique(m$study_id)[i]
  if(! is.na(working)){
    cohort <- as.character(m[which(m$study_id == working),"Cohort"])[1]
    average_timeline[[cohort]]<- c(average_timeline[[cohort]], length(which(m$study_id == working)))
  }
}
for(i in 1:length(average_timeline)){
  average_timeline[[i]] <- mean(average_timeline[[i]], na.rm = T)
}
print("Average Length of Timelines:")
print(average_timeline)

#Number of subjects with < 3 samples
average_timeline <- c(0, 0, 0)
names(average_timeline) <- c("H", "C", "D")
for(i in 1:length(unique(m$study_id))){
  working <- unique(m$study_id)[i]
  if(! is.na(working)){
    cohort <- as.character(m[which(m$study_id == working),"Cohort"])[1]
    if(length(which(m$study_id == working)) < 3){
      average_timeline[[cohort]]<- average_timeline[[cohort]] + 1
    }
  }
}
print("Number of subjects with < 3 samples:")
print(average_timeline)

