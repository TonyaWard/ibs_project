library("vegan")
library("RColorBrewer")
library("ggplot2")
library("reshape2")
library("plyr")
library('beeswarm')
library('ape')
library("grid")
library("gridExtra")
library("cowplot")
library('ROCR')
library("flux")
library('randomForest')

DATADIR <- '../data/IBS_Mayo_secondrun/'
#mapfp <- paste(DATADIR,'170118_updated.txt',sep='')
#mapfp2 <- paste(DATADIR, '170118_Sample_matrix.txt', sep='')
metadata <- paste(DATADIR, "Metadata_IBS.csv", sep='')
taxfp <- paste(DATADIR,'ibs_tax.txt',sep='')

#Sample IDs don't match the map - AT ALL.
#Reformat the names of the samples
meta <- read.delim(metadata, 
                   sep=',', 
                   head=T, 
                   comment="",
                   quote="")

for(r in 1:nrow(meta)){
  sample <- paste(meta[r,"ID_on_tube"], meta[r,"Timepoint"], sep=".T.")
  meta$SampleID[r] <- sample
}
for(r in 1:nrow(meta)){
  sample <- paste(meta[r,"study_id"], meta[r,"ID_on_tube"], sep=".")
  meta$study_vial[r] <- sample
}

#Sample IDs don't match the map - AT ALL.
#Reformat the names of the samples
x <- t(read.delim(taxfp, row=1))
sample_names <- read.table(text = rownames(x), sep = ".", colClasses = "character")
for(r in 1:nrow(sample_names)){
  sample <- paste(sample_names[r,3], sample_names[r,4], sample_names[r,5], sep=".")
  rownames(x)[r] <- sample
}
m <- meta
#Find duplicated rownames
row.names(m)[which(duplicated(row.names(m)))] #28 samples are duplicated in orginal table, fix in new table

#Keep only samples with 1000 counts
#original = 375 samples, 1005 OTUs
x <- x[rowSums(x) > 1000,] #all samples pass
m <- as.matrix(m) #This is because there were duplicates in the map
rownames(m) <- m[,"SampleID"]
m <- as.data.frame(m, stringsAsFactors = F)

keeps <- intersect(rownames(x), rownames(m)) #only 355 samples overlap between the OTU table and metadata (original)
not_in_map <- rownames(x)[which(!rownames(x) %in% rownames(m))]
cat("These samples aren't in the map:")
not_in_map
# "10.T.6" "34.T.4" "7.T.3"  "15.T.4" "10.T.3" "7.T.6"  "14.T.2" "14.T.3" "18.T.3" "14.T.4" "11.T.4" "14.T.6" "52.T.5" "38.T.4" "14.T.1" "16.T.2" "18.T.4" "14.T.5" "56.T.2" "6.T.5"
not_in_otu <- rownames(m)[which(!rownames(m) %in% rownames(x))]
cat("These samples aren't in the otu table:")
not_in_otu

#Fill in missing samples:
m2 <- as.data.frame(t(m))
for(i in 1:length(not_in_map)){
  id <- not_in_map[i]
  id2 <- paste(strsplit(id, ".",fixed=T)[[1]][1], strsplit(id, ".",fixed=T)[[1]][2], "0", sep=".")
  m2$new <- m2[,id2]
  colnames(m2)[ncol(m2)] <- id
  m2[14:nrow(m2),id] <- NA
}
m2 <- as.data.frame(t(m2))
keeps <- intersect(rownames(x), rownames(m2)) #Now 375 samples overlap between the OTU table and metadata
not_in_map <- rownames(x)[which(!rownames(x) %in% rownames(m2))]
cat("These samples aren't in the map:")
not_in_map

write.csv(m2, "../data/IBS_Mayo_secondrun/Cleaned_Metadata.csv", quote=F, row.names=F)

m <- m2[keeps,]
x <- x[keeps,]

## normalize 
m$counts <- rowSums(x)

#Remove archaea to sep table, archaea = 14
archaea <- x[,grep("Archaea", colnames(x))]
m$archaea <- rowSums(archaea)

#Remove viruses/phage to sep table = 295
viruses <- x[,grep("Viruses", colnames(x))]
m$viruses <- rowSums(viruses)

#Collapse to species--> original was 2551, bacteria only = 2441
x <- x[,grep("Bacteria", colnames(x))]
m$bacteria <- rowSums(x)
for(i in 1:ncol(x)){
    colnames(x)[i] <- strsplit(colnames(x)[i], ";", fixed=T)[[1]][7]
}
ncol(x)
x <- as.data.frame(t(x))
x <- aggregate(x, by=list(rownames(x)),sum) #now 1194 taxa
x$Group.1 <- gsub("\\[|\\]", "", x$Group.1)
x$Group.1 <- gsub("[.]", "", x$Group.1)
x$Group.1 <- gsub("s__", "", x$Group.1)
rownames(x) <- x$Group.1
x <- x[,!colnames(x) == "Group.1"]
x <- t(x)
x.raw <- x
x <- sweep(x, 1, rowSums(x), '/')

#Collapse by subject
m$ID_on_tube <- as.integer(m$ID_on_tube)
xc <- apply(x,2,function(xx) sapply(split(xx,m$ID_on_tube),mean))
rownames(xc) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

mc <- m[sapply(split(1:nrow(m),m$ID_on_tube),'[',1),,drop=TRUE]
rownames(mc) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

xc.raw <- apply(x.raw,2,function(xx) sapply(split(xx,m$ID_on_tube),sum))
rownames(xc.raw) <- sprintf('Subject_%03d',sapply(split(m$ID_on_tube,m$ID_on_tube),'[',1))

# drop rare bugs (show up in less than 10% of subjects)
# Goes from 1005 to 568 OTUs
rare.ix <- colMeans(xc > 0) < .10
x.raw <- x.raw[,!rare.ix]
xc.raw <- xc.raw[,!rare.ix]
x <- x[,!rare.ix]
xc <- xc[,!rare.ix]

m[m$Cohort == "Healthy", "Cohort"] <- "H"
mc[mc$Cohort == "Healthy", "Cohort"] <- "H"

# gather differentiation indicies
ixc.hc <- mc$Cohort == "H"
ix.hc <- m$Cohort == "H"
ixc.ibsc <- mc$Cohort == "C"
ix.ibsc <- m$Cohort == "C"
ixc.ibsd <- mc$Cohort == "D"
ix.ibsd <- m$Cohort == "D"
ixc.ibs <- mc$Cohort != "H"
ix.ibs <- m$Cohort != "H"

cols <- c("#cb1b4a", "#42aeb8", "#FDB316", "#c3c823", "#00797e", "#053058", "#aaada6", "#ae1848", "#368b90", "#ca9012", "#9ba11e", "#f1f2ec", "#d9d3df", "#348fbe", "#ff8340", "#ffAf40", "#bf503f", "#951b72", "#b75f6d")
cols2 <- colorRampPalette(cols)

mc$Cohort <- factor(mc$Cohort)
m$Cohort <- factor(m$Cohort)

#how many timelines are complete?
complete_timeline <- c(0, 0, 0)
names(complete_timeline) <- c("H", "C", "D")
for(i in 1:length(unique(m$study_id))){
  working <- unique(m$study_id)[i]
  if(length(which(m$study_id == working)) == 7){
    cohort <- as.character(m[which(m$study_id == working),"Cohort"])[1]
    complete_timeline[cohort] <- c(complete_timeline[cohort] + 1)
  }
}
print("Complete Timelines:")
print(complete_timeline)

#what are the average number of samples each person has?
average_timeline <- list(NA, NA, NA)
names(average_timeline) <- c("H", "C", "D")
for(i in 1:length(unique(m$study_id))){
  working <- unique(m$study_id)[i]
  cohort <- as.character(m[which(m$study_id == working),"Cohort"])[1]
  average_timeline[[cohort]]<- c(average_timeline[[cohort]], length(which(m$study_id == working)))
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
  cohort <- as.character(m[which(m$study_id == working),"Cohort"])[1]
  if(length(which(m$study_id == working)) < 3){
    average_timeline[[cohort]]<- average_timeline[[cohort]] + 1
  }
}
print("Number of subjects with < 3 samples:")
print(average_timeline)


# OLD Files
# m <- read.delim(mapfp, sep='\t',head=T,check=F,comment='', quote="")
# #Get the next sample map - bc having one map is too hard..?
# m2 <- read.delim(mapfp2, sep='\t',head=T,check=F,comment='',quote="")
# m2 <- m2[1:83,1:3]
# for(r in 1:nrow(m2)){
#   sample <- paste(m2[r,"Subject_ID"], m2[r,"ID_on_Tube"], sep=".")
#   m2$study_vial[r] <- sample
# }
# #M2 is 83 rows, m2 is 397, then 390 overlap
# m <- merge(m, m2, by="study_vial")

#x <- x[rownames(m),]
#x <- x[,colnames(x) != '_']
# # species only
# xsp <- x[,!grepl('_$',colnames(x))]
# colnames(xsp) <- gsub('_',' ',colnames(xsp))
# 
# # genus only
# genus <- sapply(strsplit(colnames(x),'_'),'[',1)
# xgn <- t(apply(x,1, function(xx) sapply(split(xx,genus),sum)))
# drop low-abundance samples
#xsp <- xsp[rowSums(xsp) > 1000,]
#xgn <- xgn[rownames(xsp),]


#Normalize
## xsp <- sweep(xsp, 1, rowSums(xgn), '/') # normalize species by genus-level assignments
#xsp.raw <- xsp
#xgn.raw <- xgn
#xsp <- sweep(xsp, 1, rowSums(xsp), '/')
#xgn <- sweep(xgn, 1, rowSums(xgn), '/')

# # collapse by subject
# xspc <- apply(xsp,2,function(xx) sapply(split(xx,m$Patient_no.),mean))
# rownames(xspc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
# xgnc <- apply(xgn,2,function(xx) sapply(split(xx,m$Patient_no.),mean))
# rownames(xgnc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))

# xspc.raw <- apply(xsp.raw,2,function(xx) sapply(split(xx,m$Patient_no.),sum))
# rownames(xspc.raw) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
# xgnc.raw <- apply(xgn.raw,2,function(xx) sapply(split(xx,m$Patient_no.),sum))
# rownames(xgnc.raw) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))


# rare.sp.ix <- colMeans(xspc > 0) < .25
# rare.gn.ix <- colMeans(xgnc > 0) < .25
# xsp <- xsp[,!rare.sp.ix]
# xgn <- xgn[,!rare.gn.ix]
# xspc <- xspc[,!rare.sp.ix]
# xgnc <- xgnc[,!rare.gn.ix]
# 
# xsp.raw <- xsp.raw[,!rare.sp.ix]
# xgn.raw <- xgn.raw[,!rare.gn.ix]
# xspc.raw <- xspc.raw[,!rare.sp.ix]
# xgnc.raw <- xgnc.raw[,!rare.gn.ix]

#difftest <- differentiation.test(xc, mc$Cohort)

# sink(sprintf('%s/map.txt',DATADIR)); cat('#SampleID\t'); write.table(m, sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/species.txt',DATADIR)); cat('Taxon\t'); write.table(t(xsp), sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/genus.txt',DATADIR)); cat('Taxon\t'); write.table(t(xgn), sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/map-collapsed.txt',DATADIR)); cat('#SampleID\t'); write.table(mc, sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/species-collapsed.txt',DATADIR)); cat('Taxon\t'); write.table(t(xspc), sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/genus-collapsed.txt',DATADIR)); cat('Taxon\t'); write.table(t(xgnc), sep='\t',quote=F); sink(NULL)
