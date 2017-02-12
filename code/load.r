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
mapfp <- paste(DATADIR,'170118_updated.txt',sep='')
mapfp2 <- paste(DATADIR, '170118_Sample_matrix.txt', sep='')
taxfp <- paste(DATADIR,'taxa_counts.txt',sep='')

#Sample IDs don't match the map - AT ALL.
#Reformat the names of the samples
m <- read.delim(mapfp, sep='\t',head=T,check=F,comment='', quote="")
for(r in 1:nrow(m)){
  sample <- paste(m[r,"Vial_ID"], m[r,"Timepoint"], sep=".T.")
  m$SampleID[r] <- sample
}
for(r in 1:nrow(m)){
  sample <- paste(m[r,"study_id"], m[r,"Vial_ID"], sep=".")
  m$study_vial[r] <- sample
}

#Get the next sample map - bc having one map is too hard..?
m2 <- read.delim(mapfp2, sep='\t',head=T,check=F,comment='',quote="")
m2 <- m2[1:83,1:3]
for(r in 1:nrow(m2)){
  sample <- paste(m2[r,"Subject_ID"], m2[r,"ID_on_Tube"], sep=".")
  m2$study_vial[r] <- sample
}
#M2 is 83 rows, m2 is 397, then 390 overlap
m <- merge(m, m2, by="study_vial")


#Sample IDs don't match the map - AT ALL.
#Reformat the names of the samples
x <- t(read.table(taxfp, sep='\t',head=T,row=1,check=F,comment=''))
sample_names <- read.table(text = rownames(x), sep = ".", colClasses = "character")
for(r in 1:nrow(sample_names)){
  sample <- paste(sample_names[r,3], sample_names[r,4], sample_names[r,5], sep=".")
  rownames(x)[r] <- sample
}

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
#From 375 samples to 375 (none dropped)
x <- x[rowSums(x) > 1000,]
m <- as.matrix(m) #This is because there are duplicates in the map
rownames(m) <- m[,"SampleID"]
m <- as.data.frame(m, stringsAsFactors = F)
m <- m[rownames(x),]

## normalize 
## xsp <- sweep(xsp, 1, rowSums(xgn), '/') # normalize species by genus-level assignments
#xsp.raw <- xsp
#xgn.raw <- xgn
x.raw <- x
x <- sweep(x, 1, rowSums(x), '/')
#xsp <- sweep(xsp, 1, rowSums(xsp), '/')
#xgn <- sweep(xgn, 1, rowSums(xgn), '/')

# # collapse by subject
# xspc <- apply(xsp,2,function(xx) sapply(split(xx,m$Patient_no.),mean))
# rownames(xspc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
# xgnc <- apply(xgn,2,function(xx) sapply(split(xx,m$Patient_no.),mean))
# rownames(xgnc) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
m$Vial_ID <- as.integer(m$Vial_ID)
xc <- apply(x,2,function(xx) sapply(split(xx,m$Vial_ID),mean))
rownames(xc) <- sprintf('Subject_%03d',sapply(split(m$Vial_ID,m$Vial_ID),'[',1))

mc <- m[sapply(split(1:nrow(m),m$Vial_ID),'[',1),,drop=TRUE]
rownames(mc) <- sprintf('Subject_%03d',sapply(split(m$Vial_ID,m$Vial_ID),'[',1))

# xspc.raw <- apply(xsp.raw,2,function(xx) sapply(split(xx,m$Patient_no.),sum))
# rownames(xspc.raw) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
# xgnc.raw <- apply(xgn.raw,2,function(xx) sapply(split(xx,m$Patient_no.),sum))
# rownames(xgnc.raw) <- sprintf('Subject_%03d',sapply(split(m$Patient_no.,m$Patient_no.),'[',1))
xc.raw <- apply(x.raw,2,function(xx) sapply(split(xx,m$Vial_ID),sum))
rownames(xc.raw) <- sprintf('Subject_%03d',sapply(split(m$Vial_ID,m$Vial_ID),'[',1))


# drop rare bugs (show up in less than 10% of subjects)
# Goes from 1005 to 570 OTUs
rare.ix <- colMeans(xc > 0) < .10
x.raw <- x.raw[,!rare.ix]
xc.raw <- xc.raw[,!rare.ix]
x <- x[,!rare.ix]
xc <- xc[,!rare.ix]

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

# gather differentiation indicies
ixc.hc <- mc$Cohort == "Healthy"
ix.hc <- m$Cohort == "Healthy"
ixc.ibsc <- mc$Cohort == "IBS-C"
ix.ibsc <- m$Cohort == "IBS-C"
ixc.ibsd <- mc$Cohort == "IBS-D"
ix.ibsd <- m$Cohort == "IBS-D"
ixc.ibs <- mc$Cohort != "Healthy"
ix.ibs <- m$Cohort != "Healthy"

cols <- c("#cb1b4a", "#42aeb8", "#FDB316", "#c3c823", "#00797e", "#053058", "#aaada6", "#ae1848", "#368b90", "#ca9012", "#9ba11e", "#f1f2ec", "#d9d3df", "#348fbe", "#ff8340", "#ffAf40", "#bf503f", "#951b72", "#b75f6d")
cols2 <- colorRampPalette(cols)

mc$Cohort <- factor(mc$Cohort)
m$Cohort <- factor(m$Cohort)

#difftest <- differentiation.test(xc, mc$Cohort)

# sink(sprintf('%s/map.txt',DATADIR)); cat('#SampleID\t'); write.table(m, sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/species.txt',DATADIR)); cat('Taxon\t'); write.table(t(xsp), sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/genus.txt',DATADIR)); cat('Taxon\t'); write.table(t(xgn), sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/map-collapsed.txt',DATADIR)); cat('#SampleID\t'); write.table(mc, sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/species-collapsed.txt',DATADIR)); cat('Taxon\t'); write.table(t(xspc), sep='\t',quote=F); sink(NULL)
# sink(sprintf('%s/genus-collapsed.txt',DATADIR)); cat('Taxon\t'); write.table(t(xgnc), sep='\t',quote=F); sink(NULL)
