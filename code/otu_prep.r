###Script to process the OTU tables to proper formats

#Load file paths for inputs
DATADIR <- '../data/IBS_Mayo_secondrun/'
metadata <- paste(DATADIR, "Final_Cleaned_Condensed_Metadata.csv", sep='')
taxfp <- paste(DATADIR,'shogun_output/taxatable.strain.txt',sep='')

m <- read.delim(metadata, 
                   sep=',', 
                   head=T, 
                   comment="")
rownames(m) <- m$SampleID

######################################################################
#Read in the OTU table, reformat the names of the samples, remove duplicates
x <- t(read.delim(taxfp, row=1)) #532 samples by 2783 OTUs
taxa_names <- colnames(x) #store this because it replaces ";" with "." later

#clean duplicates from the run that was messed up
rownames(x) <- gsub("_S[[:digit:]]+_R1_001.fa","",rownames(x))
rownames(x) <- gsub("_", ".", rownames(x))
x <- data.frame(x[order(rownames(x)),]) # Sort nicely by sample ID, OLDs will be second
x$samples <- rownames(x)
x$samples <- gsub(".S[[:digit:]]+.R1.001.OLD.fa","",x$samples) # remove OLD and NPM
x$samples <- gsub(".NPM", "", x$samples) #532

x$samples <- gsub("Study.ID.", "", x$samples) #now 485 samples
x2 <- x[match(unique(x$samples), x$samples),] #keep first instance of dups (new run)
rownames(x2) <- x2$samples   
x <- x2[,!colnames(x2)=="samples"]
rownames(x) <- gsub(".26.T.0", "26.T.0", rownames(x))

#Find samples not accounted for:
keeps <- intersect(rownames(x), rownames(m)) #480 samples overlap between the OTU table and metadata
not_in_map <- rownames(x)[which(!rownames(x) %in% rownames(m))]
cat("These samples aren't in the map:")
not_in_map
#[1] "Negative.Control.H12.100" "Pos.Ctrl" "Positive.Control.B01.101" "Positive.Control.G12.100" 

not_in_otu <- rownames(m)[which(!rownames(m) %in% rownames(x))]
cat("These samples aren't in the otu table:")
not_in_otu
#"33.T.1" "50.T.1" "51.T.0" "59.T.1" "63.T.4" "65.T.3" "65.T.4" "65.T.6" "58.T.0"
#Now do this for all the files!

files <- list.files("../data/IBS_Mayo_secondrun/shogun_output/")[3:9]
for(i in 1:length(files)){
  file_now <- paste("../data/IBS_Mayo_secondrun/shogun_output/", files[i], sep="")
  x <- t(read.delim(file_now, row=1))
  taxa_names <- colnames(x) #store this because it replaces ";" with "." later
  
  #clean duplicates from the run that was messed up
  rownames(x) <- gsub("_S[[:digit:]]+_R1_001.fa","",rownames(x))
  rownames(x) <- gsub("_", ".", rownames(x))
  x <- data.frame(x[order(rownames(x)),]) # Sort nicely by sample ID, OLDs will be second
  x$samples <- rownames(x)
  x$samples <- gsub(".S[[:digit:]]+.R1.001.OLD.fa","",x$samples) # remove OLD and NPM
  x$samples <- gsub(".NPM", "", x$samples) #532
  
  x$samples <- gsub("Study.ID.", "", x$samples) #now 485 samples
  x2 <- x[match(unique(x$samples), x$samples),] #keep first instance of dups (new run)
  rownames(x2) <- x2$samples   
  x <- x2[,!colnames(x2)=="samples"]
  rownames(x) <- gsub(".26.T.0", "26.T.0", rownames(x))
  colnames(x) <- taxa_names
  x <- t(x)
  
  file_name <- paste("../data/IBS_Mayo_secondrun/shogun_output/id_update/", files[i], sep="")
  cat("#OTUID\t", file=file_name)
  write.table(x, file=file_name, append=T,quote=F, row.names=T, sep="\t")
}
