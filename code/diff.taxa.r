#### Differential taxa testing ####
bT <- c(2,3,4,5,6,7)  # Levels to consider for the bugs
for (L in 1:length(bT)) {
  # Massage the taxa names
  split <- strsplit(colnames(x),";")        # Split by semicolon into levels
  taxaStrings <- sapply(split,function(x) paste(x[1:bT[L]],collapse=";"))           # Create collapsed names
  for (i in 1:7) taxaStrings <- gsub("(;[A-z]__$)?(;NA$)?","",taxaStrings,perl=T)   # Clean tips
  x.t <- x
  x.t <- data.frame(t(x.t), check.names = F)
  x.t$taxonomy <- taxaStrings
  x.t = rowsum(x.t[,-ncol(x.t)],taxaStrings) #Collapse
  x.t[is.na(x.t)] = 0                                       # pad empty columns with 0 counts
 
  # Filter out bugs that don't appear in enough samples
  select = rowSums(x.t > 0) > min(table(m$Cohort))/2 # reasonable general min
  x.t = x.t[select,]                                    # Apply drop mask
  x.t = x.t[,rownames(m)]                             # Sync with map samp order
  x.t = sweep(x.t,2,colSums(x.t),'/');                # Normalize to relative abundance

  otu.c = impRZilr(t(x.t), dl=rep(1,nrow(x.t)), maxit = 3, verbose = T, method = "lm") # zeros
  otu.c = t(cenLR(otu.c$x)$x.clr)    # Centered log-ratio transform for compositions
  colnames(otu.c) = colnames(otu.t)  # Because this gets rid of the names...
  otu.t = otu.c                      # otu.t is our active table; give it the CLR
  
  # Go through each taxon and test for significance w/group
  otu.t = as.matrix(otu.t)
  ntax = nrow(otu.t)
  CW = map$CaptiveWild %in% c("Captive","Wild") # For selecting only captive and wild
  Grp.Pvals=rep(1,ntax)
  Grp.Corrs=rep(0,ntax)
  Wld.Pvals=rep(1,ntax)
  KW.Pvals=rep(1,ntax)
  for (m.ix in 1:ntax) {  # Loop through all the rows (taxa)
    try({ # Because some correlations may be inadmissable
      ps = polyserial(otu.t[m.ix,],map$PA,ML=T,std.err = T)
      if (is.na(pchisq(ps$chisq, ps$df))) next # Invalid correlation
      Grp.Corrs[m.ix] = ps$rho             # Find intensity of correlation
      Grp.Pvals[m.ix] = 1-pchisq(ps$chisq, ps$df) # And p-value on this
    },silent=T)
    Wld.Pvals[m.ix] = wilcox.test(otu.t[m.ix,CW] ~ map$CaptiveWild[CW])$p.value
    KW.Pvals[m.ix] = kruskal.test(otu.t[m.ix,] ~ map$Lifestyle)$p.val
  }
  
  # Taxa barplots -- Top 15 most abundant (kruskal sig. + other?)
  otu.m = otu.n[-grep("Viridiplantae",rownames(otu.n)),,drop=F]
  otu.m = sweep(sqrt(otu.m),2,colSums(sqrt(otu.m)),'/')
  meanAb = apply(otu.m,1,FUN=function(x) tapply(x, map$Lifestyle, mean)) # group mean
  
  #select = apply(meanAb,2,max) > 0.01
  #KW.Pvals[!select] = 1  # coerce poor represenatives to low confidence
  #ranked = order(KW.Pvals)
  ranked = order(apply(meanAb,2,max),decreasing=T)
  otu.m = otu.m[ranked,]
  
  Taxa = rownames(otu.m)
  #lims = which(KW.Pvals[ranked]==1)
  #lim = ifelse(length(lims),min(25,lims[1]),25)
  lim = 25
  if (nrow(otu.m) > lim) Taxa[lim:nrow(otu.m)] = "Other"
  otu.m = rowsum(otu.m,Taxa)
  byAbundance = rownames(otu.m)[order(rowMeans(otu.m),decreasing=T)]
  #byAbundance = gsub(';','.',byAbundance)
  
  otu.m = data.frame(t(otu.m),check.names=F)      # flip table
  otu.m$SampleID = rownames(otu.m)  # add a column for the sample IDs
  map$SampleID = rownames(map)      # add a column for the sample IDs
  
  # The following separates taxa abundances per sample, then splices in the column of interest
  otu.m = melt(otu.m, id.vars = "SampleID", variable.name = "Taxa", value.name = "RelativeAbundance")
  otu.m = merge(otu.m, map[,c("SampleID","Lifestyle")], by="SampleID")
  otu.m$Taxa = factor(otu.m$Taxa,levels=byAbundance,ordered=T)
  
  ## Plot according to Lifestyle, sorted by abundance
  pdf(paste0("TaxaSummary_L",bT[L],".pdf"),width = 12,height=8) # Make room for legend
  plot(ggplot(otu.m, aes(x = Lifestyle, y = RelativeAbundance, fill = Taxa)) +
         geom_bar(stat ="identity", position="fill") + labs(x="Lifestyle",y="Root Relative Abundance") +
         guides(fill=guide_legend(ncol=1)) +
         scale_fill_manual(values=c("dodgerblue2","#E31A1C", # red # Kevin Wright
                                    "green4",
                                    "#6A3D9A", # purple
                                    "#FF7F00", # orange
                                    "black","gold1",
                                    "skyblue2","#FB9A99", # lt pink
                                    "palegreen2",
                                    "#CAB2D6", # lt purple
                                    "#FDBF6F", # lt orange
                                    "gray70", "khaki2",
                                    "maroon","orchid1","deeppink1","blue1","steelblue4",
                                    "darkturquoise","green1","yellow4","yellow3",
                                    "darkorange4","brown"))) 
  dev.off()
  
  # Adjust for multiple tests, sort by significance
  Grp.Pvals = p.adjust(Grp.Pvals)
  Wld.Pvals = p.adjust(Wld.Pvals)
  res = data.frame(Grp.Pvals, Grp.Corrs, Wld.Pvals,row.names=rownames(otu.t))
  res = res[order(res$Grp.Corrs),]
  
  # Add bivariate filter
  sig = 0.05
  selection = res$Grp.Pvals < sig & res$Wld.Pvals < sig & abs(res$Grp.Corrs) > 0.25
  
  # Display all significant with p < 0.05
  num_sig = sum(selection, na.rm = T) # Count how many are significant
  res = res[selection,]
  C_ix = map$CaptiveWild=="Captive"          # Stores "true" if monkey is captive, else "false"
  W_ix = map$CaptiveWild=="Wild"             # As above. Use these to select just wild/captive
  pdf(paste0("TaxaSwarms_L",bT[L],".pdf"),width = 8,height=7)
  sink(paste0("Taxa_Significance_L",bT[L],".txt"))                  # Get ready to write the significant ones
  cat("Taxon\tPolyserial_Q\tPolyserial_Cor\tCaptiveVsWild_Q\tTrendInCaptivity\n")  # Print header
  if (num_sig) for (i in 1:num_sig) {
    taxon = rownames(res)[i]
    upInCaptive = mean(otu.t[taxon,C_ix]) > mean(otu.t[taxon,W_ix]) # compare avgs
    cat(taxon,'\t',res$Grp.Pvals[i],'\t',-res$Grp.Corrs[i],'\t',res$Wld.Pvals[i],'\t',
        ifelse(upInCaptive,"UP","DOWN"),'\n',sep='')
    beeswarm(otu.t[taxon,] ~ map$PA, xlab="Lifestyle",ylab="CLR Relative Abundance",main=taxon,
             col=alpha(lscolors,0.7),cex.axis=1.1,cex.main=0.75,cex=1.1,corral="random",pch=19)
    bxplot(otu.t[taxon,] ~ map$PA, add = TRUE)
  }
  sink(NULL)
  dev.off()
  
  ####################################################################
  #heatmap of diff taxa
  
  taxa_list <- c("s__[Clostridium]_glycyrrhizinilyticum",
                 "s__[Clostridium]_innocuum",
                 "s__[Ruminococcus]_gnavus",
                 "s__[Ruminococcus]_torques",
                 "s__Adlercreutzia_equolifaciens",
                 "s__Alistipes_finegoldii",
                 "s__Alistipes_ihumii",
                 "s__Alistipes_putredinis",
                 "s__Alistipes_senegalensis",
                 "s__Alistipes_shahii",
                 "s__Anaerostipes_hadrus",
                 "s__Bifidobacterium_animalis",
                 "s__Coprococcus_sp",
                 "s__Eisenbergiella_tayi",
                 "s__Faecalitalea_cylindroides",
                 "s__Fournierella_massiliensis",
                 "s__Intestinimonas_butyriciproducens",
                 "s__Intestinimonas_massiliensis",
                 "s__Methanobrevibacter_smithii",
                 "s__Oscillibacter_sp",
                 "s__Roseburia_intestinalis",
                 "s__Ruminococcaceae_bacterium_D16",
                 "s__Ruminococcus_sp",
                 "s__Ruthenibacterium_lactatiformans",
                 "s__Streptococcus_parasanguinis",
                 "s__Subdoligranulum_variabile",)