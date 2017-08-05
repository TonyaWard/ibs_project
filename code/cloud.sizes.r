ALPHA <- .05

dotu <- as.matrix(vegdist(x))

# get all within-patient distances
m$ID_on_tube <- as.numeric(m$ID_on_tube)
patient.nos <- sort(unique(m$ID_on_tube))
wds.sp <- numeric(length(patient.nos))
names(wds.sp) <- sprintf('Subject_%03d',patient.nos)
for(i in 1:length(unique(m$ID_on_tube))){
  patient.no <- patient.nos[i]
  # within-patient distances
  patient.ix <- m$ID_on_tube == patient.no
  if(sum(patient.ix, na.rm=TRUE) <= 1){
    wds.sp[i] <- NA
  } else {
    # calculates all within-cloud distances
    tmp <- dotu[patient.ix,patient.ix]
    tmp <- tmp[upper.tri(tmp)]
    
    # calculates only distance to previous timepoint
    #tmp <- sapply(1:(sum(patient.ix)-1), function(ixx) dsp[which(patient.ix)[ixx],which(patient.ix)[ixx+1]])
    
    #total within-cloud distances
    wds.sp[i] <- mean(tmp, na.rm=TRUE)
    
    # calculates all within-cloud distances
    #tmp <- dgn[patient.ix,patient.ix]
    #tmp <- tmp[upper.tri(tmp)]
    # calculates only distance to previous timepoint
    # tmp <- sapply(1:(sum(patient.ix)-1), function(ixx) dgn[which(patient.ix)[ixx],which(patient.ix)[ixx+1]])
    #wds.gn[i] <- mean(tmp)
  }
}
# test species and genus distances
test.xs <- list(otu=wds.sp)
#test.xs <- list(species=wds.sp, genus=wds.gn)

# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
                )
# list of "reference" group names
compare.to <- list('HC v. IBS'='H',
                 'HC v. IBSD'='H',
                 'HC v. IBSC'='H',
                 'IBSC v. IBSD'='C'
                )
plot_by <- c("IBS", "Cohort", "Cohort", "Cohort")
col_list <- list(cols_ibs, cols_dh, cols_ch, cols)

# run all combinations of tests
for(i in 1:length(test.xs)){
    x.name <- names(test.xs)[i]
    test.x <- test.xs[[i]]

    for(j in 1:length(test.ixs)){
        test.name <- names(test.ixs)[j]
        compare.to.name <- compare.to[[test.name]]
        test.ix <- test.ixs[[j]] & !is.na(test.x)

        tt <- t.test(test.x[test.ix] ~ mc$Cohort[test.ix] != compare.to.name) 

        if(tt$p.value < ALPHA) {
          sink("beta_div/variability.txt", append=T)
          cat(sprintf('t-test variability, %s, %s: ', x.name, test.name))
          cat('p=',round(tt$p.value,4), ', t statistic=',round(tt$statistic,4),'\n',sep='')
          sink()
          working_mc <- subset(mc, test.ix)
          working_mc$vari <- subset(test.x, test.ix)
          plot1 <- ggplot(working_mc, aes_string(x=plot_by[j], y= "vari")) +
              geom_boxplot(outlier.shape = NA) +
              geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
              #theme(legend.position = 'bottom') + 
              labs(x="", y= "variability") +
              guides(fill=F, color=F) +
              scale_color_manual(values=col_list[[j]])
          pdf(sprintf('beta_div/variability_%s_%s.pdf',x.name, gsub(' ','_',test.name)),width=4,height=4, useDingbats=FALSE)
          print(plot1)
          dev.off()
        }
    }
}
