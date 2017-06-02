source('../code/stats.r')
source('../code/util.r')
source('../code/risk.index.r')

# source('load.r')
ALPHA <- 0.15

## test species and genus differences
#test.xs <- list(species=xspc, genus=xgnc)
test.xs <- list(otu=xc)
# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
                )
plot_by <- c("IBS", "Cohort", "Cohort", "Cohort")
col_list <- list(cols_ibs, cols_dh, cols_ch, cols)

# run all combinations of tests
for(i in 1:length(test.xs)){
    x.name <- names(test.xs)[i]
    test.x <- test.xs[[i]]

    for(j in 1:length(test.ixs)){
        test.name <- names(test.ixs)[j]
        test.ix <- test.ixs[[j]]

        cat(sprintf('Significant q-values for %s, %s:\n',x.name, test.name))
        difftest <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=TRUE)
        difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=FALSE)
        if(any(difftest$qvalues <= ALPHA)){
            signif.ix <- which(difftest$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
            pdf(sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
            for(k in signif.ix){
                if(!is.null(difftest$norm.test.pvals)){
                    norm.test <- difftest$norm.test.pvals[k]
                } else {
                    norm.test <- '0'
                }
                if(norm.test < 0.05){
                    qval <- difftest.np$qvalues[k]
                } else {
                    qval <- difftest$qvalues[k]
                }
                sink(sprintf('diff_taxa/differential_abundance_%s_%s.txt',x.name, test.name), append=T)
                cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',norm.test,'\n',sep=''))
                sink()
                # beeswarm(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                # # beeswarm(test.x[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                #         xlab='',ylab='Relative abundance', main=colnames(test.x)[k], col=cols[1:4])
                # bxplot(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
                # # bxplot(test.x[test.ix,col.ix] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
                
                working_mc <- subset(mc, test.ix)
                working_mc$vari <- data.transform(test.x)[test.ix,k]
                plot1 <- ggplot(working_mc, aes_string(x=plot_by[j], y= "vari")) +
                  geom_boxplot(outlier.shape = NA) +
                  geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
                  #theme(legend.position = 'bottom') + 
                  labs(x="", y= "Relative Abundance") +
                  ggtitle(colnames(test.x)[k]) +
                  guides(fill=F, color=F) +
                  scale_color_manual(values=col_list[[j]])
                print(plot1)
            }
            dev.off()
        } else {
            cat("None\n")
        }
    }
}
######################################################################
#Heatmap of diff taxa

## test species and genus differences
#test.xs <- list(species=xspc, genus=xgnc)
test.xs <- list(otu=xc)
# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
)
plot_by <- c("IBS", "Cohort", "Cohort", "Cohort")
col_list <- list(cols_ibs, cols_dh, cols_ch, cols)

i <- 1
x.name <- names(test.xs)[i]
test.x <- test.xs[[i]]
  
for(j in 1:length(test.ixs)){
  test.name <- names(test.ixs)[j]
  test.ix <- test.ixs[[j]]
  difftest <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=TRUE)
  difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=FALSE)
  if(any(difftest$qvalues <= ALPHA)){
    signif.ix <- which(difftest$qvalues <= ALPHA)
    signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
    #subset to keep only taxa sig.
    new_table <- data.transform(test.x)[test.ix,signif.ix]
    #new_table <- new_table[,colSums(new_table) > 0.1]
    new_map <- mc[test.ix,]
    #order according to highest in healthies
    #healthy_order <- colnames(new_table[ixc.hc,order(colSums(-new_table))])
    #new_table <- new_table[,healthy_order]
    new_table <- cbind(new_table, new_map[,c("Cohort", "IBS")])
    new_table$subject_id <- rownames(new_table)
    new_table <- melt(new_table, ids=c("Cohort", "IBS", "subject_id"))
    #new_table$variable <- factor(new_table$variable, levels=healthy_order)
    #new_table$subject_id <- factor(new_table$subject_id, levels=(new_table$subject_id)[order(new_table$Cohort)])
    ggplot(new_table, aes(x = reorder(subject_id, value), y = reorder(variable, value))) + 
      geom_tile(aes(fill = value)) + 
      scale_fill_gradient(name = 'relative abundance', low = 'white', high = 'purple') +
      facet_grid(~Cohort) +
      theme(axis.title.y = element_blank(), axis.text.x=element_blank(), axis.title.x = element_blank())
  }
}


######################################################################
#Test for pain differences
## test species and genus differences
#test.xs <- list(species=xspc, genus=xgnc)
test.xs <- list(otu=x)
m$q10_base_q4_symp <- as.numeric(droplevels(m$q10_base_q4_symp))

# different group comparisons
D.P <- m$q10_base_q4_symp == 1 & m$Cohort == "D" & !is.na(m$q10_base_q4_symp)
D.NP <- m$q10_base_q4_symp == 2 & m$Cohort == "D" & !is.na(m$q10_base_q4_symp)
C.P <-  m$q10_base_q4_symp == 1 & m$Cohort == "C" & !is.na(m$q10_base_q4_symp)
C.NP <- m$q10_base_q4_symp == 2 & m$Cohort == "C" & !is.na(m$q10_base_q4_symp)
H.P <-  m$q10_base_q4_symp == 1 & m$Cohort == "H" & !is.na(m$q10_base_q4_symp)
H.NP <- m$q10_base_q4_symp == 2 & m$Cohort == "H" & !is.na(m$q10_base_q4_symp)

test.ixs <- list('IBSD - P v. IBSD - NP'=  D.P | D.NP,
                 'IBSC - P v. IBSC - NP'= C.P | C.NP
)
m$q10_base_q4_symp <- factor(m$q10_base_q4_symp)

plot_by <- c("q10_base_q4_symp", "q10_base_q4_symp")
col_list <- list(cols_yn, cols_yn)

# run all combinations of tests
for(i in 1:length(test.xs)){
  x.name <- names(test.xs)[i]
  test.x <- test.xs[[i]]
  
  for(j in 1:length(test.ixs)){
    test.name <- names(test.ixs)[j]
    test.ix <- test.ixs[[j]]
    
    cat(sprintf('Significant q-values for %s, %s:\n',x.name, test.name))
    difftest <- differentiation.test(data.transform(test.x)[test.ix,], m$q10_base_q4_symp[test.ix], parametric=TRUE)
    difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], m$q10_base_q4_symp[test.ix], parametric=FALSE)
    if(any(na.omit(difftest$qvalues <= ALPHA))){
      signif.ix <- which(difftest$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
      pdf(sprintf('differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
      for(k in signif.ix){
        if(!is.null(difftest$norm.test.pvals)){
          norm.test <- difftest$norm.test.pvals[k]
        } else {
          norm.test <- '0'
        }
        if(norm.test < 0.05){
          qval <- difftest.np$qvalues[k]
        } else {
          qval <- difftest$qvalues[k]
        }
        
        cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',norm.test,'\n',sep=''))
        
        working_m <- subset(m, test.ix)
        working_m$vari <- data.transform(test.x)[test.ix,k]
        plot1 <- ggplot(working_m, aes_string(x=plot_by[j], y= "vari")) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
          #theme(legend.position = 'bottom') + 
          labs(x="", y= "Relative Abundance") +
          ggtitle(colnames(test.x)[k]) +
          guides(fill=F, color=F) +
          scale_color_manual(values=col_list[[j]])
        print(plot1)
      }
      dev.off()
    } else {
      cat("None\n")
    }
  }
}

######################################################################
#Test for distention differences
## test species and genus differences
#test.xs <- list(species=xspc, genus=xgnc)
test.xs <- list(otu=x)
m$q13_base_q7_symp <- as.numeric(droplevels(m$q13_base_q7_symp))

# different group comparisons
D.D <- m$q13_base_q7_symp == 1 & m$Cohort == "D" & !is.na(m$q13_base_q7_symp)
D.ND <- m$q13_base_q7_symp == 2 & m$Cohort == "D" & !is.na(m$q13_base_q7_symp)
C.D <-  m$q13_base_q7_symp == 1 & m$Cohort == "C" & !is.na(m$q13_base_q7_symp)
C.ND <- m$q13_base_q7_symp == 2 & m$Cohort == "C" & !is.na(m$q13_base_q7_symp)

test.ixs <- list('IBSD - D v. IBSD - DP'=  D.D | D.ND,
                 'IBSC - D v. IBSC - DP'= C.D | C.ND
)
m$q13_base_q7_symp <- factor(m$q13_base_q7_symp)

plot_by <- c("q13_base_q7_symp", "q13_base_q7_symp")
col_list <- list(cols_yn, cols_yn)

# run all combinations of tests
for(i in 1:length(test.xs)){
  x.name <- names(test.xs)[i]
  test.x <- test.xs[[i]]
  
  for(j in 1:length(test.ixs)){
    test.name <- names(test.ixs)[j]
    test.ix <- test.ixs[[j]]
    
    cat(sprintf('Significant q-values for %s, %s:\n',x.name, test.name))
    difftest <- differentiation.test(data.transform(test.x)[test.ix,], m$q13_base_q7_symp[test.ix], parametric=TRUE)
    difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], m$q13_base_q7_symp[test.ix], parametric=FALSE)
    if(any(na.omit(difftest$qvalues <= ALPHA))){
      signif.ix <- which(difftest$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
      pdf(sprintf('differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
      for(k in signif.ix){
        if(!is.null(difftest$norm.test.pvals)){
          norm.test <- difftest$norm.test.pvals[k]
        } else {
          norm.test <- '0'
        }
        if(norm.test < 0.05){
          qval <- difftest.np$qvalues[k]
        } else {
          qval <- difftest$qvalues[k]
        }
        
        cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',norm.test,'\n',sep=''))
        working_m <- subset(m, test.ix)
        working_m$vari <- data.transform(test.x)[test.ix,k]
        plot1 <- ggplot(working_m, aes_string(x=plot_by[j], y= "vari")) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
          #theme(legend.position = 'bottom') + 
          labs(x="", y= "Relative Abundance") +
          ggtitle(colnames(test.x)[k]) +
          guides(fill=F, color=F) +
          scale_color_manual(values=col_list[[j]])
        print(plot1)
      }
      dev.off()
    } else {
      cat("None\n")
    }
  }
}