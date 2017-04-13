library('beeswarm')
source('../code/stats.r')
source('../code/util.r')
source('../code/risk.index.r')

# source('load.r')
ALPHA <- 0.25

## test species and genus differences
#test.xs <- list(species=xspc, genus=xgnc)
test.xs <- list(otu=xc)
# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
                )

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
                beeswarm(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                # beeswarm(test.x[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                        xlab='',ylab='Relative abundance', main=colnames(test.x)[k], col=cols[1:4])
                bxplot(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
                # bxplot(test.x[test.ix,col.ix] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
            }
            dev.off()
        } else {
            cat("None\n")
        }
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
                 'IBSC - P v. IBSC - NP'= C.P | C.NP,
                 'HC - P v. HC - NP'= H.P | H.NP
)

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
    if(any(difftest$qvalues <= ALPHA)){
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
        beeswarm(data.transform(test.x)[test.ix,k] ~ droplevels(m$q10_base_q4_symp[test.ix]),
                 xlab='',ylab='Relative abundance', main=colnames(test.x)[k], col=cols[1:4])
        bxplot(data.transform(test.x)[test.ix,k] ~ droplevels(m$q10_base_q4_symp[test.ix]),add=TRUE)
        # bxplot(test.x[test.ix,col.ix] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
      }
      dev.off()
    } else {
      cat("None\n")
    }
  }
}
