#Test for Differential Taxa

# source('load.r')
ALPHA <- 0.1

test.xs <- list(otu=xc[,which(colMeans(xc) > 0.0005)])

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
        #difftest <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=TRUE)
        difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], mc[test.ix,plot_by[j]], parametric=FALSE)
        if(any(difftest.np$qvalues <= ALPHA)){
            signif.ix <- which(difftest.np$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
            res <- data.frame(Taxon=colnames(test.x)[signif.ix],
                              qvalue=round(difftest.np$qvalues[signif.ix],5),
                              pvalue=round(difftest.np$pvalues[signif.ix],5),
                              round(difftest.np$classwise.means[signif.ix,,drop=F],5))
            sink(sprintf('diff_taxa/differential_abundance_%s_%s.txt',x.name, test.name))
            #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
            write.table(res,row.names=F,sep='\t',quote=F)
            sink()
            pdf(sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
            for(k in signif.ix){
                # if(!is.null(difftest$norm.test.pvals)){
                #     norm.test <- difftest$norm.test.pvals[k]
                # } else {
                #     norm.test <- '0'
                # }
                # if(norm.test < 0.05){
                #     qval <- paste(difftest.np$qvalues[k], "np", sep= " ")
                # } else {
                #     qval <- difftest$qvalues[k]
                # }
              qval <- difftest.np$qvalues[k]
                # beeswarm(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                # # beeswarm(test.x[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),
                #         xlab='',ylab='Relative abundance', main=colnames(test.x)[k], col=cols[1:4])
                # bxplot(data.transform(test.x)[test.ix,k] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
                # # bxplot(test.x[test.ix,col.ix] ~ droplevels(mc$Cohort[test.ix]),add=TRUE)
                
                working_mc <- subset(mc, test.ix)
                working_mc$vari <- data.transform(test.x)[test.ix,k]
                #pdf(file=sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name))
                
                plot1 <- ggplot(working_mc, aes_string(x=plot_by[j], y= "vari")) +
                  geom_boxplot(outlier.shape = NA) +
                  geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
                  #theme(legend.position = 'bottom') + 
                  labs(x="", y= "Relative Abundance (sq rt)", title =  colnames(test.x)[k]) +
                  #ggtitle(colnames(test.x)[k]) +
                  theme(plot.title = element_text(size=10)) +
                  guides(fill=F, color=F) +
                  scale_y_continuous(trans='sqrt') +
                  scale_color_manual(values=col_list[[j]])
                plot(plot1)
            }
            dev.off()
            cat("see file\n")
        } else {
            cat("None\n")
        }
    }
}

####################################################################
#heatmap of diff taxa

taxa_list <- c("s__[Clostridium]_glycyrrhizinilyticum",
               "s__[Clostridium]_innocuum",
               "s__Alistipes_obesi",
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
               "s__Coprococcus_sp._HPP0048",
               "s__Eisenbergiella_tayi",
               "s__Faecalitalea_cylindroides",
               "s__Fournierella_massiliensis",
               "s__Intestinimonas_butyriciproducens",
               "s__Intestinimonas_massiliensis",
               "s__Methanobrevibacter_smithii",
               "s__Oscillibacter_sp._ER4",
               "s__Roseburia_intestinalis",
               "s__Ruminococcaceae_bacterium_D16",
               "s__Ruminococcus_sp._AT10",
               "s__Ruthenibacterium_lactatiformans",
               "s__Streptococcus_parasanguinis",
               "s__Subdoligranulum_variabile")

taxa_sub <- x[,taxa_list]
taxa_sub <- melt(taxa_sub)
colnames(taxa_sub) <- c("SampleID", "Taxa", "Relative_Abundance")
taxa_sub <- merge(taxa_sub, m[,c("Cohort", "SampleID")], by="SampleID")
taxa_sub$Cohort <- factor(taxa_sub$Cohort, levels=c("C", "H", "D"))
taxa_sub$Taxa <- gsub("s__", "", taxa_sub$Taxa)
taxa_sub$Taxa <- factor(taxa_sub$Taxa, levels=rev(c(
  "[Clostridium]_innocuum",
  "Ruthenibacterium_lactatiformans",
  "Alistipes_finegoldii",
  "Ruminococcus_sp._AT10",
  "[Ruminococcus]_torques",
  "Streptococcus_parasanguinis",
  "Alistipes_obesi",
  "Fournierella_massiliensis",
  "Intestinimonas_massiliensis",
  "Alistipes_senegalensis",
  "Alistipes_shahii",
  "Coprococcus_sp._HPP0048",
  "[Ruminococcus]_gnavus",
  "Oscillibacter_sp._ER4",
  "Intestinimonas_butyriciproducens",
  "[Clostridium]_glycyrrhizinilyticum",
  "Eisenbergiella_tayi",
  "Subdoligranulum_variabile",
  "Alistipes_ihumii",
  "Anaerostipes_hadrus",
  "Adlercreutzia_equolifaciens",
  "Alistipes_putredinis",
  "Methanobrevibacter_smithii",
  "Roseburia_intestinalis",
  "Ruminococcaceae_bacterium_D16",
  "Bifidobacterium_animalis",
  "Faecalitalea_cylindroides"
)))
taxa_sub$Relative_Abundance2 <- cut(taxa_sub$Relative_Abundance, c(0,0.1,0.2,0.3,0.4))
pdf("diff_taxa/heat_map.pdf", height = 6, width=5)
ggplot(taxa_sub, aes(x=SampleID, y=Taxa))+
  geom_tile(aes(fill=Cohort, alpha=Relative_Abundance)) +
  facet_grid(. ~ Cohort, scales="free", space="free") +
  scale_fill_manual(values=cols[c(1,3,2)]) +
  theme(axis.text.x =element_blank()) +
  guides(alpha=F, fill=F)
dev.off()

######################################################################
#Test for flare differences
test.xs <- list(otu=x[,which(colMeans(x) > 0.0005)])

# different group comparisons
D.flare <- m$Cohort == "D" & !is.na(m$Flare)
D.Nflare <- m$Cohort == "D" & is.na(m$Flare)
C.flare <-  m$Cohort == "C" & !is.na(m$Flare)
C.Nflare <- m$Cohort == "C" & is.na(m$Flare)

test.ixs <- list('IBSDflare v. IBSDNF'=  D.flare | D.Nflare,
                 'IBSCflare v. IBSCNF'= C.flare | C.Nflare
)
m$Flare <- as.character(m$Flare)
m[is.na(m$Flare),"Flare"] <- "Normal"
m$Flare <- factor(m$Flare)

plot_by <- c("Flare", "Flare")
col_list <- list(cols_yn, cols_yn)

# run all combinations of tests
for(i in 1:length(test.xs)){
  x.name <- names(test.xs)[i]
  test.x <- test.xs[[i]]
  
  for(j in 1:length(test.ixs)){
    test.name <- names(test.ixs)[j]
    test.ix <- test.ixs[[j]]
    cat(sprintf('Significant q-values for %s, %s:\n',x.name, test.name))
    difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], m[test.ix,plot_by[j]], parametric=FALSE)
    if(any(difftest.np$qvalues <= ALPHA)){
      signif.ix <- which(difftest.np$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
      res <- data.frame(Taxon=colnames(test.x)[signif.ix],
                        qvalue=round(difftest.np$qvalues[signif.ix],5),
                        pvalue=round(difftest.np$pvalues[signif.ix],5),
                        round(difftest.np$classwise.means[signif.ix,,drop=F],5))
      sink(sprintf('diff_taxa/differential_abundance_%s_%s.txt',x.name, test.name))
      #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
      write.table(res,row.names=F,sep='\t',quote=F)
      sink()
      pdf(sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
      for(k in signif.ix){
        qval <- difftest.np$qvalues[k]
        working_m <- subset(m, test.ix)
        working_m$vari <- data.transform(test.x)[test.ix,k]
        #pdf(file=sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name))
        plot1 <- ggplot(working_m, aes_string(x=plot_by[j], y= "vari")) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
          #theme(legend.position = 'bottom') + 
          labs(x="", y= "Relative Abundance (sq rt)", title =  colnames(test.x)[k]) +
          theme(plot.title = element_text(size=10)) +
          guides(fill=F, color=F) +
          scale_y_continuous(trans='sqrt') +
          scale_color_manual(values=col_list[[j]])
        plot(plot1)
      }
      dev.off()
      cat("see file\n")
    } else {
      cat("None\n")
    }
  }
}

######################################################################
#Test for pain differences
test.xs <- list(otu=x[,which(colMeans(x) > 0.0005)])
m$q10_base_q4_symp <- as.numeric(droplevels(m$q10_base_q4_symp))

# different group comparisons
D.P <- m$q10_base_q4_symp == 1 & m$Cohort == "D" & !is.na(m$q10_base_q4_symp)
D.NP <- m$q10_base_q4_symp == 2 & m$Cohort == "D" & !is.na(m$q10_base_q4_symp)
C.P <-  m$q10_base_q4_symp == 1 & m$Cohort == "C" & !is.na(m$q10_base_q4_symp)
C.NP <- m$q10_base_q4_symp == 2 & m$Cohort == "C" & !is.na(m$q10_base_q4_symp)
H.P <-  m$q10_base_q4_symp == 1 & m$Cohort == "H" & !is.na(m$q10_base_q4_symp)
H.NP <- m$q10_base_q4_symp == 2 & m$Cohort == "H" & !is.na(m$q10_base_q4_symp)

test.ixs <- list('IBSDP v. IBSDNP'=  D.P | D.NP,
                 'IBSCP v. IBSCNP'= C.P | C.NP
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
    difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], m[test.ix,plot_by[j]], parametric=FALSE)
    if(any(difftest.np$qvalues <= ALPHA)){
      signif.ix <- which(difftest.np$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
      res <- data.frame(Taxon=colnames(test.x)[signif.ix],
                        qvalue=round(difftest.np$qvalues[signif.ix],5),
                        pvalue=round(difftest.np$pvalues[signif.ix],5),
                        round(difftest.np$classwise.means[signif.ix,,drop=F],5))
      sink(sprintf('diff_taxa/differential_abundance_%s_%s.txt',x.name, test.name))
      #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
      write.table(res,row.names=F,sep='\t',quote=F)
      sink()
      pdf(sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
      for(k in signif.ix){
        qval <- difftest.np$qvalues[k]
        working_m <- subset(m, test.ix)
        working_m$vari <- data.transform(test.x)[test.ix,k]
        #pdf(file=sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name))
        plot1 <- ggplot(working_m, aes_string(x=plot_by[j], y= "vari")) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
          #theme(legend.position = 'bottom') + 
          labs(x="", y= "Relative Abundance (sq rt)", title =  colnames(test.x)[k]) +
          scale_x_discrete(labels=c("Pain", "No Pain")) +
          theme(plot.title = element_text(size=10)) +
          guides(fill=F, color=F) +
          scale_y_continuous(trans='sqrt') +
          scale_color_manual(values=col_list[[j]])
        plot(plot1)
      }
      dev.off()
      cat("see file\n")
    } else {
      cat("None\n")
    }
  }
}


######################################################################
#Test for distention differences
test.xs <- list(otu=x[,which(colMeans(x) > 0.0005)])
#m$q13_base_q7_symp <- as.numeric(droplevels(m$q13_base_q7_symp))

# different group comparisons
D.D <- m$q13_base_q7_symp == 1 & m$Cohort == "D" & !is.na(m$q13_base_q7_symp)
D.ND <- m$q13_base_q7_symp == 2 & m$Cohort == "D" & !is.na(m$q13_base_q7_symp)
C.D <-  m$q13_base_q7_symp == 1 & m$Cohort == "C" & !is.na(m$q13_base_q7_symp)
C.ND <- m$q13_base_q7_symp == 2 & m$Cohort == "C" & !is.na(m$q13_base_q7_symp)

test.ixs <- list('IBSDD v. IBSDND'=  D.D | D.ND,
                 'IBSCD v. IBSCND'= C.D | C.ND
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
    difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], m[test.ix,plot_by[j]], parametric=FALSE)
    if(any(difftest.np$qvalues <= ALPHA)){
      signif.ix <- which(difftest.np$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
      res <- data.frame(Taxon=colnames(test.x)[signif.ix],
                        qvalue=round(difftest.np$qvalues[signif.ix],5),
                        pvalue=round(difftest.np$pvalues[signif.ix],5),
                        round(difftest.np$classwise.means[signif.ix,,drop=F],5))
      sink(sprintf('diff_taxa/differential_abundance_%s_%s.txt',x.name, test.name))
      #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
      write.table(res,row.names=F,sep='\t',quote=F)
      sink()
      pdf(sprintf('diff_taxa/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
      for(k in signif.ix){
        qval <- difftest.np$qvalues[k]
        working_m <- subset(m, test.ix)
        working_m$vari <- data.transform(test.x)[test.ix,k]
        plot1 <- ggplot(working_m, aes_string(x=plot_by[j], y= "vari")) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color=plot_by[j])) +
          #theme(legend.position = 'bottom') + 
          labs(x="", y= "Relative Abundance (sq rt)", title =  colnames(test.x)[k]) +
          scale_x_discrete(labels=c("Distention", "No Distention")) +
          theme(plot.title = element_text(size=10)) +
          guides(fill=F, color=F) +
          scale_y_continuous(trans='sqrt') +
          scale_color_manual(values=col_list[[j]])
        plot(plot1)
      }
      dev.off()
      cat("see file\n")
    } else {
      cat("None\n")
    }
  }
}
