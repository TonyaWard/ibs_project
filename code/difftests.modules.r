#Test for differentiation across groups

ALPHA <- 0.1
test.xs <- list(otu=modulesc)

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
        difftest.np$qvalues[is.na(difftest.np$qvalues)] <- 1
        if(any(difftest.np$qvalues <= ALPHA)){
            signif.ix <- which(difftest.np$qvalues <= ALPHA)
            signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
            res <- data.frame(Module=colnames(test.x)[signif.ix],
                              qvalue=round(difftest.np$qvalues[signif.ix],5),
                              pvalue=round(difftest.np$pvalues[signif.ix],5),
                              round(difftest.np$classwise.means[signif.ix,,drop=F],5))
            sink(sprintf('modules/diff_modules/differential_abundance_%s_%s.txt',x.name, test.name))
            write.table(res,row.names=F,sep='\t',quote=F)
            sink()
            pdf(sprintf('modules/diff_modules/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
            for(k in signif.ix){
              qval <- difftest.np$qvalues[k]
                working_mc <- subset(mc, test.ix)
                working_mc$vari <- data.transform(test.x)[test.ix,k]
       
                plot1 <- ggplot(working_mc, aes_string(x=plot_by[j], y= "vari")) +
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

# ######################################################################
# #Heatmap of diff taxa
# 
# ## test species and genus differences
# #test.xs <- list(species=xspc, genus=xgnc)
# test.xs <- list(otu=xc[,which(colMeans(xc) > 0.00001)])
# # different group comparisons
# test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
#                  'HC v. IBSD'=ixc.hc | ixc.ibsd,
#                  'HC v. IBSC'=ixc.hc | ixc.ibsc,
#                  'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
# )
# plot_by <- c("IBS", "Cohort", "Cohort", "Cohort")
# col_list <- list(cols_ibs, cols_dh, cols_ch, cols)
# 
# i <- 1
# x.name <- names(test.xs)[i]
# test.x <- test.xs[[i]]
#   
# for(j in 1:length(test.ixs)){
#   test.name <- names(test.ixs)[j]
#   test.ix <- test.ixs[[j]]
#   difftest <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=TRUE)
#   difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], mc$Cohort[test.ix], parametric=FALSE)
#   if(any(difftest$qvalues <= ALPHA)){
#     signif.ix <- which(difftest$qvalues <= ALPHA)
#     signif.ix <- signif.ix[order(difftest$pvalues[signif.ix])]
#     #subset to keep only taxa sig.
#     new_table <- data.transform(test.x)[test.ix,signif.ix]
#     #new_table <- new_table[,colSums(new_table) > 0.1]
#     new_map <- mc[test.ix,]
#     #order by pval
#     taxa_order <- names(signif.ix)
#     #or by taxon abundance
#     #taxa_order <- colnames(new_table[1:(ncol(new_table)-3)])[order(colSums(new_table[,1:(ncol(new_table)-3)]))] #set the taxa order
#     new_table <- new_table[,taxa_order]
#     sample_order <- rownames(new_table[order(rowSums(new_table)),])
#     new_table <- cbind(new_table, new_map[,c("Cohort", "IBS")])
#     new_table$subject_id <- rownames(new_table)
#     new_table2 <- melt(new_table, ids=c("Cohort", "IBS", "subject_id"))
#     #new_table$variable <- factor(new_table$variable, levels=healthy_order)
#     #new_table2$subject_id <- factor(new_table2$subject_id, levels=sample_order)
#     plot1 <- ggplot(new_table2, aes(x = subject_id, y = variable)) + 
#       geom_tile(aes(fill = value)) + 
#       scale_fill_gradient(name = 'relative abundance', low = 'white', high = 'purple') +
#       facet_grid(~Cohort, scales="free") +
#       theme(axis.title.y = element_blank(), axis.text.x=element_blank(),axis.title.x = element_blank())
#     pdf_name <- paste("diff_taxa/",test.name, "heatmap.pdf", sep="")
#     pdf(pdf_name,width=12,height=8)
#     plot(plot1)
#     dev.off()
#   }
# }

######################################################################
#Test for flare differences
test.xs <- list(otu=modules)

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
    difftest.np$qvalues[is.na(difftest.np$qvalues)] <- 1
    if(any(difftest.np$qvalues <= ALPHA)){
      signif.ix <- which(difftest.np$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
      res <- data.frame(Module=colnames(test.x)[signif.ix],
                        qvalue=round(difftest.np$qvalues[signif.ix],5),
                        pvalue=round(difftest.np$pvalues[signif.ix],5),
                        round(difftest.np$classwise.means[signif.ix,,drop=F],5))
      sink(sprintf('modules/diff_modules/differential_abundance_%s_%s.txt',x.name, test.name))
      #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
      write.table(res,row.names=F,sep='\t',quote=F)
      sink()
      pdf(sprintf('modules/diff_modules/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
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
test.xs <- list(otu=modules)

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
    difftest.np$qvalues[is.na(difftest.np$qvalues)] <- 1
    if(any(difftest.np$qvalues <= ALPHA)){
      signif.ix <- which(difftest.np$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
      res <- data.frame(Module=colnames(test.x)[signif.ix],
                        qvalue=round(difftest.np$qvalues[signif.ix],5),
                        pvalue=round(difftest.np$pvalues[signif.ix],5),
                        round(difftest.np$classwise.means[signif.ix,,drop=F],5))
      sink(sprintf('modules/diff_modules/differential_abundance_%s_%s.txt',x.name, test.name))
      #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
      write.table(res,row.names=F,sep='\t',quote=F)
      sink()
      pdf(sprintf('modules/diff_modules/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
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
test.xs <- list(otu=modules)

# different group comparisons
D.D <- m$q13_base_q7_symp == 1 & m$Cohort == "D" & !is.na(m$q13_base_q7_symp)
D.ND <- m$q13_base_q7_symp == 2 & m$Cohort == "D" & !is.na(m$q13_base_q7_symp)
C.D <-  m$q13_base_q7_symp == 1 & m$Cohort == "C" & !is.na(m$q13_base_q7_symp)
C.ND <- m$q13_base_q7_symp == 2 & m$Cohort == "C" & !is.na(m$q13_base_q7_symp)

test.ixs <- list('IBSDD v. IBSDND'=  D.D | D.ND,
                 'IBSCD v. IBSCND'= C.D | C.ND
)

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
    difftest.np <- differentiation.test(data.transform(test.x)[test.ix,], m$q13_base_q7_symp[test.ix], parametric=FALSE)
    difftest.np$qvalues[is.na(difftest.np$qvalues)] <- 1
    if(any(difftest.np$qvalues <= ALPHA)){
      signif.ix <- which(difftest.np$qvalues <= ALPHA)
      signif.ix <- signif.ix[order(difftest.np$pvalues[signif.ix])]
      res <- data.frame(Module=colnames(test.x)[signif.ix],
                        qvalue=round(difftest.np$qvalues[signif.ix],5),
                        pvalue=round(difftest.np$pvalues[signif.ix],5),
                        round(difftest.np$classwise.means[signif.ix,,drop=F],5))
      sink(sprintf('modules/diff_modules/differential_abundance_%s_%s.txt',x.name, test.name))
      #cat(paste('q=',qval,' taxon: ',colnames(test.x)[k],' ks.test pval=',difftest.np$pvalues[k],'\n',sep=''))
      write.table(res,row.names=F,sep='\t',quote=F)
      sink()
      pdf(sprintf('modules/diff_modules/differential_abundance_%s_%s.pdf',x.name, test.name),width=4,height=4)
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
