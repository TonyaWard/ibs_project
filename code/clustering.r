ALPHA <- 0.05

# test species and genus differences
#test.xs <- list(species=xspc, genus=xgnc)
test.xs <- list(OTU=xc)
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

        # ANOSIM test for clustering
        aa <- anosim(as.dist(as.matrix(vegdist(test.x))[test.ix,test.ix]), droplevels(mc$Cohort[test.ix]))
        if(aa$signif <= ALPHA){
          sink("beta_div/cluster.txt", append=T)
          cat(sprintf('Significant ANOSIM clustering for %s, %s: ',x.name, test.name))
          cat(paste('R=',aa$statistic,', p=', aa$signif,'\n',sep=''))
          sink()
        }

        # PERMANOVA test for clustering
        aa <- adonis(as.dist(as.matrix(vegdist(test.x))[test.ix,test.ix]) ~ droplevels(mc$Cohort[test.ix]))
        pval <- aa$aov.tab[1,'Pr(>F)']
        if(pval <= ALPHA){
          sink("beta_div/cluster.txt", append=T)
          cat(sprintf('Significant ADONIS clustering for %s, %s: ',x.name, test.name))
          cat(paste('Partial R2=',round(aa$aov.tab[1,'R2'],3),', p=', pval,'\n',sep=''))
          sink()
          pc <- cmdscale(vegdist(test.x))
          

          pdf(sprintf('beta_div/clustering_%s_%s.pdf',x.name, gsub(' ','_',test.name)),width=4,height=4)
          
          #getting the convex hull of each unique point set
          working_table <- cbind(mc, pc)
          working_table <- subset(working_table, test.ix)
          colnames(working_table)[(ncol(working_table)-1):ncol(working_table)] <- c("PC1", "PC2")
          plot1 <- ggplot(working_table, aes_string(y ="PC1", x="PC2", colour= plot_by[j], fill = plot_by[j])) +
            geom_point(size = 3) + 
            stat_chull(alpha = 0.1, geom = "polygon") +
            scale_fill_manual(values=col_list[[j]]) +
            scale_color_manual(values=col_list[[j]])
          print(plot1)
          dev.off()
        }
    }
}