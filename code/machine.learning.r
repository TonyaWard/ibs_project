#Machine Learning

#Subset table
test.xc <- xc[,which(colMeans(xc) > 0.0005)]

# different group comparisons
test.ixs <- list('HC v. IBS'=ixc.hc | ixc.ibs,
                 'HC v. IBSD'=ixc.hc | ixc.ibsd,
                 'HC v. IBSC'=ixc.hc | ixc.ibsc,
                 'IBSC v. IBSD'=ixc.ibsc | ixc.ibsd
)
plot_by <- c("IBS", "Cohort", "Cohort", "Cohort")
col_list <- list(cols_ibs, cols_dh, cols_ch, cols)

#Run all comparisons
for(i in 1:length(test.ixs)){
  ix <- test.ixs[[i]] #set test
  res <- randomForest(test.xc[ix,],droplevels(mc[ix,plot_by[i]]),ntree=2000) #predict
  sink(paste("machine_learning/", names(test.ixs)[i], ".txt", sep=""))
  print(res) #print outputs
  print(table(mc[ix,plot_by[i]])/sum(ix))
  print(res$importance[order(res$importance[,1],decreasing=TRUE),1,drop=F][1:10,,drop=F]) #print top microbes
  sink()
  #plot top bugs
  microbes <- data.frame(res$importance[order(res$importance[,1],decreasing=TRUE),1,drop=F][1:10,,drop=F])
  colnames(microbes) <- "Importance"
  rownames(microbes) <- gsub("s__", "", rownames(microbes))
  microbes$names <- factor(rownames(microbes), levels = rownames(microbes))
  
  pdf(paste('machine_learning/RF_importance', names(test.ixs)[i], '.pdf', sep=""), width=6, height=3)
  plot1 <- ggplot(microbes, aes(x=names, y= Importance)) +
    geom_bar(stat="identity", fill="#053058") +
    coord_flip() +
    labs(x="", y="mean decrease accuracy")
  plot(plot1)
  dev.off()
}
