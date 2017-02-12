#Calculate Beta Div
beta_div <- as.matrix(vegdist(x, method = "bray"))

#PCOA of beta div
PCOA <- pcoa(beta_div)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}
  
PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
mapping2 <- m
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]
  
plot1 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.75)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3)) +
  theme(legend.title=element_blank())
plot2 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC2", y = "PC3", color = "Cohort", alpha=0.75)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3)) +
  theme(legend.title=element_blank())
plot3 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC3", y = "PC4", color = "Cohort", alpha=0.75)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3)) +
  theme(legend.title=element_blank())
PCOAs <- plot_grid(plot1, plot2, plot3, ncol=2, nrow=2)
save_plot("PCOAs.pdf", PCOAs, ncol=2)

#Compare beta div distances
#make stats file
file_name <- paste("Beta_Stats.txt", sep='')
sink(file_name)
sink()

#Set up pvals to get 
pvals<- c()
distances <- data.frame(nrow=1, ncol=2)
colnames(distances) <- c("value", "L1")
distances2 <- distances

test.ixs <- list(ix.hc, ix.ibsc, ix.ibsd)
names(test.ixs) <- c("Healthy", "IBS-C", "IBS-D")

for(n in 1:(length(test.ixs)-1)){
  for(i in (n+1):length(test.ixs)){
    test.x <- rownames(m[test.ixs[[n]],])
    test.y <- rownames(m[test.ixs[[i]],])
    full_set <- c(test.x, test.y)
    ktest <- c()
    set1_within <- c()
    set2_within <- c()
    between_sets <- c()
      
    if(length(test.x) > 2 && length(test.y) > 2){
      for(c in 1:length(colnames(beta_div))){
        for(r in (c+1):length(rownames(beta_div))){
          if(rownames(beta_div)[r] %in% test.x && colnames(beta_div)[c] %in% test.x){
            set1_within <- c(set1_within, beta_div[r,c])
          } else {
            if(rownames(beta_div)[r] %in% test.x && colnames(beta_div)[c] %in% test.y){
              between_sets <- c(between_sets, beta_div[r,c])
            } else {
              if(rownames(beta_div)[r] %in% test.y && colnames(beta_div)[c] %in% test.x){
                between_sets <- c(between_sets, beta_div[r,c])
              } else {
                if(rownames(beta_div)[r] %in% test.y && colnames(beta_div)[c] %in% test.y){
                  set2_within <- c(set2_within, beta_div[r,c])
                }
              }
            }
          }
        }
      }
    }
    sets_test <- list(set2_within, between_sets, set1_within)
    ktest <- kruskal.test(sets_test)
    names(sets_test) <- c(names(test.ixs[i]), paste(names(test.ixs[i]), names(test.ixs[n]), sep='-'), names(test.ixs[n]))
    pvals <- c(pvals,ktest$p.value)
      
    distances <- rbind(distances,melt(sets_test))
      
    #print stats to screen
    cat(sprintf('\n%s,%s:\n',names(test.ixs[n]), names(test.ixs[i])))
    print(ktest)
      
    #write stats to file
    sink(file_name, append =TRUE)
    cat(sprintf('\n%s, %s:\n',names(test.ixs[n]), names(test.ixs[i])))
    print(ktest)
    sink()
      
    #assign pdf name for plot
    name1 <- paste(names(test.ixs[n]), names(test.ixs[i]), sep="-")
    name1 <- paste(name1, ".pdf", sep='')
    pdf(name1, height=4,width=6)
      
    #make beta div box plots
    title <- sprintf('%s, %s',names(test.ixs[n]), names(test.ixs[i]))
    boxplot(sets_test,
            xlab='',ylab='Bray Curtis', main=title,
            col=cols[1:4])
    dev.off()
  }
}