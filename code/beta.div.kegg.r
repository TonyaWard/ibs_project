#Calculate Beta Div
beta_div <- as.matrix(vegdist(kegg, method = "bray"))
write.csv(beta_div, file="KEGG/beta_div/IBS_BrayCurtis.csv", quote=F, row.names=T)

######################################################################
#Adonis to check differences in centroid
beta_dist = as.dist(beta_div)
ad = adonis(beta_div ~ m[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, m$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]

######################################################################
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
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)
plot2 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC2", y = "PC3", color = "Cohort", alpha=0.65)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())
plot3 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC3", y = "PC4", color = "Cohort", alpha=0.65)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())
PCOAs <- plot_grid(plot1, plot2, plot3, ncol=2, nrow=2)
pdf("KEGG/beta_div/PCOAs.pdf", height=8, width=9)
plot(PCOAs)
dev.off()


######################################################################
#Compare beta div distances
#make stats file
file_name <- paste("KEGG/beta_div/Beta_Stats.txt", sep='')
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
    name1 <- paste("KEGG/beta_div/", name1, ".pdf", sep='')
    pdf(name1, height=4,width=6)
      
    #make beta div box plots
    title <- sprintf('%s, %s',names(test.ixs[n]), names(test.ixs[i]))
    boxplot(sets_test,
            xlab='',ylab='Bray Curtis', main=title,
            col=cols[1:4])
    dev.off()
  }
}


######################################################################
#Beta Div for each cohort separate, with flares noted
plotD <- ggplot(PCOA[PCOA$Cohort=="D",]) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.35)) + 
  scale_color_manual(values=cols[2]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="D",],  aes_string(x = "PC1", y = "PC2", color = "Cohort"),size=4, stroke=1.5, pch=2)


plotC <- ggplot(PCOA[PCOA$Cohort=="C",]) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.35)) + 
  scale_color_manual(values=cols[1]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="C",],  aes_string(x = "PC1", y = "PC2", color = "Cohort"),size=4, stroke=1.5, pch=2)

together <- plot_grid(plotD, plotC, ncol=2)
pdf("KEGG/beta_div/cohort_flares.pdf", width=10, height=4)
plot(together)
dev.off()

plotC_patient <- ggplot(PCOA[PCOA$Cohort=="C",]) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "ID_on_tube", alpha=0.35)) + 
  scale_color_manual(values=cols2(23)) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="C",],  aes_string(x = "PC1", y = "PC2", color = "ID_on_tube"),size=4, stroke=1.5, pch=2)

plotD_patient <- ggplot(PCOA[PCOA$Cohort=="D",]) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "ID_on_tube", alpha=0.35)) + 
  scale_color_manual(values=cols2(29)) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="D",],  aes_string(x = "PC1", y = "PC2", color = "ID_on_tube"),size=4, stroke=1.5, pch=2)

#or with bursts
ggC <- merge(PCOA[PCOA$Cohort=="C",],aggregate(cbind(mean.x=PC1,mean.y=PC2)~ID_on_tube,PCOA[PCOA$Cohort=="C",],mean),by="ID_on_tube")
burst_C <- ggplot(ggC, aes(x=PC1,y=PC2,color=factor(ID_on_tube)))+
  geom_point(size=3, alpha=0.6)+
  scale_color_manual(values=cols2(29)) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  geom_point(aes(x=mean.x,y=mean.y),size=3, shape=46)+
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2))+
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="C",],  alpha=0.7, aes_string(x = "PC1", y = "PC2", color = "ID_on_tube"),size=4, stroke=1.5, pch=2)

ggD <- merge(PCOA[PCOA$Cohort=="D",],aggregate(cbind(mean.x=PC1,mean.y=PC2)~ID_on_tube,PCOA[PCOA$Cohort=="D",],mean),by="ID_on_tube")
burst_D <- ggplot(ggD, aes(x=PC1,y=PC2,color=factor(ID_on_tube)))+
  geom_point(size=3, alpha=0.6)+
  scale_color_manual(values=cols2(29)) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  geom_point(aes(x=mean.x,y=mean.y),size=3, shape=46)+
  geom_segment(aes(x=mean.x, y=mean.y, xend=PC1, yend=PC2))+
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="D",],  alpha=0.7, aes_string(x = "PC1", y = "PC2", color = "ID_on_tube"),size=4, stroke=1.5, pch=2)

together <- plot_grid(plotD_patient, plotC_patient, burst_D, burst_C, ncol=2, nrow=2)
pdf("KEGG/beta_div/patients_flares.pdf", width=10, height=9)
plot(together)
dev.off()


######################################################################
#PCOA of beta div collapsed
#Calculate Beta Div
beta_div <- as.matrix(vegdist(keggc, method = "bray"))
write.csv(beta_div, file="KEGG/beta_div/IBS_BrayCurtis.collapsed.csv", quote=F, row.names=T)

######################################################################
#Adonis to check differences in centroid
beta_dist = as.dist(beta_div)
ad = adonis(beta_div ~ mc[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, mc$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]

######################################################################

PCOA <- pcoa(beta_div)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}

PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
mapping2 <- mc
mapping2$SampleID <- rownames(mapping2)
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]

plot1 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)
plot2 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC2", y = "PC3", color = "Cohort", alpha=0.65)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())
plot3 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC3", y = "PC4", color = "Cohort", alpha=0.65)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())
PCOAs <- plot_grid(plot1, plot2, plot3, ncol=2, nrow=2)
pdf("KEGG/beta_div/PCOAs_Collapsed.pdf", height=8, width=9)
plot(PCOAs)
dev.off()

