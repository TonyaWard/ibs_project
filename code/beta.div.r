#Calculate Beta Div
beta_div <- as.matrix(vegdist(x[,which(colMeans(x) > 0.0005)], method = "bray"))
write.csv(beta_div, file="beta_div/IBS_BrayCurtis.csv", quote=F, row.names=T)
# eps = 0.5
# otu_table = x.raw*(1 - rowSums(x.raw==0)*eps/rowSums(x.raw))
# otu_table[otu_table==0]=eps
# otu_table = sweep(otu_table,1,rowSums(otu_table),'/')
# ls = log(otu_table)
# otu_table = ls - rowMeans(ls)
# otu_table = otu_table[!is.nan(rowSums(otu_table)),]
# 
# beta_div <- as.matrix(vegdist(otu_table, method="euclidean"))
# 
# adaptive <- adaptivegpca(t(otu_table), beta_div, k=2)
# 
# this <- data.frame(adaptive$V)
# this$Cohort <- m$Cohort
# ggplot(this) +
#   geom_point(aes(x = Axis1, y = Axis2, color=Cohort)) +
#   xlab("Axis 1") + ylab("Axis 2")

beta_div_bio <- as.matrix(vegdist(x_bio[,which(colMeans(x_bio) > 0.0005)], method = "bray"))
write.csv(beta_div_bio, file="beta_div/IBS_BrayCurtis_Biopsy.csv", quote=F, row.names=T)

beta_div_all <- as.matrix(vegdist(x_all[,which(colMeans(x_all) > 0.0005)], method = "bray"))
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
  #geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare),], size=4, aes(x=PC1, y=PC2, color=Cohort), pch=2)
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
pdf("beta_div/PCOAs.pdf", height=8, width=9)
plot(PCOAs)
dev.off()


######################################################################
#Compare beta div distances
#make stats file
file_name <- paste("beta_div/Beta_Stats.txt", sep='')
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
    name1 <- paste("beta_div/", name1, ".pdf", sep='')
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
pdf("beta_div/cohort_flares.pdf", width=10, height=4)
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
pdf("beta_div/patients_flares.pdf", width=10, height=9)
plot(together)
dev.off()

#Or elipse
#or with bursts
el_C <- ggplot(ggC, aes(x=PC1,y=PC2,color=factor(ID_on_tube)))+
  geom_point(size=3, alpha=0.6)+
  scale_color_manual(values=cols2(29)) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  geom_polygon(alpha=0.15, fill=NA, linetype="dashed") +
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="C",],  alpha=0.7, aes_string(x = "PC1", y = "PC2", color = "ID_on_tube"),size=4, stroke=1.5, pch=2)

el_D <- ggplot(ggD, aes(x=PC1,y=PC2,color=factor(ID_on_tube)))+
  geom_point(size=3, alpha=0.6)+
  scale_color_manual(values=cols2(29)) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  geom_polygon(aes(alpha=0.05), fill=NA, linetype="dashed") +
  geom_point(data= PCOA[PCOA$Flare == "Flare" & !is.na(PCOA$Flare) & PCOA$Cohort =="D",],  alpha=0.7, aes_string(x = "PC1", y = "PC2", color = "ID_on_tube"),size=4, stroke=1.5, pch=2)

together2 <- plot_grid(el_D, el_C, ncol=2)
pdf("beta_div/patients_flares_grouped.pdf", width=10, height=5)
plot(together2)
dev.off()


######################################################################
#PCOA of beta div collapsed
#Calculate Beta Div
beta_div <- as.matrix(vegdist(xc[,which(colMeans(xc) > 0.0005)], method = "bray"))
write.csv(beta_div, file="beta_div/IBS_BrayCurtis.collapsed.csv", quote=F, row.names=T)

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
pdf("beta_div/PCOAs_Collapsed.pdf", height=8, width=9)
plot(PCOAs)
dev.off()


ggplot(PCOA, aes_string(x = "PC1", y = "PC2")) +
  geom_point(size = 3, aes_string(color = "study_id", alpha=0.65)) + 
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) + 
  geom_text(aes(label=study_id),hjust=0, vjust=0)

######################################################################
##Check distance of Flare to collapsed space:
Flares <- rownames(m[!is.na(m$Flare),])
xc_plusF <- rbind(xc, x[Flares,,drop=T])
map_cols <- intersect(colnames(mc), colnames(m))
mc_plusF <- rbind(mc[,map_cols], m[Flares,map_cols])
beta_div <- as.matrix(vegdist(xc_plusF[,which(colMeans(xc_plusF) > 0.0005)], method = "bray"))
write.csv(beta_div, file="beta_div/IBS_BrayCurtis.collapsed_plusF.csv", quote=F, row.names=T)

paired <- mc_plusF$ID_on_tube[duplicated(mc_plusF$ID)]
paired_sams <- rownames(mc_plusF[mc_plusF$ID_on_tube %in% paired,])
paired_div <- beta_div[paired_sams, paired_sams]
paired_m <- mc_plusF[rownames(paired_div),]
paired_m$Cohort <- as.character(paired_m$Cohort)
within_C <- c()
within_D <- c()
for(i in 1:length(Flares)){
  Flare <- Flares[i]
  ID <- paired_m[Flare, "ID_on_tube"]
  match <- rownames(paired_m[is.na(paired_m$Flare) & paired_m$ID_on_tube == ID,])
  if(paired_m[Flares[i],"Cohort"]=="D"){
    within_D <- c(within_D, paired_div[Flare, match])
  } else {
    within_C <- c(within_C, paired_div[Flare,match])
  }
}

paired_table <- data.frame(matrix(NA, ncol=2, nrow=6))
paired_table$D <- within_D
paired_table$C <- within_C
paired_table <- paired_table[,3:4]
paired_table <- melt(paired_table)
pdf("beta_div/paired_flare_dist.pdf", height=3, width=3)
ggplot(paired_table, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape = NA, aes(color=variable)) +
  geom_jitter(position=position_jitter(0.1), size=3, aes(color=variable)) +
  labs(x="", y= "Flare Distance") +
  guides(fill=F, color=F) +
  scale_color_manual(values=cols[2:1])
dev.off()
test_out <- t.test(paired_table$value ~ paired_table$variable)
sink("beta_div/paired_flare_collapse.txt")
test_out
sink()

#make PCoA
PCOA_bio <- data.frame(pcoa(paired_div)$vectors)
for(i in 1:ncol(PCOA_bio)){
  colnames(PCOA_bio)[i] <- paste("PC",i, sep="")
}

PCOA_bio$SampleID <- rownames(PCOA_bio)

mapping2 <- paired_m
mapping2$SampleID <- rownames(mapping2)
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA_bio <- merge(PCOA_bio, mapping2, by="SampleID")
#PCOA_bio$PC1 <- as.numeric(levels(PCOA_bio$PC1))[PCOA_bio$PC1]
#PCOA_bio$PC2 <- as.numeric(levels(PCOA_bio$PC2))[PCOA_bio$PC2]
#PCOA_bio$PC3 <- as.numeric(levels(PCOA_bio$PC3))[PCOA_bio$PC3]
#PCOA_bio$PC4 <- as.numeric(levels(PCOA_bio$PC4))[PCOA_bio$PC4]
PCOA_bio$Flare[is.na(PCOA_bio$Flare)] <- "Normal"
plot1 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", shape="Flare"), alpha=0.85) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  geom_line(aes_string(x = "PC1", y = "PC2", group = "study_id"), color="grey")
pdf("beta_div/flare_paired_PCOA.pdf", heigh=3, width=4.5)
plot1
dev.off()
###BIOPSY
######################################################################
#Adonis to check differences in centroid
beta_dist_bio = as.dist(beta_div_bio)
ad = adonis(beta_div_bio ~ m_bio[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out_bio <- betadisper(beta_dist_bio, m_bio$Cohort)
p_val_disp <- permutest(beta_out_bio)$tab[1, 6]

######################################################################
#PCOA of beta div
PCOA_bio <- pcoa(beta_div_bio)$vectors
for(i in 1:ncol(PCOA_bio)){
  colnames(PCOA_bio)[i] <- paste("PC",i, sep="")
}

PCOA_bio <- cbind(PCOA_bio, rownames(PCOA_bio))
colnames(PCOA_bio)[ncol(PCOA_bio)] <- "SampleID"
mapping2 <- m_bio
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA_bio <- merge(PCOA_bio, mapping2, by="SampleID")
PCOA_bio$PC1 <- as.numeric(levels(PCOA_bio$PC1))[PCOA_bio$PC1]
PCOA_bio$PC2 <- as.numeric(levels(PCOA_bio$PC2))[PCOA_bio$PC2]
PCOA_bio$PC3 <- as.numeric(levels(PCOA_bio$PC3))[PCOA_bio$PC3]
PCOA_bio$PC4 <- as.numeric(levels(PCOA_bio$PC4))[PCOA_bio$PC4]

plot1 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  #geom_text(aes(x = PC1, y = PC2, label=study_id),hjust=0, vjust=0) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)
plot2 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC2", y = "PC3", color = "Cohort", alpha=0.65)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())
plot3 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC3", y = "PC4", color = "Cohort", alpha=0.65)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())
PCOAs_bio <- plot_grid(plot1, plot2, plot3, ncol=2, nrow=2)
pdf("beta_div/PCOAs_bio.pdf", height=8, width=9)
plot(PCOAs_bio)
dev.off()

#By patient
plot_p <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort"), alpha=0.65) + 
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank()) +
  scale_color_manual(values=cols[1:3]) +
  geom_line(aes_string(x = "PC1", y = "PC2", group = "study_id"), alpha=0.45, color= "darkgrey")
pdf("beta_div/PCOAs_bio_patient.pdf", height=3.5, width=4)
plot(plot_p)
dev.off()

#Calculate average biopsy distance:
file_name <- paste("beta_div/Beta_Stats_Bio.txt", sep='')
sink(file_name)
sink()

bio_list <- unique(m_bio$study_id)

bio_H <- c()
bio_C <- c()
bio_D <- c()

for(z in 1:length(bio_list)){
  samples <- rownames(m_bio[m_bio$study_id == bio_list[z],])
  type <- as.character(m_bio[samples[1],"Cohort"])
  if(length(samples) > 1){
    if(type == "H"){
      bio_H <- c(bio_H, beta_div_bio[samples[1],samples[2]])
    } else {
      if(type == "C"){
        bio_C <- c(bio_C, beta_div_bio[samples[1],samples[2]])
      } else {
        bio_D <-  c(bio_D, beta_div_bio[samples[1],samples[2]])
      }
    }
  } else {
    print("only one sample for:")
    cat(samples)
    print(type)
  }
}

names(bio_C) <- rep("C", length(bio_C))
names(bio_D) <- rep("D", length(bio_D))
names(bio_H) <- rep("H", length(bio_H))

bio_C <- data.frame(bio_C)
bio_C$Cohort <- rep("C", nrow(bio_C))
bio_H <- data.frame(bio_H)
bio_H$Cohort <- rep("H", nrow(bio_H))
bio_D <- data.frame(bio_D)
bio_D$Cohort <- rep("D", nrow(bio_D))
colnames(bio_C)[1] <- "dist"
colnames(bio_D)[1] <- "dist"
colnames(bio_H)[1] <- "dist"
together_table <- rbind(bio_C, bio_D)
together_table <- rbind(together_table, bio_H)

test_out <- aov(together_table[,"dist"] ~ factor(together_table[,"Cohort"]))
anova_out <- summary(test_out)
tukey_out <- TukeyHSD(test_out)
sink(file_name)
cat("anova results for biopsies:")
anova_out
tukey_out
sink()
hc_p <- round(tukey_out$`factor(together_table[, "Cohort"])`[2,4], digits = 3)

pdf("beta_div/beta_dist_bio.pdf", height=3.5, width=4)
ggplot(together_table, aes(x=Cohort, y=dist)) +
  geom_boxplot(outlier.shape = NA, aes(color=Cohort)) +
  geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes(color=Cohort)) +
  guides(fill=F, color=F) +
  labs(y= "Bray Curtis Distance") +
  scale_color_manual(values=cols[1:3]) + 
  geom_text(aes(label= paste("H-C pval:", hc_p), x="H", y=0.75))
dev.off()


#Try beta for all samples
######################################################################
#Adonis to check differences in centroid
beta_dist = as.dist(beta_div_all)
ad = adonis(beta_div_all ~ factor(m_all[,"Cohort"]), data=m_all, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, m_all$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]

######################################################################
#PCOA of beta div
PCOA_all <- pcoa(beta_div_all)$vectors
for(i in 1:ncol(PCOA_all)){
  colnames(PCOA_all)[i] <- paste("PC",i, sep="")
}

PCOA_all <- cbind(PCOA_all, rownames(PCOA_all))
colnames(PCOA_all)[ncol(PCOA_all)] <- "SampleID"
mapping2 <- m_all
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA_all <- merge(PCOA_all, mapping2, by="SampleID")
PCOA_all$PC1 <- as.numeric(levels(PCOA_all$PC1))[PCOA_all$PC1]
PCOA_all$PC2 <- as.numeric(levels(PCOA_all$PC2))[PCOA_all$PC2]
PCOA_all$PC3 <- as.numeric(levels(PCOA_all$PC3))[PCOA_all$PC3]
PCOA_all$PC4 <- as.numeric(levels(PCOA_all$PC4))[PCOA_all$PC4]

plot1 <- ggplot(PCOA_all) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "sType", alpha=0.25)) + 
  scale_color_manual(values=cols[c(5,4,7)]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())

plot2 <- ggplot(PCOA_all) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank())

plot3 <- ggplot(PCOA_all) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "study_id", alpha=0.25)) + 
  scale_color_manual(values=cols2(length(unique(m_all$study_id)))) +
  guides(color=F, alpha=F) +
  theme(legend.title=element_blank())

PCOAs <- plot_grid(plot1, plot2, plot3, ncol=2, nrow=2)
pdf("beta_div/PCOAs_allTypes.pdf", height=8, width=9)
plot(PCOAs)
dev.off()

######################################################################
#try combining collapsed with biopsy
just_stool <- x_all[rownames(m_all[m_all$sType == "Stool",]),]
stool_m <- m_all[rownames(just_stool),]
just_stool <- apply(just_stool,2,function(xx) sapply(split(xx,stool_m$study_id),mean))
rownames(just_stool) <- sprintf('Subject_%03d',sapply(split(stool_m$study_id,stool_m$study_id),'[',1))
stool_m <- stool_m[sapply(split(1:nrow(stool_m),stool_m$study_id),'[',1),,drop=TRUE]
rownames(stool_m) <- sprintf('Subject_%03d',sapply(split(stool_m$study_id,stool_m$study_id),'[',1))

just_stool <- rbind(just_stool, x_all[rownames(m_all[m_all$sType == "Biopsy",]),])
stool_m <- rbind(stool_m, m_all[rownames(m_all[m_all$sType == "Biopsy",]),])

beta_div <- as.matrix(vegdist(just_stool[,which(colMeans(just_stool) > 0.0005)], method = "bray"))

paired <- stool_m$study_id[duplicated(stool_m$study_id)]
paired_sams <- rownames(stool_m[stool_m$study_id %in% paired,])
paired_div <- beta_div[paired_sams, paired_sams]
paired_m <- stool_m[rownames(paired_div),]
paired_m$Cohort <- as.character(paired_m$Cohort)

within_C <- c()
within_D <- c()
within_H <- c()

Bios <- unique(stool_m[stool_m$sType == "Biopsy","study_id"])
for(i in 1:length(Bios)){
  if(nrow(stool_m[stool_m$study_id == Bios[i],]) == 3){
    bio1 <- rownames(stool_m[stool_m$study_id == Bios[i] & stool_m$Timepoint == 1,])
    bio2 <- rownames(stool_m[stool_m$study_id == Bios[i] & stool_m$Timepoint == 2,])
    match <- rownames(stool_m[stool_m$sType == "Stool" & stool_m$study_id == Bios[i],])
    if(stool_m[bio1,"Cohort"]=="D"){
      within_D <- c(within_D, paired_div[bio1, match])
      within_D <- c(within_D, paired_div[bio2, match])
    } else {
      if(stool_m[bio1,"Cohort"]=="C"){
        within_C <- c(within_C, paired_div[bio1,match])
        within_C <- c(within_C, paired_div[bio2,match])
      }
      if(stool_m[bio1,"Cohort"]=="H"){
        within_H <- c(within_H, paired_div[bio1,match])
        within_H <- c(within_H, paired_div[bio2,match])
      }
    }
  }
}

paired_table <- data.frame(matrix(NA, ncol=3, nrow=22))
within_C <- c(within_C, rep(NA, 22-length(within_C)))
within_H <- c(within_H, rep(NA, 22-length(within_H)))
paired_table$D <- within_D
paired_table$C <- within_C
paired_table$H <- within_H
paired_table <- paired_table[,4:6]
paired_table <- melt(paired_table)
pdf("beta_div/paired_biopsy_dist.pdf", height=3, width=3)
ggplot(paired_table, aes(x=variable, y=value)) +
  geom_boxplot(outlier.shape = NA, aes(color=variable)) +
  geom_jitter(position=position_jitter(0.1), size=3, aes(color=variable), alpha=0.8) +
  labs(x="", y= "Biopsy Distance") +
  guides(fill=F, color=F) +
  scale_color_manual(values=cols[c(2,1,3)])
dev.off()
test_out <- aov(paired_table$value ~ paired_table$variable)
sink("beta_div/paired_biopsy_collapse.txt")
summary(test_out)
TukeyHSD(test_out)
sink()

#make PCoA
PCOA_bio <- data.frame(pcoa(paired_div)$vectors)
for(i in 1:ncol(PCOA_bio)){
  colnames(PCOA_bio)[i] <- paste("PC",i, sep="")
}

PCOA_bio$SampleID <- rownames(PCOA_bio)

mapping2 <- stool_m
mapping2$SampleID <- rownames(mapping2)
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA_bio <- merge(PCOA_bio, mapping2, by="SampleID")
#PCOA_bio$PC1 <- as.numeric(levels(PCOA_bio$PC1))[PCOA_bio$PC1]
#PCOA_bio$PC2 <- as.numeric(levels(PCOA_bio$PC2))[PCOA_bio$PC2]
#PCOA_bio$PC3 <- as.numeric(levels(PCOA_bio$PC3))[PCOA_bio$PC3]
#PCOA_bio$PC4 <- as.numeric(levels(PCOA_bio$PC4))[PCOA_bio$PC4]

plot1 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", shape="sType"), alpha=0.75) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  geom_line(aes_string(x = "PC1", y = "PC2", group = "study_id"), color="grey", linetype=2)
plot2 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", shape="sType"), alpha=0.75) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
pdf("beta_div/biopsy_paired_PCOA.pdf", height=3, width=8)
plot_grid(plot1, plot2, ncol=2)
dev.off()
