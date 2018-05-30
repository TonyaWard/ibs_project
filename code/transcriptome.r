# Read in transcriptome data and plot

trans_1 <- read.table(trans1_fp, sep = "\t", head=T, stringsAsFactors = F)
trans_2 <- read.table(trans2_fp, sep = "\t", head=T, stringsAsFactors = F)

trans_1[which(duplicated(trans_1$GeneName)), "GeneName"] <- paste(trans_1$GeneName[which(duplicated(trans_1$GeneName))], "_2", sep="")

rownames(trans_1) <- trans_1$GeneName
rownames(trans_2) <- trans_2$GeneName

trans_1 <- trans_1[,! colnames(trans_1) %in% c("GeneName", "GeneID", "Chromosome", "Start", "Stop", "Chr")]
trans_2 <- trans_2[,! colnames(trans_2) %in% c("GeneName", "GeneID", "Chromosome", "Start", "Stop", "Chr")]

colnames(trans_1) <- gsub("X", "", colnames(trans_1))
colnames(trans_2) <- gsub("X", "", colnames(trans_2))

#Line up these with biopsy map

t1_biopsy <- m_bio[m_bio$Timepoint == "1" & m_bio$study_id %in% colnames(trans_1),]
t2_biopsy <- m_bio[m_bio$Timepoint == "2" & m_bio$study_id %in% colnames(trans_2),]

#the are in the same order, so make the IDs the SampleID from biopsy
colnames(trans_1) <- rownames(t1_biopsy)
colnames(trans_2) <- rownames(t2_biopsy)

trans_1 <- data.frame(t(trans_1))
trans_2 <- data.frame(t(trans_2))

#See if they are correlated with subtype - timepoint 1
bio_otu <- x_bio[rownames(trans_1),]

#None are sig.
correlations <- associate(bio_otu, trans_1, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

Ds <- rownames(trans_1)[t1_biopsy$Cohort == "D"]
Cs <-rownames(trans_1)[t1_biopsy$Cohort == "C"]
Hs <- rownames(trans_1)[t1_biopsy$Cohort == "H"]

#None are sig.
correlation.table.D <- associate(bio_otu[Ds,], trans_1[Ds,], method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
correlation.table.C <- associate(bio_otu[Cs,], trans_1[Cs,], method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
correlation.table.H <- associate(bio_otu[Hs,], trans_1[Hs,], method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

#Try PCoA of transcripts and microbiome
beta_trans1 <- as.matrix(vegdist(trans_1))

beta_dist = as.dist(beta_trans1)
ad = adonis(beta_trans1 ~ t1_biopsy[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, t1_biopsy$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]

######################################################################
#PCOA of beta div
PCOA <- pcoa(beta_trans1)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}

PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
mapping2 <- t1_biopsy
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]

trans1_pcoa <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)

# Trans2
beta_trans2 <- as.matrix(vegdist(trans_2))

beta_dist = as.dist(beta_trans2)
ad = adonis(beta_trans2 ~ t2_biopsy[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, t2_biopsy$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]

######################################################################
#PCOA of beta div
PCOA <- pcoa(beta_trans2)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}

PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
mapping2 <- t2_biopsy
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]

trans2_pcoa <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)

# compare this to the bacterial mb
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
beta_div_bio1 <- beta_div_bio[rownames(trans_1), rownames(trans_1)]
PCOA_bio1 <- pcoa(beta_div_bio1)$vectors
for(i in 1:ncol(PCOA_bio1)){
  colnames(PCOA_bio1)[i] <- paste("PC",i, sep="")
}

PCOA_bio1 <- cbind(PCOA_bio1, rownames(PCOA_bio1))
colnames(PCOA_bio1)[ncol(PCOA_bio1)] <- "SampleID"
mapping2 <- m_bio[rownames(trans_1),]
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA_bio <- merge(PCOA_bio1, mapping2, by="SampleID")
PCOA_bio$PC1 <- as.numeric(levels(PCOA_bio$PC1))[PCOA_bio$PC1]
PCOA_bio$PC2 <- as.numeric(levels(PCOA_bio$PC2))[PCOA_bio$PC2]
PCOA_bio$PC3 <- as.numeric(levels(PCOA_bio$PC3))[PCOA_bio$PC3]
PCOA_bio$PC4 <- as.numeric(levels(PCOA_bio$PC4))[PCOA_bio$PC4]

PCOA_bio1 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  #geom_text(aes(x = PC1, y = PC2, label=study_id),hjust=0, vjust=0) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)


#PCOA of beta div
beta_div_bio2 <- beta_div_bio[rownames(trans_2), rownames(trans_2)]
PCOA_bio2 <- pcoa(beta_div_bio2)$vectors
for(i in 1:ncol(PCOA_bio2)){
  colnames(PCOA_bio2)[i] <- paste("PC",i, sep="")
}

PCOA_bio2 <- cbind(PCOA_bio2, rownames(PCOA_bio2))
colnames(PCOA_bio2)[ncol(PCOA_bio2)] <- "SampleID"
mapping2 <- m_bio[rownames(trans_2),]
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA_bio <- merge(PCOA_bio2, mapping2, by="SampleID")
PCOA_bio$PC1 <- as.numeric(levels(PCOA_bio$PC1))[PCOA_bio$PC1]
PCOA_bio$PC2 <- as.numeric(levels(PCOA_bio$PC2))[PCOA_bio$PC2]
PCOA_bio$PC3 <- as.numeric(levels(PCOA_bio$PC3))[PCOA_bio$PC3]
PCOA_bio$PC4 <- as.numeric(levels(PCOA_bio$PC4))[PCOA_bio$PC4]

PCOA_bio2 <- ggplot(PCOA_bio) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  #geom_text(aes(x = PC1, y = PC2, label=study_id),hjust=0, vjust=0) +
  annotate("text", x=-0.3, y=-0.2, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.3, y=-0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)


#Procrustes of trans and mb - time 1
#Generate PCoAs
beta_div_bio1 <- beta_div_bio[rownames(trans_1), rownames(trans_1)]
PCOA_bio1 <- pcoa(beta_div_bio1)$vectors
PCOA_trans1 <- pcoa(beta_trans1)$vectors

pro <- procrustes(PCOA_bio1, PCOA_trans1)
pro <- procrustes(PCOA_bio1[,1:3], PCOA_trans1[,1:3])
summary(pro)
pro2 <- protest(PCOA_bio1[,1:3], PCOA_trans1[,1:3], permutations = 999)

bio_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)

colnames(trans_pro) <- colnames(bio_pro)

trans_pro$type <- "transcript"
trans_pro <- cbind(trans_pro, t1_biopsy)
bio_pro$type <- "microbiome"
bio_pro <- cbind(bio_pro, t1_biopsy)

print(paste( "P=", pro2$signif))
print(paste( "M=", pro2$ss))
p_val <- pro2$signif
m_val <- pro2$ss

PCOA <- rbind(bio_pro, trans_pro)
PCOA$SampleID <- rownames(PCOA)
colnames(PCOA) <- gsub("Axis.", "PC", colnames(PCOA))

pro_bio_trans <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", shape="type")) +
  geom_line(aes(x= PC1, y=PC2, group=study_id, color=Cohort)) +
  scale_shape_manual(values=c(16, 15)) +
  scale_color_manual(values=cols) +
  annotate("text", x=0.2, y=0.25, label= paste("P=", p_val), size=2) +
  annotate("text", x=0.2, y=0.20, label= paste("M2=", round(m_val, digits=3)), size=2)

#Procrustes of trans and mb - time 1
#Generate PCoAs
beta_div_bio2 <- beta_div_bio[rownames(trans_2), rownames(trans_2)]
PCOA_bio2 <- pcoa(beta_div_bio2)$vectors
PCOA_trans2 <- pcoa(beta_trans2)$vectors

pro <- procrustes(PCOA_bio2, PCOA_trans2)
pro <- procrustes(PCOA_bio2[,1:3], PCOA_trans2[,1:3])
summary(pro)
pro2 <- protest(PCOA_bio2[,1:3], PCOA_trans2[,1:3], permutations = 999)

bio_pro <- data.frame(pro$X)
trans_pro <- data.frame(pro$Yrot)

colnames(trans_pro) <- colnames(bio_pro)

trans_pro$type <- "transcript"
trans_pro <- cbind(trans_pro, t2_biopsy)
bio_pro$type <- "microbiome"
bio_pro <- cbind(bio_pro, t2_biopsy)

print(paste( "P=", pro2$signif))
print(paste( "M=", pro2$ss))
p_val <- pro2$signif
m_val <- pro2$ss

PCOA <- rbind(bio_pro, trans_pro)
PCOA$SampleID <- rownames(PCOA)
colnames(PCOA) <- gsub("Axis.", "PC", colnames(PCOA))

pro_bio_trans2 <- ggplot(PCOA) +
  geom_point(size = 3, aes_string(x = "PC1", y = "PC2", color = "Cohort", shape="type")) +
  geom_line(aes(x= PC1, y=PC2, group=study_id, color=Cohort)) +
  scale_shape_manual(values=c(16, 15)) +
  scale_color_manual(values=cols) +
  annotate("text", x=0.2, y=0.25, label= paste("P=", p_val), size=2) +
  annotate("text", x=0.2, y=0.20, label= paste("M2=", round(m_val, digits=3)), size=2)
