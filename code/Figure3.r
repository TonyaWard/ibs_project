# Figure 3 - Biopsy

####
# PCoA
####
beta_div_bio <- as.matrix(vegdist(x_bio[,which(colMeans(x_bio) > 0.0005)], method = "bray"))
hist(colMeans(beta_div_bio))
#beta_div_bio <- beta_div_bio[colMeans(beta_div_bio)<0.85,colMeans(beta_div_bio)<0.85]


#Adonis to check differences in centroid
beta_dist_bio = as.dist(beta_div_bio)
ad = adonis(beta_div_bio ~ m_bio[colnames(beta_div_bio),"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]
#Run Stats for diff. dispersion
beta_out_bio <- betadisper(beta_dist_bio, m_bio[colnames(beta_div_bio),"Cohort"])
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


OTU_PCoA_bio <- ggplot(PCOA_bio) +
  geom_point(size = 2.5, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.35)) + 
  geom_line(aes_string(x="PC1", y="PC2", group="ID.on.tube"), color="grey", alpha=0.35) +
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  scale_fill_manual(values=cols[1:3]) +
  annotate("text", x=-0.3, y=0.35, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.4, y=0.45, label= paste("R2=", round(r_sq, digits=3)), size=2.5)

####
# Test for distances between 1/2
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

biopsy_dist_plot <- ggplot(together_table, aes(x=Cohort, y=dist)) +
  geom_jitter(position=position_jitter(0.1), size=2.5, alpha=0.65, aes(color=Cohort)) +
  geom_boxplot(outlier.shape = NA, fill=NA) +
  guides(fill=F, color=F) +
  labs(y= "Bray Curtis Distance") +
  scale_color_manual(values=cols[1:3]) + 
  geom_text(aes(label= paste("H-C pval:", hc_p), x="H", y=0.75))

right_side <- plot_grid(biopsy_dist_plot, NULL, ncol=1, rel_heights = c(2,1))
Figure3 <- plot_grid(OTU_PCoA_bio, right_side, rel_widths = c(2,1))
