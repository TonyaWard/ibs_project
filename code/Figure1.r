# Figure 1 - Taxa and Kegg Beta Diversity PCoA, Alpha-Time, Beta-Time
# and within person variability

######
# Filter the data
######

# Keep only women & no flares
w_samples <- m$SampleID[as.character(m$Gender) == "F"]
w_sample_bio <- m_bio$SampleID[as.character(m_bio$Gender)=="F"]


# Keep only women & no flares
m <- m[w_samples,]
flares <- m[m$Flare=="Flare" & ! is.na(m$Flare),"SampleID"]
xF<- x[flares,]
mF <- m[flares,]
x <- x[w_samples,]
kegg <- kegg[w_samples,]
x_bio <- x_bio[w_sample_bio,]
m_bio <- m_bio[w_sample_bio,]


noF <- m[is.na(m$Flare), "SampleID"]
m_noF <- m[noF,]
kegg_noF <- kegg[noF,]
x_noF <- x[noF,]

#Calculate Alpha Div for all
alpha_div <- as.data.frame(diversity(x_noF, index="shannon"))
colnames(alpha_div) <- "shannon"
alpha_div$simpson <- diversity(x_noF, index="simpson")
alpha_div$obs_species <- rowSums(x_noF > 0)
alpha_div <- alpha_div[rownames(m_noF),]
m_noF$shannon <- alpha_div$shannon
m_noF$observed_species <- alpha_div$obs_species
m_noF$simpson <- alpha_div$simpson

#Alpha Div for collapsed
alpha_divC <- as.data.frame(diversity(xc, index="shannon"))
colnames(alpha_divC) <- "shannon"
alpha_divC$simpson <- diversity(xc, index="simpson")
alpha_divC$obs_species <- rowSums(xc > 0)
alpha_divC <- alpha_divC[rownames(mc),]
mc$shannon <- alpha_divC$shannon
mc$observed_species <- alpha_divC$obs_species
mc$simpson <- alpha_divC$simpson

#Alpha Dive for flares
alpha_divF <- as.data.frame(diversity(xF, index="shannon"))
colnames(alpha_divF) <- "shannon"
alpha_divF$simpson <- diversity(xF, index="simpson")
alpha_divF$obs_species <- rowSums(xF > 0)
alpha_divF <- alpha_divF[rownames(mF),]
mF$shannon <- alpha_divF$shannon
mF$observed_species <- alpha_divF$obs_species
mF$simpson <- alpha_divF$simpson

#Alpha Div for biopsies
alpha_divB <- as.data.frame(diversity(x_bio, index="shannon"))
colnames(alpha_divB) <- "shannon"
alpha_divB$simpson <- diversity(x_bio, index="simpson")
alpha_divB$obs_species <- rowSums(x_bio > 0)
alpha_divB <- alpha_divB[rownames(m_bio),]
m_bio$shannon <- alpha_divB$shannon
m_bio$observed_species <- alpha_divB$obs_species
m_bio$simpson <- alpha_divB$simpson

##################################################


######
# PCoA of taxa
######

# Normalize the OTU table
#CLR_otutable <- x
#CLR_otutable[x == 0] <- 0.65 #Convert any 0 to 0.65 to allow for CLR transform

#CLR_otutable <- cenLR(CLR_otutable)$x.clr  #transform


#Beta div calc.
#beta_div <- as.matrix(vegdist(CLR_otutable, method = "euclidean"))
#beta_div <- as.matrix(vegdist(x_noF, method = "bray"))
beta_div <- as.matrix(vegdist(x_noF[,which(colMeans(x_noF) > 0.0005)], method = "bray"))

#Adonis to check differences in centroid
beta_dist = as.dist(beta_div)
ad = adonis(beta_div ~ m_noF[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]

#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, m_noF$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]
PCOA <- pcoa(beta_div)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}

PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
mapping2 <- m_noF
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]

range(PCOA$PC1)
range(PCOA$PC2)

centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$Cohort,PCOA,mean)
colnames(centroids) <- c("Cohort", "PC1", "PC2")

#pairwise adonis for differences:
pair_ad_list <- list()
for(t in 1:2){
  for(u in (t+1):3){
    t1 <- unique(m$Cohort)[[t]]
    t2 <- unique(m$Cohort)[[u]]
    sams <- rownames(m_noF[m_noF$Cohort == t1 | m_noF$Cohort == t2,])
    beta_dist = beta_div[sams,sams]
    ad = adonis(beta_dist ~ m[sams,"Cohort"], data=m, permutations=999)
    p_val2 <- ad$aov.tab[1,6]
    r_sq2 <- ad$aov.tab[1,5]
    pair_ad_list[[paste(t1, "_", t2, sep="")]] <- c("pval", p_val2, "r_sq", r_sq2)
  }
}

OTU_PCoA <- ggplot(PCOA) +
  geom_point(size = 2.5, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  stat_chull(aes(x=PC1,y=PC2,color=Cohort, fill=Cohort),alpha = 0.05, geom = "polygon") +
  scale_fill_manual(values=cols[1:3]) +
  annotate("text", x=-0.3, y=0.35, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.4, y=0.45, label= paste("R2=", round(r_sq, digits=3)), size=2.5)
  #geom_point(data=centroids, aes(x=PC1, y=PC2, fill=Cohort), 
  #           shape = 23, size=3, alpha=0.8, show.legend = F)


######
# Use Collapsed Data
#####
mc$SampleID_OG <- mc$SampleID
mc$SampleID <- rownames(mc)
w_collapsed <- mc[mc$Gender=="F", "SampleID"]

mc <- mc[w_collapsed,]
xc <- xc[w_collapsed,]
keggc <- keggc[w_collapsed,]

# H vs D Clustering
#Beta div calc.
beta_div_c <- as.matrix(vegdist(xc[,which(colMeans(xc) > 0.0005)], method = "bray"))

#Adonis to check differences in centroid
beta_dist = as.dist(beta_div_c)
ad = adonis(beta_div_c ~ mc[,"Cohort"], data=m, permutations=999)
p_val <- ad$aov.tab[1,6]
r_sq <- ad$aov.tab[1,5]

#Run Stats for diff. dispersion
beta_out <- betadisper(beta_dist, mc$Cohort)
p_val_disp <- permutest(beta_out)$tab[1, 6]
PCOA <- pcoa(beta_div_c)$vectors
for(i in 1:ncol(PCOA)){
  colnames(PCOA)[i] <- paste("PC",i, sep="")
}

PCOA <- cbind(PCOA, rownames(PCOA))
colnames(PCOA)[ncol(PCOA)] <- "SampleID"
mapping2 <- mc
mapping2 <- data.frame(lapply(mapping2, as.character), stringsAsFactors=FALSE)
PCOA <- merge(PCOA, mapping2, by="SampleID")
PCOA$PC1 <- as.numeric(levels(PCOA$PC1))[PCOA$PC1]
PCOA$PC2 <- as.numeric(levels(PCOA$PC2))[PCOA$PC2]
PCOA$PC3 <- as.numeric(levels(PCOA$PC3))[PCOA$PC3]
PCOA$PC4 <- as.numeric(levels(PCOA$PC4))[PCOA$PC4]

range(PCOA$PC1)
range(PCOA$PC2)

centroids <- aggregate(cbind(PCOA$PC1,PCOA$PC2) ~ PCOA$Cohort,PCOA,mean)
colnames(centroids) <- c("Cohort", "PC1", "PC2")

#pairwise adonis for differences:
pair_ad_list <- list()
for(t in 1:2){
  for(u in (t+1):3){
    t1 <- unique(m$Cohort)[[t]]
    t2 <- unique(m$Cohort)[[u]]
    sams <- rownames(mc[mc$Cohort == t1 | mc$Cohort == t2,])
    beta_dist = beta_div_c[sams,sams]
    ad = adonis(beta_dist ~ mc[sams,"Cohort"], data=mc, permutations=999)
    p_val2 <- ad$aov.tab[1,6]
    r_sq2 <- ad$aov.tab[1,5]
    pair_ad_list[[paste(t1, "_", t2, sep="")]] <- c("pval", p_val2, "r_sq", r_sq2)
  }
}

OTU_C_PCoA <- ggplot(PCOA) +
  geom_point(size = 2.5, aes_string(x = "PC1", y = "PC2", color = "Cohort", alpha=0.25)) + 
  scale_color_manual(values=cols[1:3]) +
  guides(color=guide_legend(nrow=3), alpha=F) +
  theme(legend.title=element_blank()) +
  stat_chull(aes(x=PC1,y=PC2,color=Cohort, fill=Cohort),alpha = 0.05, geom = "polygon") +
  scale_fill_manual(values=cols[1:3]) +
  annotate("text", x=-0.23, y=0.35, label= paste("P=", p_val), size=2.5) +
  annotate("text", x=-0.23, y=0.25, label= paste("R2=", round(r_sq, digits=3)), size=2.5)
#geom_point(data=centroids, aes(x=PC1, y=PC2, fill=Cohort), 
#           shape = 23, size=3, alpha=0.8, show.legend = F)

#####
# Plot between vs witin 
####

# get all within-cohort distances
distance_table <- data.frame("Subject"=NA, "Variability"=NA, "Comparison"=NA)
mc$Cohort <- as.character(mc$Cohort)
people <- sort(unique(mc$SampleID))
for(i in 1:length(people)){
  
  cohort <- mc[mc$SampleID == people[i],"Cohort"]
  withins <- mc[!mc$SampleID == people[i] & mc$Cohort == cohort,"SampleID"]
  within_name <- paste("within_", cohort, sep ="")
  dist_within <- mean(beta_div_c[withins, people[i]])
  distance_table <- rbind(distance_table, c(people[i], dist_within, within_name))
  
  non <- unique(mc$Cohort)[!unique(mc$Cohort) == cohort]
  for(j in 1:length(non)){
    name <- paste(cohort, "_", non[j], sep ="")
    betweens <-  mc[!mc$SampleID == people[i] & mc$Cohort == non[j],"SampleID"]
    dist_btwn <- mean(beta_div_c[betweens,people[i]])
    distance_table <- rbind(distance_table, c(people[i], dist_btwn, name))
  }
  
  
}
distance_table <- distance_table[2:nrow(distance_table),]
distance_table$Variability <- as.numeric(as.character(distance_table$Variability))
distance_table[distance_table$Comparison == "D_C", "Comparison"] <- "C_D" 
distance_table[distance_table$Comparison == "D_H", "Comparison"] <- "H_D" 
distance_table[distance_table$Comparison == "C_H", "Comparison"] <- "H_C" 

keep_dists <- c("within_C", "H_C", "within_H", "H_D", "within_D")
distance_table <- distance_table[distance_table$Comparison %in% keep_dists,]
distance_table$Comparison <- factor(distance_table$Comparison, 
                                    levels = c("within_C", "H_C", "within_H",
                                               "H_D", "within_D"))
# Pairwise tests
pval_table <- data.frame("Test_Name"=NA, "Pval"=NA)
for(i in 1:(length(unique(distance_table$Comparison))-1)){
  for(j in (i+1):length(unique(distance_table$Comparison))){
    name1 <- unique(as.character(distance_table$Comparison))[i]
    name2 <- unique(as.character(distance_table$Comparison))[j]
    test_name <- paste(name1,"_vs_",name2, sep="")
    test.p <- t.test(distance_table[as.character(distance_table$Comparison)==name1,"Variability"],
                    distance_table[as.character(distance_table$Comparison)==name2,"Variability"])$p.val
    pval_table <- rbind(pval_table, c(test_name, test.p))
    }
}

Distance_Collapsed <- ggplot(distance_table, aes(x=Comparison, y=Variability)) +
  geom_jitter(width=0.09, aes(color=Comparison), alpha=0.65, size=2.5) +
  geom_boxplot(fill=NA, outlier.color = NA) +
  scale_color_manual(values=c(cols[1], cols[4], cols[3], cols[6], cols[2])) +
  labs(y = "Mean Distance") +
  guides(color=F)

##### 
# Plot Diversity by Time
#####
# get all within-patient distances
vari_table <- data.frame("Subject"=NA, "Timepoint"=NA, "Cohort"=NA, "Distance"=NA)
m_noF$ID_on_tube <- as.numeric(m_noF$ID_on_tube)
patient.nos <- sort(unique(m_noF$ID_on_tube))
for(i in 1:length(unique(m_noF$ID_on_tube))){
  patient.no <- patient.nos[i]
  sub_map <- m_noF[m_noF$ID_on_tube == patient.no,]
  sub_map2 <- sub_map[with(sub_map,order(sub_map$Timepoint)),]
  
  cohort <- sub_map$Cohort[1] 
  
  # calculates only distance to previous timepoint
  if(nrow(sub_map2)> 1){
    for(j in 1:(nrow(sub_map2)-1)){
      timepoint <- sub_map2[j,"Timepoint"]
      distance <- beta_div[sub_map2[j,"SampleID"], sub_map2[j+1, "SampleID"]]
      vari_table <- rbind(vari_table, c(patient.no, timepoint, cohort,distance))
    }
  }
}
vari_table <- vari_table[2:nrow(vari_table),]
vari_table$Distance <- as.numeric(as.character(vari_table$Distance))

Variability <- ggplot(vari_table, aes(x=Timepoint, y=Distance, group=Cohort, color=Cohort)) +
  geom_jitter(width=0.2, alpha=0.3, size=1.5, color="grey") +
  geom_smooth(alpha=0.25) +
  scale_color_manual(values=cols[1:3]) +
  labs(y = "Taxonomic Variability (Bray Curtis Dissimilarity)")



Figure1 <- plot_grid(OTU_PCoA, Variability,OTU_C_PCoA, Distance_Collapsed)
