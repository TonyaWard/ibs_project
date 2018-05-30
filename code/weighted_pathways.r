# Weighted analysis for Butyrate

BUT <- read.table("../data/IBS_Mayo_secondrun/Weighted_KEGG", 
                  sep="\t", 
                  header=T)

# Gaussian function where:
# Sigma^2 = 0.2
# Mu = 0
BUT$weight <- (1/0.4472136*(sqrt(2*pi)))*exp(1)^(-0.5*((BUT$Level - 0) / 0.4472136)^2)

#Take the KEGG table and multiply by weights
but_keeps <- intersect(BUT$KEGG_ID, colnames(kegg)) #only 11 of 26 are in here
but_kegg <- kegg[, but_keeps]
but_keggC <- keggc[,but_keeps]
rownames(BUT) <- BUT$KEGG_ID
BUT <- BUT[but_keeps,]

but_kegg_weight <- sweep(but_kegg,MARGIN=2,BUT$weight,`*`)
but_kegg_weight <- cbind(but_kegg_weight, Total = rowSums(but_kegg_weight))

but_kegg_weightC <- sweep(but_keggC,MARGIN=2,BUT$weight,`*`)
but_kegg_weightC <- cbind(but_kegg_weightC, Total = rowSums(but_kegg_weightC))
but_kegg_weightC <- data.frame(but_kegg_weightC)

# Test for Collapsed data
ktest <- kruskal.test(but_kegg_weightC[,'Total'] ~ mc[rownames(but_kegg_weightC),"Cohort"])
but_kegg_weightC <-  cbind(but_kegg_weightC, Cohort = mc[rownames(but_kegg_weightC),"Cohort"])

Bplot <- ggplot(but_kegg_weightC) + 
  geom_boxplot(aes(x=Cohort, y=Total), outlier.colour = NA) +
  geom_jitter(aes(x=Cohort, y=Total, color=Cohort), width=0.08, show.legend = F) +
  scale_color_manual(values=cols)
pdf("Weighted_Butyrate.pdf", width=4, height=2.5)
Bplot
dev.off()

# Test for Flare vs Non Flare
# Use D first
IBSD_all <- rownames(m[m$Cohort=="D",])
D_overlap <- intersect(IBSD_all, rownames(but_kegg_weight))

D_but <- data.frame(but_kegg_weight[D_overlap,])
D_but <- cbind(D_but, Flare = m[D_overlap, "Flare"])
D_but$Flare <- as.character(D_but$Flare)
D_but$Flare[is.na(D_but$Flare)] <- "Normal"

ktest <- kruskal.test(D_but[,'Total'] ~ factor(D_but[,"Flare"]))

Bplot2 <- ggplot(D_but) + 
  geom_boxplot(aes(x=Flare, y=Total), outlier.colour = NA) +
  geom_jitter(aes(x=Flare, y=Total, color=Flare), width=0.08, show.legend = F) +
  scale_color_manual(values=c("#c3c823", cols[2]))
pdf("Weighted_Butyrate_DFlares.pdf", width=4, height=2.5)
Bplot2
dev.off()

### Same but C
IBSD_all <- rownames(m[m$Cohort=="C",])
D_overlap <- intersect(IBSD_all, rownames(but_kegg_weight))

D_but <- data.frame(but_kegg_weight[D_overlap,])
D_but <- cbind(D_but, Flare = m[D_overlap, "Flare"])
D_but$Flare <- as.character(D_but$Flare)
D_but$Flare[is.na(D_but$Flare)] <- "Normal"

ktest <- kruskal.test(D_but[,'Total'] ~ factor(D_but[,"Flare"]))

Bplot2 <- ggplot(D_but) + 
  geom_boxplot(aes(x=Flare, y=Total), outlier.colour = NA) +
  geom_jitter(aes(x=Flare, y=Total, color=Flare), width=0.08, show.legend = F) +
  scale_color_manual(values=c("#c3c823", cols[1]))
pdf("Weighted_Butyrate_CFlares.pdf", width=4, height=2.5)
Bplot2
dev.off()

