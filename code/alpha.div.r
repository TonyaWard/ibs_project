#Calculate Alpha Div
alpha_div <- as.data.frame(diversity(x, index="shannon"))
colnames(alpha_div) <- "shannon"
alpha_div$simpson <- diversity(x, index="simpson")
alpha_div$obs_species <- rowSums(x > 0)
alpha_div <- alpha_div[rownames(m),]
m$shannon <- alpha_div$shannon
m$observed_species <- alpha_div$obs_species
m$simpson <- alpha_div$simpson

alpha_divC <- as.data.frame(diversity(xc, index="shannon"))
colnames(alpha_divC) <- "shannon"
alpha_divC$simpson <- diversity(xc, index="simpson")
alpha_divC$obs_species <- rowSums(xc > 0)
alpha_divC <- alpha_divC[rownames(mc),]
mc$shannon <- alpha_divC$shannon
mc$observed_species <- alpha_divC$obs_species
mc$simpson <- alpha_divC$simpson

#test alpha diversity across cohorts
#mc$Cohort <- factor(mc$Cohort, levels=c("H", "C", "D"))

ktest <- kruskal.test(alpha_divC[,'shannon'], mc[,"Cohort"])
ktest2 <- kruskal.test(alpha_divC[,'simpson'], mc[,"Cohort"])
ktest3 <- kruskal.test(alpha_divC[,'obs_species'], mc[,"Cohort"])

#write stats to file
file_name <- "alpha_div.txt"
sink(file_name, append =TRUE)
cat("Shannon kruskal test:\n")
print(ktest)
cat("\nSimpson kruskal test:\n")
print(ktest2)
cat("\nObserved species kruskal test:\n")
print(ktest3)
sink()

#assign pdf name for plot
file_path <- "alpha_div.pdf"
pdf(file_path, height=4,width=6)
boxplot(alpha_divC[,"shannon"] ~ mc[,"Cohort"],
        xlab='',ylab="shannon", col = cols[1:3])
boxplot(alpha_divC[,"simpson"] ~ mc[,"Cohort"],
        xlab='',ylab="simpson", col = cols[1:3])
boxplot(alpha_divC[,"obs_species"] ~ mc[,"Cohort"],
        xlab='',ylab="obs_species", col = cols[1:3])
dev.off()

