#Calculate Alpha Div for all
alpha_div <- as.data.frame(diversity(kegg, index="shannon"))
colnames(alpha_div) <- "shannon"
alpha_div$simpson <- diversity(kegg, index="simpson")
alpha_div$obs_species <- rowSums(kegg > 0)
alpha_div <- alpha_div[rownames(m),]
m$shannon_kegg <- alpha_div$shannon
m$observed_species_kegg <- alpha_div$obs_species
m$simpson_kegg <- alpha_div$simpson

#Alpha Div for collapsed
alpha_divC <- as.data.frame(diversity(keggc, index="shannon"))
colnames(alpha_divC) <- "shannon"
alpha_divC$simpson <- diversity(keggc, index="simpson")
alpha_divC$obs_species <- rowSums(keggc > 0)
alpha_divC <- alpha_divC[rownames(mc),]
mc$shannon_kegg <- alpha_divC$shannon
mc$obs_spec_kegg <- alpha_divC$obs_species
mc$simpson_kegg <- alpha_divC$simpson

######################################################################
#test alpha diversity across cohorts using collapsed
#using a shapiro test, they aren't normal
alphas <- c("shannon", "simpson", "obs_species")
for(i in 1:3){
  print(shapiro.test((alpha_divC[,alphas[i]])))
}

#Use non-parametric tests:
ktest <- kruskal.test(alpha_divC[,'shannon'] ~ mc[,"Cohort"])
ktest2 <- kruskal.test(alpha_divC[,'simpson'], mc[,"Cohort"])
ktest3 <- kruskal.test(alpha_divC[,'obs_species'], mc[,"Cohort"])

#Test pairwise for shannon:
shannon_p <- c()
test.ixs <- list(ix.hc, ix.ibsc, ix.ibsd)
names(test.ixs) <- c("Healthy", "IBS-C", "IBS-D")
for(n in 1:(length(test.ixs)-1)){
  for(i in (n+1):length(test.ixs)){
    shannon_p <- c(shannon_p, wilcox.test(alpha_divC[test.ixs[[n]], "shannon"], alpha_divC[test.ixs[[i]], "shannon"])$p.value)
    names(shannon_p)[length(shannon_p)] <- paste(names(test.ixs)[n], " vs ", names(test.ixs)[i], sep="")
  }
}

#write stats to file
file_name <- "KEGG/alpha_div/alpha_div.txt"
sink(file_name, append =TRUE)
cat("Shannon kruskal test:\n")
print(ktest)
cat("Pairwsie Shannon Test:\n")
print(shannon_p)
cat("\nSimpson kruskal test:\n")
print(ktest2)
cat("\nObserved species kruskal test:\n")
print(ktest3)
sink()

#Note - MC doesn't have flares
#assign pdf name for plot
alpha_plot <- cbind(alpha_divC, mc[,"Cohort"])
colnames(alpha_plot)[4] <- "Cohort"

file_path <- "KEGG/alpha_div/alpha_div.pdf"
pdf(file_path, height=4,width=6)
for(i in 1:3){
  plot1 <- ggplot(alpha_plot, aes_string(y=colnames(alpha_plot)[i], x="Cohort")) +
    geom_boxplot(outlier.shape = NA, aes(color=Cohort)) +
    geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes(color=Cohort)) +
    labs(x="", y= colnames(alpha_plot)[i]) +
    guides(fill=F, color=F) +
    scale_color_manual(values=cols[1:3])
  plot(plot1)
}
dev.off()

######################################################################
#test alpha diversity changes over time, use shannon
set.seed(30)

working_alpha <- melt(m, id.vars = c("SampleID", "Cohort", "Timepoint", "ID_on_tube"), measure.vars = c("shannon_kegg"))
working_alpha <- droplevels(working_alpha[! working_alpha$Timepoint == "Flare",])
working_alpha$Timepoint <- as.numeric(droplevels(working_alpha$Timepoint))
working_alpha$ID_on_tube <- factor(working_alpha$ID_on_tube)
  
##Permuation based test to see if different from random
working_alpha1 <- droplevels(working_alpha[working_alpha$Cohort =="C",])
alpha <- working_alpha1$value
subject <- working_alpha1$ID_on_tube
time <- working_alpha1$Timepoint
obs1 <- -mean(sapply(split(1:nrow(working_alpha1), subject), 
                     function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
mc.stats1 <- -replicate(999,mean(sapply(split(1:nrow(working_alpha1), subject), 
                                        function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
pval1 <- mean(c(obs1,mc.stats1) >= obs1)

working_alpha2 <- droplevels(working_alpha[working_alpha$Cohort =="D",])
alpha <- working_alpha2$value
subject <- working_alpha2$ID_on_tube
time <- working_alpha2$Timepoint
obs2 <- -mean(sapply(split(1:nrow(working_alpha2), subject), 
                     function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
mc.stats2 <- -replicate(999,mean(sapply(split(1:nrow(working_alpha2), subject), 
                                        function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
pval2 <- mean(c(obs2,mc.stats2) >= obs2)

working_alpha3 <- droplevels(working_alpha[working_alpha$Cohort =="H",])
alpha <- working_alpha3$value
subject <- working_alpha3$ID_on_tube
time <- working_alpha3$Timepoint
obs3 <- -mean(sapply(split(1:nrow(working_alpha3), subject), 
                     function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], time[ixx], method='spear')$statistic))
mc.stats3 <- -replicate(999,mean(sapply(split(1:nrow(working_alpha3), subject), 
                                        function(ixx) if(length(ixx) < 3) 0 else cor.test(alpha[ixx], sample(time[ixx]), method='spear')$statistic)))
pval3 <- mean(c(obs3,mc.stats3) >= obs3)
  
figure <- ggplot(working_alpha, aes_string(x="Timepoint", y="value", color="Cohort", group="Cohort")) +
  geom_jitter(width = 0.25, size=4, alpha=0.65) +
  geom_smooth(method=lm, linetype = "dashed", se=FALSE) +
  scale_color_manual(values=c("#cb1b4a", "#42aeb8", "#FDB316")) +
  theme(legend.title=element_blank()) +
  annotate("text", x=6, y=7.05, label= paste("P=", round(pval1, digits=3)), size=2, col="#cb1b4a") +
  annotate("text", x=6, y=7.15, label= paste("P=", round(pval2, digits=3)), size=2, col="#42aeb8") +
  annotate("text", x=6, y=7.25, label= paste("P=", round(pval3, digits=3)), size=2, col="#FDB316") +
  labs(x="Timepoint", y="shannon")

file_path <- "KEGG/alpha_div/alpha_div_time.pdf"
pdf(file_path, height=4,width=6)
plot(figure)
dev.off()

#drop people with < 3 samples
working_alpha2 <- with(working_alpha,by(value,ID_on_tube,function(xx)sum(xx > 0)))

#get the flare samples 
figure2 <- ggplot(working_alpha[working_alpha$ID_on_tube %in% names(which(working_alpha2>=3)),]) +
  geom_line(aes(x=Timepoint, y=value, color=Cohort, group=ID_on_tube), alpha = .25, stat = "smooth", method = "loess") +
  scale_color_manual(values=c("#cb1b4a", "#42aeb8", "#FDB316")) +
  theme(legend.title=element_blank()) +
  scale_y_continuous(limits = c(min(working_alpha$value), max(working_alpha$value)))+
  geom_smooth(method = "loess", se=FALSE,  aes_string(x="Timepoint", y="value", color="Cohort", group="Cohort")) +
  annotate("text", x=6, y=4.35, label= paste("P=", round(pval1, digits=3)), size=2, col="#cb1b4a") +
  annotate("text", x=6, y=4.25, label= paste("P=", round(pval2, digits=3)), size=2, col="#42aeb8") +
  annotate("text", x=6, y=4.15, label= paste("P=", round(pval3, digits=3)), size=2, col="#FDB316") +
  labs(x="Timpoint", y="shannon") +
  geom_jitter(size=1.5, data=m[m$Flare == "Flare" & !is.na(m$Flare),], aes(y=shannon_kegg, x= Flare_timepoint.x, color= Cohort), width = 0.25)

file_path <- "KEGG/alpha_div/alpha_div_time2.pdf"
pdf(file_path, height=4,width=6)
plot(figure2)
dev.off()

######################################################################
#Test difference of within-person variance across treatments
working_alpha1 <- droplevels(working_alpha[working_alpha$Cohort =="C",])
alpha <- working_alpha1$value
subject <- working_alpha1$ID_on_tube
time <- working_alpha1$Timepoint
var_c <- sapply(split(1:nrow(working_alpha1), subject),
                function(ixx) if(length(ixx) < 3) NA else var(alpha[ixx]))
new1 <- t(data.frame(as.list(var_c), check.names = F))
new1 <- cbind(new1, rep("C", nrow(new1)))
colnames(new1) <- c("variance", "Cohort")

working_alpha2 <- droplevels(working_alpha[working_alpha$Cohort =="D",])
alpha <- working_alpha2$value
subject <- working_alpha2$ID_on_tube
time <- working_alpha2$Timepoint
var_d <- sapply(split(1:nrow(working_alpha2), subject),
                function(ixx) if(length(ixx) < 3) NA else var(alpha[ixx]))
new2 <- t(data.frame(as.list(var_d), check.names = F))
new2 <- cbind(new2, rep("D", nrow(new2)))
colnames(new2) <- c("variance", "Cohort")

working_alpha3 <- droplevels(working_alpha[working_alpha$Cohort =="H",])
alpha <- working_alpha3$value
subject <- working_alpha3$ID_on_tube
time <- working_alpha3$Timepoint
var_h <- sapply(split(1:nrow(working_alpha3), subject),
                function(ixx) if(length(ixx) < 3) NA else var(alpha[ixx]))
new3 <- t(data.frame(as.list(var_h), check.names = F))
new3 <- cbind(new3, rep("H", nrow(new3)))
colnames(new3) <- c("variance", "Cohort")

total <- data.frame(rbind(new1, new2, new3))
total <- total[!is.na(total$variance),]
total$variance <- as.numeric(as.character(total$variance))

file_path <- "KEGG/alpha_div/Alpha_Variance.pdf"
pdf(file_path, height=4,width=6)
plot1 <- ggplot(total, aes(y=variance, x=Cohort)) +
    geom_boxplot(outlier.shape = NA, aes(color=Cohort)) +
    geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes(color=Cohort)) +
    labs(x="", y= "variance") +
    guides(fill=F, color=F) +
    scale_color_manual(values=cols[1:3])
plot(plot1)
dev.off()

######################################################################
#Plot flare kegg alphas
# Just plot the alpha as a boxplot + the point for the flare sample:
file_path <- "KEGG/alpha_div/alpha_div_flares_patients.pdf"
flare_patients <- m[m$Timepoint == "Flare", "study_id"]
working_alpha <- m[m$study_id %in% flare_patients,]
working_alpha$study_id <- factor(working_alpha$study_id)
pdf(file_path, height=4,width=6)
plot1 <- ggplot() +
  geom_boxplot(data = working_alpha[is.na(working_alpha$Flare),], 
               outlier.shape = NA, aes(y=shannon_kegg, x=study_id, color=Cohort)) +
  geom_jitter(data = working_alpha[is.na(working_alpha$Flare),], 
              position=position_jitter(0.1), size=3, alpha=0.75, aes(y=shannon_kegg, x=study_id, color=Cohort)) +
  guides(fill=F, color=F) +
  scale_color_manual(values=cols[1:3]) +
  geom_jitter(size=3, shape=17, color= "#c3c823",
              data=working_alpha[working_alpha$Flare == "Flare" & !is.na(working_alpha$Flare),], 
              width = 0.25, aes(y=shannon_kegg, x=study_id)) + 
  theme(axis.text.x=element_text(angle=90,hjust=1))+
  facet_grid(.~ Cohort, space="free", scales="free")
plot(plot1)
dev.off()
#Test for paired difference:
c_outs <- c()
d_outs <- c()
all_outs <- list()
for(i in 1:length(flare_patients)){
  working_ID <- flare_patients[i]
  working_table <- m[m$study_id == working_ID,]
  flare_sample <- working_table[working_table$Flare == "Flare" & !is.na(working_table$Flare), "shannon_kegg"]
  others <-  working_table[is.na(working_table$Flare), "shannon_kegg"]
  if(length(others) < 2){
    print("only one non-flare")
  } else {
    outcome <- t.test(others, mu=flare_sample, alternative="two.sided", conf.level=0.99)$p.value
  }
  if("D" %in% working_table$Cohort){
    d_outs <- c(d_outs, outcome)
  } else {
    c_outs <- c(c_outs, outcome)
  }
  all_outs[[i]] <- outcome
  names(all_outs)[i] <- as.character(working_ID)
}
sink("KEGG/alpha_div/paired_flare_stats.txt")
print(all_outs)
sink()

