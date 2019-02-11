####
# Figure 2  Flares
####
d_alpha_table <- data.frame("Subject"=NA, "mean_noF" = NA, "sd_noF"=NA, "f_alpha"=NA, "norm_alpha" = NA, "Cohort"=NA)

##USE:
# round(flare - mean(alpha for person)/sd(alpha/person), 2)
people <- unique(m_noF$ID_on_tube)
for(i in 1:length(people)){
  ID <- people[i]
  if(ID %in% mF$ID_on_tube){
    mean_noF <- mean(m_noF[m_noF$ID_on_tube == ID, "shannon"])
    sd_noF <- sd(m_noF[m_noF$ID_on_tube == ID, "shannon"])
    f_alpha <-  mF[mF$ID_on_tube == ID,"shannon"]
    norm_alpha <-  round((f_alpha - mean_noF)/sd_noF,2)
    cohort <- mc[mc$ID_on_tube==ID,"Cohort"]
    d_alpha_table <- rbind(d_alpha_table, c(ID, mean_noF, sd_noF, f_alpha, norm_alpha, cohort))
  }
}
d_alpha_table <- d_alpha_table[complete.cases(d_alpha_table),]
d_alpha_table$norm_alpha <- as.numeric(as.character(d_alpha_table$norm_alpha))
d_alpha_table <- d_alpha_table[order(d_alpha_table$Cohort, d_alpha_table$norm_alpha),]
d_alpha_table$Subject <- factor(d_alpha_table$Subject, levels =d_alpha_table$Subject)

alphaF_plot <- ggplot(d_alpha_table, aes(x=Subject, y=norm_alpha, label=norm_alpha)) +
  geom_bar(stat='identity', aes(fill=Cohort)) +
  scale_fill_manual(name="Cohort", values = cols[1:2] )+
  coord_flip() +
  labs(y= "Normalized change in alpha diversity") +
  theme(axis.text.y = element_blank(),text = element_text(size=7),
        axis.text.x  = element_text(size=7))
  

####
# Find diff taxa in IBSD flares vs non and C vs non
test.otu.features<-function(otu, response, sig.level)
{
  pvals <- apply(otu, 2, function(feature) 
    (kruskal.test(feature~response, data.frame(feature=feature, response=response)))$p.value)
  adj.pvals <- p.adjust(pvals, "fdr")
  
  diff.features <- names(adj.pvals)[adj.pvals <= sig.level & !is.na(adj.pvals)]
  list(features=diff.features, pvals=adj.pvals)
}


working_table <- data.frame(t(x_noF), check.names = F)
working_table <- working_table[rowMeans(working_table)> 0.01,]
working_table2 <- data.frame(t(xF), check.names = F)
working_table2 <- working_table2[rowMeans(working_table2)> 0.01,]


#Set up tests to run
C <- m_noF[m_noF$Cohort == "C", "SampleID"]
D <- m_noF[m_noF$Cohort == "D", "SampleID"]
DF <- mF[mF$Cohort == "D", "SampleID"]
CF <-  mF[mF$Cohort == "C", "SampleID"]

test.ixs <- list(C,D)
test.ixs2 <- list(CF,DF)
names(test.ixs) <- c("C","D")
names(test.ixs2) <- c("CF", "DF")

pvals<- c()
ALPHA <- 0.05
diff_list <- list()

working_table$OTU <- rownames(working_table)
working_table2$OTU <- rownames(working_table2)

for(n in 1:length(test.ixs)){
    set1  <- c(test.ixs[n][[1]], "OTU")
    set2 <- c(test.ixs2[n][[1]], "OTU")
    df <- merge(working_table[,set1], working_table2[,set2], by="OTU", all=TRUE)
    rownames(df) <- df$OTU
    df <- df[,2:ncol(df)]
    map <- rbind(m_noF[test.ixs[n][[1]],], mF[test.ixs2[n][[1]],])
    rownames(map) <- map$SampleID
  
    #keep taxa and the samples you are testing
    test_table <- t(df)
    test_table[is.na(test_table)] <- 0
    map_test <- map[rownames(test_table),]
    map_test$Flare <- as.character(map_test$Flare)
    map_test[is.na(map_test$Flare),"Flare"] <- "Non_Flare"
    difftest <- test.otu.features(test_table, response=map_test$Flare, sig.level = 0.10)
    #difftest <- differentiation.test(test_table, map_test$Treatment2, parametric=FALSE)
      
    if(any(difftest$pvals <= ALPHA)){
      signif.ix <- which(difftest$pvals <= ALPHA)
      signif.ix <- signif.ix[order(difftest$pvals[signif.ix])]
      #sink("Diff_Taxa_Flares.txt", append=T)
      this_compare <- paste(names(test.ixs)[n], " vs ", names(test.ixs2)[n])
      #cat(this_compare)
      #cat("\n")
      #print(difftest$pvals[signif.ix])
      #cat("\n")
      #sink()
      diff_list[[this_compare]] <- names(signif.ix)
    } 
} 

diff_list$`C  vs  CF`

df <- merge(working_table, working_table2, by="OTU", all=TRUE)
rownames(df) <- df$OTU
df <- df[,2:ncol(df)]
df <- df[,c(C,D,CF,DF)]
map <- rbind(m_noF, mF)
rownames(map) <- map$SampleID
map <- map[c(C,D,DF,CF),]

#keep taxa and the samples you are testing
test_table2 <- t(df)
test_table2[is.na(test_table2)] <- 0
map_test <- map[rownames(test_table2),]
map_test$Flare <- as.character(map_test$Flare)
map_test[is.na(map_test$Flare) & map_test$Cohort == "C","Flare"] <- "C"
map_test[is.na(map_test$Flare) & map_test$Cohort == "D","Flare"] <- "D"

keep_taxa <- unique(diff_list$`D  vs  DF`, diff_list$`C  vs  CF`)
test_table2 <- test_table2[,keep_taxa]
row.order <- hclust(dist(test_table2))$order # clustering
col.order <- hclust(dist(t(test_table2)))$order
test_table2 <- test_table2[row.order, col.order]
plot_table <- data.frame(test_table2)
plot_table$SampleID <- rownames(plot_table)
plot_table <- melt(plot_table, id=c("SampleID"))
plot_table <- merge(plot_table, map_test, by="SampleID")

plot_table$variable <- factor(plot_table$variable, levels=colnames(test_table2), ordered=TRUE)
plot_table <- plot_table[!is.na(plot_table$variable),]


taxa_F_plot <- ggplot(plot_table, aes(x=SampleID, y=variable)) +
  geom_tile(aes(fill = value)) + 
  scale_fill_gradient(low = "white", high = "steelblue") +
  facet_grid(.~Flare, scales="free") +
  theme_grey() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Sample", y="Taxa") +
  theme(axis.text.x = element_blank() ,
        axis.text.y = element_text(size=7), 
        text=element_text(size=7), legend.position = "bottom")






####
# Plot taxa per person
####
flare_ppl <- unique(plot_table$ID_on_tube[plot_table$Flare == "Flare" & ! is.na(plot_table$Flare)])
new_table <- plot_table[plot_table$ID_on_tube %in% flare_ppl,]
taxa <- unique(new_table$variable)


new_table <- data.frame(lapply(new_table, as.character), stringsAsFactors=FALSE)

taxa_time_table <- new_table[1,]
taxa_time_table[1,] <- NA
taxa_time_table$new_timeline <- NA
for(j in 1:length(taxa)){
  this_taxon <- taxa[j]
  for(i in 1:length(unique(new_table$ID_on_tube))){
    thisID <- unique(new_table$ID_on_tube)[i]
    subtable <- new_table[new_table$ID_on_tube == thisID & new_table$variable == this_taxon,]
    subtable[subtable$Flare =="Flare", "Timepoint"] <- as.numeric(subtable[subtable$Flare =="Flare", "Flare_timepoint.x"])
    subtable$new_timeline <- as.numeric(subtable$Timepoint) - as.numeric(subtable[subtable$Flare == "Flare", "Timepoint"])
    subtable <- subtable[subtable$new_timeline > -1.0,]
    taxa_time_table <- rbind(taxa_time_table, subtable)
  }
}

taxa_time_table$value <- as.numeric(taxa_time_table$value)
taxa_time_table <- taxa_time_table[!is.na(taxa_time_table$variable),]
taxa
#Ups:
UPS <- c( "s__Bacteroides_uniformis",
         "s__Roseburia_hominis",
         "s__Ruminococcus_faecis",
         "s__Ruthenibacterium_lactatiformans",
         "s__Alistipes_shahii",
         "s__Escherichia_coli")
DOWNS<- c("s__Bacteroides_coprocola",
         "s__Parabacteroides_distasonis",
         "s__Alistipes_putredinis",
         "s__Bacteroides_thetaiotaomicron",
         "s__Akkermansia_muciniphila",
         "s__Anaerostipes_hadrus",
         "s__Ruminococcus_bicirculans")

taxa_time_table$joined <- paste(taxa_time_table$ID_on_tube, taxa_time_table$variable)

Up_taxa <-ggplot(taxa_time_table[taxa_time_table$variable %in% UPS,]) +
  #geom_line(aes(x=new_timeline, y=value, group=joined, color=ID_on_tube),
  #          alpha=0.1) +
  geom_smooth(aes(x=new_timeline, y=value, group=variable, color=variable), alpha=0.1, size=0) +
  stat_smooth(aes(x=new_timeline, y=value, group=variable, color=variable), geom="line", size=1.3, alpha=0.65) +
  #geom_smooth(aes(x=new_timeline, y=value, group=variable, color=variable), 
  #            alpha=0.1)+
  scale_color_manual(values = cols[c(4,5,6,7,8,17,11)]) +
  theme( axis.text.x  = element_text(size=7),
        axis.text.y = element_text(size=7), text=element_text(size=7),
        legend.key.size = unit(0.5, "cm"))

down_taxa <- ggplot(taxa_time_table[taxa_time_table$variable %in% DOWNS,]) +
  #geom_line(aes(x=new_timeline, y=value, group=joined, color=ID_on_tube),
  #          alpha=0.1) +
  geom_smooth(aes(x=new_timeline, y=value, group=variable, color=variable), alpha=0.1, size=0) +
  stat_smooth(aes(x=new_timeline, y=value, group=variable, color=variable), geom="line", size=1.3, alpha=0.65) +
  #geom_smooth(aes(x=new_timeline, y=value, group=variable, color=variable), 
  #            alpha=0.1) +
  scale_color_manual(values = cols[c(4,5,6,7,8,17,20)]) +
  theme(axis.text.x  = element_text(size=7),
        axis.text.y = element_text(size=7), text=element_text(size=7),
        legend.key.size = unit(0.5, "cm"))
  



left_side <- plot_grid(alphaF_plot, Up_taxa, down_taxa, ncol=1)
Figure2 <- plot_grid(left_side, taxa_F_plot)
