#Plot the metadata to show sample completeness
m_all$label <- paste(m_all$Timepoint, m_all$sType, sep="")
m_all$label <- factor(m_all$label, levels=c("1Biopsy", "0Stool","1Stool","2Stool","3Stool","4Stool","5Stool","6Stool", "2Biopsy", "FlareFlare"))
m_all$Flare <- as.character(m_all$Timepoint)
m_all$Flare[is.na(m_all$Flare)] <- "Normal"
m_all$Flare <- factor(m_all$Flare, levels = c("Normal", "Flare"))
m_all$SB <- as.character(m_all$sType)
m_all$SB[m_all$SB == "Flare"] <- "Stool"

pdf("Sample_overview.pdf", width=5, height=4)
ggplot(m_all) +
  geom_tile(aes(x=label, y=factor(study_id), fill=Cohort), color="white") +
  scale_x_discrete(labels=c("Biopsy","0", "1", "2", "3", "4", "5", "6", "Biopsy", "Flare")) +
  facet_grid(Cohort ~ Flare, scales="free", space="free") +
  scale_fill_manual(values = cols) +
  theme(axis.text.x=element_text(angle=90,hjust=1),
        axis.text.y =element_blank(),
        axis.ticks.y =element_blank(),
        strip.text.x = element_text(size = 8),
        strip.text.y = element_text(size = 8)) +
  labs(y="", x="Timeline") +
  guides(fill=F, color=F, alpha=F)
dev.off()
