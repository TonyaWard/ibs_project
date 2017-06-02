ix <- ixc.hc | ixc.ibsd; disease.case <- "D"
#ri <- cv.risk.index(sqrt(xgnc[ix,]), droplevels(mc$Cohort[ix]) == disease.case, 
#                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
ri <- cv.risk.index(sqrt(xc[ix,]), droplevels(mc$Cohort[ix]) == disease.case, 
                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
sink("risk/risk.txt", append=T)
cat('t-test: hold-out risk index IBS-D vs. healthy:\n')
cat('p=',t.test(ri$risk.index ~ droplevels(mc$Cohort[ix]) != "H")$p.value, '\n', sep='')
sink()

pdf('risk/risk_index_IBSD_v_Healthy.pdf',width=4,height=4)

working_mc <- subset(mc, ix)
working_mc$vari <- ri$risk.index
plot1 <- ggplot(working_mc, aes_string(x="Cohort", y= "vari")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color="Cohort")) +
  #theme(legend.position = 'bottom') + 
  labs(x="", y= "Risk Index") +
  guides(fill=F, color=F) +
  scale_color_manual(values=cols_dh)
print(plot1)
dev.off()

######################################################################
ix <- ixc.hc | ixc.ibsc; disease.case <- "C"
#ri <- cv.risk.index(sqrt(xgnc[ix,]), droplevels(mc$Cohort[ix]) == disease.case, 
#                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
ri <- cv.risk.index(sqrt(xc[ix,]), droplevels(mc$Cohort[ix]) == disease.case, 
                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
sink("risk/risk.txt", append=T)
cat('t-test: hold-out risk index IBS-C vs. healthy:\n')
cat('p=',t.test(ri$risk.index ~ droplevels(mc$Cohort[ix]) != "H")$p.value, '\n', sep='')
sink()

pdf('risk/risk_index_IBSC_v_Healthy.pdf',width=4,height=4)
working_mc <- subset(mc, ix)
working_mc$vari <- ri$risk.index
plot1 <- ggplot(working_mc, aes_string(x="Cohort", y= "vari")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color="Cohort")) +
  #theme(legend.position = 'bottom') + 
  labs(x="", y= "Risk Index") +
  guides(fill=F, color=F) +
  scale_color_manual(values=cols_ch)
print(plot1)
dev.off()

######################################################################
ix <- ixc.hc | ixc.ibs
#ri <- cv.risk.index(sqrt(xgnc[ix,]), droplevels(mc$Cohort[ix]) != "Healthy", 
#                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
ri <- cv.risk.index(sqrt(xc[ix,]), droplevels(mc$Cohort[ix]) != "H", 
                    alpha=.25, correct.pvalues=FALSE, nfolds=-1)
sink("risk/risk.txt", append=T)
cat('t-test: hold-out genus-level risk index IBS vs. healthy:\n')
cat('p=',t.test(ri$risk.index ~ droplevels(mc$Cohort[ix]) != "H")$p.value, '\n', sep='')
sink()

pdf('risk/risk_index_IBS_v_Healthy.pdf',width=4,height=4)
working_mc <- subset(mc, ix)
working_mc$vari <- ri$risk.index
plot1 <- ggplot(working_mc, aes_string(x="IBS", y= "vari")) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(position=position_jitter(0.1), size=3, alpha=0.75, aes_string(color="IBS")) +
  #theme(legend.position = 'bottom') + 
  labs(x="", y= "Risk Index") +
  guides(fill=F, color=F) +
  scale_color_manual(values=cols_ibs)
print(plot1)
dev.off()
