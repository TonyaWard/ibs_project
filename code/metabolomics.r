#Load in metabonomics
library(microbiome)
metabo <- read.table("../data/IBS_Mayo_secondrun/metabonomics_v1.txt",
                     sep="\t",
                     header=T,
                     check.names=F,
                     stringsAsFactors = F,
                     row=1)
keeps <- intersect(rownames(metabo), rownames(m))

metabo <- metabo[keeps,c(1:7)]
met_m <- m[keeps,]
met_otu <- x[keeps,]

Ds <- rownames(met_m)[met_m$Cohort == "D"]
Cs <- rownames(met_m)[met_m$Cohort == "C"]
Hs <- rownames(met_m)[met_m$Cohort == "H"]

correlations <- associate(met_otu[Ds,], metabo[Ds,], method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
correlation.table.D <- associate(met_otu[Ds,], metabo[Ds,], method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
correlation.table.C <- associate(met_otu[Cs,], metabo[Cs,], method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)
correlation.table.H <- associate(met_otu[Hs,], metabo[Hs,], method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

cor.D <- heat(correlation.table.D, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.005, star= "p.adj", colours = c("orange", "white", "purple")) 
cor.C <- heat(correlation.table.C, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.005, star= "p.adj", colours = c("orange", "white", "purple")) 
cor.H <- heat(correlation.table.H, "X1", "X2", fill = "Correlation", p.adj.threshold = 0.005, star= "p.adj", colours = c("orange", "white", "purple")) 

overlap_taxa_disease <- intersect(correlation.table.C$X1, correlation.table.D$X1)
overlap_taxa_3 <- intersect(correlation.table.H$X1, overlap_taxa_disease)

correlation.table.C$Type <- "IBS-C"
correlation.table.D$Type <- "IBS-D"
correlation.table.H$Type <- "IBS-H"

together_table <- rbind(correlation.table.C, correlation.table.D)
together_table <- rbind(together_table, correlation.table.H)  

together_table_D <- together_table[together_table$X1 %in% overlap_taxa_disease,]
together_table_3 <- together_table[together_table$X1 %in% overlap_taxa_3,]
keeps <- as.character(together_table_3$X1[which(duplicated(together_table_3$X1))])

heat2(together_table_3[together_table_3$X1 %in% keeps,], "X1", "X2", fill= "Correlation", colours = c("orange", "white", "purple"))
heat2(together_table_3, "X1", "X2", fill= "Correlation")

View(heat)

together <-plot_grid(cor.D, cor.C, cor.H, nrow=3, labels=c("IBS-D", "IBS-C", "IBS-H"))
pdf("metabolite_bacteria_correlations.pdf", height=20, width=10.5)
together
dev.off()


heat2 <- function (df, Xvar = names(df)[[1]], Yvar = names(df)[[2]], 
          fill = names(df)[[3]], LB = names(df)[[5]], star = NULL, p.adj.threshold = 1, 
          association.threshold = 0, step = 0.2, colours = c("darkblue", 
                                                             "blue", "white", "red", "darkred"), limits = NULL, legend.text = "", 
          order.rows = TRUE, order.cols = TRUE, text.size = 10, filter.significant = TRUE, 
          star.size = NULL, plot.values = FALSE) 
{
  if (is.null(limits)) {
    maxval <- max(abs(df[[fill]]))
    if (maxval <= 1) {
      limits <- c(-1, 1)
    }
    else {
      xmaxval <- ceiling(maxval)
      limits <- c(-maxval, maxval)
    }
  }
  df[df[[fill]] < limits[[1]], fill] <- limits[[1]]
  df[df[[fill]] > limits[[2]], fill] <- limits[[2]]
  if (nrow(df) == 0) {
    warning("Input data frame is empty.")
    return(NULL)
  }
  if (filter.significant & !is.null(star)) {
    keep.X <- as.character(unique(df[((df[[star]] < p.adj.threshold) & 
                                        (abs(df[[fill]]) > association.threshold)), Xvar]))
    keep.Y <- as.character(unique(df[((df[[star]] < p.adj.threshold) & 
                                        (abs(df[[fill]]) > association.threshold)), Yvar]))
    df <- df[((df[[Xvar]] %in% keep.X) & (df[[Yvar]] %in% 
                                            keep.Y)), ]
  }
  theme_set(theme_bw(text.size))
  if (any(c("XXXX", "YYYY", "ffff") %in% names(df))) {
    stop("XXXX, YYYY, ffff are not allowed in df")
  }
  df[[Xvar]] <- factor(df[[Xvar]])
  df[[Yvar]] <- factor(df[[Yvar]])
  if (is.logical(order.rows) || is.logical(order.cols)) {
    rnams <- unique(as.character(df[[Xvar]]))
    cnams <- unique(as.character(df[[Yvar]]))
    mat <- matrix(0, nrow = length(rnams), ncol = length(cnams))
    rownames(mat) <- rnams
    colnames(mat) <- cnams
    for (i in 1:nrow(df)) {
      mat[as.character(df[i, Xvar]), as.character(df[i, 
                                                     Yvar])] <- df[i, fill]
    }
    mat <- t(mat)
    cind <- 1:ncol(mat)
    rind <- 1:nrow(mat)
  }
  if (is.logical(order.rows)) {
    if (order.rows) {
      if (nrow(mat) > 1 && ncol(mat) > 1) {
        rind <- hclust(as.dist(1 - cor(t(mat), use = "pairwise.complete.obs")))$order
      }
      if (nrow(mat) > 1 && ncol(mat) == 1) {
        rind <- order(mat[, 1])
      }
      order.rows <- rownames(mat)[rind]
    }
    else {
      order.rows <- rownames(mat)[rind]
    }
  }
  if (is.logical(order.cols)) {
    if (order.cols) {
      if (ncol(mat) > 1 && nrow(mat) > 1) {
        cind <- hclust(as.dist(1 - cor(mat, use = "pairwise.complete.obs")))$order
      }
      else {
        cind <- order(mat[1, ])
      }
      order.cols <- colnames(mat)[cind]
    }
    else {
      order.cols <- colnames(mat)[cind]
    }
  }
  df[[Yvar]] <- factor(df[[Yvar]], levels = order.rows)
  df[[Xvar]] <- factor(df[[Xvar]], levels = order.cols)
  XXXX <- YYYY <- ffff <- NULL
  df[["XXXX"]] <- df[[Xvar]]
  df[["YYYY"]] <- df[[Yvar]]
  df[["ffff"]] <- df[[fill]]
  p <- ggplot(df, aes(x = XXXX, y = YYYY, fill = ffff))
  p <- p + geom_tile()
  p <- p + scale_fill_gradientn(legend.text, breaks = seq(from = min(limits), 
                                                          to = max(limits), by = step), colours = colours, limits = limits)
  p <- p + xlab("") + ylab("")
  p <- p + theme(axis.text.x = element_text(angle = 90))
  p <- p + facet_grid(df[[LB]] ~ ., scales="free", space = "free")
  if (!is.null(star)) {
    inds <- which((df[[star]] < p.adj.threshold) & (abs(df[[fill]]) > 
                                                      association.threshold))
    if (!is.null(star) & length(inds) > 0) {
      df.sub <- df[inds, ]
      if (is.null(star.size)) {
        star.size <- max(1, floor(text.size/2))
      }
      p <- p + geom_text(data = df.sub, aes(x = XXXX, 
                                            y = YYYY, label = "+"), col = "white", size = star.size)
    }
  }
  if (plot.values) {
    p <- p + geom_text(aes(label = round(ffff, 2)), size = 3)
  }
  p
}



