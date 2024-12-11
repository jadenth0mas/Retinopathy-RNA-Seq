# Pubhbio 5280 Final Project
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221521
# Jaden Thomas
#setwd("/Users/jadenthomas/Desktop/OSU/OSU Senior/Fall/PUBHBIO 5280/Final Project")


library(tidyverse)
library(Biobase)
library(BiocManager)
library(GEOquery)
library(DESeq2)
library(ggrepel)
library(gt)
library(cowplot)
library(clusterProfiler)
library(glmnet)
library(preprocessCore)
library(org.Hs.eg.db)
library(pROC)
library(randomForest)
library(cluster)
library(Rtsne)
library(umap)
library(factoextra)
library(gplots)
library(made4)
library(dendextend)
library(ggthemes)
library(webshot2)

#dir()
#load("pubhbio-final-plots.RData")

# Download Data
pheno <- getGEO(GEO="GSE221521")[[1]]
sup <- getGEOSuppFiles("GSE221521")
gznames <- file.path("GSE221521", basename(rownames(sup)))
count.data <- read.table(gzfile(gznames[1]), header=T, sep="\t", quote = "", row.names=1)

# Get Genes
genes <- data.frame(count.data$gene_name)
rownames(genes) <- rownames(count.data)


# Get Raw Counts
expression <- count.data %>% dplyr::select(gene_name, locus, contains("_count")) %>% 
  dplyr::select(!c(locus, gene_name))

# Separate to get group and ID
ph <- pData(pheno) %>% separate(title, into=c("cell", "group_id"), sep=",", remove=T) %>%
  separate(group_id, into=c("Group", "id"), sep=" group ", remove=T) %>%
  dplyr::select(c(Group, id, geo_accession, growth_protocol_ch1, extract_protocol_ch1, 
                  extract_protocol_ch1.1, data_processing, data_processing.1, 
                  data_processing.2, data_processing.3, data_processing.4, 
                  data_processing.5, data_processing.6))



# Change rowname to IDs found in count data
rownames(ph) <- NULL
ph <- ph %>% column_to_rownames(var="id")

# Remove _count from colnames to match with ph
l <- colnames(expression)
clean_l <- gsub("_count", "", l)
colnames(expression) <- clean_l

# Fix Group having extra space
ph <- ph %>% dplyr::select(Group) %>% mutate(Group=gsub(" ", "", Group, fixed=T))
ph$Group <- as.factor(ph$Group)

# **Pre-processing**


## Filtering
expression_filter <- expression[rowMeans(expression)>10,]


## Get DESeqDataset
dds <- DESeqDataSetFromMatrix(countData=as.matrix(expression_filter), colData=ph, design=~Group)

## Get Negative Binomial Distribution Differential Expression
set.seed(2024)
dds <- DESeq(dds)

# VST
ex_vst <- vst(assay(dds))
temp <- data.frame(rowname=rownames(ex_vst)) %>% left_join(rownames_to_column(genes, "rowname"))
rownames(ex_vst) <- temp$count.data.gene_name

# Quantile Normalization
ex_qn <- normalize.quantiles(ex_vst)
rownames(ex_qn) <- rownames(ex_vst)
colnames(ex_qn) <- colnames(ex_vst)

# 16697 original genes and 193 original patients
dim(ex_qn)


# Method 1 **Differential Expression**

## Do Pairwise comparisons
results.diabetes.vs.control <- results(dds, contrast=c("Group", "DM", "Control"))
results.diabetes.vs.retinopathy <- results(dds, contrast=c("Group", "DR", "DM"))
results.retinopathy.vs.control <- results(dds, contrast=c("Group", "DR", "Control"))

## Get summaries
summary(results.diabetes.vs.control, alpha=0.05)
summary(results.diabetes.vs.retinopathy, alpha=0.05)
summary(results.retinopathy.vs.control, alpha=0.05)


# Not included in paper
lfc_drdm <- lfcShrink(dds, type="ashr", contrast=c("Group", "DR", "DM"))
lfc_drc <- lfcShrink(dds, type="ashr", contrast=c("Group", "DR", "Control"))

# Not included in paper, used for initial EDA and analysis
par(mfrow=c(1, 2))
plotMA(lfc_drdm, alpha=0.05, main="DR vs DM")
plotMA(lfc_drc, alpha=0.05, main="DR vs C")
plotMA(results.diabetes.vs.retinopathy, colSig="blue", colNonSig="grey70", alpha=0.05, main="DR vs DM")
plotMA(results.retinopathy.vs.control, colSig="blue", colNonSig="grey70", alpha=0.05, main="DR vs C")
par(mfrow=c(1, 1))

# Get top over and underexpressed genes for DE
d.v.dr.results %>% filter(padj<0.05 & log2FoldChange< -1) %>% arrange(padj)
d.v.dr.results %>% filter(padj<0.05 & log2FoldChange>1) %>% arrange(padj)


dr.v.c.results %>% filter(padj<0.05 & log2FoldChange>1) %>% arrange(padj)
dr.v.c.results %>% filter(padj<0.05 & log2FoldChange< -1) %>% arrange(padj)



# Convert to df to use ggplot
d.v.c.results <- results.diabetes.vs.control %>% as.data.frame() %>% rownames_to_column(var="rowname") %>% left_join(rownames_to_column(genes, "rowname"))
d.v.dr.results <- results.diabetes.vs.retinopathy %>% as.data.frame() %>% rownames_to_column(var="rowname") %>% left_join(rownames_to_column(genes, "rowname"))
dr.v.c.results <- results.retinopathy.vs.control %>% as.data.frame() %>% rownames_to_column(var="rowname") %>% left_join(rownames_to_column(genes, "rowname"))

# Get differentially expressed genes for each
gene_to_label <- d.v.c.results[abs(d.v.c.results$log2FoldChange) > 2 & d.v.c.results$padj<0.05,]$count.data.gene_name
gene_to_label <- gene_to_label[!is.na(gene_to_label)]

dr_gene_to_label <- d.v.dr.results[abs(d.v.dr.results$log2FoldChange) > 2 & d.v.dr.results$padj<0.05,]$count.data.gene_name
dr_gene_to_label <- dr_gene_to_label[!is.na(dr_gene_to_label)]

dr_v_gene_to_label <- dr.v.c.results[abs(dr.v.c.results$log2FoldChange) > 2 & dr.v.c.results$padj<0.05,]$count.data.gene_name
dr_v_gene_to_label <- dr_v_gene_to_label[!is.na(dr_v_gene_to_label)]


# Table 1
table_1 <- table(ph) %>% as_tibble()
table_1_gt <- gt(table_1) %>% tab_header(
  title = "Distribution of Response Variable"
)

#gtsave(data=table_1_gt, filename="table1.png")



# DM vs Control, not included in analysis
dm_v_c_volcano_plot <- ggplot(d.v.c.results, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj<0.05 & abs(log2FoldChange)>2), alpha=0.9) +
  theme_minimal() +
  labs(x=expression(log[2](FoldChange)), y=expression(-log[10](FDR)), title="Diabetic Mellitus Vs Control") +
  geom_label_repel(
    data = subset(d.v.c.results, count.data.gene_name %in% gene_to_label),
    aes(label = count.data.gene_name),
    box.padding = 0.5,   # Adjust padding
    point.padding = 0.3, # Adjust point padding
    segment.color = 'gray50',
    max.overlaps = 10
  )  +
  scale_color_manual(values=c('grey60', 'blue')) +
  theme(legend.position="none") +
  geom_hline(yintercept=-log10(0.05), color="black", alpha=0.5, lty=2) +
  geom_vline(xintercept=c(-2, 2), color="black", alpha=0.5, lty=2)


# DM vs DR
dm_vs_dr_volcano_plot <- ggplot(d.v.dr.results, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj<0.05 & abs(log2FoldChange)>2), alpha=0.9) +
  theme_minimal() +
  labs(x=expression(log[2](FoldChange)), y=expression(-log[10](FDR)), title="Diabetic Retinopathy vs. Diabetic Mellitus") +
  geom_label_repel(
    data = subset(d.v.dr.results, count.data.gene_name %in% dr_gene_to_label),
    aes(label = count.data.gene_name),
    box.padding = 0.5,   # Adjust padding
    point.padding = 0.3, # Adjust point padding
    segment.color = 'gray50',
    max.overlaps = 15
  ) +
  scale_color_manual(values=c('grey60', 'blue')) +
  theme(legend.position="none") +
  geom_hline(yintercept=-log10(0.05), color="black", alpha=0.5, lty=2) +
  geom_vline(xintercept=c(-2, 2), color="black", alpha=0.5, lty=2)


# DR vs Control

dr_v_c_volcano_plot <- ggplot(dr.v.c.results, aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(aes(color=padj<0.05 & abs(log2FoldChange)>2), alpha=0.9) +
  theme_minimal() +
  labs(x=expression(log[2](FoldChange)), y=expression(-log[10](FDR)), title="Diabetic Retinopathy Vs Control") +
  geom_label_repel(
    data = subset(dr.v.c.results, count.data.gene_name %in% dr_v_gene_to_label),
    aes(label = count.data.gene_name),
    box.padding = 0.5,   # Adjust padding
    point.padding = 0.3, # Adjust point padding
    segment.color = 'gray50',
    max.overlaps = 15
  )  +
  scale_color_manual(values=c('grey60', 'blue')) +
  theme(legend.position="none") +
  geom_hline(yintercept=-log10(0.05), color="black", alpha=0.5, lty=2) +
  geom_vline(xintercept=c(-2, 2), color="black", alpha=0.5, lty=2)


# Top of figure 1
fig_1_top <- plot_grid(dm_vs_dr_volcano_plot,
                       dr_v_c_volcano_plot,
                       ncol=2,
                       labels = c('A', 'B'), label_size = 12)
fig_1_top

# Gene ontology analysis
dr_go <- d.v.dr.results[d.v.dr.results$padj<0.05,]$count.data.gene_name
dr_go <- dr_go[!is.na(dr_go)]


ego_dr_v_dm <- enrichGO(gene=dr_go,
                OrgDb=org.Hs.eg.db,
                keyType="SYMBOL",
                ont="BP",
                pAdjustMethod="BH",
                qvalueCutoff=0.05)



drdm_dotplot <- dotplot(ego_dr_v_dm, font.size=10)

dr_v_c_go <- dr.v.c.results[dr.v.c.results$padj<0.05,]$count.data.gene_name
dr_v_c_go <- dr_v_c_go[!is.na(dr_v_c_go)]


ego_dr_v_c <- enrichGO(gene=dr_v_c_go,
                        OrgDb=org.Hs.eg.db,
                        keyType="SYMBOL",
                        ont="BP",
                        pAdjustMethod="BH",
                        qvalueCutoff=0.05)

drc_dotplot <- dotplot(ego_dr_v_c, font.size=10)

fig_1_bottom <- plot_grid(drdm_dotplot, drc_dotplot, nrow=1, labels=c("C", "D"), label_size = 12)

# Assemble figure 1
figure1 <- plot_grid(fig_1_top, fig_1_bottom, nrow=2)


# Method 2 - Lasso
set.seed(2024)

lasso_data <- filter(ph, Group %in% c("DR", "DM")) %>% mutate(Group=factor(Group, levels=c("DM", "DR")))
x <- t(ex_vst)
x <- x[row.names(x) %in% row.names(lasso_data),] %>% scale()

train <- sample(c(T, F), nrow(x), replace=T)
x.train <- x[train,]
x.test <- x[!train,]
y.train <- lasso_data[train,]
y.test <- lasso_data[!train,]

set.seed(2024)
cv_lasso <- cv.glmnet(x.train, y.train, alpha=1, family="binomial")

# Figure 2A
plot(cv_lasso, xlab="Log(Lambda)", font.lab=1, cex.lab=1.9, cex.axis=1.1)
fig2_a <- recordPlot()

# Final lasso model
lasso.fit <- glmnet(x.train, y.train, alpha=1, family="binomial", lambda=cv_lasso$lambda.min)

y.pred <- predict(lasso.fit, x.test, type="response")
y.pred.class <- ifelse(y.pred>.5, "DR", "DM")

# Get nonzero coefs
lasso.coef <- coef(lasso.fit)
nonzero_coefs <- rownames(lasso.coef)[abs(lasso.coef[,1])>0]
nonzero_coefs <- nonzero_coefs[nonzero_coefs!="(Intercept)"]
lasso_nonzero_coefs <- data.frame(Gene=nonzero_coefs) %>% mutate("Lasso Beta"=lasso.coef[abs(lasso.coef[,1])>0][-1])

# Number of nonzero coefs
length(nonzero_coefs)

# Table 2
table2 <- lasso_nonzero_coefs %>% arrange(desc(abs(`Lasso Beta`))) %>% head(n=10) %>% gt()
#gtsave(data=table2, filename="table2.png")


# Lasso test error
(1-mean(y.test==y.pred.class))*100

# ROC and AUC
lasso_roc <- roc(y.test, c(y.pred))
lasso_auc <- auc(lasso_roc)

lasso_auc


# Method 3 Random Forest


set.seed(2024)
rf_dr <- randomForest(x=x.train, y=y.train, importance=T, ntree=500)


# Random forest OOB error rate, DM and DR error rate
rf_dr$err.rate[nrow(rf_dr$err.rate),]
1-0.2328767


# Figure 2B
varImpPlot(rf_dr, type=1, main="Random Forest Variable Importance",
 font.main=1, cex.main=2, cex=0.8, cex.lab=1.7, cex.axis=.9, cex.sub=.8) # Change rownames before t to get gene names

fig2_b <- recordPlot()


# Not included in analysis, OOB, DR, and DM error rate across trees
plot(rf_dr, col=c("darkred", "black", "darkgreen"), main="Random Forest Error",
     cex.main=1, cex.lab=.8, cex.axis=.9, lty=1)
legend("bottomright", cex=.9, legend=c("OOB", "DM", "DR"), col=c("darkred", "black", "darkgreen"), lwd=2)

fig2_c <- recordPlot()

# Predictive probabilities
rf_preds <- predict(rf_dr, x.test, type="prob")[,2]
rf.pred.class <- ifelse(rf_preds>.5, "DR", "DM")

# RF error rate
(1-mean(y.test==rf.pred.class))*100


# Confusion matrices
lasso_confusion <- table(Predicted=y.pred.class, Actual=y.test)
rf_confusion <- table(Predicted=rf.pred.class, Actual=y.test)

lasso_confusion <- as.data.frame(as.table(lasso_confusion))
rf_confusion <- as.data.frame(as.table(rf_confusion))

# Table 3
lasso_cf_gt <- lasso_confusion %>% gt() %>%
  tab_header("Lasso Confusion Matrix") %>%
  cols_label(Predicted="Predicted Group",
             Actual="Actual Group",
             Freq="Count") %>%
  fmt_number(columns=c(Freq),
            decimals=0)

# Table 4
rf_cf_gt <- rf_confusion %>% gt() %>%
  tab_header("Random Forest Confusion Matrix") %>%
  cols_label(Predicted="Predicted Group",
             Actual="Actual Group",
             Freq="Count") %>%
  fmt_number(columns=c(Freq),
             decimals=0)

#gtsave(data=rf_cf_gt, filename="rf_cmatrix.png")
#gtsave(data=lasso_cf_gt, filename="lasso_cmatrix.png")

# RF AUC and ROC
rf_roc <- roc(y.test, c(rf_preds))
rf_auc <- auc(rf_roc)

rf_auc

# Plot LASSO ROC curve, Figure 2C
plot(lasso_roc, col = "blue", font.main=1, cex.main=2, cex.axis=1.2, cex.lab=1.8, main = "ROC Curves", lwd = 2)
plot(rf_roc, col = "red", add = TRUE, cex.axis=1.2, lwd = 2, cex.main=2, cex.lab=1.8)
legend("bottomright", cex=.6, legend = c("LASSO", "Random Forest"), col = c("blue", "red"), lwd = 2)

fig2_d <- recordPlot()

# Not used, fig2_c not in analysis
figure2 <- plot_grid(fig2_a, fig2_b, fig2_c, fig2_d, rel_widths=c(2, 1), rel_heights=c(1.5, 1), scale=c(.9, 1, 1, .9),  nrow=2, ncol=2, labels=c("A", "B", "C", "D"), label_size = 12)
figure2

figure2_top <- plot_grid(fig2_a, labels=c("A"), label_size = 15)
figure2_bottom <- plot_grid(fig2_b, fig2_d, nrow=1, scale=c(.9, .8), align="v", rel_widths=c(2, 1),labels=c("B", "C"), label_size=15)

# Figure 2
plot_grid(figure2_top, figure2_bottom, nrow=2, scale=c(.8, 1), rel_heights=c(1.5, 1))

# Method 3 - Sample Clustering, not included in analysis (DE, LASSO, RF, Clustering)
cluster_data <- filter(lasso_data, Group %in% c("DR"))
x_c <- t(ex_qn)
x_c <- x_c[row.names(x_c) %in% row.names(lasso_data),]
x_cluster <- x_c[rownames(x_c) %in% rownames(cluster_data),]


euclidean <- dist(x_cluster, method="euclidean")
hcl.ward <- hclust(euclidean, method="ward.D2")
htree.ward <- as.dendrogram(hcl.ward)

h <- cutree(hcl.ward, h=100)
table(h)
plot(htree.ward, main="Ward Hierarchical Clustering",
     cex.main=1, cex.lab=.8, cex.axis=.9)

htree.ward %>% set("branches_k_color", k=2) %>%
  set("labels_cex", .8) %>% 
  plot(main="k=2 clusters", font.main=1)

fig3_a <- recordPlot()

set.seed(2024)
z1 <- fviz_nbclust(x_cluster, FUN=kmeans, method="silhouette")
z1


# Optimal K 2
set.seed(2024)
kmeans.2 <- kmeans(x_cluster, centers=2)
table(kmeans.2$cluster)

z2 <- fviz_cluster(kmeans.2, data=x_cluster) +
  z1_theme
z2

plot_grid(z1, z2, nrow=1)

z1_theme <- z1$theme

hcl.2 <- cutree(hcl.ward, k=2)

table(hcl.2)

means <- rowMeans(ex_qn)
sd <- apply(ex_qn, 1, sd)
coef.var <- sd/means

ex_qn.filter <- ex_qn[coef.var > 0.10,]

z3 <- heatplot(ex_qn.filter)

fig3_d <- recordPlot()

figure3_top <- plot_grid(fig3_a, labels=c("A"), label_size=12)
figure3_bottom <- plot_grid(z1, z2, labels=c("B", "C"), label_size=12, nrow=1)
figure3 <- plot_grid(figure3_top, figure3_bottom, nrow=2)


figure3

# Method 4 - Dimension Reduction
set.seed(2024)
tsne.dr <- Rtsne(x_c)
pca.dr <- prcomp(x_c)

summary(pca.dr)

umap.dr <- umap(x_c)

# PCA scree, Fig 3A
screeplot(pca.dr, oma=c(0, 0, 0, 0), cex.main=1, cex.lab=.8, cex.axis=.9, 
          type="lines", main="PCA Scree Plot", font.main=1) +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))
fig4_a <- recordPlot()

# Fig 3B
plot(tsne.dr$Y, col=as.factor(lasso_data$Group),
     oma=c(0, 0, 0, 0),
     cex.main=1, cex.lab=.8, cex.axis=.9,
     font.main=1,
     pch=15,
     xlab="First Component",
     ylab="Second Component",
     main="T-SNE") +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

fig4_b <- recordPlot()

# Fig 3C
plot(pca.dr$x[,1:2], col=as.factor(lasso_data$Group),
     cex.main=1, cex.lab=.8, cex.axis=.9,
     font.main=1,
     oma=c(0, 0, 0, 0),
     pch=15,
     xlab="PC1 (13.45%)",
     ylab="PC2 (12.97%)",
     main="PCA") +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

fig4_c <- recordPlot()

# Fig 3D
plot(umap.dr$layout, col=as.factor(lasso_data$Group),
     cex.main=1, cex.lab=.8, cex.axis=.9,
     font.main=1,
     pch=15,
     oma=c(0, 0, 0, 0),
     xlab="First Component",
     ylab="Second Component",
     main="UMAP") +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

fig4_d <- recordPlot()

#Figure 3
figure4 <- plot_grid(fig4_a, fig4_b, fig4_c, fig4_d,
                     labels=c("A", "B", "C", "D"),
                     label_size = 12,
                     nrow=2, ncol=2,
                     scale=c(.8,.8,.8,.8)
                     )

#save.image(file="pubhbio-final-plots.RData")

figure4
