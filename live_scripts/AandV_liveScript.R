library(AnnotationDbi)
library(AnnotationHub)
library(ensembldb)
library(DESeq2)
library(tidyverse)

ddsObj.interaction <- readRDS("RNAseq/RObjects/DESeqDataSet.interaction.rds")
results.interaction.11 <- readRDS("RNAseq/RObjects/DESeqResults.interaction_d11.rds")
results.interaction.33 <- readRDS("RNAseq/RObjects/DESeqResults.interaction_d33.rds")

# Annotation

ah <- AnnotationHub()
ah
ah[1]

MouseEnsDb <- query(ah, c("EnsDb", "Mus musculus", "102"))[[1]]
MouseEnsDb

annotations <- genes(MouseEnsDb, return.type = "data.frame")

colnames(annotations)

annot <- annotations %>%
  select(gene_id, gene_name, entrezid) %>%
  filter(gene_id %in% rownames(results.interaction.11))
head(annot)

length(unique(annot$entrezid))

ensemblAnnot <- readRDS("RNAseq/RObjects/Ensembl_annotations.rds")
colnames(ensemblAnnot)

annot.interaction.11 <- as.data.frame(results.interaction.11) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)

# Visualisation

hist(annot.interaction.11$pvalue)

ddsShrink.11 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.11,
                          type = "ashr")

ShrinkTab.11 <-as.data.frame(ddsShrink.11) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)

# MA plots

plotMA(results.interaction.11, alpha = 0.05)
plotMA(ddsShrink.11, alpha = 0.05)

# Volcano plot

volcanoTab.11 <- ShrinkTab.11 %>%
  mutate(`-log10(pvalue)` = -log10(pvalue))

ggplot(volcanoTab.11, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = pvalue < 0.05), size = 1) +
  geom_text(data=~top_n(.x, 1, wt=-FDR), aes(label=Symbol))

# Exercise

ddsShrink.33 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.33,
                          type = "ashr")

ShrinkTab.33 <- as.data.frame(ddsShrink.33) %>%
  rownames_to_column("GeneID") %>%
  left_join(ensemblAnnot, "GeneID") %>%
  rename(logFC=log2FoldChange, FDR=padj)

volcanoTab.33 <- ShrinkTab.33 %>%
  mutate(`-log10(pvalue)` = -log10(pvalue))

ggplot(volcanoTab.33, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = pvalue < 0.05), size = 1) +
  geom_text(data=~top_n(.x, 1, wt=-FDR), aes(label=Symbol))

# Venn

library(ggvenn)

vennDat <- tibble(Geneid=rownames(results.interaction.11)) %>% 
  mutate(Upregulated_11 = results.interaction.11$padj < 0.05 & !is.na(results.interaction.11$padj) & results.interaction.11$log2FoldChange > 0) %>% 
  mutate(Downregulated_11 = results.interaction.11$padj < 0.05 & !is.na(results.interaction.11$padj) & results.interaction.11$log2FoldChange < 0) %>%
  mutate(Upregulated_33 = results.interaction.33$padj < 0.05 & !is.na(results.interaction.33$padj) & results.interaction.33$log2FoldChange > 0) %>%
  mutate(Downregulated_33 = results.interaction.33$padj < 0.05 & !is.na(results.interaction.33$padj) & results.interaction.33$log2FoldChange < 0) 

ggvenn(vennDat, set_name_size = 4)

# heatmaps

library(ComplexHeatmap)
library(circlize)

sigGenes <- ShrinkTab.11 %>%
  top_n(300, wt=-FDR) %>%
  pull("GeneID")

plotDat <- vst(ddsObj.interaction)[sigGenes,] %>%
  assay()


z.mat <- t(scale(t(plotDat), center = TRUE, scale = TRUE))

myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2,0,2), myPalette)

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE)

ha1 <- HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status","TimePoint")])

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        split = 3,
        rect_gp = gpar(col = "lightgrey", lwd = 0.3),
        top_annotation = ha1)

ha1 = HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")],
                        col = list(Status = c("Uninfected" = "hotpink", "Infected" = "purple"), 
                                   TimePoint = c("d11" = "lightblue", "d33" = "darkblue")))

Heatmap(z.mat, name = "z-score",
        col = myRamp,            
        show_row_name = FALSE,
        split=3,
        rect_gp = gpar(col = "lightgrey", lwd=0.3),
        top_annotation = ha1)