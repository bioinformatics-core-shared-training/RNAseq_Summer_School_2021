library(DESeq2)
library(tidyverse)

# Load in the data

txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv",
                       col_types = "cccc")

all(colnames(txi$counts)==sampleinfo$SampleName)

# Create DESeq2 object

## Create the model

simple.model <- as.formula(~ Status)

model.matrix(simple.model, data = sampleinfo)

sampleinfo <- mutate(sampleinfo, Status = fct_relevel(Status, "Uninfected"))

model.matrix(simple.model, data = sampleinfo)

### Build the DESeq2 object

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colDat = sampleinfo,
                                       design = simple.model)

## Filter out unexpressed

keep <- rowSums(counts(ddsObj.raw)) > 5
table(keep)

ddsObj.filt <- ddsObj.raw[keep, ]

# The DESeq2 workflow

## estimate Size Factors

ddsObj <- estimateSizeFactors(ddsObj.filt)

normalizationFactors(ddsObj.filt)
head(normalizationFactors(ddsObj))

### MA plots for SRR7657882

logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)
limma::plotMA(logcounts, array = 5, ylim=c(-5, 5))
abline(h=0, col="red")


logcounts <- log2(counts(ddsObj, normalized = TRUE) + 1)
limma::plotMA(logcounts, array = 5, ylim=c(-5, 5))
abline(h=0, col="red")

## Estimate Dispersion

ddsObj <- estimateDispersions(ddsObj)

plotDispEsts(ddsObj)

## Fit negative binomial model and run Wald test

ddsObj <- nbinomWaldTest(ddsObj)

# The DESeq command

ddsObj <- DESeq(ddsObj.filt)

## Generate a results table

results.simple <- results(ddsObj, alpha = 0.05)
results.simple

### How many genes are differentially expressed at padj < 0.05

sum(results.simple$padj < 0.05, na.rm = TRUE)

### How many are upregulated?

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = TRUE)


### How many are downregulated?

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = TRUE)

### Exercise 1

additive.model <- as.formula(~ TimePoint + Status)
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj <- DESeq(ddsObj.filt)

results.additive <- results(ddsObj, alpha=0.05)
results.additive

sum(results.additive$padj < 0.05, na.rm = TRUE)

model.matrix(additive.model, data = sampleinfo)

resultsNames(ddsObj)

results.InfectedvUninfected <- results.additive
rm(results.additive)

## Pull top 100 genes by adjusted p-values

topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>%
  rownames_to_column("GeneID") %>% 
  top_n(100, wt = -padj)

# Exercise 2

resultsNames(ddsObj)

results.d33vd11 <- results(ddsObj, 
                           name = "TimePoint_d33_vs_d11",
                           alpha = 0.05)
results.d33vd11

sum(results.d33vd11$padj < 0.05, na.rm = TRUE)

## Should we use the interaction model

vstcounts <- vst(ddsObj.raw, blind = TRUE)
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))

# Comparing two models

ddsObj.LRT <- DESeq(ddsObj, test = "LRT", reduced = simple.model)

results.Additive_v_Simple <- results(ddsObj.LRT, alpha=0.05)
results.Additive_v_Simple

sum(results.Additive_v_Simple$padj < 0.05, na.rm = TRUE)

# Running the interaction model


interaction.model <- as.formula(~ TimePoint * Status)
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)
ddsObj.filt <- ddsObj.raw[keep, ]

ddsObj.interaction <- DESeq(ddsObj.filt)

# test interaction model v addititive

ddsObj.LRT <- DESeq(ddsObj.interaction, 
                    test = "LRT",
                    reduced = additive.model)

results.Interaction_v_Additive <- results(ddsObj.LRT)
table(results.Interaction_v_Additive$padj < 0.05)


ddsObj.LRT <- DESeq(ddsObj.interaction, 
                    test = "LRT",
                    reduced = simple.model)

results.Interaction_v_Simple <- results(ddsObj.LRT)
table(results.Interaction_v_Simple$padj < 0.05)

# Extract specific contrasts from an interactions

resultsNames(ddsObj.interaction)

# At day 11 Infected v Uninfected
results.interaction.11 <- results(ddsObj.interaction,
                                  name = "Status_Infected_vs_Uninfected",
                                  alpha = 0.05)

results.interaction.11

# A day 33 Infected v Uninfected
results.interaction.33 <- results(ddsObj.interaction,
                                  contrast = list(c("Status_Infected_vs_Uninfected",
                                                     "TimePointd33.StatusInfected")),
                                  alpha = 0.05)

sum(results.interaction.11$padj < 0.05, na.rm = TRUE)
sum(results.interaction.33$padj < 0.05, na.rm = TRUE)

## Extract d33 v d11 for Uninfected

results.d33_v_d11_uninfected <- results(ddsObj.interaction,
                                        name = "TimePoint_d33_vs_d11",
                                        alpha = 0.05)

table(results.d33_v_d11_uninfected$padj < 0.05)

## Extract d33 v d11 for Infected

results.d33_v_d11_infected <- results(ddsObj.interaction,
                                        contrast = list(c("TimePoint_d33_vs_d11",
                                                          "TimePointd33.StatusInfected")),
                                        alpha = 0.05)

table(results.d33_v_d11_infected$padj < 0.05)

saveRDS(ddsObj.interaction, "results/DESeqDataSet.interaction.rds")
saveRDS(results.interaction.11, "results/DESeqResults.interaction_d11.rds")
saveRDS(results.interaction.33, "results/DESeqResults.interaction_d33.rds")
