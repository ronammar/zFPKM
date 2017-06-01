## ---- reproducible_settings, echo=FALSE----------------------------------
rm(list=ls())
set.seed(8675309)
options(stringsAsFactors=FALSE)
library(printr)

knitr::opts_chunk$set(message=FALSE, fig.align="center", fig.width=8, fig.height=6)

## ------------------------------------------------------------------------
library(dplyr)
library(GEOquery)
library(stringr)
library(tidyr)

getSpecificGEOSupp <- function(url) {
  temp <- tempfile()
  download.file(url, temp)
  out <- read.csv(gzfile(temp), row.names=1, check.names=FALSE)
  out <- select(out, -MGI_Symbol)
  return(out)
}

gse94802_fpkm <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_normalized_FPKM.csv.gz"
gse94802_counts <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE94nnn/GSE94802/suppl/GSE94802_Minkina_etal_raw_counts.csv.gz"

fpkm <- getSpecificGEOSupp(gse94802_fpkm)
counts <- getSpecificGEOSupp(gse94802_counts)

esetlist <- getGEO("gse94802")
doe <- pData(esetlist[[1]]) %>%
  mutate(condition=ifelse(str_detect(title, regex("control", ignore_case=TRUE)), 
                          "control", "mutant"),
         sample_id=(str_match(title, "rep\\d_(.+)")[, 2]))
rownames(doe) <- doe$sample_id

rm(esetlist, gse94802_fpkm, gse94802_counts)  # clear namespace

## ------------------------------------------------------------------------
count(doe, condition) %>% as.data.frame()

## ------------------------------------------------------------------------
library(zFPKM)
zfpkm <- zFPKMTransformDF(fpkm)

## ------------------------------------------------------------------------
activeGenes <- zfpkm %>%
  mutate(gene=rownames(zfpkm)) %>%
  gather(sample_id, zfpkm, -gene) %>%
  left_join(select(doe, sample_id, condition), by="sample_id") %>%
  group_by(gene, condition) %>%
  summarize(median_zfpkm=median(zfpkm)) %>%
  ungroup() %>%
  mutate(active=(median_zfpkm > -3)) %>%
  filter(active) %>%
  select(gene) %>%
  distinct()

counts <- counts[activeGenes$gene, ]

## ------------------------------------------------------------------------
library(limma)
library(edgeR)

# Generate normalized log2CPM from counts AFTER we filter for protein-coding
# genes that are detectably expressed.
dge <- DGEList(counts=counts)
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + condition, data=doe)
vq <- voomWithQualityWeights(dge, design, plot=TRUE)
fit <- lmFit(vq, design)
contrastMatrix <- makeContrasts(conditioncontrol - conditionmutant, levels=design)
fit <- contrasts.fit(fit, contrastMatrix)
fit <- eBayes(fit, robust=TRUE)
deGenes <- topTable(fit, number=Inf)

## ------------------------------------------------------------------------
sessionInfo()

