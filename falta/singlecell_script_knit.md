---
title: "single_nuclei_knit"
author: "Julia Faura"
date: "30/05/2022"
output: html_document
---

#### Install and load packages

```{r}
#install.packages("Seurat", lib="/home/julia/R/x86_64-pc-linux-gnu-library/4.1", force=TRUE)
#BiocManager::install("scater")
#install.packages("ggvenn")
#install.packages("dplyr")
library(Seurat)
library(scater)
library(ggvenn)
library(dplyr)
```

### Load data

```{r}
setwd("/home/julia/single_cell_RNAseq/Human nuclei_10x trial_#3 (EMC_2013-098)/cellranger_output/outs/filtered_feature_bc_matrix")
data <- Read10X("/home/julia/single_cell_RNAseq/Human nuclei_10x trial_#3 (EMC_2013-098)/cellranger_output/outs/filtered_feature_bc_matrix")
```

### Clean data

What percentage of counts is equal to 0?
The percentage of zero counts is equal to the total number of zeros divided by the total number of counts

```{r}
total.zeros <- sum(data==0) 
total.counts <- nrow(data) * ncol(data)
(total.zeros/total.counts) * 100
```

How many genes are expressed in each cell?

- A gene is considered expressed when count > 0
- Cells are the columns in the count matrix hence MARGIN=2
- sum(x > 0) calculates number of genes with non-zero counts

```{r}
cellCounts <- apply(data,2,function(x) sum(x > 0))
cellCounts[1:5]
dim(data)
```

We will clean the huge matrix to make it easier to process:

Cells with less than 200 expressed genes are removed.

```{r}
data <- data[,cellCounts>=200]
dim(data)
```

For each gene: in how many cells is it expressed?

- A gene is considered expressed if count > 0
- Genes are the rows in the count matrix hence MARGIN=1
- sum(x > 0) calculates number of cells with non-zero counts

Genes that are expressed in less than 3 cells are removed.

```{r}
geneCounts <- apply(data,1,function(x) sum(x > 0))
geneCounts[1:5]
data <- data[geneCounts>=3,]
dim(data)
```

```{r}
rm(total.counts)
rm(total.zeros)
rm(cellCounts)
rm(geneCounts)
saveRDS(data,file="data_clean.rds")
```

### Convert to SingleCellExperiment object

Scater needs the data in a Single Cell Experiment object, where the counts are stored in the assays slot.

```{r}
sce <- SingleCellExperiment(assays=list(counts=data)) 
rm(data)
saveRDS(sce,file="sce_first.rds")

```

### Quality control of cells user scater

Identify mitochondrial transcripts
Gene names are used as row names in the count matrix. To find the mitochondrial genes we need to look for row names that start with mt-

grepl() checks if words contain a pattern and returns booleans:

- ^ in the pattern argument means starts with. Check out this tutorial on regular expressions in R.
- ignore.case=TRUE: in mouse names of mitochondrial genes start with mt, in human with MT

The function returns a vector with Booleans:

- TRUE if the name of a row (=gene) starts with mt or MT
- FALSE if the name of a row (=gene) does not start with mt or MT

Mitochondrial genes can be used to identify unhealthy cells. When cells are exposed to stress, their cell membrane becomes leaky. When this happens:

- transcripts go though the tears in the membrane
- large mitochondria stay inside the cell
- mitochondria are the last organelles to degrade

Enrichment of mitochondrial transcripts is therefore a clear indication of cell stress.

```{r}
is.mito <- grepl("^mt-",rownames(sce),ignore.case=TRUE) 

```

#### Calculate QC metrics per cell

Use addPerCellQC() from scater.

Arguments of this function:

- input: count matrix
- subsets: a named list (objects in the list have labels) with booleans or names of genes that can be used as controls like mitochondrial genes

```{r}
sce <- addPerCellQC(sce,subsets=list(Mt=is.mito))
names(colData(sce))[c(1,2,5)] <- c("nUMI","nGene","mito")
```

Six quality metrics are calculated for every cell (=row):

- detected: number of expressed genes (count > 0) (nGene)
- sum: number of transcripts (total count or library size) (nUMI)
- subsets_Mt_percent: percentage counts of mitochondrial genes (mito)

We are going to rename these columns:

```{r}
colData(sce)$nUMI.out.low <- isOutlier(colData(sce)$nUMI,nmads=3,type="lower",log=TRUE) 
colData(sce)$nUMI.out.high <- isOutlier(colData(sce)$nUMI,nmads=3,type="higher",log=TRUE)
sum(colData(sce)$nUMI.out.low | colData(sce)$nUMI.out.high)
```

### Detect outliers cells using scater

Outlier cells are:

- cells with low number of genes (or UMIs) expressed
- cells with high percentage of mitochondrial transcripts

isOutlier() identifies outlier cells for each of these metrics:

- cells with extreme log-library size
- cells with extreme percentage of mitochondrial genes

Extreme is defined by a certain number of MADs from the median.

This is a trial and error process:

- you choose a number of MADs
- create plots
- inspect the plots
- change the threshold for number of MADs
- repeat until it looks ok
- Use isOutlier() from scater.

Arguments of this function:

- input: values for the metric you want to use
- nmads: threshold for number of MADs. The lower you set this, the more outliers will be found.
- type: find outliers at both tails (default) or only at lower end (= too few) or only at higher end (= too many)
- log: do you want to take log10 of metric values before computing MAD?

#### UMI counts per cell

```{r}
colData(sce)$nUMI.out.low <- isOutlier(colData(sce)$nUMI,nmads=3,type="lower",log=TRUE) 
colData(sce)$nUMI.out.high <- isOutlier(colData(sce)$nUMI,nmads=3,type="higher",log=TRUE)
sum(colData(sce)$nUMI.out.low | colData(sce)$nUMI.out.high)
```

Histogram:

```{r}
metaData <- as.data.frame(colData(sce))
ggplot(metaData,aes(nUMI)) + 
  geom_histogram(binwidth=100) + 
  xlab("Count depth (total UMI count)") +
  ylab("Frequency") +
  ggtitle("Histogram of total UMI count per cell") + 
  theme_bw()
```


#### Number of expressed genes per cell

```{r}
colData(sce)$nGene.out.low <- isOutlier(colData(sce)$nGene,nmads=3,type="lower",log=TRUE) 
colData(sce)$nGene.out.high <- isOutlier(colData(sce)$nGene,nmads=3,type="higher",log=TRUE) 
sum(colData(sce)$nGene.out.low | colData(sce)$nGene.out.high)
```

Histogram:

```{r}
cut.nGene <- 2^(median(log2(metaData$nGene))-3*mad(log2(metaData$nGene),na.rm=TRUE))
ggplot(metaData,aes(nGene)) + 
  geom_histogram(binwidth=20) +
  xlab("Number of Genes") +
  ylab("Frequency") +
  ggtitle("Histogram of number of genes per cell") + 
  geom_vline(xintercept=cut.nGene,color="red") +
  theme_bw()
```


#### Mitochondrial count percentages

```{r}
colData(sce)$mito.out.high <- isOutlier(colData(sce)$mito,nmads=3,type="higher") 
sum(colData(sce)$mito.out.high)
```

Histogram:

```{r}
cut.mito <- median(metaData$mito) + 3*mad(metaData$mito,na.rm=TRUE)
ggplot(metaData,aes(mito)) + 
  geom_histogram(binwidth=0.1) +
  xlab("% Mitochondrial counts") +
  ylab("Frequency") +
  ggtitle("Histogram of % mitogenes per cell") +
  geom_vline(xintercept=cut.mito,color="red") +
  theme_bw()
```

#### Violin plots

```{r}
ggplot(metaData,aes("",nUMI)) +
  geom_jitter(height=0,width=0.3,aes(color=nUMI.out.low)) +
   geom_violin(fill="gray80",alpha=0.5) +
  scale_color_manual(values=c("#00BFC4","#F8766D")) + 
  ggtitle("Total UMI counts per cell") +
  theme_classic()
```

#### Violin plots after filtering

```{r}
non.out <- !(metaData$nUMI.out.low | metaData$nUMI.out.high | metaData$nGene.out.low | metaData$nGene.out.high | metaData$mito.out.high)
metaData.filtered <- metaData[non.out,]
dim(metaData.filtered)
```


#### Venn diagram

So we have detected outliers cells with:

- few expressed genes
- a lot of mitochondrial transcripts

We are going to check if these two groups of cells overlap by making a Venn diagram.

```{r}
ggvenn(list(nGene=rownames(metaData[metaData$nGene.out.low,]),
            mito=rownames(metaData[metaData$mito.out.high,])),
       fill_color=c("green","orange"),show_percentage=FALSE) 
```

#### Remove outliers

Cells with quality values < threshold on number of MADs thresholds are considered outliers. They are removed under the assumption that they correspond to low-quality cells.

```{r}
keep <- !(metaData$mito.out.high) 
sce <- sce[,keep] 
dim(sce)
```

```{r}
rm(metaData.filtered)
rm(metaData)
rm(keep)
rm(non.out)
saveRDS(sce,file="sceUni.rds")
```

### Multivariate outlier detection by PCA

We will identify outliers using a principal component analysis (PCA) since it allows to include multiple QC metrics. Hence it performs multivariate outlier detection.

First we define which metrics we want to take into account.

PCA will identify cells with substantially different QC metrics as outliers. These outliers are marked and can be deleted from the dataset.

PCA is done using runColDataPCA() from the scater package:

- object input is a SCE object
- variables which quality metrics to use for the PCA?
- outliers by default it will not identify outlier cells based on the PCA, if you do want to identify outliers set to TRUE

```{r}
varsToUse <- c("nUMI","nGene","mito")
sce <- runColDataPCA(sce,variables=varsToUse,
                     outliers=TRUE)
```

This adds a column called outlier to the colData, which is accessible as colData(sce)$outlier. It contains Booleans:

- TRUE for outlier cells
- FALSE for ok cells
- Number of detected outliers:

```{r}
sum(colData(sce)$outlier)

```

#### Draw PCA plot and color outliers

The plot can be created using plotReducedDim() of the scater package. Arguments are:

- object input is a sce object
- dimred more specifically the results of the PCA calculations which are stored in the PCA_coldata object of the list generated by reducedDims(sce)
- colour_by the name of the column in coldata that contains the values you want to use to color the points, in this case: is the cell (point) an outlier or not?

```{r}
reducedDims(sce)
plotReducedDim(sce,dimred="PCA_coldata",colour_by="outlier")

```

#### Violin plots

```{r}
metaData <- as.data.frame(colData(sce))
ggplot(metaData,aes("",nUMI)) +
  geom_jitter(height=0,width=0.3,aes(color=outlier)) +
   geom_violin(fill="gray80",alpha=0.5) +
  scale_color_manual(values=c("#00BFC4","#F8766D")) + 
  ggtitle("Total UMI counts per cell") +
  theme_classic()
```

```{r}
ggplot(metaData,aes("",nGene)) +
  geom_jitter(height=0,width=0.3,aes(color=outlier)) +
   geom_violin(fill="gray80",alpha=0.5) +
  scale_color_manual(values=c("#00BFC4","#F8766D")) + 
  ggtitle("Number of genes per cell") +
  theme_classic()
```

#### Violin plots after filtering

```{r}
non.out <- !(metaData$outlier)
metaData.filtered <- metaData[non.out,]
dim(metaData.filtered)
```

```{r}
ggplot(metaData.filtered,aes("",nUMI)) +
  geom_jitter(height=0,width=0.3,aes(color=outlier)) +
   geom_violin(fill="gray80",alpha=0.5) +
  scale_color_manual(values=c("#00BFC4","#F8766D")) + 
  ggtitle("Total UMI counts per cell") +
  theme_classic()
```

```{r}
ggplot(metaData.filtered,aes("",nGene)) +
  geom_jitter(height=0,width=0.3,aes(color=outlier)) +
   geom_violin(fill="gray80",alpha=0.5) +
  scale_color_manual(values=c("#00BFC4","#F8766D")) + 
  ggtitle("Number of genes per cell") +
  theme_classic()
```

```{r}
keep <- !(metaData$outlier) 
sce <- sce[,keep] 
dim(sce)
```

```{r}
saveRDS(sce,file="sce.rds")
rm(list=ls())
```

```{r}
sce <- readRDS("sce.rds")
counts <- counts(sce)
metaData <- colData(sce) 
rm(sce)
seuratObj <- CreateSeuratObject(counts=counts)
seuratObj <- AddMetaData(seuratObj,as.data.frame(metaData)) 
rm(counts)
```

```{r}
seuratObj@assays
```

```{r}
seuratObj <- NormalizeData(object=seuratObj)
seuratObj@assays$RNA@data[1:5,1:10]
gene3.raw <- seuratObj@assays$RNA@counts[3,]
gene3.norm <- seuratObj@assays$RNA@data[3,]
sum(gene3.raw!=gene3.norm)
```

```{r}
ggplot(seuratObj@meta.data,aes(nUMI)) + geom_histogram(bins=60)

```

```{r}
counts.norm <- as.matrix(seuratObj@assays$RNA@data)
nUMILN <- data.frame(nUMI=colSums(counts.norm))
ggplot(nUMILN,aes(nUMI)) + geom_histogram(bins=60)
```

```{r}
seuratObj <- FindVariableFeatures(object=seuratObj)
length(VariableFeatures(seuratObj)) 
head(VariableFeatures(seuratObj)) 
head(HVFInfo(seuratObj))
```

```{r}
VariableFeaturePlot(object=seuratObj)
```

```{r}
seuratObj <- ScaleData(object=seuratObj)
```

```{r}
seuratObj@assays$RNA@scale.data[1:5,1:4]
seuratObj <- RunPCA(object=seuratObj,
                    features=VariableFeatures(seuratObj),
                    nfeatures.print=10)
```

```{r}
seuratObj@reductions$pca@cell.embeddings[1:5,1:5]
seuratObj@reductions$pca@feature.loadings[1:5,1:5]
DimPlot(object=seuratObj,reduction="pca")
```

```{r}
DimHeatmap(seuratObj,dims=1:10,cells=500,balanced=TRUE)
ElbowPlot(seuratObj)
```

```{r}
seuratObj <- FindNeighbors(object=seuratObj,dims=1:7)
names(seuratObj@graphs)
```

```{r}
resToUse <- 0.5 
seuratObj <- FindClusters(object=seuratObj,
                          resolution=resToUse,graph.name="RNA_snn") 
head(seuratObj@meta.data) 

```

```{r}
seuratObj <- RunTSNE(object=seuratObj,dims=1:8,
                     check_duplicates=FALSE)
DimPlot(seuratObj,reduction="tsne",label=TRUE,label.size=8,pt.size=2) + NoLegend()

```

```{r}
seuratObj <- RunUMAP(seuratObj,dims=1:7)
DimPlot(object=seuratObj)
pdf("UMAP_shortread.pdf")
DimPlot(object=seuratObj)
dev.off()
```

Astrocytes

```{r}
pdf("AQP4.pdf")
FeaturePlot(seuratObj, features = "AQP4", min.cutoff = "q9", label = TRUE, repel = TRUE)
dev.off()
FeaturePlot(seuratObj, features = "GFAP", min.cutoff = "q9", label = TRUE, repel = TRUE)
```

Excitatory neurons

```{r}
pdf("NRGN.pdf")
FeaturePlot(seuratObj, features = "NRGN", min.cutoff = "q9", label = TRUE, repel = TRUE)
dev.off()
```

Inhibitory neurons

```{r}
pdf("GAD1.pdf")
FeaturePlot(seuratObj, features = "GAD1", min.cutoff = "q9", label = TRUE, repel = TRUE)
dev.off()
```

Neurons

```{r}
FeaturePlot(seuratObj, features = "RBFOX3", min.cutoff = "q9", label = TRUE, repel = TRUE)
```

Oligodendrocytes

```{r}
pdf("MBP.pdf")
FeaturePlot(seuratObj, features = "MBP", min.cutoff = "q9",  label = TRUE, repel = TRUE)
dev.off()
```

Microglia

```{r}
pdf("CSF1R.pdf")
FeaturePlot(seuratObj, features = "CSF1R", min.cutoff = "q9",  label = TRUE, repel = TRUE)
dev.off()
FeaturePlot(seuratObj, features = "CD74", min.cutoff = "q9",  label = TRUE, repel = TRUE)
```

OPC

```{r}
pdf("VCAN.pdf")
FeaturePlot(seuratObj, features = "VCAN", min.cutoff = "q9",  label = TRUE, repel = TRUE)
dev.off()
```

Endothelial cells

```{r}
FeaturePlot(seuratObj, features = "FLT1", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "PECAM1", min.cutoff = "q9", label = TRUE, repel = TRUE)
```

```{r}
FeaturePlot(seuratObj, features = "GRIA2", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "NPY", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "PTN", min.cutoff = "q9", label = TRUE, repel = TRUE)

FeaturePlot(seuratObj, features = "PCP4", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "NEFM", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "CRYM", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "RRM2", min.cutoff = "q9", label = TRUE, repel = TRUE)
FeaturePlot(seuratObj, features = "DEPTOR", min.cutoff = "q9", label = TRUE, repel = TRUE)

```

