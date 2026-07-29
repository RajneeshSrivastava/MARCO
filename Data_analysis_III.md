## Project: *Staphylococcus aureus Biofilm Infection Impairs Efferocytosis in Wound-Macrophages via Induction of a MARCO-high MerTK-low Subset*

##### Upload libraries
```
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(hdf5r)
library(future)
library(sctransform)

plan("multisession", workers = 4)
options(future.globals.maxSize = 15000 * 1024^2)
options(stringsAsFactors = FALSE)
```
#### Single-cell data analysis using Seurat
##### Upload processed single-cell RNA-seq data obtained from GSE165816 
##### Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165816
###### Create seurat objects
```
DFU.dir = "path/GSE165816_dataset/"
DFU.list=list("G6", "G9", "G33", "G34", "G39", "G7", "G8", "G49", "G42", "G45", "G23", "G4", "G2", "G15")

for (file in DFU.list)
              {DFU_data <- read.csv(file = paste0(DFU.dir, file,".csv.gz"), header = TRUE, row.names = 1)
               DFU_obj <- CreateSeuratObject(counts = DFU_data,
                                              min.cells=3,
                                              min.features = 200,
                                              project = file)
                                              assign(file, DFU_obj)
               }

sample.list=list(G6,G9,G33,G34,G39,G7,G8,G49,G42,G45,G23,G4,G2,G15)
```
###### Quality Filtering and normalization using SCTransform
```
for (i in 1:length(sample.list)) {
  sample.list[[i]][["percent.mt"]] <- PercentageFeatureSet(sample.list[[i]], pattern = "^MT-")
  sample.list[[i]] <- subset(sample.list[[i]], 
                             subset = nFeature_RNA > 200 & 
                               nFeature_RNA < 5000 & percent.mt < 15 & 
                               nCount_RNA < 25000 & nCount_RNA > 2000)
#SCT transformation
  sample.list[[i]] <- SCTransform(sample.list[[i]],
                                vars.to.regress = "percent.mt", 
                                verbose = FALSE)
                                 }
```
###### Data integration
```
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- lapply(X = sample.list, FUN = function(x) {
                                x <- RunPCA(x, features = sample.features)
                                                         })
sample.list <- PrepSCTIntegration(object.list = sample.list, 
                                  anchor.features = sample.features,
                                  verbose = FALSE)           
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, 
                                         normalization.method = "SCT", 
                                         anchor.features = sample.features, 
                                         reference = c(1, 2), 
                                         reduction = "rpca", 
                                         verbose = TRUE)
sample.integrated <- IntegrateData(anchorset = sample.anchors, 
                                   normalization.method = "SCT",
                                   verbose = TRUE)
```
###### Dimensionality reduction and clustering
```
sample.integrated <- RunPCA(object = sample.integrated, verbose = FALSE) 
sample.integrated <- RunUMAP(sample.integrated, dims = 1:30)
sample.integrated <- FindNeighbors(sample.integrated)

```
###### tweak-in to assign groups; HDHU and NHDHU
```
meta=read.table("metadata.txt",sep="\t", header=T) #This file includes metadata for each samples, labled with "HG" and "NHG" groups

GSM=sample.integrated@meta.data
GSM$group=1
for(i in 1:nrow(GSM)){
  for (j in 1:nrow(meta)){
    if (GSM[i,1] == meta[j,1]){
      GSM[i,7] = meta[j,2] }
  }
}

sample.integrated@meta.data=GSM
```
###### Find clusters with defined resolution
```
sample25 <- FindClusters(sample.integrated, resolution = 0.25)
sample25$group <- factor(sample25$group, levels = c("NHG", "HG"))

```
##### Figure S3A-B
```
DimPlot(sample25)
#Find all markers
sample25_markers=FindAllMarkers(sample25,only.pos = T,
                                logfc.threshold = 0.30, min.pct = 0.10, recorrect_umi=FALSE)
#write.csv(sample25_markers,"All_markers_sample25.csv")
##### Explore the table and identify cell types using signatures from published articles and PanglaoDb (https://panglaodb.se/))
```
##### Figure S3C
###### Isolate myeloid cluster and do re-clustering analysis
```
#DefaultAssay(sample25)="integrated"
M=subset(sample25,subset=seurat_clusters==5)
M_re <- RunPCA(object = M, verbose = FALSE)
M_re <- RunUMAP(M_re, dims = 1:30)
M_re <- FindNeighbors(M_re)
M_re20 <- FindClusters(M_re, resolution = 0.20)
DimPlot(M_re20)
table(M_re20$seurat_clusters,M_re20$group)
```
##### Figure S3D
###### Identify top 5 signatures for subsets
```
DefaultAssay(MK_re20)="SCT"
M_re20_markers=FindAllMarkers(M_re20,only.pos = T,logfc.threshold = 0.30, min.pct = 0.10,recorrect_umi=FALSE)
write.csv(M_re20_markers,"M_re20_markers_SCT.csv")

top5 <- M_re20_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5genes=unique(top5$gene)
DotPlot(M_re20, features = top10genes) + 
    theme(axis.text.x=element_text(size=8.5)) + 
    theme(axis.text.y=element_text(size=10.0))+
    scale_colour_gradient2(low = "grey", mid = "red", high = "#000000")
```
##### Figure S3E
###### MARCO cell abundance and expression profile in "HG" and "NHG" groups
```
FeaturePlot(M_re20, features = "MARCO",order=T,split.by="group") & scale_colour_gradientn(colours =c("lightgrey","orange","red","black"))
```
##### Figure S3F
###### MARCO cell abundance and expression profile in "HG" and "NHG" groups
```
VlnPlot(M_re20, features=c("TYRO3", "AXL", "MERTK"),stack=T,flip=T,split.by="group")
```
### Thank you
