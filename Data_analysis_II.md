#### This file contains the script for stepwise analysis of single-cell RNA-seq data, utilised for figures []

```
setwd('path_to_project')

library(Seurat)
library(SeuratWrappers)
library(ggplot2)

object_rex = 'path_to_sample/REX1'
object_sar = 'path_to_sample/SAR1'
object_control = 'path_to_sample/CON_control'

object_rex = Read10X(object_rex)
object_sar = Read10X(object_sar)
object_control = Read10X(object_control)

object_rex = CreateSeuratObject(counts = object_rex, min.cells = 0, min.features = 0, project = "REX1")
object_sar = CreateSeuratObject(counts = object_sar, min.cells = 0, min.features = 0, project = "SAR1")
object_control = CreateSeuratObject(counts = object_control, min.cells = 0, min.features = 0, project = "Control")

object_merged = merge(x = object_rex, y = c(object_sar,object_control))
object_merged$orig.ident = factor(x = object_merged$orig.ident, levels = c('Control','SAR1','REX1'))
Idents(object_merged) = 'orig.ident'

object_merged[["percent.mt"]] <- PercentageFeatureSet(object_merged, pattern = "^MT-")

VlnPlot(object_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
plot1 <- FeatureScatter(object_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.3)
plot2 <- FeatureScatter(object_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.3)
CombinePlots(plots = list(plot1, plot2))

object_merged = subset(object_merged, subset = nFeature_RNA > 200 & percent.mt < 25 & nCount_RNA < 30000 & nCount_RNA > 500)

VlnPlot(object_merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
plot1 <- FeatureScatter(object_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.3)
plot2 <- FeatureScatter(object_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",pt.size = 0.3)
CombinePlots(plots = list(plot1, plot2))

object_merged <- NormalizeData(object_merged)
object_merged <- FindVariableFeatures(object_merged)
object_merged <- ScaleData(object_merged)
object_merged <- RunPCA(object_merged)
ElbowPlot(object_merged, ndims = 40)

# Integrate samples
object_merged <- IntegrateLayers(
  object = object_merged, method = FastMNNIntegration,
  new.reduction = "integrated.mnn",
  verbose = FALSE
)
object_merged = JoinLayers(object_merged)

object_merged <- RunUMAP(object_merged, reduction = "integrated.mnn", dims = 1:25, reduction.name = "umap.mnn")
DimPlot(object_merged,split.by = 'orig.ident',ncol = 3)
DimPlot(object_merged,ncol = 3)

object_merged <- FindNeighbors(object_merged, reduction = "integrated.mnn", dims = 1:25)
object_merged <- FindClusters(object_merged, resolution = 0.1, cluster.name = "mnn_clusters_0.1")
object_merged <- FindClusters(object_merged, resolution = 0.2, cluster.name = "mnn_clusters_0.2")
object_merged <- FindClusters(object_merged, resolution = 0.3, cluster.name = "mnn_clusters_0.3")
object_merged <- FindClusters(object_merged, resolution = 0.4, cluster.name = "mnn_clusters_0.4")
object_merged <- FindClusters(object_merged, resolution = 0.5, cluster.name = "mnn_clusters_0.5")
object_merged <- FindClusters(object_merged, resolution = 0.6, cluster.name = "mnn_clusters_0.6")

plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.1',label = T)
jpeg(paste0("DimPlot (clusters 0.1).jpg"),res=600,width=3500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.2',label = T)
jpeg(paste0("DimPlot (clusters 0.2).jpg"),res=600,width=3500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.3',label = T)
jpeg(paste0("DimPlot (clusters 0.3).jpg"),res=600,width=3500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.4',label = T)
jpeg(paste0("DimPlot (clusters 0.4).jpg"),res=600,width=3500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.5',label = T)
jpeg(paste0("DimPlot (clusters 0.5).jpg"),res=600,width=3500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.6',label = T)
jpeg(paste0("DimPlot (clusters 0.6).jpg"),res=600,width=3500,height=2800);print(plot);dev.off()

plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.1',label = T,split.by = 'orig.ident',pt.size = 0.1)
jpeg(paste0("DimPlot splitted (clusters 0.1).jpg"),res=600,width=7500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.2',label = T,split.by = 'orig.ident',pt.size = 0.1)
jpeg(paste0("DimPlot splitted (clusters 0.2).jpg"),res=600,width=7500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.3',label = T,split.by = 'orig.ident',pt.size = 0.1)
jpeg(paste0("DimPlot splitted (clusters 0.3).jpg"),res=600,width=7500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.4',label = T,split.by = 'orig.ident',pt.size = 0.1)
jpeg(paste0("DimPlot splitted (clusters 0.4).jpg"),res=600,width=7500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.5',label = T,split.by = 'orig.ident',pt.size = 0.1)
jpeg(paste0("DimPlot splitted (clusters 0.5).jpg"),res=600,width=7500,height=2800);print(plot);dev.off()
plot = DimPlot(object_merged,group.by = 'mnn_clusters_0.6',label = T,split.by = 'orig.ident',pt.size = 0.1)
jpeg(paste0("DimPlot splitted (clusters 0.6).jpg"),res=600,width=7500,height=2800);print(plot);dev.off()

table(object_merged$mnn_clusters_0.1)
table(object_merged$mnn_clusters_0.2)
table(object_merged$mnn_clusters_0.3)
table(object_merged$mnn_clusters_0.4)
table(object_merged$mnn_clusters_0.5)
table(object_merged$mnn_clusters_0.6)

table(object_merged$mnn_clusters_0.1,object_merged$orig.ident)
table(object_merged$mnn_clusters_0.2,object_merged$orig.ident)
table(object_merged$mnn_clusters_0.3,object_merged$orig.ident)
table(object_merged$mnn_clusters_0.4,object_merged$orig.ident)
table(object_merged$mnn_clusters_0.5,object_merged$orig.ident)
table(object_merged$mnn_clusters_0.6,object_merged$orig.ident)


temp_0.1 =data.frame(table(object_merged$mnn_clusters_0.1,object_merged$orig.ident));colnames(temp_0.1) = c('cluster','condition','freq')
temp_0.2 =data.frame(table(object_merged$mnn_clusters_0.2,object_merged$orig.ident));colnames(temp_0.2) = c('cluster','condition','freq')
temp_0.3 =data.frame(table(object_merged$mnn_clusters_0.3,object_merged$orig.ident));colnames(temp_0.3) = c('cluster','condition','freq')
temp_0.4 =data.frame(table(object_merged$mnn_clusters_0.4,object_merged$orig.ident));colnames(temp_0.4) = c('cluster','condition','freq')
temp_0.5 =data.frame(table(object_merged$mnn_clusters_0.5,object_merged$orig.ident));colnames(temp_0.5) = c('cluster','condition','freq')
temp_0.6 =data.frame(table(object_merged$mnn_clusters_0.6,object_merged$orig.ident));colnames(temp_0.6) = c('cluster','condition','freq')

temp_0.3 = subset(temp_0.3, subset = cluster !=5)

plot = ggplot(temp_0.1, aes(fill=cluster, y=freq, x=condition)) +
  geom_bar(position="fill", stat="identity",color = 'black')+
  theme_minimal()
jpeg(paste0("Barplot (clusters 0.1).jpg"),res=600,width=3000,height=2800);print(plot);dev.off()
plot = ggplot(temp_0.2, aes(fill=cluster, y=freq, x=condition)) +
  geom_bar(position="fill", stat="identity",color = 'black')+
  theme_minimal()
jpeg(paste0("Barplot (clusters 0.2).jpg"),res=600,width=3000,height=2800);print(plot);dev.off()
plot = ggplot(temp_0.3, aes(fill=cluster, y=freq, x=condition)) +
  geom_bar(position="fill", stat="identity",color = 'black')+
  theme_minimal()
jpeg(paste0("Barplot (clusters 0.3).jpg"),res=600,width=3000,height=2800);print(plot);dev.off()
plot = ggplot(temp_0.4, aes(fill=cluster, y=freq, x=condition)) +
  geom_bar(position="fill", stat="identity",color = 'black')+
  theme_minimal()
jpeg(paste0("Barplot (clusters 0.4).jpg"),res=600,width=3000,height=2800);print(plot);dev.off()
plot = ggplot(temp_0.5, aes(fill=cluster, y=freq, x=condition)) +
  geom_bar(position="fill", stat="identity",color = 'black')+
  theme_minimal()
jpeg(paste0("Barplot (clusters 0.5).jpg"),res=600,width=3000,height=2800);print(plot);dev.off()
plot = ggplot(temp_0.6, aes(fill=cluster, y=freq, x=condition)) +
  geom_bar(position="fill", stat="identity",color = 'black')+
  theme_minimal()
jpeg(paste0("Barplot (clusters 0.6).jpg"),res=600,width=3000,height=2800);print(plot);dev.off()

Idents(object_merged) = 'mnn_clusters_0.3'

DimPlot(object_merged,label = T)|VlnPlot(object_merged, 'MARCO')

VlnPlot(object_merged, 'IL7R')

table(object_merged$orig.ident)
table(object_merged$mnn_clusters_0.3)

Idents(object_merged) = 'orig.ident'
VlnPlot(object_merged, 'MARCO')|VlnPlot(object_merged, 'MERTK')|VlnPlot(object_merged, 'CD163')

Idents(object_merged) = 'mnn_clusters_0.3'
VlnPlot(object_merged, 'MARCO', idents = c(3,4))|VlnPlot(object_merged, 'MERTK', idents = c(3,4))|VlnPlot(object_merged, 'CD163', idents = c(3,4))
FeaturePlot(object_merged, c('MARCO','MERTK','CD163'),split.by = 'orig.ident')
FeaturePlot(object_merged, c('MARCO','MERTK','CD163'))

marker = FindMarkers(object_merged, ident.1 = 4, ident.2 = 3, Features = 'CD163')
View(marker)
write.csv(marker, '4 vs 3.csv')


markers = FindAllMarkers(object_merged, features =  'MARCO')

Idents(object_merged) = 'mnn_clusters_0.3'
markers = FindAllMarkers(object_merged, features =  'MARCO')

object_merged$cluster_condition = paste0(object_merged$mnn_clusters_0.3,'_', object_merged$orig.ident)
Idents(object_merged) = 'cluster_condition'
markers = FindAllMarkers(object_merged, features =  'MARCO')

cluster0 = subset(object_merged, subset = mnn_clusters_0.3 == 0)
DimPlot(cluster0, group.by = 'orig.ident')
markers = FindAllMarkers(cluster0, features =  'MARCO')


saveRDS(object_merged,'object_merged (before cluster 5 removal).rds')

Idents(object_merged) = 'orig.ident'

markers = FindMarkers(object_merged, ident.1 = 'REX1',ident.2 = 'SAR1',logfc.threshold = 0, min.pct=0.2)
write.csv(markers, 'REX1 vs SAR1.csv')

VlnPlot(object_merged,'FABP5')

# Refine dataset after removing the small cluster
Idents(object_merged) = 'mnn_clusters_0.3'
DimPlot(object_merged,label = T)
table(object_merged$mnn_clusters_0.3)
#    0     1     2     3     4     5
# 10417  8398  8155  3820  3518    92
object_merged = subset(object_merged, subset = mnn_clusters_0.3 != 5)

DimPlot(object_merged,label = T,split.by = 'orig.ident', pt.size = 0.2)

markers = FindAllMarkers(object_merged, logfc.threshold = 0, min.pct = 0.2)
write.csv(markers, 'Markers per cluster (res = 0.3) min.pct = 0.2 logThreshold = 0.csv')

markers_4_vs_3 = FindMarkers(object_merged, ident.1 = 4, ident.2 = 3, logfc.threshold = 0, min.pct = 0.2)
write.csv(markers_4_vs_3, 'markers_4_vs_3 (min.pct = 0.2).csv')

Idents(object_merged) = 'orig.ident'
markers_REX_vs_SAR = FindMarkers(object_merged, ident.1 = 'REX1',ident.2 = 'SAR1',logfc.threshold = 0, min.pct = 0.2)
write.csv(markers_REX_vs_SAR, 'markers_REX_vs_SAR (min.pct = 0.2).csv')

markers_REX_vs_control = FindMarkers(object_merged, ident.1 = 'REX1',ident.2 = 'Control',logfc.threshold = 0, min.pct = 0.2)
write.csv(markers_REX_vs_control, 'markers_REX_vs_control (min.pct = 0.2).csv')

markers_SAR_vs_control = FindMarkers(object_merged, ident.1 = 'SAR1',ident.2 = 'Control',logfc.threshold = 0, min.pct = 0.2)
write.csv(markers_SAR_vs_control, 'markers_SAR_vs_control (min.pct = 0.2).csv')


DimPlot(object_merged,label = T,split.by = 'orig.ident', pt.size = 0.2)

gene = 'TNFSF10' # Increase with REX1 and decrease with SAR1
gene = 'PDK4' # Increase with SAR1 and decrease with REX1
Idents(object_merged) = 'orig.ident';plot1 = VlnPlot(object_merged,gene)
Idents(object_merged) = 'mnn_clusters_0.3';plot2 = VlnPlot(object_merged,gene)
plot1|plot2

FeaturePlot(object_merged, gene,split.by = 'orig.ident',pt.size = 0.2)

Idents(object_merged) = 'orig.ident';plot1 = VlnPlot(object_merged,gene)
Idents(object_merged) = 'mnn_clusters_0.3';plot2 = VlnPlot(object_merged,gene,split.by = 'orig.ident')
plot1|plot2


FeaturePlot(object_merged, c('CD200R1','IL7R','TNFRSF4'),split.by = 'orig.ident',pt.size = 0.2)
gene = c('CD200R1','IL7R','TNFRSF4')
Idents(object_merged) = 'orig.ident';plot1 = VlnPlot(object_merged,gene)
Idents(object_merged) = 'mnn_clusters_0.3';plot2 = VlnPlot(object_merged,gene)
plot1|plot2

FeaturePlot(object_merged, c('CD36'),split.by = 'orig.ident',pt.size = 0.2)
VlnPlot(object_merged,'CLEC2D')

for(i in genes)
{
  Idents(object_merged) = 'orig.ident';plot1 = VlnPlot(object_merged,i)
  Idents(object_merged) = 'mnn_clusters_0.3';plot2 = VlnPlot(object_merged,i)
  plot3 = FeaturePlot(object_merged, c(i),split.by = 'orig.ident',pt.size = 0.1)

  jpeg(paste0(i," (VlnPlot by sample).jpg"),res=500,width=2700,height=2800);print(plot1);dev.off()
  jpeg(paste0(i," (VlnPlot by cluster).jpg"),res=500,width=3000,height=2800);print(plot2);dev.off()
  jpeg(paste0(i," (UMAP).jpg"),res=600,width=6400,height=2200);print(plot3);dev.off()
}



gene = 'TNFSF10' # Increase with REX1 and decrease with SAR1
gene = 'PDK4' # Increase with SAR1 and decrease with REX1

gene = 'IGSF6' # Increase with REX1 and decrease with SAR1
gene = 'CD300A' # Increase with SAR1 and decrease with REX1

Idents(object_merged) = 'orig.ident';plot1 = VlnPlot(object_merged,gene)
Idents(object_merged) = 'mnn_clusters_0.3';plot2 = VlnPlot(object_merged,gene)
plot1|plot2

FeaturePlot(object_merged, gene,split.by = 'orig.ident',pt.size = 0.2)

saveRDS(object_merged, 'object (after removing the smallest cluster).rds')

Idents(object_merged) = 'mnn_clusters_0.3'
Idents(object_merged) = 'orig.ident'

plot = DotPlot(object_merged, features = unique(genes))+theme(axis.text.x = element_text(angle = 45, hjust=1))
jpeg(paste0("Top 15 transmembrane receptors per sample.jpg"),res=600,width=9300,height=2300);print(plot);dev.off()

# Compare each sample vs other samples
Idents(object_merged) = 'orig.ident'
markers = FindAllMarkers(object_merged,logfc.threshold = 0, min.pct = 0.2)
write.csv(markers, 'Sample markers (each sample vs other samples).csv')


genes = c('MARCO','MERTK','CD163')


for(i in genes)
{
  Idents(object_merged) = 'orig.ident';plot1 = VlnPlot(object_merged,i)
  Idents(object_merged) = 'mnn_clusters_0.3';plot2 = VlnPlot(object_merged,i)
  plot3 = FeaturePlot(object_merged, c(i),split.by = 'orig.ident',pt.size = 0.1)

  jpeg(paste0(i," (VlnPlot by sample).jpg"),res=500,width=2700,height=2800);print(plot1);dev.off()
  jpeg(paste0(i," (VlnPlot by cluster).jpg"),res=500,width=3000,height=2800);print(plot2);dev.off()
  jpeg(paste0(i," (UMAP).jpg"),res=600,width=6400,height=2200);print(plot3);dev.off()
}



object_macro_inf = readRDS('object (after removing the smallest cluster).rds')
Idents(object_macro_inf) = 'orig.ident'

i = 'MARCO'
i = 'MERTK'
i = 'CD163'

plot = VlnPlot(object_macro_inf, i, pt.size = 0.01)
plot = VlnPlot(object_macro_inf, i, pt.size = 0.01, idents = c('SAR1','REX1'))
jpeg(paste0(i," (VlnPlot by sample) 2.jpg"),res=500,width=2700,height=2800);print(plot);dev.off()

FeaturePlot(object_macro_inf, gene)
FeaturePlot(object_macro_inf, gene,split.by = 'orig.ident')
```
