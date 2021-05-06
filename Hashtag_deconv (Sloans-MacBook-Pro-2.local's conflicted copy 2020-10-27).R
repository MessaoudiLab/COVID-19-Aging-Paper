library(Seurat)
library(cowplot)
library(dplyr) 
library(Matrix) 
library(gdata)
library(DESeq2)
library(ggplot2)

#Aged d1
aged_d1_1.data <- Read10X(data.dir = "../cellranger1/Aged_D1_aged_d1.combined_counts_1/outs/filtered_feature_bc_matrix/")
aged_d1_2.data <- Read10X(data.dir = "../cellranger1/Aged_D1_aged_d1.combined_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_aged_d1_1 = CreateSeuratObject(counts = aged_d1_1.data$`Gene Expression`, project = "aged_d1_1")
seurat_object_aged_d1_2 = CreateSeuratObject(counts = aged_d1_2.data$`Gene Expression`, project = "aged_d1_2")
seurat_object_aged_d1_1[["ADT"]] <- CreateAssayObject(counts = aged_d1_1.data$`Antibody Capture`)
seurat_object_aged_d1_2[["ADT"]] <- CreateAssayObject(counts = aged_d1_2.data$`Antibody Capture`)
aged_d1.combined <- merge(seurat_object_aged_d1_1, y = seurat_object_aged_d1_2, add.cell.ids = c("r1", "r2"))

aged_d1.combined <- PercentageFeatureSet(aged_d1.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("aged_d1_qc_violins.pdf")
VlnPlot(aged_d1.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

aged_d1.combined <- subset(aged_d1.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
aged_d1.combined$sample <- "aged_d1"



aged_d1.combined <- SCTransform(aged_d1.combined, vars.to.regress = "percent.mt", verbose = FALSE)

aged_d1.combined <- NormalizeData(aged_d1.combined, assay = "ADT", normalization.method = "CLR")
aged_d1.combined.hashtag <- HTODemux(aged_d1.combined, assay = "ADT", positive.quantile = 0.99)
table(aged_d1.combined.hashtag$ADT_classification.global)
Idents(aged_d1.combined.hashtag) <- "ADT_classification.global"
pdf("vln_adt.pdf")
VlnPlot(aged_d1.combined.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

# First, we will remove negative cells from the object
aged_d1.combined.hashtag.subset <- subset(aged_d1.combined.hashtag, idents = "Negative", invert = TRUE)

# Calculate a distance matrix using HTO
hto.dist.mtx <- as.matrix(dist(t(GetAssayData(object = aged_d1.combined.hashtag.subset, assay = "HTO"))))

# Calculate tSNE embeddings with a distance matrix
aged_d1.combined.hashtag.subset <- RunTSNE(aged_d1.combined.hashtag.subset, distance.matrix = hto.dist.mtx, perplexity = 100)
pdf("tsne_doublet.pdf")
DimPlot(aged_d1.combined.hashtag.subset)
dev.off()
# Extract the singlets
aged_d1.combined.singlet <- subset(aged_d1.combined.hashtag, idents = "Singlet")
# Select the top 1000 most variable features
aged_d1.combined.singlet <- FindVariableFeatures(aged_d1.combined.singlet, selection.method = "mean.var.plot")

# Scaling RNA data, we only scale the variable features here for efficiency
aged_d1.combined.singlet <- ScaleData(aged_d1.combined.singlet, features = VariableFeatures(aged_d1.combined.singlet))

# Run PCA
aged_d1.combined.singlet <- RunPCA(aged_d1.combined.singlet, features = VariableFeatures(aged_d1.combined.singlet))
# We select the top 10 PCs for clustering and tSNE based on PCElbowPlot
aged_d1.combined.singlet <- FindNeighbors(aged_d1.combined.singlet, reduction = "pca", dims = 1:10)
aged_d1.combined.singlet <- FindClusters(aged_d1.combined.singlet, resolution = 0.6, verbose = FALSE)
aged_d1.combined.singlet <- RunUMAP(aged_d1.combined.singlet, reduction = "pca", dims = 1:20)

pdf("umap.pdf")
DimPlot(aged_d1.combined.singlet, reduction= "umap", label = TRUE) + NoLegend()
dev.off()

pdf("umap_adt.pdf")
DimPlot(aged_d1.combined.singlet, reduction= "umap", label = TRUE, group.by = "ADT_classification")
dev.off()

saveRDS(aged_d1.combined, "aged_d1_merged_raw.rds") 



#Aged d8
aged_d8_1.data <- Read10X(data.dir = "../cellranger1/Aged_D8_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
aged_d8_2.data <- Read10X(data.dir = "../cellranger1/Aged_D8_PBMC_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_aged_d8_1 = CreateSeuratObject(counts = aged_d8_1.data$`Gene Expression`, project = "aged_d8_1")
seurat_object_aged_d8_2 = CreateSeuratObject(counts = aged_d8_2.data$`Gene Expression`, project = "aged_d8_2")
seurat_object_aged_d8_1[["ADT"]] <- CreateAssayObject(counts = aged_d8_1.data$`Antibody Capture`)
seurat_object_aged_d8_2[["ADT"]] <- CreateAssayObject(counts = aged_d8_2.data$`Antibody Capture`)
aged_d8.combined <- merge(seurat_object_aged_d8_1, y = seurat_object_aged_d8_2, add.cell.ids = c("r1", "r2"))

aged_d8.combined <- PercentageFeatureSet(aged_d8.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("aged_d8_qc_violins.pdf")
VlnPlot(aged_d8.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

aged_d8.combined <- subset(aged_d8.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
aged_d8.combined$sample <- "aged_d8"

#Aged d11
aged_d11_1.data <- Read10X(data.dir = "../cellranger1/Aged_D11_PBMC_counts/outs/filtered_feature_bc_matrix/")
seurat_object_aged_d11_1 = CreateSeuratObject(counts = aged_d11_1.data$`Gene Expression`, project = "aged_d11_1")
seurat_object_aged_d11_1[["ADT"]] <- CreateAssayObject(counts = aged_d11_1.data$`Antibody Capture`)
aged_d11.combined <- seurat_object_aged_d11_1
aged_d11.combined <- PercentageFeatureSet(aged_d11.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("aged_d11_qc_violins.pdf")
VlnPlot(aged_d11.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

aged_d11.combined <- subset(aged_d11.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
aged_d11.combined$sample <- "aged_d11"

#Aged controls
aged_control_1.data <- Read10X(data.dir = "../cellranger1/Aged_control_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
aged_control_2.data <- Read10X(data.dir = "../cellranger1/Aged_control_PBMC_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_aged_control_1 = CreateSeuratObject(counts = aged_control_1.data$`Gene Expression`, project = "aged_control_1")
seurat_object_aged_control_2 = CreateSeuratObject(counts = aged_control_2.data$`Gene Expression`, project = "aged_control_2")
seurat_object_aged_control_1[["ADT"]] <- CreateAssayObject(counts = aged_control_1.data$`Antibody Capture`)
seurat_object_aged_control_2[["ADT"]] <- CreateAssayObject(counts = aged_control_2.data$`Antibody Capture`)
aged_control.combined <- merge(seurat_object_aged_control_1, y = seurat_object_aged_control_2, add.cell.ids = c("r1", "r2"))

aged_control.combined <- PercentageFeatureSet(aged_control.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("aged_control_qc_violins.pdf")
VlnPlot(aged_control.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

aged_control.combined <- subset(aged_control.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
aged_control.combined$sample <- "aged_control"

#Aged mild
aged_mild_1.data <- Read10X(data.dir = "../cellranger1/Aged_mild_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
aged_mild_2.data <- Read10X(data.dir = "../cellranger1/Aged_mild_PBMC_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_aged_mild_1 = CreateSeuratObject(counts = aged_mild_1.data$`Gene Expression`, project = "aged_mild_1")
seurat_object_aged_mild_2 = CreateSeuratObject(counts = aged_mild_2.data$`Gene Expression`, project = "aged_mild_2")
seurat_object_aged_mild_1[["ADT"]] <- CreateAssayObject(counts = aged_mild_1.data$`Antibody Capture`)
seurat_object_aged_mild_2[["ADT"]] <- CreateAssayObject(counts = aged_mild_2.data$`Antibody Capture`)
aged_mild.combined <- merge(seurat_object_aged_mild_1, y = seurat_object_aged_mild_2, add.cell.ids = c("r1", "r2"))

aged_mild.combined <- PercentageFeatureSet(aged_mild.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("aged_mild_qc_violins.pdf")
VlnPlot(aged_mild.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

aged_mild.combined <- subset(aged_mild.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
aged_mild.combined$sample <- "aged_mild"



#Young controls
young_control_1.data <- Read10X(data.dir = "../cellranger1/young_control_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
seurat_object_young_control_1 = CreateSeuratObject(counts = young_control_1.data$`Gene Expression`, project = "young_control_1")
seurat_object_young_control_1[["ADT"]] <- CreateAssayObject(counts = young_control_1.data$`Antibody Capture`)
young_control.combined <- seurat_object_young_control_1
young_control.combined <- PercentageFeatureSet(young_control.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("young_control_qc_violins.pdf")
VlnPlot(young_control.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

young_control.combined <- subset(young_control.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
young_control.combined$sample <- "young_control"

#Young d1
young_d1_1.data <- Read10X(data.dir = "../cellranger1/Young_D1_PBMC_counts_1/outs/filtered_feature_bc_matrix_Remdes_PBMC_5/")
young_d1_2.data <- Read10X(data.dir = "../cellranger1/Young_D1_PBMC_counts_2_nofb/outs/filtered_feature_bc_matrix/")
seurat_object_young_d1_1 = CreateSeuratObject(counts = young_d1_1.data$`Gene Expression`, project = "young_d1_1")
seurat_object_young_d1_2 = CreateSeuratObject(counts = young_d1_2.data, project = "young_d1_2")
seurat_object_young_d1_1[["ADT"]] <- CreateAssayObject(counts = young_d1_1.data$`Antibody Capture`)
young_d1_fb <- seurat_object_young_d1_1
young_d1_fb <- PercentageFeatureSet(young_d1_fb, pattern = "^MT-", col.name = "percent.mt")
pdf("young_d1_fb_qc_violins.pdf")
VlnPlot(young_d1_fb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
young_d1_fb <- subset(young_d1_fb, subset = nFeature_RNA > 200 & percent.mt < 25)
young_d1_fb$sample <- "young_d1"
young_d1_nofb <- seurat_object_young_d1_2
young_d1_nofb <- PercentageFeatureSet(young_d1_nofb, pattern = "^MT-", col.name = "percent.mt")
pdf("young_d1_nofb_qc_violins.pdf")
VlnPlot(young_d1_nofb, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

young_d1_nofb <- subset(young_d1_nofb, subset = nFeature_RNA > 200 & percent.mt < 25)
young_d1_nofb$sample <- "young_d1"


#Young d8
young_d8_1.data <- Read10X(data.dir = "../cellranger1/Young_D8_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
young_d8_2.data <- Read10X(data.dir = "../cellranger1/Young_D8_PBMC_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_young_d8_1 = CreateSeuratObject(counts = young_d8_1.data$`Gene Expression`, project = "young_d8_1")
seurat_object_young_d8_2 = CreateSeuratObject(counts = young_d8_2.data$`Gene Expression`, project = "young_d8_2")
seurat_object_young_d8_1[["ADT"]] <- CreateAssayObject(counts = young_d8_1.data$`Antibody Capture`)
seurat_object_young_d8_2[["ADT"]] <- CreateAssayObject(counts = young_d8_2.data$`Antibody Capture`)
young_d8.combined <- merge(seurat_object_young_d8_1, y = seurat_object_young_d8_2, add.cell.ids = c("r1", "r2"))

young_d8.combined <- PercentageFeatureSet(young_d8.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("young_d8_qc_violins.pdf")
VlnPlot(young_d8.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

young_d8.combined <- subset(young_d8.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
young_d8.combined$sample <- "young_d8"


#Young d8
young_d11_1.data <- Read10X(data.dir = "../cellranger1/Young_D11_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
young_d11_2.data <- Read10X(data.dir = "../cellranger1/Young_D11_PBMC_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_young_d11_1 = CreateSeuratObject(counts = young_d11_1.data$`Gene Expression`, project = "young_d11_1")
seurat_object_young_d11_2 = CreateSeuratObject(counts = young_d11_2.data$`Gene Expression`, project = "young_d11_2")
seurat_object_young_d11_1[["ADT"]] <- CreateAssayObject(counts = young_d11_1.data$`Antibody Capture`)
seurat_object_young_d11_2[["ADT"]] <- CreateAssayObject(counts = young_d11_2.data$`Antibody Capture`)
young_d11.combined <- merge(seurat_object_young_d11_1, y = seurat_object_young_d11_2, add.cell.ids = c("r1", "r2"))

young_d11.combined <- PercentageFeatureSet(young_d11.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("young_d11_qc_violins.pdf")
VlnPlot(young_d11.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

young_d11.combined <- subset(young_d11.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
young_d11.combined$sample <- "young_d11"

#Young mild
young_mild_1.data <- Read10X(data.dir = "../cellranger1/Young_mild_PBMC_counts_1/outs/filtered_feature_bc_matrix/")
young_mild_2.data <- Read10X(data.dir = "../cellranger1/Young_mild_PBMC_counts_2/outs/filtered_feature_bc_matrix/")
seurat_object_young_mild_1 = CreateSeuratObject(counts = young_mild_1.data$`Gene Expression`, project = "young_mild_1")
seurat_object_young_mild_2 = CreateSeuratObject(counts = young_mild_2.data$`Gene Expression`, project = "young_mild_2")
seurat_object_young_mild_1[["ADT"]] <- CreateAssayObject(counts = young_mild_1.data$`Antibody Capture`)
seurat_object_young_mild_2[["ADT"]] <- CreateAssayObject(counts = young_mild_2.data$`Antibody Capture`)
young_mild.combined <- merge(seurat_object_young_mild_1, y = seurat_object_young_mild_2, add.cell.ids = c("r1", "r2"))

young_mild.combined <- PercentageFeatureSet(young_mild.combined, pattern = "^MT-", col.name = "percent.mt")
pdf("young_mild_qc_violins.pdf")
VlnPlot(young_mild.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

young_mild.combined <- subset(young_mild.combined, subset = nFeature_RNA > 200 & percent.mt < 25)
young_mild.combined$sample <- "young_mild"

saveRDS(young_mild.combined, "young_mild_combined.rds")
saveRDS(aged_mild.combined, "aged_mild_combined.rds")
saveRDS(young_control.combined, "young_control_combined.rds")
saveRDS(aged_control.combined, "aged_control_combined.rds")
saveRDS(young_d1_nofb, "young_d1_nofb.rds")
saveRDS(young_d1_fb, "young_d1_fb.rds")
saveRDS(aged_d1.combined, "aged_d1_combined.rds")
saveRDS(young_d8.combined, "young_d8_combined.rds")
saveRDS(aged_d8.combined, "aged_d8_combined.rds")
saveRDS(young_d11.combined, "young_d11_combined.rds")
saveRDS(aged_d11.combined, "aged_d11_combined.rds")

young.combined <- merge(young_control.combined, y=c(young_mild.combined, young_d1_nofb, young_d1_fb, young_d8.combined, young_d11.combined), add.cell.ids = c("HC", "MILD", "D1_nofb", "D1_fb", "D8", "D11"), project = "young")
aged.combined <- merge(aged_control.combined, y=c(aged_mild.combined, aged_d1.combined, aged_d8.combined, aged_d11.combined), add.cell.ids = c("HC", "MILD", "D1", "D8", "D11"), project = "AGED")
all.combined <- merge(young.combined, y= aged.combined, add.cell.ids = c("young", "aged"), project = "ALL_PBMC")
saveRDS(all.combined, "all.combined.rds")

all.list <- SplitObject(all.combined, split.by = "sample")
for (i in names(all.list)) {
  all.list[[i]] <- SCTransform(all.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
all.features <- SelectIntegrationFeatures(object.list = all.list, nfeatures = 3000)
all.list <- PrepSCTIntegration(object.list = all.list, anchor.features = all.features)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(all.list) == "young_control")

all.anchors <- FindIntegrationAnchors(object.list = all.list, normalization.method = "SCT", 
                                       anchor.features = all.features, reference = reference_dataset)
all.integrated <- IntegrateData(anchorset = all.anchors, normalization.method = "SCT")

all.integrated <- RunPCA(object = all.integrated, verbose = FALSE)
all.integrated <- RunUMAP(object = all.integrated, dims = 1:30)
all.integrated <- FindNeighbors(all.integrated, reduction = "pca", dims = 1:30)
all.integrated <- FindClusters(all.integrated, resolution = 0.5)
pdf("umap_all.pdf")
a <- DimPlot(all.integrated, reduction= "umap", label = TRUE) + NoLegend()
b <- DimPlot(all.integrated, reduction= "umap", label = FALSE, group.by = "COVID")
a+b
dev.off()

all.integrated <- NormalizeData(all.integrated, assay = "ADT", normalization.method = "CLR")
all.integrated.hashtag <- HTODemux(all.integrated, assay = "ADT", positive.quantile = 0.99)
table(all.integrated.hashtag$ADT_classification.global)
Idents(all.integrated.hashtag) <- "ADT_classification.global"
pdf("vln_adt.pdf")
VlnPlot(all.integrated.hashtag, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
dev.off()

all.integrated.hashtag$sample.ADT = paste(all.integrated.hashtag$sample,all.integrated.hashtag$ADT_classification)

DefaultAssay(all.integrated.hashtag) <- "RNA"
all.integrated.hashtag <- NormalizeData(all.integrated.hashtag,  verbose = FALSE)

pdf("umap_all_patient.pdf")
DimPlot(all.integrated.hashtag, reduction= "umap", group.by = "patient")
dev.off()

#remove doublets
all.integrated.hashtag <- subset(all.integrated.hashtag, idents = "Doublet", invert = TRUE)

Idents(all.integrated.hashtag) <- "seurat_clusters"

pdf("umap_all.pdf")
DimPlot(all.integrated.hashtag, reduction= "umap", label=TRUE)
dev.off()

pdf("umap_all_group.pdf")
DimPlot(all.integrated.hashtag, reduction= "umap", label=TRUE, split.by = "sample")
dev.off()

t <- table(all.integrated.hashtag@meta.data$covsev_death, all.integrated.hashtag@meta.data$celltype)
write.table(t, "all.integrated.hashtag_distribution_covsevdeath_12042020.xls", sep="\t", quote=F, col.names=NA)

combined.markers <- FindAllMarkers(all.integrated.hashtag, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(combined.markers, "./all.integrated.hashtag_combined_markers_1.xls", sep="\t", col.names=NA)
DefaultAssay(all.integrated.hashtag) <- "RNA"
x <- DimPlot(all.integrated.hashtag, reduction= "umap", label=TRUE) + NoLegend()
y <- FeaturePlot(all.integrated.hashtag, features = c("GZMH"), cols = c("grey90", "darkmagenta"), min.cutoff = "q9")
pdf("delete.pdf")
x+y
dev.off()

FeaturePlot(all.integrated.hashtag, features = c("CD3G", "CD3E", "CD4","CD8A", "IL7R", "CD27", "SELL", "CCR7", "CD28"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
FeaturePlot(all.integrated.hashtag, features = c("FAS", "GZMB", "GZMK", "GNLY", "PRF1", "S100A4", "FOXP3", "IL2RA", "TBX21"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
FeaturePlot(all.integrated.hashtag, features = c("GATA3", "RORC", "MKI67", "FCGR3A", "NCAM1", "NKG7", "CD14", "HLA-DRA", "PPBP"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
FeaturePlot(all.integrated.hashtag, features = c("MZB1", "ALAS2", "CD34"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
FeaturePlot(all.integrated.hashtag, features = c("SLC4A10", "TRDV2", "LEF1", "GATA2", "KLRF1", "TRAV1-2"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
#paper
FeaturePlot(all.integrated.hashtag, features = c("CD3G", "CD8A", "MKI67", "CD79A"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
FeaturePlot(all.integrated.hashtag, features = c("FCER1A", "FCGR3A", "KLRF1", "CD14"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
FeaturePlot(all.integrated.hashtag, features = c("CD34", "MZB1", "PPBP", "ALAS2"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))

all.integrated.hashtag <- RenameIdents(all.integrated.hashtag, `0`= "CD4 naive", `1`= "Monocytes CD14", `2`= "CD4 EM", `3`= "CD8 memory", `4`= "NK CD16", `5`= "B cells", `6`= "CD8 memory", 
                               `7`="CD4 memory", `8`="NKT", `9`= "Monocytes Activated", `10`= "CD8 naive", `11`= "Treg/prolif T", `12`= "CD8 viral", `13`= "Monocytes CD16", `14`= "Erythroid", 
                               `15`= "Platelets", `16`="NK CD56", `17`="DC", `18`="Erythroid", `19`= "Progenitor CD34", `20`= "Plasmablasts")
all.integrated.hashtag$celltype <- Idents(all.integrated.hashtag)

Idents(all.integrated.hashtag) <- "seurat_clusters"
all.integrated.hashtag <- RenameIdents(all.integrated.hashtag, `0`= "CD4 naive", `1`= "Monocytes CD14", `2`= "CD4 memory", `3`= "CD8 memory", `4`= "NK", `5`= "B cells", `6`= "CD8 memory", 
                                       `7`="CD4 memory", `8`="NKT", `9`= "Monocytes CD14", `10`= "CD8 naive", `11`= "Treg/prolif T", `12`= "CD8 memory", `13`= "Monocytes CD16", `14`= "Erythroid", 
                                       `15`= "Platelets", `16`="NK", `17`="DC", `18`="Erythroid", `19`= "Progenitor CD34", `20`= "Plasmablasts")
all.integrated.hashtag$celltype <- Idents(all.integrated.hashtag)

#Paper
Idents(all.integrated.hashtag) <- "seurat_clusters"
all.integrated.hashtag <- RenameIdents(all.integrated.hashtag, `0`= "CD4 naive", `1`= "Monocytes CD14", `2`= "CD4 memory", `3`= "CD8 memory", `4`= "NK", `5`= "B cells", `6`= "CD8 memory", 
                                       `7`="CD4 memory", `8`="NKT", `9`= "Monocytes CD14", `10`= "CD8 naive", `11`= "Treg/prolif T", `12`= "CD8 memory", `13`= "Monocytes CD16", `14`= "Erythroid", 
                                       `15`= "Platelets", `16`="NK", `17`="DC", `18`="Erythroid", `19`= "Progenitor CD34", `20`= "Plasmablasts")
all.integrated.hashtag$celltype_paper <- Idents(all.integrated.hashtag)


pdf("dimplot_all_celltype.pdf")
DimPlot(all.integrated.hashtag, reduction= "umap", label=TRUE, group.by = "celltype") + NoLegend()
dev.off()

a <- DimPlot(all.integrated.hashtag, reduction= "umap", label = TRUE,  pt.size=0.1) + NoLegend()
pdf("umap_all_patient_grouped.pdf")
b <- DimPlot(all.integrated.hashtag, reduction= "umap", label = FALSE, shuffle=TRUE, group.by = "patient_grouped", pt.size=0.01, cols= DiscretePalette(30, "alphabet"))
a+b
dev.off()

#marker dotplot

DefaultAssay(all.integrated.hashtag) <- "RNA"
my_levels <- c("B cells", "Plasmablasts", "CD4 naive", "CD4 memory", "CD8 naive", "CD8 memory", "Treg/prolif T", "NK", "NKT", "DC", "Monocytes CD14", "Monocytes CD16", "Progenitor CD34", "Erythroid", "Platelets")
factor(Idents(all.integrated.hashtag), levels= rev(my_levels))
Idents(all.integrated.hashtag) <- factor(Idents(all.integrated.hashtag), levels= rev(my_levels))
markers.to.plot <- c("CD3G", "CD8A", "GNLY", "GZMK", "CD27", "FOXP3", "MKI67", "NCAM1", "KLRF1", "NKG7", "CD79A", "MZB1", "CD14", "FCGR3A", "FCER1A", "KIT", "CD34", "ALAS2", "PPBP")
markers.to.plot <- c("KIT", "CD34", "FCER1A", "LYZ", "FCN1", "MKI67", "GZMH", "KLRF1", "NKG7", "TRGV9", "CD4", "CCR7", "CD8A", "CD3E", "IL7R", "CD79A", "MZB1", "ALAS2", "PPBP")
pdf("dotplot_all_paper.pdf")
DotPlot(all.integrated.hashtag, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

markers.to.plot <- c("CD79A", "MZB1", "CD3E",  "CCR7", "IL7R", "CD8A", "GZMH", "MKI67", "KLRF1", "NKG7", "NCAM1", "FCER1A", "LYZ", "FCN1" ,"KIT", "CD34", "ALAS2", "PPBP")
pdf("dotplot_all_paper.pdf", width=10, height=10)
DotPlot(all.integrated.hashtag, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()



#Subsets
#TNK
all.integrated.hashtag.TNK <- SubsetData(object = all.integrated.hashtag, ident.use = c("CD4 EM", "CD4 memory", "CD4 naive", "CD8 memory", "CD8 naive", "CD8 viral", "NK CD16", "NK CD56", "NKT", "Treg/prolif T"))
allTNK.list <- SplitObject(all.integrated.hashtag.TNK, split.by = "sample")
for (i in names(allTNK.list)) {
  allTNK.list[[i]] <- SCTransform(allTNK.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
allTNK.features <- SelectIntegrationFeatures(object.list = allTNK.list, nfeatures = 3000)
allTNK.list <- PrepSCTIntegration(object.list = allTNK.list, anchor.features = allTNK.features)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(allTNK.list) == "young_control")

allTNK.anchors <- FindIntegrationAnchors(object.list = allTNK.list, normalization.method = "SCT", 
                                      anchor.features = allTNK.features, reference = reference_dataset)
allTNK.integrated <- IntegrateData(anchorset = allTNK.anchors, normalization.method = "SCT")

allTNK.integrated <- RunPCA(object = allTNK.integrated, verbose = FALSE)
allTNK.integrated <- RunUMAP(object = allTNK.integrated, dims = 1:30)
allTNK.integrated <- FindNeighbors(allTNK.integrated, reduction = "pca", dims = 1:30)
allTNK.integrated <- FindClusters(allTNK.integrated, resolution = 0.5)
pdf("umap_allTNK.pdf")
a <- DimPlot(allTNK.integrated, reduction= "umap", label = TRUE) + NoLegend()
b <- DimPlot(allTNK.integrated, reduction= "umap", label = FALSE, group.by = "COVID")
a+b
dev.off()

combined.markers <- FindAllMarkers(allTNK.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(combined.markers, "./all.integrated.hashtag_combined_markers_TNK.xls", sep="\t", col.names=NA)

#slingshot subsets
allCD4.integrated <- SubsetData(object = allTNK.integrated, ident.use = c("0", "1", "4", "6", "10"))
allCD8.integrated <- SubsetData(object = allTNK.integrated, ident.use = c("5", "2", "9", "7", "11"))

Idents(allTNK.integrated) <- "seurat_clusters"  
allTNK.integrated <- RenameIdents(allTNK.integrated, `0`= "CD4 1", `1`= "CD4 2", `2`= "CD8 1", 
                                  `3`= "NK CD16", `4`= "CD4 3", `5`= "CD8 2", `6`= "CD4 4",`7`="CD8 3", `8`="NK 3", `9`="CD8 naive", 
                                  `10`="Treg", `11`="CD8 4", `12`="NK CD56", `13`="MAIT/iNKT", `14`="T prolif")
allTNK.integrated$celltype_subset <- Idents(allTNK.integrated)                                      
Idents(allTNK.integrated) <- "celltype_subset" 

t <- table(allTNK.integrated@meta.data$patient, allTNK.integrated@meta.data$celltype_subset)
write.table(t, "all.integrated.hashtag_distribution_TNK_DPS_bins.xls", sep="\t", quote=F, col.names=NA)

DefaultAssay(allTNK.integrated) <- "RNA"
allTNK.integrated <- NormalizeData(allTNK.integrated, verbose = FALSE)

my_levels <- c("NK 3", "NK CD56", "NK CD16", "T prolif", "MAIT/iNKT","CD8 2", "CD8 1", "CD8 3", "CD8 4", "CD8 naive", "Treg", "CD4 1", "CD4 4","CD4 3", "CD4 2")
Idents(allTNK.integrated) <- factor(Idents(allTNK.integrated), levels= my_levels)


my_levels <- c("NK 3", "NK CD56", "NK CD16")
Idents(NK) <- factor(Idents(NK), levels= my_levels)
markers.to.plot <- c("NKG7", "FCGR3A", "PRF1", "CCL3","GZMH", "KLRD1", "XCL1", "NCAM1", "GZMK", "GZMB", "CD8A", "GNLY", "IFNG", "PLCG2", "CXCR3")
pdf("dotplot_NKs.pdf", width=7, height=5)
DotPlot(NK, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

my_levels <- c("T prolif", "MAIT/iNKT","CD8 2", "CD8 1", "CD8 3", "CD8 4", "CD8 naive", "Treg", "CD4 4", "CD4 1","CD4 3", "CD4 2")
Idents(T) <- factor(Idents(T), levels= my_levels)
markers.to.plot <- c("CD3E", "IL7R", "CCR7", "LEF1", "LTB", "GPR183", "PTPRC", "RORA", "S100A8", "S100A9", "FOXP3", "CTLA4", "CD8A", "IFIT2", "IFIT3", "OASL", "IFNG", "TNF", "CCL3", "CCL4", "CD69", "S100A4", "KLRG1", "GZMH", "CXCR4", "NKG7", "GNLY", "PRF1", "GZMB", "GZMK", "CXCR3", "SLC4A10",  "MKI67")
pdf("dotplot_Tcells.pdf", width=15, height=5)
DotPlot(T, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()



markers.to.plot <- c("LEF1", "ATM", "SELL", "KLF2", "ITGA6", "IL7R", "CD52", "S100A4", "TGFB3", "AQP3", "NLRP3", "ITGB7", "IFIT3", "IFIT2", "STAT1", "MX1", "IRF7", "ISG15", "IFITM3", "OAS2", "JAK2", "SOCS1", "TRIM21")
pdf("dotplot_TNK_cd44.pdf", width=15, height=5)
DotPlot(allTNK.integrated, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

#NKT genes
pdf("NKT1_feature.pdf")
FeaturePlot(allTNK.integrated, features = c("FCGR3A", "NCAM1", "B3GAT1","LAG3", "PDCD1", "HAVCR2", "CCL3", "CCL4", "IFNG"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
dev.off()

pdf("NKT3_feature.pdf")
FeaturePlot(allTNK.integrated, features = c("BLK", "CCR3", "IL17RB", "RORA", "DAO1"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
dev.off()


features <- c("FCGR3A", "NCAM1", "B3GAT1","LAG3", "PDCD1", "HAVCR2", "CCL3", "CCL4", "IFNG")
pdf("stacked_violin_NKT.pdf", width=5, height=10)
StackedVlnPlot(obj = NK, features = features)
dev.off()

my_levels <- c("HD_young", "MILD_young", "SEVERE_young", "HD_aged", "MILD_aged", "SEVERE_aged")
Idents(NK16) <- factor(Idents(NK16), levels= my_levels)

pdf("violin_GZMB_young.pdf")
VlnPlot(NK56, features = "GZMB", group.by = "COVID", idents = c("HD_young", "MILD_young", "SEVERE_young"),
                 pt.size = 0, combine = FALSE)
dev.off()
pdf("violin_GZMB_aged.pdf")
VlnPlot(NK56, features = "GZMB", group.by = "COVID", idents = c("HD_aged", "MILD_aged", "SEVERE_aged"),
        pt.size = 0, combine = FALSE)
dev.off()


markers.to.plot <- c("CD3E", "FCGR3A", "NCAM1", "KLRC2", "B3GAT1", "HLA-DRA", "CD38", "SELL", "CCL3", "CCL4", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "HAVCR2", "LAG3", "PDCD1", "NKTR", "S100A9", "CD6", "PLCG2", "NLRC5")
pdf("dotplot_NK_COVID.pdf", width=15, height=5)
DotPlot(NK, features = markers.to.plot, split.by="COVID", dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

markers.to.plot <- c("CD3E", "FCGR3A", "NCAM1", "KLRC2", "B3GAT1", "HLA-DRA", "CD38", "SELL", "CCL3", "CCL4", "CCL5", "GZMA", "GZMB", "PRF1", "MKI67", "HAVCR2", "LAG3", "PDCD1", "NKTR", "S100A9", "CD6", "PLCG2", "NLRC5")
pdf("dotplot_NK_COVID.pdf", width=15, height=5)
DotPlot(NK, features = markers.to.plot, split.by="COVID", dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

markers.to.plot <- c("FCGR3A", "NCAM1", "B3GAT1","LAG3", "PDCD1", "HAVCR2", "CCL3", "CCL4", "IFNG", "CD3E", "CD3D", "CD3G", "IFI44L", "ISG20", "NKG7", "TRDC", "TRGV9", "KLRC1", "KLRC2", "GZMH", "GZMK", "LTB")
pdf("NKT_genes.pdf", width=15, height=5)
DotPlot(NK, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

pdf("ridge_GZMB_dps.pdf")
RidgePlot(NK56, features = "GZMB", idents=c("young_2", "young_3", "young_4", "aged_1", "aged_2", "aged_3", "aged_4"))
dev.off()

pdf("violin_CCL4.pdf")
VlnPlot(NK3, features = "CCL4", group.by = "COVID",
        pt.size = 0, combine = FALSE)
dev.off()


features<- c("GZMB", "TGFB1", "ICAM1", "BIRC2", "PRF1", "PIK3R1")
features<- c("XCL2", "XCL1", "GZMK", "CCR7", "IL7R", "TNFSF10", "FCGR3A", "NCAM1", "KLRC1", "CD2", "LTB", "SELL")
features <- c("XIST", "PLCG2", "IFI44L", "IFIT2", "NKG7")
pdf("stacked_violin_nk56_young4.pdf", width=5, height = 20)
StackedVlnPlot(obj = NK56, features = features, idents=c("HD_young", "MILD_young", "SEVERE_young"))
dev.off()


features<- c("GZMB", "TGFB1", "ICAM1", "S100A6", "CD82", "FURIN", "ISG20", "VIM", "CD52", "BIRC2", "PRF1", "PIK3R1")
features <- c("XIST", "PLCG2", "IFI44L", "IFIT2", "NKG7")
pdf("stacked_violin_nk56_aged4.pdf", width=5, height = 20)
StackedVlnPlot(obj = NK56, features = features, idents=c("HD_aged", "MILD_aged", "SEVERE_aged"))
dev.off()

features <- c("FOSB", "IFNG", "CCL3", "NKG7", "KLRB1", "B2M", "IFI44L", "PLCG2", "XAF1", "HLA-C", "HLA-B")
pdf("stacked_violin_nk3_young.pdf", width=5, height = 25)
StackedVlnPlot(obj = NK3, features = features, idents=c("HD_young", "MILD_young", "SEVERE_young"))
dev.off()

features <- c("FOSB", "IFNG", "CCL3", "NKG7", "KLRB1", "B2M", "IFI44L", "PLCG2", "XAF1", "HLA-C", "HLA-B")
pdf("stacked_violin_nk3_aged.pdf", width=5, height = 25)
StackedVlnPlot(obj = NK3, features = features, idents=c("HD_aged", "MILD_aged", "SEVERE_aged"))
dev.off()

markers.to.plot <- c("SAMSN1", "GZMB", "ICAM1", "S100A6", "CD82", "FURIN", "ISG20", "KLF2", "CXCR4", "IL2RB", "DUSP1", "JUN", "REL", "NFKBIA")
pdf("dotplot_NK56_COVID.pdf", width=15, height=5)
DotPlot(nk56, features = markers.to.plot, split.by="COVID", dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()




#DEG CD4
DefaultAssay(CD4_1_subset.integrated) <- "RNA"
CD4_1_subset.integrated <- SubsetData(object = cd4, ident.use = "0")
Idents(CD4_1_subset.integrated) <- "COVID"
deg1 <- FindMarkers(CD4_1_subset.integrated, ident.1="MILD_young", ident.2="HD_young", method="MAST")
write.table(deg1, "CD4_1_subset_mildyoung-hdyoung_mast.xls", sep="\t", quote=F, col.names=NA)
deg2 <- FindMarkers(CD4_1_subset.integrated, ident.1="SEVERE_young", ident.2="HD_young", method="MAST")
write.table(deg2, "CD4_1_subset_severeyoung-hdyoung_mast.xls", sep="\t", quote=F, col.names=NA)
deg3 <- FindMarkers(CD4_1_subset.integrated, ident.1="MILD_aged", ident.2="HD_aged", method="MAST")
write.table(deg3, "CD4_1_subset_mildaged-hdaged_mast.xls", sep="\t", quote=F, col.names=NA)
deg4 <- FindMarkers(CD4_1_subset.integrated, ident.1="SEVERE_aged", ident.2="HD_aged", method="MAST")
write.table(deg4, "CD4_1_subset_severeaged-hdaged_mast.xls", sep="\t", quote=F, col.names=NA)

Idents(CD4_1_subset.integrated) <- "covsev_death"
deg5 <- FindMarkers(CD4_1_subset.integrated, ident.1="SEVERE_aged", ident.2="HD_aged", method="MAST")
write.table(deg5, "CD4_1_subset_severeaged-hdaged_mast_covsevdeath.xls", sep="\t", quote=F, col.names=NA)
deg6 <- FindMarkers(CD4_1_subset.integrated, ident.1="DEATH_aged", ident.2="HD_aged", method="MAST")
write.table(deg6, "CD4_1_subset_deathaged-hdaged_mast_covsevdeath.xls", sep="\t", quote=F, col.names=NA)


Idents(cd4) <- "seurat_clusters"
CD4_4_subset.integrated <- SubsetData(object = cd4, ident.use = "6")
DefaultAssay(CD4_4_subset.integrated) <- "RNA"
Idents(CD4_4_subset.integrated) <- "COVID"
deg1 <- FindMarkers(CD4_4_subset.integrated, ident.1="MILD_young", ident.2="HD_young", method="MAST")
write.table(deg1, "CD4_4_subset_mildyoung-hdyoung_mast.xls", sep="\t", quote=F, col.names=NA)
deg2 <- FindMarkers(CD4_4_subset.integrated, ident.1="SEVERE_young", ident.2="HD_young", method="MAST")
write.table(deg2, "CD4_4_subset_severeyoung-hdyoung_mast.xls", sep="\t", quote=F, col.names=NA)
deg3 <- FindMarkers(CD4_4_subset.integrated, ident.1="MILD_aged", ident.2="HD_aged", method="MAST")
write.table(deg3, "CD4_4_subset_mildaged-hdaged_mast.xls", sep="\t", quote=F, col.names=NA)
deg4 <- FindMarkers(CD4_4_subset.integrated, ident.1="SEVERE_aged", ident.2="HD_aged", method="MAST")
write.table(deg4, "CD4_4_subset_severeaged-hdaged_mast.xls", sep="\t", quote=F, col.names=NA)

Idents(CD4_4_subset.integrated) <- "covsev_death"
deg5 <- FindMarkers(CD4_4_subset.integrated, ident.1="SEVERE_aged", ident.2="HD_aged", method="MAST")
write.table(deg5, "CD4_4_subset_severeaged-hdaged_mast_covsevdeath.xls", sep="\t", quote=F, col.names=NA)
deg6 <- FindMarkers(CD4_4_subset.integrated, ident.1="DEATH_aged", ident.2="HD_aged", method="MAST")
write.table(deg6, "CD4_4_subset_deathaged-hdaged_mast_covsevdeath.xls", sep="\t", quote=F, col.names=NA)

#DEG CD8
DefaultAssay(cd81) <- "RNA"
cd81 <- NormalizeData(cd81, verbose = FALSE)

features <- c("B2M", "BAX", "BCL11B", "BST2", "CCL5", "CD160", "CD2", "CD3E", "CD3G", "CD5", "CD6")
features <- c("CD7", "CD74", "CD8A", "CD8B", "CDC42", "CTNNB1", "CYLD", "FOXP1", "HLA-A", "HLA-DPA1", "HLA-DPB1")
features <- c("HLA-E", "HLA-F", "HSPD1", "ID2", "IFNG", "IKZF1", "IKZF3", "IL7R", "IRF1", "ITGA4", "ITGB1")
features <- c("ITGB2", "JMJD6", "LAG3", "LGALS1", "MAP3K8", "MIF", "MSN", "MYH9", "NDFIP1", "NFKBIZ", "PAK2")
features <- c("PTGER4", "PTPRC", "RAC1", "RASGRP1", "RHOH", "RIF1", "RSAD2", "RUNX3", "SATB1", "SELENOK", "SOD1")
features <- c("STAT3", "TCF7", "TNFAIP3", "TNFRSF1B", "TNFSF9", "TRBC1", "TRBC2", "XBP1", "ZBTB1", "ZFP36L1", "ZNF683")
features <- c("S100A10", "CD8A", "ITGAM", "TNFSF9", "CCL5", "FOS", "JUNB", "CD69", "ISG20", "IFITM1")
features <- c("LAG3", "IL7R", "TNF", "IFNG", "HLA-DRA", "IRF9", "STAT3", "CCL4L2", "XCL2", "TRIM22", "IL32")
pdf("stacked_violin_cd8_1_both2.pdf", width=5, height = 20)
StackedVlnPlot(obj = cd81, features = features)
dev.off()



#myeloid
all.integrated.hashtag.Mye <- SubsetData(object = all.integrated.hashtag, ident.use = c("Monocytes CD14", "Monocytes CD16", "DC", "Progenitor CD34"))
all.integrated.hashtag.Mye <- SubsetData(object = all.integrated.hashtag.Mye, ident.remove="5")
Monos_dc.integrated <- SubsetData(object = all.integrated.hashtag.Mye, ident.remove="7")
Monos_only_integrated <- SubsetData(object = all.integrated.hashtag.Mye, ident.remove=c("6", "7", "8"))
allMye.list <- SplitObject(Monos_dc.integrated, split.by = "sample")
for (i in names(allMye.list)) {
  allMye.list[[i]] <- SCTransform(allMye.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
allMye.features <- SelectIntegrationFeatures(object.list = allMye.list, nfeatures = 3000)
allMye.list <- PrepSCTIntegration(object.list = allMye.list, anchor.features = allMye.features)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(allMye.list) == "young_control")

allMye.anchors <- FindIntegrationAnchors(object.list = allMye.list, normalization.method = "SCT", 
                                         anchor.features = allMye.features, reference = reference_dataset, k.filter = 100)
allMye.integrated <- IntegrateData(anchorset = allMye.anchors, normalization.method = "SCT")

allMye.integrated <- RunPCA(object = allMye.integrated, verbose = FALSE)
allMye.integrated <- RunUMAP(object = allMye.integrated, dims = 1:30)
allMye.integrated <- FindNeighbors(allMye.integrated, reduction = "pca", dims = 1:30)
allMye.integrated <- FindClusters(allMye.integrated, resolution = 0.5)
pdf("umap_monos_dc.pdf")
a <- DimPlot(allMye.integrated, reduction= "umap", label = TRUE) + NoLegend()
b <- DimPlot(allMye.integrated, reduction= "umap", label = FALSE, group.by = "covsev_death")
a+b
dev.off()

DefaultAssay(allMye.integrated) <- "RNA"
allMye.integrated <- NormalizeData(allMye.integrated, verbose = FALSE)

pdf("feature_mye_dc.pdf")
FeaturePlot(allMye.integrated, features = c("CD14", "FCGR3A", "FCER1A", "LILRA4", "CLEC4C", "CD1C", "LYZ"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
dev.off()

combined.markers <- FindAllMarkers(allMye.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(combined.markers, "./all.integrated.hashtag_combined_markers_Mono_dc.xls", sep="\t", col.names=NA)

Idents(allMye.integrated) <- "seurat_clusters"  
allMye.integrated <- RenameIdents(allMye.integrated, `0`= "MS1 CD14", `1`= "MS2 CD14", `2`= "Intermediate", 
                                `3`= "Non-classical", `4`= "MS3 CD14", `5`= "MS4 CD14", `6`= "mDC/cDC",`7`="pDC")

allMye.integrated$celltype_subset <- Idents(allMye.integrated)                                      
Idents(allMye.integrated) <- "celltype_subset" 

t <- table(allMye.integrated@meta.data$COVID, allMye.integrated@meta.data$celltype_subset)
write.table(t, "all.integrated.hashtag_distribution_Mono_dc_COVID.xls", sep="\t", quote=F, col.names=NA)

my_levels <- c("pDC", "mDC/cDC", "Non-classical", "Intermediate", "MS4 CD14", "MS3 CD14","MS2 CD14", "MS1 CD14")
Idents(allMye.integrated) <- factor(Idents(allMye.integrated), levels= my_levels)

markers.to.plot <- c("CD14", "S100A8", "VCAN", "CXCL8", "IL1B", "HLA-DRA", "ISG20", "ISG15", "HERC5", "FCGR3A", "FCER1A", "LILRA4", "CLEC4C")
pdf("dotplot_mono_dc.pdf")
DotPlot(allMye.integrated, features = markers.to.plot, dot.scale = 8, cols= "RdYlBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()


pdf("CD64_violin_mono.pdf")
VlnPlot(allMye.integrated, features = "FCGR1A", group.by = "COVID", 
        pt.size = 0, combine = FALSE)
dev.off()



#DC pathway gene bubbles
markers.to.plot <- c("PLAC8", "ISG20" , "MX1", "ISG15", "IFITM1","IFI6", "IFI27",  "STAT1", "IFITM3", "IFITM2", "IFI44L", "TNFAIP3", "IRF1")
markers.to.plot <- c("AREG", "EREG", "SPRY2","HIF1A", "EZR",  "MYADM", "KLF9", "KLF2", "RGCC", "RNF213", "HBB", "S100A8", "TIMP1", "ZFP36L2", "ATP2B1", "DUSP1", "FOSB", "DDIT4")
pdf("dotplot_mdc_thromb.pdf", width=9, height=7)
DotPlot(dc, features = markers.to.plot, dot.scale = 8, cols= "RdBu", split.by="COVID")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

pdf("FCER1A_violin_dc.pdf")
VlnPlot(dc, features = "FCER1A", split.by = "COVID", group.by = "celltype_subset", 
                 pt.size = 0, combine = FALSE)
dev.off()


pdf("CD11b_violin_mono.pdf")
VlnPlot(dc, features = "ITGAM", group.by = "COVID", 
        pt.size = 0, combine = FALSE)
dev.off()



#subset objects to the same cell number
all.list <- SplitObject(monos, split.by = "COVID")
set.seed(seed = 1)
cells.subset.1 <- base::sample(x = colnames(all.list$HD_young), size = 442, replace = F)
HD_young.integrated_subset <- SubsetData(object = all.list$HD_young, cell = cells.subset.1)

merged_subsets <- merge(HD_young.integrated_subset, y=c(HD_aged.integrated_subset, MILD_aged.integrated_subset, MILD_young.integrated_subset, SEVERE_young.integrated_subset, SEVERE_aged.integrated_subset))




#bcells
all.integrated.hashtag.Bpla <- SubsetData(object = all.integrated.hashtag, ident.use = c("B cells", "Plasmablasts"))
all.integrated.hashtag.Bpla <- SubsetData(object = allB.integrated, ident.remove = "7")
allB.list <- SplitObject(all.integrated.hashtag.Bpla, split.by = "sample")
for (i in names(allB.list)) {
  allB.list[[i]] <- SCTransform(allB.list[[i]], vars.to.regress = "percent.mt", verbose = FALSE)
}
allB.features <- SelectIntegrationFeatures(object.list = allB.list, nfeatures = 3000)
allB.list <- PrepSCTIntegration(object.list = allB.list, anchor.features = allB.features)

# This command returns dataset 5.  We can also specify multiple refs. (i.e. c(5,6))
reference_dataset <- which(names(allB.list) == "young_control")

allB.anchors <- FindIntegrationAnchors(object.list = allB.list, normalization.method = "SCT", 
                                         anchor.features = allB.features, reference = reference_dataset)
allB.integrated <- IntegrateData(anchorset = allB.anchors, normalization.method = "SCT")

allB.integrated <- RunPCA(object = allB.integrated, verbose = FALSE)
allB.integrated <- RunUMAP(object = allB.integrated, dims = 1:30)
allB.integrated <- FindNeighbors(allB.integrated, reduction = "pca", dims = 1:30)
allB.integrated <- FindClusters(allB.integrated, resolution = 0.5)
pdf("umap_B_pla.pdf")
a <- DimPlot(allB.integrated, reduction= "umap", label = TRUE) + NoLegend()
b <- DimPlot(allB.integrated, reduction= "umap", label = FALSE, group.by = "covsev_death")
a+b
dev.off()

DefaultAssay(allB.integrated) <- "RNA"
allB.integrated <- NormalizeData(allB.integrated, verbose = FALSE)

pdf("feature_B1.pdf")
FeaturePlot(allB.integrated, features = c("CD79A", "MS4A1", "MZB1", "CD38", "IGHD", "CD27", "NEIL1", "MKI67", "CD74"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
dev.off()

pdf("feature_B2.pdf")
FeaturePlot(allB.integrated, features = c("IGHA1", "IGHG1", "IGHM", "FAS", "IGHD"), min.cutoff = "q9", cols = c("grey90", "darkmagenta"))
dev.off()

combined.markers <- FindAllMarkers(allB.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(combined.markers, "./all.integrated.hashtag_combined_markers_Bpla.xls", sep="\t", col.names=NA)

pdf("violin_B1.pdf")
VlnPlot(allB.integrated, features = c("CD79A", "MS4A1", "MZB1", "CD38", "IGHD", "CD27", "NEIL1", "MKI67", "CD74"), pt.size=0.1)
dev.off()

pdf("violin_B2.pdf")
VlnPlot(allB.integrated, features = c("IGHA1", "IGHG1", "IGHM", "FAS", "IGHD"), pt.size=0.1)
dev.off()


Idents(allB.integrated) <- "seurat_clusters"  
allB.integrated <- RenameIdents(allB.integrated, `0`= "naive B", `1`= "intermediate B", `2`= "memory B", 
                                       `3`= "naive B", `4`= "intermediate B", `5`= "intermediate B", `6`= "GC B",`7`="intermediate B", `8`="plasma B")

                                    
allB.integrated$celltype_subset <- Idents(allB.integrated)                                      
Idents(allB.integrated) <- "celltype_subset"        

my_levels <- c("naive B", "GC B", "intermediate B", "memory B", "plasma B")
Idents(allB.integrated) <- factor(Idents(allB.integrated), levels= rev(my_levels))

markers.to.plot <- c("CD79A", "MS4A1", "CD74", "IGHD", "NEIL1", "FAS", "MZB1", "CD38", "MKI67")
markers.to.plot <- c( "ITGA6","IRF4", "SSR4", "XBP1", "JCHAIN","IGHG3", "IGHM", "IGHD","MKI67","MZB1", "CD38", "SDC1")
pdf("dotplot_plasma.pdf")
DotPlot(plasma, features = markers.to.plot, dot.scale = 8, cols= "RdBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()

t <- table(allB.integrated@meta.data$COVID, allB.integrated@meta.data$celltype_subset)
write.table(t, "all.integrated.hashtag_distribution_B_COVID.xls", sep="\t", quote=F, col.names=NA)

features <- c("CD40",	"NFKB1",	"NFKB2", "TNFAIP3",	"CFLAR",	"EZH2",	"EIF2AK2",	"SEMA7A",	"TRAF4",	"NINJ1",	"PLEK",	"CD70",	"IGLC2",	"HBA2",	"HSP90B1")

pdf("violin_bmemdeg_aged.pdf", width=5, height = 20)
StackedVlnPlot(obj = Bmem, features = features, idents=c("HD_aged", "MILD_aged", "SEVERE_aged"))
dev.off()
pdf("violin_bmemdeg_young.pdf", width=5, height = 20)
StackedVlnPlot(obj = Bmem, features = features, idents=c("HD_young", "MILD_young", "SEVERE_young"))
dev.off()

#DEG analysis
allBmem_subset.integrated <- SubsetData(object = allB.integrated, ident.use = c("GC B", "intermediate B", "memory B"))
Idents(allBmem_subset.integrated) <- "COVID"
deg1 <- FindMarkers(allBmem_subset.integrated, ident.1="MILD_young", ident.2="HD_young", method="MAST")
write.table(deg1, "allBmem_subset_mildyoung-hdyoung_mast.xls", sep="\t", quote=F, col.names=NA)
deg2 <- FindMarkers(allBmem_subset.integrated, ident.1="SEVERE_young", ident.2="HD_young", method="MAST")
write.table(deg2, "allBmem_subset_severeyoung-hdyoung_mast.xls", sep="\t", quote=F, col.names=NA)
deg3 <- FindMarkers(allBmem_subset.integrated, ident.1="MILD_aged", ident.2="HD_aged", method="MAST")
write.table(deg3, "allBmem_subset_mildaged-hdaged_mast.xls", sep="\t", quote=F, col.names=NA)
deg4 <- FindMarkers(allBmem_subset.integrated, ident.1="SEVERE_aged", ident.2="HD_aged", method="MAST")
write.table(deg4, "allBmem_subset_severeaged-hdaged_mast.xls", sep="\t", quote=F, col.names=NA)

Idents(allBmem_subset.integrated) <- "covsev_death"
deg5 <- FindMarkers(allBmem_subset.integrated, ident.1="SEVERE_aged", ident.2="HD_aged", method="MAST")
write.table(deg5, "allBmem_subset_severeaged-hdaged_mast_covsevdeath.xls", sep="\t", quote=F, col.names=NA)
deg6 <- FindMarkers(allBmem_subset.integrated, ident.1="DEATH_aged", ident.2="HD_aged", method="MAST")
write.table(deg6, "allBmem_subset_deathaged-hdaged_mast_covsevdeath.xls", sep="\t", quote=F, col.names=NA)

DefaultAssay(b) <- "RNA"
b <- NormalizeData(b, verbose = FALSE)  

markers.to.plot <- c("IGHM", "IGHG3", "HLA-DRB1", "CD40", "TRAF4", "NFKB1", "NFKB2", "REL", "IFI44L", "HBA2", "HBB")
pdf("dotplot_Bmem_deg.pdf", width=8, height=3)
DotPlot(b, features = markers.to.plot, dot.scale = 8, cols= "RdBu")+ theme(axis.text.x = element_text(size = rel(.5), angle = 45))
dev.off()



                                     