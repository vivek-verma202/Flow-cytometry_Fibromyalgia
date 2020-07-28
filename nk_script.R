# ----------------------------------------------------------- #
#                 GRAHAM
cd /project/6007297/vivek22/FM_flow_cytometry/NK
salloc --time=2:59:0 --mem=150G --account=def-ldiatc --verbose
module load gcc/8.3.0 r/4.0.0
R --no-save
# ----------------------------------------------------------- #
# compensated in Flowjo
## wsp: S:\FM_FLOW_CYTOMETRY\clean_fcs\1_NK{1..4}.wsp (compensation)
## S:\FM_FLOW_CYTOMETRY\clean_fcs\NK_file_selection.wsp
# 01. flowClean
library(flowClean)
f <- list.files(
  path = "./NK/data/fcs/01_compensated",
  pattern = "\\.fcs$", full.names = T
)
for(i in seq_along(1:length(f))){
  ff <- read.FCS(f[i], truncate_max_range = F)
  of <- gsub("./NK/data/fcs/01_compensated/",
             "./NK/data/fcs/tmp/",f[i])
  fc <- clean(fF = ff,
              vectMarkers = c(7:21),
              filePrefixWithDir = "c_",
              ext = "fcs",
              binSize = 0.01,
              nCellCutoff = 500,
              cutoff = "median",
              fcMax = 1.2,
              announce = T,
              diagnostic = F
  )
  write.FCS(fc,of)
}
## good, alive, dump-ve singlets: "./data/fcs/02_cleaned"
## wsp: S:\FM_FLOW_CYTOMETRY\clean_fcs\2_NK_gating.wsp

# 02. Normalized using swiftReg (individual, NDCR only)
## 80000000 cells (concat) used for reference
## normalized fcs: "./NK/data/fcs/03_normalized"

# 03. make flowSet, metadata, panel, and SCE
library("readxl")
library("flowCore")
library("CATALYST")
library("diffcyt")
library("cowplot")
library("ggplot2")
theme_set(theme_bw())
files <- list.files(
  path = "./data/fcs/03_normalized", full.names = T
)
PBMC_fs <- read.flowSet(files,
                        transformation = F,
                        truncate_max_range = F
)
file_name <- as.character(pData(PBMC_fs)$name)
# https://github.com/HelenaLC/CATALYST/issues/103
for (i in 1:length(file_name)) {
  keyword(PBMC_fs@frames[[file_name[i]]])[["$CYT"]] <- "FACS"
}
ID <- gsub("\\.fcs$", "", file_name)
df1 <- data.frame(file_name, ID)
pheno <- read_excel("./data/FM_pheno.xlsx")
PBMC_md <- merge(df1, pheno, by = "ID", all.x = T)
PBMC_md <- PBMC_md[, c(1:4)]
PBMC_md$date <- NULL
for (i in 1:nrow(PBMC_md)) {
  PBMC_md$date[i] <-
    PBMC_fs@frames[[PBMC_md$file_name[i]]]@description[["$DATE"]]
}
names(PBMC_md) <- c("sample_id", "file_name", "condition", "age", "date")
PBMC_md$condition <- factor(PBMC_md$condition,
                            levels = c("Control", "Case")
)
PBMC_md$date <- factor(PBMC_md$date)
PBMC_md <- PBMC_md[, c(2, 1, 3:5)]
PBMC_md$patient_id <- PBMC_md$sample_id
colnames(PBMC_fs) <- gsub("FJComp-", "", gsub("-A$", "", colnames(PBMC_fs)))
fcs_colname <- colnames(PBMC_fs)
antigen <- c(
  "TIGIT", "CD16", "CD57", "CD226", "CD3", "CD56", "CD107a",
  "CD335", "CD159c", "CD158e", "CD314", "CD96", "CD8a", "CD159a"
)
marker_class <- c(
  "state", "type", "type", "state", "type", "type", "state",
  "state", "state", "state", "state", "state", "type", "state"
)
PBMC_panel <- data.frame(fcs_colname, antigen, marker_class,
                         stringsAsFactors = F
)
PBMC_panel <- PBMC_panel[c(5, 13, 6, 2, 3, 8, 14, 9, 11, 10, 7, 1, 12, 4), ]
row.names(PBMC_panel) <- NULL
PBMC_panel
PBMC_md
PBMC_fs
sce <- prepData(PBMC_fs, PBMC_panel, PBMC_md,
                md_cols = list(
                  file = "file_name", id = "sample_id",
                  factors = c("condition", "age", "date")),
                transform = T, cofactor = 150
)
sce@metadata[["experiment_info"]][["age"]] <- as.numeric(as.character(sce@metadata[["experiment_info"]][["age"]]))
## QC
sce # 57576092 cells
summary(sce@metadata[["experiment_info"]][["age"]]) # median = 56
summary(sce@metadata[["experiment_info"]][["condition"]]) # 45 cntrl, 42 FM
summary(sce@metadata[["experiment_info"]][["date"]]) # 13 | 18 | 15 | 41
rm(list = (setdiff(ls(), "sce"))); gc()

# 04. perform flowSOM clustering
sce <- cluster(sce,features = NULL,xdim = 20,
               ydim = 20, maxK = 50,
               verbose = T, seed = 1
); gc()
saveRDS(sce,"./data/nk_sce_clust_all.RDS")
tiff("./plots/meta50_abun.tiff", units="in", width=8, height=8, res=300)
plotAbundances(sce, k = "meta50", by = "cluster_id",
               group_by = "condition") +
  scale_color_manual(values = c("#4BB269", "#BC3A3A")) +
  scale_fill_manual(values = c("#4BB269", "#BC3A3A"))
dev.off()

# 05. explore clusters
tiff("./plots/meta50_hm_all.tiff", units="in", width=10, height=20, res=300)
plotExprHeatmap(sce, by = "cluster_id", k = "meta50",
                scale = "first", col_clust = F, row_clust = F)
dev.off()
tiff("./plots/meta50_exp_norm.tiff", units="in", width=10, height=20, res=300)
plotClusterExprs(sce, k = "meta50", features = "type")
dev.off()
gc()

## add UMAP & TSNE dimensions
set.seed(1)
sce <- runDR(sce, dr = "UMAP", cells = 1000, features = NULL); gc()
sce <- runDR(sce, dr = "TSNE", cells = 1000, features = NULL); gc()
saveRDS(sce,"./data/nk_sce_clust_all.RDS")
tiff("./plots/meta50_dim_red.tiff", units="in", width=20, height=10, res=300)
plot_grid(
  plotDR(sce, dr = "TSNE", color_by = "meta50"),
  plotDR(sce, dr = "UMAP", color_by = "meta50"))
dev.off()
tiff("./plots/pbmds_both.tiff", units="in", width=15, height=10, res=300)
pbMDS(sce, by = "both", k = "meta50",
      shape_by = "condition", size_by = F)
dev.off();gc()
tiff("./plots/pca_clrDR.tiff", units="in", width=15, height=10, res=300)
clrDR(sce, by = "cluster_id", dr = "PCA",
      k = "meta50", arrows = F, size_by = F)
dev.off();gc()
tiff("./plots/mds_clrDR.tiff", units="in", width=15, height=10, res=300)
clrDR(sce, by = "cluster_id", dr = "MDS",
      k = "meta50", arrows = F, size_by = F)
dev.off();gc()
tiff("./plots/umap_clrDR.tiff", units="in", width=15, height=10, res=300)
clrDR(sce, by = "cluster_id", dr = "UMAP",
      k = "meta50", arrows = F, size_by = F)
dev.off();gc()
tiff("./plots/tsne_clrDR.tiff", units="in", width=15, height=10, res=300)
clrDR(sce, by = "cluster_id", dr = "TSNE",
      k = "meta50", arrows = F, size_by = F)
dev.off();gc()
#-------------------------------------------------------------------------
#                            sce <- readRDS("./data/nk_sce_clust_all.RDS")
# 06. Merge clusters
old_cluster <- c(1:50)
new_cluster <- c("undefined","CD57+CD16+NK",
                 "CD8+CD56-CD57-T","CD8+CD56-CD57-T",
                 "CD8+CD56-CD57-T","CD8-CD56-CD57-T",
                 "CD8-CD56-CD57-T","CD8-CD56-CD57+T",
                 "undefined","CD8+CD56-CD57+T",
                 "CD8+CD56-CD57-T","CD8-CD56-CD57-T",
                 "CD8+CD56-CD57-T","undefined",
                 "CD8+CD56-CD57+T","CD8+CD56-CD57-T",
                 "CD8+CD56-CD57-T","CD8-CD56-CD57-T",
                 "CD8+CD56-CD57-T","CD8+CD56-CD57+T",
                 "CD8+CD56-CD57-T","CD8+CD56-CD57-T",
                 "CD8+CD56-CD57-T","CD8+CD56-CD57-T",
                 "CD8+CD56-CD57+T","CD8+CD56-CD57+T",
                 "CD8+CD56-CD57+T","CD8+CD56-CD57+T",
                 "CD8+CD56-CD57+T","CD8+CD56-CD57+T",
                 "CD8-CD56-CD57+T","CD8+CD56+CD57+T",
                 "undefined","CD8+CD56+CD57+T",
                 "CD8-CD56+CD57+T","CD8-CD56+CD57+T",
                 "CD57+CD16+NK","CD57+CD16+NK",
                 "CD8+CD56+CD57+T","CD57+CD16+NK",
                 "CD57-CD16+NK","CD57+CD16+NK",
                 "CD57-CD16+NK","CD56_bri_NK",
                 "undefined","undefined",
                 "CD56_bri_NK","CD56_bri_NK",
                 "CD56-CD16+NK","undefined"
)
m1 <- data.frame(old_cluster,new_cluster)
sce <- mergeClusters(sce, k = "meta50", table = m1, id = "m1", overwrite = T)
rm(list = (setdiff(ls(), "sce"))); gc()

## explore merger:
tiff("./plots/m1_dim_red.tiff", units="in", width=20, height=10, res=300)
plot_grid(
  plotDR(sce, dr = "TSNE", color_by = "meta50"),
  plotDR(sce, dr = "TSNE", color_by = "m1"),
  plotDR(sce, dr = "UMAP", color_by = "meta50"),
  plotDR(sce, dr = "UMAP", color_by = "m1")
)
dev.off()
tiff("./plots/m1_abun.tiff", units="in", width=5, height=10, res=600)
p <- plotAbundances(sce, k = "m1", by = "cluster_id",
                    group_by = "condition")
p + scale_colour_manual(values = c("#5FB84F","#8A0200")) +
  scale_fill_manual(values = c("#5FB84F","#8A0200"))
dev.off()
tiff("./plots/m1_counts.tiff", units="in", width=5, height=5, res=300)
plotCounts(sce,
           group_by = "condition")
dev.off()


saveRDS(sce,"./data/nk_sce_clust_all.RDS")
#-------------------------------------------------------------------------
#                            sce <- readRDS("./data/nk_sce_clust_all.RDS")
# 07. run Diffcyt
design <- createDesignMatrix(ei(sce),
                             cols_design = c("condition","age","date")
)
contrast <- createContrast(c(0,1,0,0,0,0))
## QC
nrow(contrast) == ncol(design)
data.frame(parameters = colnames(design), contrast)

## DA
res_DA <- diffcyt(sce,
                  design = design,
                  contrast = contrast,
                  analysis_type = "DA",
                  clustering_to_use = "m1",
                  verbose = T
)
saveRDS(res_DA,"./data/res_DA_m1.RDS")
res_DS <- diffcyt(sce,
                  design = design,
                  contrast = contrast,
                  analysis_type = "DS",
                  clustering_to_use = "m1",
                  verbose = T
)
saveRDS(res_DS,"./data/res_DS_m1.RDS")
as.data.frame(topTable(res_DA, order = T, all = T,
                       order_by = "p_val", show_all_cols = T))
as.data.frame(topTable(res_DS, order = T, all = T,
                       order_by = "p_val", show_all_cols = T))



res_DA <- readRDS("C:/Users/vverma3/Desktop/Repos/FM_flow_cytometry/res_DA_meta50.RDS")
df_da <- as.data.frame(topTable(res_DA, order = T, all = T,
                                order_by = "p_val", show_all_cols = T))

res_DS <- readRDS("C:/Users/vverma3/Desktop/Repos/FM_flow_cytometry/res_DS_meta50.RDS")
df_ds <- as.data.frame(topTable(res_DS, order = T, all = T,
                                order_by = "p_val", show_all_cols = T))


ds1 <- df_ds[df_ds$cluster_id %in% c(2,9,37,38,40:50) & df_ds$p_val<0.05,]











p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot(fill = c("#8DE887","#8A0200",
                          "#8DE887","#8A0200",
                          "#8DE887","#8A0200",
                          "#8DE887"))



tiff("./plots/m2_counts.tiff", units="in", width=5, height=5, res=300)
plotCounts(sce,
           group_by = "condition",
           color_by = "m2")
dev.off()


## DS
res_DS <- diffcyt(sce,
                  design = design,
                  contrast = contrast,
                  analysis_type = "DS",
                  clustering_to_use = "m3",
                  plot = T,
                  verbose = T
)
as.data.frame(topTable(res_DS, order = T,
                       order_by = "p_val", show_all_cols = T))
write.csv(as.data.frame(topTable(res_DS, all = T, order = T,
                                 order_by = "p_val",
                                 show_all_cols = T,
                                 show_counts = F)),
          "res_ds_m3.csv",row.names = F)

## plots
tiff("./plots/m2_abun.tiff", units="in", width=8, height=8, res=300)
plotAbundances(sce, k = "m2", by = "cluster_id",
               group_by = "condition") +
  scale_color_manual(values = c("#4BB269", "#BC3A3A"))
dev.off()

plotClusterHeatmap(sce,
                   hm2 = "state_markers",  m = NULL,  k = "m1",
                   cluster_anno = F, scale = F, draw_freqs = T)


res_DS <- diffcyt(sce,
                  design = design,
                  contrast = contrast,
                  analysis_type = "DS",
                  clustering_to_use = "m1",
                  verbose = T
)
topTable(res_DS, order = T,order_by = "p_val")


saveRDS(res_DS,"res_DS.RDS")

as.data.frame(topTable(res_DS, order = T,
                       order_by = "p_val", show_all_cols = T))




