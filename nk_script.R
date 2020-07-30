
## %######################################################%##
#                                                           #
####                        SETUP                        ####
#                                                           #
## %######################################################%##

cd /project/6007297/vivek22/FM_flow_cytometry/NK
# salloc --time=2:59:0 --mem=150G --account=def-ldiatc
# module load gcc/8.3.0 r/4.0.0
R - -no - save
library("flowClean")
library("readxl")
library("flowCore")
library("CATALYST")
library("diffcyt")
library("cowplot")
library("ggplot2")
theme_set(theme_bw())

## %######################################################%##
#                                                           #
####                   01. flowClean                     ####
#                                                           #
## %######################################################%##

# compensated in Flowjo
## wsp: S:\FM_FLOW_CYTOMETRY\clean_fcs\1_NK{1..4}.wsp (compensation)
## S:\FM_FLOW_CYTOMETRY\clean_fcs\NK_file_selection.wsp
f <- list.files(
  path = "./NK/data/fcs/01_compensated",
  pattern = "\\.fcs$", full.names = T
)
for (i in seq_along(1:length(f))) {
  ff <- read.FCS(f[i], truncate_max_range = F)
  of <- gsub(
    "./NK/data/fcs/01_compensated/",
    "./NK/data/fcs/tmp/", f[i]
  )
  fc <- clean(
    fF = ff,
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
  write.FCS(fc, of)
}
## good, alive, dump-ve singlets: "./data/fcs/02_cleaned"
## wsp: S:\FM_FLOW_CYTOMETRY\clean_fcs\2_NK_gating.wsp

##%######################################################%##
#                                                          #
####                 02. Normalization                  ####
#                                                          #
##%######################################################%##

##  Normalized using swiftReg (individual, NDCR only)
## 80000000 cells (concat) used for reference
## normalized fcs: "./NK/data/fcs/03_normalized"
## configuration details: nk_swift.config
## metadata: ./data/swift_metadata.xlsx

## %######################################################%##
#                                                           #
####                  03. Input prep                     ####
#                                                           #
## %######################################################%##

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
ID           <- gsub("\\.fcs$", "", file_name)
df1          <- data.frame(file_name, ID)
pheno        <- read_excel("./data/FM_pheno.xlsx")
PBMC_md      <- merge(df1, pheno, by = "ID", all.x = T)
PBMC_md      <- PBMC_md[, c(1:4)]
PBMC_md$date <- NULL
for (i in 1:nrow(PBMC_md)) {
  PBMC_md$date[i] <-
    PBMC_fs@frames[[PBMC_md$file_name[i]]]@description[["$DATE"]]
}
names(PBMC_md)    <- c("sample_id", "file_name", "condition", "age", "date")
PBMC_md$condition <- factor(PBMC_md$condition,
                            levels = c("Control", "Case")
)
PBMC_md$date       <- factor(PBMC_md$date)
PBMC_md            <- PBMC_md[, c(2, 1, 3:5)]
PBMC_md$patient_id <- PBMC_md$sample_id
colnames(PBMC_fs)  <- gsub("FJComp-", "", gsub("-A$", "", colnames(PBMC_fs)))
fcs_colname        <- colnames(PBMC_fs)
antigen            <- c(
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
  md_cols      = list(
    file       = "file_name", id = "sample_id",
    factors    = c("condition", "age", "date")
  ), transform = T, cofactor = 150
)
sce@metadata[["experiment_info"]][["age"]] <- as.numeric(as.character(sce@metadata[["experiment_info"]][["age"]]))
## QC
sce # 57576092 cells
summary(sce@metadata[["experiment_info"]][["age"]]) # median = 56
summary(sce@metadata[["experiment_info"]][["condition"]]) # 45 cntrl, 42 FM
summary(sce@metadata[["experiment_info"]][["date"]]) # 13 | 18 | 15 | 41
rm(list = (setdiff(ls(), "sce")))
gc()

## %######################################################%##
#                                                           #
####                   04. Clustering                    ####
#                                                           #
## %######################################################%##

sce <- cluster(sce,
  features = NULL, xdim = 20,
  ydim     = 20, maxK = 50,
  verbose  = T, seed = 1
)
gc()
saveRDS(sce, "./data/nk_sce_clust_all.RDS")
tiff("./plots/meta50_abun.tiff", units = "in",
     width = 8, height = 8, res = 300)
plotAbundances(sce,
  k = "meta50", by = "cluster_id",
  group_by = "condition"
) +
  scale_color_manual(values = c("#4BB269", "#BC3A3A")) +
  scale_fill_manual(values = c("#4BB269", "#BC3A3A"))
dev.off()
## explore clusters
tiff("./plots/meta50_hm_all.tiff", units = "in",
     width = 10, height = 20, res = 300)
plotExprHeatmap(sce,
  by = "cluster_id", k = "meta50",
  scale = "first", col_clust = F, row_clust = F
)
dev.off()
tiff("./plots/meta50_exp_norm.tiff", units = "in",
     width = 10, height = 20, res = 300)
plotClusterExprs(sce, k = "meta50", features = "type")
dev.off()
gc()
## add UMAP & TSNE dimensions
set.seed(1)
sce <- runDR(sce, dr = "UMAP", cells = 1000, features = NULL);gc()
sce <- runDR(sce, dr = "TSNE", cells = 1000, features = NULL);gc()
saveRDS(sce, "./data/nk_sce_clust_all.RDS")
tiff("./plots/meta50_dim_red.tiff", units = "in",
     width = 20, height = 10, res = 300)
plot_grid(
  plotDR(sce, dr = "TSNE", color_by = "meta50"),
  plotDR(sce, dr = "UMAP", color_by = "meta50")
)
dev.off();gc()
tiff("./plots/pbmds_both.tiff", units = "in",
     width = 15, height = 10, res = 300)
pbMDS(sce,
  by = "both", k = "meta50",
  shape_by = "condition", size_by = F
)
dev.off();gc()
tiff("./plots/pca_clrDR.tiff", units = "in",
     width = 15, height = 10, res = 300)
clrDR(sce,
  by = "cluster_id", dr = "PCA",
  k = "meta50", arrows = F, size_by = F
)
dev.off();gc()
tiff("./plots/mds_clrDR.tiff", units = "in",
     width = 15, height = 10, res = 300)
clrDR(sce,
  by = "cluster_id", dr = "MDS",
  k = "meta50", arrows = F, size_by = F
)
dev.off();gc()
tiff("./plots/umap_clrDR.tiff", units = "in",
     width = 15, height = 10, res = 300)
clrDR(sce,
  by = "cluster_id", dr = "UMAP",
  k = "meta50", arrows = F, size_by = F
)
dev.off();gc()
tiff("./plots/tsne_clrDR.tiff", units = "in",
     width = 15, height = 10, res = 300)
clrDR(sce,
  by = "cluster_id", dr = "TSNE",
  k = "meta50", arrows = F, size_by = F
)
dev.off();gc()

## %######################################################%##
#                                                           #
####                     05. Merger                      ####
#                                                           #
## %######################################################%##

## merge details: ./data/merge.xlsx
old_cluster <- c(1:50)
new_cluster <- c(
  "CD8+CD56+T"  ,  "CD57+NK",
  "CD8+CD56-T"  ,  "CD8+CD56-T",
  "CD8+CD56-T"  ,  "CD8-CD56-T",
  "CD8-CD56-T"  ,  "CD8-CD56+T",
  "undefined"   ,  "CD8+CD56+T",
  "CD8+CD56-T"  ,  "CD8-CD56-T",
  "CD8+CD56-T"  ,  "CD8-CD56+T",
  "CD8+CD56+T"  ,  "CD8+CD56-T",
  "CD8+CD56-T"  ,  "CD8-CD56-T",
  "CD8+CD56-T"  ,  "CD8+CD56+T",
  "CD8+CD56-T"  ,  "CD8+CD56-T",
  "CD8+CD56-T"  ,  "CD8+CD56-T",
  "CD8+CD56+T"  ,  "CD8+CD56+T",
  "CD8+CD56+T"  ,  "CD8+CD56+T",
  "CD8+CD56+T"  ,  "CD8+CD56+T",
  "CD8-CD56+T"  ,  "CD8+CD56+T",
  "CD8-CD56-T"  ,  "CD8+CD56+T",
  "CD8-CD56+T"  ,  "CD8-CD56+T",
  "CD57+NK"     ,  "CD57+NK",
  "CD8+CD56+T"  ,  "CD57+NK",
  "CD57-NK"     ,  "CD57+NK",
  "CD57-NK"     ,  "CD57-NK",
  "undefined"   ,  "undefined",
  "CD56_bri_NK" ,  "CD57-NK",
  "CD56-CD16+NK",  "undefined"
)
m1  <- data.frame(old_cluster, new_cluster)
sce <- mergeClusters(sce, k = "meta50", table = m1, id = "m1", overwrite = T)
rm(list = (setdiff(ls(), "sce")));gc()
## explore merger:
tiff("./plots/m1_dim_red.tiff", units = "in",
     width = 20, height = 10, res = 300)
plot_grid(
  plotDR(sce, dr = "TSNE", color_by = "meta50"),
  plotDR(sce, dr = "TSNE", color_by = "m1"),
  plotDR(sce, dr = "UMAP", color_by = "meta50"),
  plotDR(sce, dr = "UMAP", color_by = "m1")
)
dev.off();gc()
tiff("./plots/m1_abun_all_samples.tiff", units = "in",
     width = 15, height = 8, res = 600)
plotAbundances(sce, k = "m1", by = "sample_id", group_by = "condition")
dev.off();gc()
p <- plotAbundances(sce,
  k = "m1", by = "cluster_id",
  group_by = "condition"
)
p + scale_colour_manual(values = c("#5FB84F", "#8A0200")) +
  scale_fill_manual(values = c("#5FB84F", "#8A0200"))
dev.off();gc()
saveRDS(sce, "./data/nk_sce_clust_all.RDS")

## %######################################################%##
#                                                           #
####                    06. Diffcyt                      ####
#                                                           #
## %######################################################%##

design   <- createDesignMatrix(ei(sce),
                               cols_design = c("condition", "age", "date")
)
contrast <- createContrast(c(0, 1, 0, 0, 0, 0))
## QC
nrow(contrast) == ncol(design)
data.frame(parameters = colnames(design), contrast)
## DA
res_DA <- diffcyt(sce,
  design            = design,
  contrast          = contrast,
  analysis_type     = "DA",
  clustering_to_use = "m1",
  verbose           = T
)
saveRDS(res_DA, "./data/clean_res_DA_m1.RDS");gc()
res_DS <- diffcyt(sce,
  design            = design,
  contrast          = contrast,
  analysis_type     = "DS",
  clustering_to_use = "m1",
  verbose           = T
)
saveRDS(res_DS, "./data/clean_res_DS_m1.RDS");gc()
da <- as.data.frame(topTable(res_DA,
  order    = T, all = T,
  order_by = "p_val", show_all_cols = T
))
ds <- as.data.frame(topTable(res_DS,
  order    = T, all = T,
  order_by = "p_val", show_all_cols = T
))
da <- da[da$cluster_id != "undefined",-6]
ds <- ds[ds$cluster_id != "undefined",-8]
da <- da[order(da$p_val),]
ds <- ds[order(ds$p_val),]

da$fdr <- p.adjust(da$p_val, method = "fdr")
ds$fdr <- p.adjust(ds$p_val, method = "fdr")
write.csv(da,"./data/da_res.csv",row.names = F)
write.csv(ds,"./data/ds_res.csv",row.names = F)

##%######################################################%##
#                                                          #
####                 07. NK activation                  ####
#                                                          #
##%######################################################%##

# cytoqc to fix channels
library(flowCore)
library(flowWorkspace)
library(cytoqc)
files <- list.files(
  path       = "./NKA/data/fcs/01_compensated",
  full.names = T
)
cqc_data <- cqc_load_fcs(files)
res      <- cqc_check(cqc_data, type = "channel")
res
table(res$channel,res$group_id)
res1 <- cqc_match(res, ref = 1)
res1 <- cqc_match_update(res1, map = c("FJComp-V500-A"="FJComp-BV510-A"),
                         group_id = 2)
res1 <- cqc_match_remove(res1,"FJComp-V500-A",group_id = 2)
res1 <- cqc_match_update(res1, map = c("FJComp-Qdot 605-A"="FJComp-BV605-A"),
                         group_id = 2)
res1 <- cqc_match_remove(res1,"FJComp-Qdot 605-A",group_id = 2)
res1 <- cqc_match_update(res1, map = c("FJComp-QDOT-705-A"="FJComp-BV711-A"),
                         group_id = 2)
res1 <- cqc_match_remove(res1,"FJComp-QDOT-705-A",group_id = 2)
cqc_fix(res1)
cqc_data <- cqc_get_data(res)
cqc_data
cqc_write_fcs(cqc_data,
              "./NKA/data/fcs/01_cleaned")
# fix markers:
rm(list=ls());gc()
files <- list.files(
  path       = "./NKA/data/fcs/01_cleaned",
  full.names = T
)
cs        <- load_cytoset_from_fcs(files)
file_name <- as.character(pData(cs)$name)
for (i in seq_along(1:length(file_name))) {
  write.FCS(cs[[i]],paste0("./NKA/data/fcs/02_cleaned/",
                           file_name[i]))
}

## %######################################################%##
#                                                           #
####                   00.begin here                     ####
#                                                           #
## %######################################################%##

files <- list.files(
  path       = "./NKA/data/fcs/02_cleaned",
  full.names = T
)
# prepare adnka_fs
adnka    <- files[grepl("./NKA/data/fcs/02_cleaned/ADNKA*",
                     files)]
adnka_fs <- read.flowSet(adnka,
                         transformation     = F,
                         truncate_max_range = F
)
adnka_fn <- as.character(pData(adnka_fs)$name)
for (i in seq_along(1:length(adnka_fn))) {
  keyword(adnka_fs@frames[[adnka_fn[i]]])[["$CYT"]] <- "FACS"
}
nka      <- files[grepl("./NKA/data/fcs/02_cleaned/NKA*",
                        files)]
nka_fs   <- read.flowSet(nka,
                         transformation     = F,
                         truncate_max_range = F
)
nka_fn <- as.character(pData(nka_fs)$name)
for (i in seq_along(1:length(nka_fn))) {
  keyword(nka_fs@frames[[nka_fn[i]]])[["$CYT"]] <- "FACS"
}
# adnka_fn      <- c("ADNKA_A_E40_STIM.fcs","ADNKA_B_F9_UNSTIM.fcs")
ID            <- gsub("_.*","",gsub("^ADNKA_._", "",adnka_fn));ID
state         <- gsub(".*[_]([^.]+)[.].*", "\\1",
                 gsub("^ADNKA_._", "",adnka_fn));state
df1           <- data.frame(file_name, ID, state)
pheno         <- read_excel("./data/FM_pheno.xlsx")
pheno         <- read_excel("./data/FM_pheno.xlsx")
adnka_md      <- merge(df1, pheno, by = "ID", all.x = T)
adnka_md      <- adnka_md[, c(1:5)]
adnka_md$date <- NULL
for (i in 1:nrow(adnka_md)) {
  adnka_md$date[i] <-
    adnka_fs@frames[[adnka_md$file_name[i]]]@description[["$DATE"]]
}
names(adnka_md) <- c("sample_id", "file_name", "condition", "age", "date")
adnka_md$condition <- factor(adnka_md$condition,
                            levels = c("Control", "Case")
)
adnka_md$date <- factor(adnka_md$date)
adnka_md <- adnka_md[, c(2, 1, 3:5)]
adnka_md$patient_id <- adnka_md$sample_id
colnames(adnka_fs) <- gsub("FJComp-", "", gsub("-A$", "", colnames(adnka_fs)))
fcs_colname <- colnames(adnka_fs)
# <<---------------------- code set till this point
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
                  factors = c("condition", "age", "date")
                ),
                transform = T, cofactor = 150
)
sce@metadata[["experiment_info"]][["age"]] <- as.numeric(as.character(sce@metadata[["experiment_info"]][["age"]]))
## QC
sce # 57576092 cells
summary(sce@metadata[["experiment_info"]][["age"]]) # median = 56
summary(sce@metadata[["experiment_info"]][["condition"]]) # 45 cntrl, 42 FM
summary(sce@metadata[["experiment_info"]][["date"]]) # 13 | 18 | 15 | 41
rm(list = (setdiff(ls(), "sce")))
gc()






























## %######################################################%##
#                                                           #
####                            EXTRA                    ####
#                                                           #
## %######################################################%##
cd /project/6007297/vivek22/FM_flow_cytometry/NK
salloc --time=2:59:0 --mem=150G --account=def-ldiatc
module load gcc/8.3.0 r/4.0.0
R - -no - save
library("flowClean")
library("readxl")
library("flowCore")
library("CATALYST")
library("diffcyt")
library("cowplot")
library("ggplot2")
theme_set(theme_bw())














p <- ggplot(mpg, aes(class, hwy))
p + geom_boxplot(fill = c(
  "#8DE887", "#8A0200",
  "#8DE887", "#8A0200",
  "#8DE887", "#8A0200",
  "#8DE887"
))

tiff("./plots/m2_counts.tiff", units = "in", width = 5, height = 5, res = 300)
plotCounts(sce,
  group_by = "condition",
  color_by = "m2"
)
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
as.data.frame(topTable(res_DS,
  order = T,
  order_by = "p_val", show_all_cols = T
))
write.csv(as.data.frame(topTable(res_DS,
  all = T, order = T,
  order_by = "p_val",
  show_all_cols = T,
  show_counts = F
)),
"res_ds_m3.csv",
row.names = F
)

## plots
tiff("./plots/m2_abun.tiff", units = "in", width = 8, height = 8, res = 300)
plotAbundances(sce,
  k = "m2", by = "cluster_id",
  group_by = "condition"
) +
  scale_color_manual(values = c("#4BB269", "#BC3A3A"))
dev.off()

plotClusterHeatmap(sce,
  hm2 = "state_markers", m = NULL, k = "m1",
  cluster_anno = F, scale = F, draw_freqs = T
)


res_DS <- diffcyt(sce,
  design = design,
  contrast = contrast,
  analysis_type = "DS",
  clustering_to_use = "m1",
  verbose = T
)
topTable(res_DS, order = T, order_by = "p_val")


saveRDS(res_DS, "res_DS.RDS")

as.data.frame(topTable(res_DS,
  order = T,
  order_by = "p_val", show_all_cols = T
))
