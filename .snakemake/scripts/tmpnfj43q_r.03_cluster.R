
######## snakemake preamble start (automatically inserted, do not edit) ########
library(methods)
Snakemake <- setClass(
    "Snakemake",
    slots = c(
        input = "list",
        output = "list",
        params = "list",
        wildcards = "list",
        threads = "numeric",
        log = "list",
        resources = "list",
        config = "list",
        rule = "character",
        bench_iteration = "numeric",
        scriptdir = "character",
        source = "function"
    )
)
snakemake <- Snakemake(
    input = list('result_out/02_scRNA_Qc2Sct/Adult_PeritonealCavity/seurat_SCT.Rdata', "seurat_path" = 'result_out/02_scRNA_Qc2Sct/Adult_PeritonealCavity/seurat_SCT.Rdata'),
    output = list('result_out/03_cluster/Adult_PeritonealCavity/seurat_cluster.Rdata', "seurat_rdata" = 'result_out/03_cluster/Adult_PeritonealCavity/seurat_cluster.Rdata'),
    params = list('orig.ident', 'result_out/03_cluster/Adult_PeritonealCavity', 'figure_out/03_cluster/Adult_PeritonealCavity', 0.3, "orig_ident" = 'orig.ident', "result_out_dir" = 'result_out/03_cluster/Adult_PeritonealCavity', "figure_out_dir" = 'figure_out/03_cluster/Adult_PeritonealCavity', "cluster_resolution" = 0.3),
    wildcards = list('Adult_PeritonealCavity', "tissue_all" = 'Adult_PeritonealCavity'),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("tissue" = 'Adult_PeritonealCavity', "input01" = list("rawdatadir" = 'data/{tissue}'), "params01" = list("projectname" = 'B1_{tissue}', "samples" = c('HC_1', 'HC_2'), "species" = 'human'), "output01" = list("integrated_pdf" = 'figure_out/01_rawdata_merge/{tissue}/integrated.pdf', "seurat_rdata" = 'result_out/01_rawdata_merge/{tissue}/seurat_raw.Rdata'), "tissue_non10x" = c('Adult_BoneMarrow', 'Adult_PeripheralBlood', 'Adult_Spleen'), "input01_non" = list("rawdatadir" = 'data/{tissue_non10x}'), "params01_non" = list("projectname_prefix" = 'B1_', "species" = 'mouse'), "output01_non" = list("integrated_pdf" = 'figure_out/01_rawdata_merge/{tissue_non10x}/integrated.pdf', "seurat_rdata" = 'result_out/01_rawdata_merge/{tissue_non10x}/seurat_raw.Rdata'), "tissue_all" = c('Adult_BoneMarrow', 'Adult_PeripheralBlood', 'Adult_Spleen', 'Adult_PeritonealCavity'), "input02" = list("seurat_rdata" = 'result_out/01_rawdata_merge/{tissue_all}/seurat_raw.Rdata'), "params02" = list("species" = 'mouse', "nFeature_RNA_min" = 100L, "nFeature_RNA_max" = 4000L, "mt_max" = 15L, "group_var" = 'orig.ident', "output_pdf_dir" = 'figure_out/02_scRNA_Qc2Sct/{tissue_all}'), "output02" = list("seurat_rdata" = 'result_out/02_scRNA_Qc2Sct/{tissue_all}/seurat_SCT.Rdata'), "input03" = list("seurat_path" = 'result_out/02_scRNA_Qc2Sct/{tissue_all}/seurat_SCT.Rdata'), "params03" = list("orig_ident" = 'orig.ident', "result_out_dir" = 'result_out/03_cluster/{tissue_all}', "figure_out_dir" = 'figure_out/03_cluster/{tissue_all}', "cluster_resolution" = 0.3), "output03" = list("seurat_rdata" = 'result_out/03_cluster/{tissue_all}/seurat_cluster.Rdata')),
    rule = 'cluster',
    bench_iteration = as.numeric(NA),
    scriptdir = '/home/test/workshop/B1/scDB_mine/pipeline',
    source = function(...){
        old_wd <- getwd()
        on.exit(setwd(old_wd), add = TRUE)

        is_url <- grepl("^https?://", snakemake@scriptdir)
        file <- ifelse(is_url, file.path(snakemake@scriptdir, ...), ...)
        if (!is_url) setwd(snakemake@scriptdir)
        source(file)
    }
)


######## snakemake preamble end #########
#' 02_0_cluster.R
#' Single-cell RNA-seq Clustering Pipeline
#' 
#' @description 该脚本用于单细胞数据聚类分析，包括PCA降维、维度选择、Harmony整合、聚类和标记基因识别。
#' @param config 包含以下参数的列表：
#'               - seurat_rdata: 输入Seurat对象路径
#'               - orig_ident: 样本标识顺序
#'               - Condition: 实验条件顺序
#'               - cluster_resolution: 聚类分辨率
#' @param outputs 输出文件列表：
#'               - figures: 结果图路径
#'               - results: 结果数据路径
#' 
#' @return 处理后的Seurat对象并保存结果

# 加载依赖包
.libPaths("~/R/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(future)
  library(ggrepel)
  library(stringr)
  library(dplyr)
  library(patchwork)
})

#------------------- 核心函数 -------------------#
load_seurat_data <- function(seurat_path) {
  message("[1/10] Loading Seurat object from: ", seurat_path)
  load(seurat_path)
  message("✓ Data loaded. Cells: ", ncol(seurat_filtered), 
          ", Genes: ", nrow(seurat_filtered))
  return(seurat_filtered)
}

run_pca_analysis <- function(seurat_filtered) {
  message("\n[2/10] Running PCA analysis")
  seurat_obj <- RunPCA(seurat_filtered, verbose = FALSE, seed.use = 10)
  message("✓ PCA completed. Dimensions: ", 
          length(seurat_obj@reductions$pca))
  return(seurat_obj)
}

visualize_pca_results <- function(seurat_obj, figure_out_dir) {
  message("\n[3/10] Visualizing PCA results")
  if(!dir.exists(figure_out_dir)) dir.create(figure_out_dir)
  # PCA散点图
  p1 <- DimPlot(seurat_obj, reduction = "pca", raster = FALSE)
  figure1 <- file.path(figure_out_dir,"pca_scatter.pdf")
  ggsave(figure1, plot = p1, width = 7.5, height = 6.5,create.dir = TRUE)
  
  # PCA热图
  PCAnum <- length(seurat_obj@reductions$pca)
  figure2 <- file.path(figure_out_dir,"pca_heatmap.pdf")
  pdf(figure2, width = 15, height = 5)
  for (i in 1:ceiling(PCAnum/5)) {
    dims <- if (5*i <= PCAnum) (5*i-4):(5*i) else (5*i-4):PCAnum
    DimHeatmap(seurat_obj, dims = dims, nfeatures = 20, ncol = 5,
               cells = 1000, balanced = TRUE)
  }
  dev.off()
  message("✓ PCA visualizations saved")
}

determine_optimal_pcs <- function(seurat_obj, figure_out_dir) {
  message("\n[4/10] Determining optimal PCs")
  if(!dir.exists(figure_out_dir)) dir.create(figure_out_dir)
  # 计算主成分贡献
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  
  # 确定关键截点
  co1 <- which(cumu > 80 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct)-1)] - pct[2:length(pct)]) > 0.1), 
             decreasing = TRUE)[1] + 1
  pcs <- min(co1, co2, na.rm = TRUE)
  
  # 可视化肘部图
  plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
  p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + 
    geom_vline(xintercept = 80, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5], na.rm = TRUE), color = "grey") +
    labs(title = "PCA Elbow Plot", 
         x = "Cumulative Variance (%)", 
         y = "Individual Variance (%)") +
    theme_bw()
  
  figure3 <- file.path(figure_out_dir,"pca_elbow.pdf")
  ggsave(figure3, plot = p, width = 7, height = 5,create.dir = TRUE)
  message("✓ Optimal PCs determined: ", pcs)
  return(pcs)
}

run_harmony_clustering <- function(seurat_obj, pcs, orig_ident) {
  message("\n[5/10] Running Harmony integration and clustering")
  
  # Harmony整合
  if(length(unique(seurat_obj$orig.ident))>1){
    seurat_obj <- RunHarmony(seurat_obj, "orig.ident", verbose = FALSE)
  
    # 降维
    seurat_obj <- RunUMAP(seurat_obj, 
                          reduction = "harmony",
                          dims = 1:pcs,
                          seed.use = 10)

    seurat_obj <- RunTSNE(seurat_obj, 
                          reduction = "harmony",
                          dims = 1:pcs,
                          seed.use = 10)

    # 聚类
    seurat_obj <- FindNeighbors(seurat_obj,
                                reduction = "harmony",
                                dims = 1:pcs)

    seurat_obj <- FindClusters(seurat_obj, 
                               resolution = seq(0.1, 1.0, 0.2),
                               random.seed = 10)

    # 设置样本顺序
    seurat_obj$orig.ident <- factor(seurat_obj$orig.ident, 
                                    levels = orig_ident)
  
    message("✓ Harmony integration and clustering completed")
  } else {
        # 降维
    seurat_obj <- RunUMAP(seurat_obj, 
                          dims = 1:pcs,
                          seed.use = 10)

    seurat_obj <- RunTSNE(seurat_obj, 
                          dims = 1:pcs,
                          seed.use = 10)

    # 聚类
    seurat_obj <- FindNeighbors(seurat_obj,
                                dims = 1:pcs)

    seurat_obj <- FindClusters(seurat_obj, 
                               resolution = seq(0.1, 1.0, 0.2),
                               random.seed = 10)
    message("✓ clustering completed")
  }
  return(seurat_obj)
}

visualize_umap_results <- function(seurat_obj,figure_out_dir) {
  message("\n[6/10] Visualizing UMAP results")
  # 按聚类分面UMAP
  p_cluster <- DimPlot(seurat_obj, 
                       reduction = "umap",
                       label = TRUE,
                       label.size = 3,
                       repel = TRUE,
                       raster = FALSE) +
    ggtitle("By Cluster")
  figure4 <- file.path(figure_out_dir,"umap_cluster.pdf")
  ggsave(figure4, plot = p_cluster, width = 8, height = 10,create.dir = TRUE)  
  message("✓ UMAP visualizations saved")
}

identify_markers <- function(seurat_obj,result_out_dir) {
  message("\n[7/10] Identifying cluster markers")
  
  # 缩放数据
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  
  # 寻找标记基因
  markers <- FindAllMarkers(seurat_obj, 
                            only.pos = TRUE,
                            min.pct = 0.25, 
                            logfc.threshold = 0.25)
  
  # 保存完整标记列表
  result4 <- file.path(result_out_dir,"DEG.csv")
  write.csv(markers, result4, row.names = FALSE)
  
  # 提取top10标记基因
  top10 <- markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
  
  result5 <- file.path(result_out_dir,"DEG_top10.csv")
  write.csv(top10, result5, row.names = FALSE)
  message("✓ Cluster markers identified. Top markers per cluster saved")
  return(markers)
}

visualize_markers <- function(seurat_obj, markers,figure_out_dir) {
  message("\n[8/10] Visualizing marker genes")
  
  # 高变基因可视化
  top10_first <- head(VariableFeatures(seurat_obj), 10)
  plot1 <- VariableFeaturePlot(seurat_obj)
  plot2 <- LabelPoints(plot = plot1, points = top10_first, repel = TRUE)
  figure6 <- file.path(figure_out_dir,"HVG.pdf")
  ggsave(figure6, plot = plot2, width = 8, height = 6,create.dir = TRUE)
  
  # 标记基因热图
  top10_genes <- markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    pull(gene)
  
  p_heatmap <- DoHeatmap(subset(seurat_obj, downsample = 1000), 
                         features = top10_genes,
                         size = 3) + 
    theme(axis.text.y = element_text(size = 6))
  figure7 <- file.path(figure_out_dir,"heatmap_top10gene.pdf")  
  ggsave(figure7, plot = p_heatmap, width = 12, height = 15,create.dir = TRUE)
  message("✓ Marker visualizations saved")
}

save_cluster_results <- function(seurat_obj, 
                                 cluster_resolution, 
                                 result_out_dir,
                                 seurat_rdata
                                 ) {
  message("\n[9/10] Saving clustering results")
  
  # 保存细胞聚类信息
  cluster_col <- paste0("SCT_snn_res.", cluster_resolution)
  cell_cluster <- data.frame(
    cell_ID = names(seurat_obj$seurat_clusters),
    cell_cluster = seurat_obj[[cluster_col]][, 1]
  )
  results1 <- file.path(result_out_dir,"cell_cluster.csv")
  write.csv(cell_cluster, results1, row.names = FALSE)
  
  # 保存UMAP坐标
  embed_umap <- Embeddings(seurat_obj, "umap")
  results2 <- file.path(result_out_dir,"embed_umap.csv")  
  write.csv(embed_umap, results2)
  
  # 保存Seurat对象
  saveRDS(seurat_obj, file = seurat_rdata)
  
  # 保存每个cluster的标记基因列表
  result4 <- file.path(result_out_dir,"DEG.csv")
  markers <- read.csv(result4, header = TRUE)
  cluster_marker_row <- markers %>%
    group_by(cluster) %>%
    summarise(markers = paste(gene, collapse = ","))
  
  results6 <- file.path(result_out_dir,"cluster_marker.csv")
  write.csv(cluster_marker_row, results6, row.names = FALSE)
  
  # 保存最终数据
  results7 <- file.path(result_out_dir,"seu_marker.RData")
  save(seurat_obj, markers, file = results7)
  
  message("✓ All results saved")
}

#------------------- 主流程函数 -------------------#
run_cluster_pipeline <- function(seurat_path, 
                                 orig_ident="orig.ident", 
                                 result_out_dir, 
                                 cluster_resolution,
                                 figure_out_dir,
                                 seurat_rdata
                                 ) {
  dir.create(dirname(result_out_dir), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(figure_out_dir), recursive = TRUE, showWarnings = FALSE)

  # 1. 加载数据
  seuratdata <- load_seurat_data(seurat_path)
  
  # 2-4. PCA分析
  seuratdata <- run_pca_analysis(seuratdata)
  visualize_pca_results(seuratdata, figure_out_dir)
  best.num <- determine_optimal_pcs(seuratdata, figure_out_dir)
  
  # 5-6. Harmony整合与聚类
  seuratdata <- run_harmony_clustering(seuratdata, best.num, orig_ident)
  visualize_umap_results(seuratdata,figure_out_dir)
  
  # 7-8. 标记基因分析
  markers <- identify_markers(seuratdata,result_out_dir)
  visualize_markers(seuratdata, markers, figure_out_dir)
  
  # 9. 保存结果
  save_cluster_results(seuratdata, cluster_resolution, result_out_dir,seurat_rdata)
  
  message("\n[10/10] Clustering pipeline completed successfully!")
  return(seuratdata)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 配置并行处理
  options(future.globals.maxSize = 500000 * 1024^2)
  plan("multisession", workers = 10)
  
  # 从snakemake获取输入参数
  seurat_path <- as.character(snakemake@input$seurat_path)
  orig_ident <- as.character(snakemake@params$orig_ident)
  result_out_dir <- as.character(snakemake@params$result_out_dir)
  figure_out_dir <- as.character(snakemake@params$figure_out_dir)
  cluster_resolution <- as.numeric(snakemake@params$cluster_resolution)
  seurat_rdata <- as.character(snakemake@output$seurat_rdata)
  
  # 执行主流程
  run_cluster_pipeline(
    seurat_path=seurat_path, 
    orig_ident=orig_ident, 
    result_out_dir=result_out_dir, 
    cluster_resolution=cluster_resolution,
    figure_out_dir=figure_out_dir,
    seurat_rdata=seurat_rdata
  )
}
