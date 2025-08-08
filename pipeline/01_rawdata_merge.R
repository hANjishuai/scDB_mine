#' 00_rawdatamerge.R
#' Single-cell RNA-seq Raw Data Merging Pipeline
#' 
#' @description 该脚本用于合并多个样本的原始数据，进行整合分析并去除双细胞。
#' @param config 包含以下参数的列表：
#'               - rawdatadir: 原始数据目录
#'               - projectname: 项目名称
#'               - samples: 样本名称向量
#'               - species: 物种 ("human"/"mouse")
#' @param outputs 输出文件列表：
#'               - integrated_pdf: 整合结果PDF路径
#'               - seurat_rdata: Seurat对象保存路径
#' @return 处理后的Seurat对象并保存结果

# 加载依赖包
.libPaths("/home/test/R/library")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(harmony)
  library(dplyr)
})

#------------------- 核心函数 -------------------#
load_and_merge_data <- function(rawdatadir, samples, projectname) {
  message("\n[1/8] Loading and merging samples from: ", rawdatadir)
    
  sceList <- lapply(seq_along(samples), function(i) {
    pro <- samples[i]
    folder <- file.path(rawdatadir, pro)
    CreateSeuratObject(
      counts = Read10X(folder), 
      project = pro
    )
  })
  
  # 合并样本
  merged_seurat <- merge(
    x = sceList[[1]],
    y = sceList[-1],
    project = projectname
  )
  message("✓ Data merged. Total cells: ", ncol(merged_seurat), 
          ", Features: ", nrow(merged_seurat))
  return(merged_seurat)
}

add_qc_metrics <- function(seurat_obj, species) {
  message("\n[2/8] Calculating QC metrics for species: ", species)
  
  # 计算线粒体基因百分比
  mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  message("✓ QC metrics added: percent.mt")
  return(seurat_obj)
}

run_unintegrated_analysis <- function(seurat_obj) {
  message("\n[3/8] Running unintegrated analysis")
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunHarmony(seurat_obj, "orig.ident")
  seurat_obj <- RunUMAP(
    seurat_obj, 
    dims = 1:30, 
    reduction = "harmony", 
    reduction.name = "umap.unintegrated"
  )
  
  message("✓ Unintegrated analysis completed")
  return(seurat_obj)
}

run_integration <- function(seurat_obj) {
  message("\n[4/8] Running data integration")
  
  seurat_obj <- IntegrateLayers(
    object = seurat_obj, 
    method = CCAIntegration,
    orig.reduction = "harmony", 
    new.reduction = "integrated.cca",
    verbose = FALSE
  )
  seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "integrated.cca", dims = 1:30)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction = "integrated.cca")
  
  message("✓ Data integration completed")
  return(seurat_obj)
}

generate_integration_plots <- function(seurat_obj, output_pdf) {
  message("\n[5/8] Generating integration plots")
  
  umap_ui <- DimPlot(seurat_obj, 
                     reduction = "umap.unintegrated", 
                     group.by = "orig.ident") +
    ggtitle("Unintegrated")
  
  umap_i <- DimPlot(seurat_obj, 
                    reduction = "umap", 
                    group.by = "orig.ident") +
    ggtitle("Integrated")
  
  p <- umap_ui + umap_i
  ggsave(output_pdf, plot = p, width = 12, height = 5)
  message("✓ Integration plots saved to: ", output_pdf)
}

remove_doublets <- function(seurat_obj) {
  message("\n[6/8] Removing doublets using scDblFinder")
  # 使用scDblFinder检测双细胞
  set.seed(123)  # 设置随机种子保证结果可重复

  # 转换为SingleCellExperiment对象
  message("  > Converting to SingleCellExperiment object...")
  sce <- as.SingleCellExperiment(seurat_obj)

  # 使用scDblFinder检测双细胞
  message("  > Running scDblFinder...")
  sce <- scDblFinder(sce)

  # 将结果添加回Seurat对象
  seurat_obj$scDblFinder.class <- sce$scDblFinder.class
  seurat_obj$scDblFinder.score <- sce$scDblFinder.score

  # 提取双细胞信息
  n_doublets <- sum(seurat_obj$scDblFinder.class == "doublet")
  doublet_rate <- n_doublets / ncol(seurat_obj)

  message("  > Detected doublets: ", n_doublets, " (", round(doublet_rate * 100, 1), "%)")

  # 可视化双细胞检测结果
  if (!is.null(seurat_obj@reductions$umap)) {
    print(
      DimPlot(seurat_obj, group.by = "scDblFinder.class", 
              reduction = "umap", pt.size = 0.5) +
        ggtitle("Doublet Detection") +
        theme(plot.title = element_text(hjust = 0.5)) +
        scale_color_manual(values = c("doublet" = "red", "singlet" = "gray"))
    )
  }

  # 过滤双细胞
  message("  > Removing doublets...")
  seurat_obj <- subset(seurat_obj, subset = scDblFinder.class == "singlet")
  
  message("✓ Doublets removed: ", n_doublets, " (", 
          round(doublet_rate * 100, 1), "%)")
  message("✓ Cells after filtering: ", ncol(seurat_obj))


  # 移除临时添加的列
  seurat_obj@meta.data$scDblFinder.class <- NULL
  seurat_obj@meta.data$scDblFinder.score <- NULL
  
  return(seurat_obj)  
}

add_metadata <- function(seurat_obj) {
  message("\n[7/8] Adding metadata")
  
  # 添加条件信息
  seurat_obj$condition <- gsub("_\\d$", "", seurat_obj$orig.ident)
  message("✓ Added condition metadata: ", 
          paste(unique(seurat_obj$condition), collapse = ", "))
  
  return(seurat_obj)
}

#------------------- 主流程函数 -------------------#
process_raw_merge_pipeline <- function(rawdatadir,
                                       samples,
                                       projectname,
                                       species,
                                       integrated_pdf,
                                       seurat_rdata) {

  
  # 1. 加载和合并数据
  seurat_obj <- load_and_merge_data(
    rawdatadir, 
    samples, 
    projectname
  )
  
  # 2. 添加QC指标
  seurat_obj <- add_qc_metrics(seurat_obj, species)
  
  # 3. 未整合分析
  seurat_obj <- run_unintegrated_analysis(seurat_obj)
  
  # 4. 数据整合
  seurat_obj <- run_integration(seurat_obj)
  
  # 5. 生成整合图
  generate_integration_plots(seurat_obj, integrated_pdf)
  
  # 6. 去除双细胞
  seurat_obj <- remove_doublets(seurat_obj)
  
  # 7. 添加元数据
  seurat_obj <- add_metadata(seurat_obj)
  
  # 8. 保存结果
  message("\n[8/8] Saving results to: ", seurat_rdata)
  save(seurat_obj, file = seurat_rdata)
  message("\n✓ Pipeline completed successfully!")
  
  return(seurat_obj)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 从snakemake获取配置
  
  rawdatadir = snakemake@params$rawdatadir
  projectname = snakemake@params$projectname
  samples = snakemake@params$samples
  species = snakemake@params$species
  integrated_pdf = snakemake@params$integrated_pdf
  seurat_rdata = snakemake@params$seurat_rdata
  
  # 确保输出目录存在
  dir.create(dirname(integrated_pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(seurat_rdata), recursive = TRUE, showWarnings = FALSE)
  
  # 执行主流程
  process_raw_merge_pipeline(rawdatadir,
                             samples,
                             projectname,
                             species,
                             integrated_pdf,
                             seurat_rdata)
}

