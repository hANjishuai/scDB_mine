#' 01_scRNA.QC.R
#' Single-cell RNA-seq Quality Control Pipeline
#' 
#' @description 该脚本用于单细胞数据质量控制，包括QC指标计算、可视化、过滤和标准化。
#' @param config 包含以下参数的列表：
#'               - seurat_rdata: 输入Seurat对象路径
#'               - species: 物种 ("human"/"mouse")
#'               - nFeature_RNA_min: 最小基因数
#'               - nFeature_RNA_max: 最大基因数
#'               - percent_mt_max: 最大线粒体百分比
#'               - group: 分组变量
#' @param outputs 输出文件列表：
#'               - figure1-4: QC图路径
#'               - result1: 处理后数据路径
#' 
#' @return 处理后的Seurat对象并保存结果

# 加载依赖包
.libPaths("/home/test/R/library")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(glmGamPoi)
})

#------------------- 核心函数 -------------------#
load_seurat_data <- function(seurat_obj) {
  message("[1/7] Loading Seurat object from: ", seurat_obj)
  load(seurat_obj)

  message("✓ Data loaded. Cells: ", ncol(seurat_obj), 
          ", Genes: ", nrow(seurat_obj))
  return(seurat_obj)
}

calculate_qc_metrics <- function(seurat_obj, species) {
  message("\n[2/7] Calculating QC metrics for species: ", species)
  
  # 计算线粒体基因百分比
  mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  # 计算核糖体基因百分比
  if (species == "mouse") {
    rb_pattern <- "^Rp[sl]"
    rb_genes <- grep(rb_pattern, rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
    
    if (length(rb_genes) > 0) {
      C <- GetAssayData(seurat_obj, layer = "counts")
      percent.ribo <- Matrix::colSums(C[rb_genes, ]) / Matrix::colSums(C) * 100
      seurat_obj[["percent.ribo"]] <- percent.ribo
    } else {
      warning("No ribosomal genes found with pattern: ", rb_pattern)
      seurat_obj[["percent.ribo"]] <- 0
    }
  } else { # human
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(
      seurat_obj, pattern = "^RP[SL]")
  }
  
  message("✓ QC metrics calculated: percent.mt, percent.ribo")
  return(seurat_obj)
}

plot_qc_metrics <- function(seurat_obj, group_var="orig.ident", output_pdf_dir) {
  message("\n[3/7] Generating QC plots")
  if(dir.exists(output_pdf_dir)) dir.create(output_pdf_dir,recursive= T)
  # 质控前violin图
  p1 <- VlnPlot(seurat_obj, alpha = 0.5,
                features = c("percent.ribo", "percent.mt", 
                             "nFeature_RNA", "nCount_RNA"),
                group.by = group_var,
                ncol = 4, pt.size = 0.00001) +
    ggtitle("Pre-QC Metrics")
  
  # 质控前散点图
  plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- plot1 + plot2 + plot_annotation(title = "Pre-QC Feature Relationships")
  
  # 保存图形
  preQC_violin <- file.path(output_pdf_dir, "preQC_violin.pdf")
  preQC_scatter <- file.path(output_pdf_dir, "preQC_scatter.pdf")
  ggsave(preQC_violin, plot = p1, width = 12, height = 6,create.dir = TRUE)
  ggsave(preQC_scatter, plot = p2, width = 10, height = 5,create.dir = TRUE)
  
  message("✓ Pre-QC plots saved: ", output_pdf_dir)
  return(list(violin = p1, scatter = p2))
}

filter_low_quality_cells <- function(seurat_obj, nFeature_min, nFeature_max, mt_max) {
  message("\n[4/7] Filtering low-quality cells")
  
  pre_filter_cells <- ncol(seurat_obj)
  
  # 应用过滤条件
  seurat_filtered <- subset(seurat_obj,
                           subset = nFeature_RNA > nFeature_min & 
                                    nFeature_RNA < nFeature_max & 
                                    percent.mt < mt_max)
  
  post_filter_cells <- ncol(seurat_filtered)
  removed_cells <- pre_filter_cells - post_filter_cells
  
  message("✓ Cells filtered: ", removed_cells, " (", 
          round(removed_cells/pre_filter_cells*100, 1), "%)")
  message("✓ Cells remaining: ", post_filter_cells)
  
  return(seurat_filtered)
}

plot_post_qc <- function(seurat_filtered, group_var, output_pdf_dir) {
  message("\n[5/7] Generating post-QC plots")
  if(dir.exists(output_pdf_dir)) dir.create(output_pdf_dir,recursive= T)
  # 质控后violin图
  p1 <- VlnPlot(seurat_filtered,
                features = c("percent.ribo", "percent.mt", 
                             "nFeature_RNA", "nCount_RNA"),
                group.by = group_var,
                ncol = 4, pt.size = 0.01) +
    ggtitle("Post-QC Metrics")
  
  # 质控后散点图
  plot1 <- FeatureScatter(seurat_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p2 <- plot1 + plot2 + plot_annotation(title = "Post-QC Feature Relationships")
  
  # 保存图形
  postQC_violin <- file.path(output_pdf_dir, "postQC_violin.pdf")
  postQC_scatter <- file.path(output_pdf_dir, "postQC_scatter.pdf")
  ggsave(postQC_violin, plot = p1, width = 12, height = 6,create.dir = TRUE)
  ggsave(postQC_scatter, plot = p2, width = 10, height = 5,create.dir = TRUE)
  
  message("✓ Post-QC plots saved: ", output_pdf_dir)
}

normalize_data <- function(seurat_filtered) {
  message("\n[6/7] Normalizing data with SCTransform")
  
  # 设置并行处理
  options(future.globals.maxSize = 500000 * 1024^2)
  
  # 执行标准化
  seurat_filtered <- SCTransform(seurat_filtered,
                            vars.to.regress = "percent.mt",
                            verbose = FALSE)
  
  message("✓ Normalization completed. Variable features: ", 
          length(VariableFeatures(seurat_filtered)))
  return(seurat_filtered)
}

save_results <- function(seurat_filtered, output_rdata) {
  if(dir.exists(dirname(output_rdata))) dir.create(dirname(output_rdata),recursive= T)
  message("\n[7/7] Saving results to: ", output_rdata)
  save(seurat_filtered, file = output_rdata)
  message("✓ Results saved successfully!")
}

#------------------- 主流程函数 -------------------#
run_qc_pipeline <- function(seurat_obj,
                            species,
                            group_var,
                            nFeature_RNA_min,
                            nFeature_RNA_max,
                            mt_max,
                            output_rdata,
                            output_pdf_dir) {
  
  # 1. 加载数据
  seurat_obj <- load_seurat_data(seurat_obj)
  
  # 2. 计算QC指标
  seurat_obj <- calculate_qc_metrics(seurat_obj, species)
  
  # 3. 质控前可视化
  plot_qc_metrics(seurat_obj, group_var, output_pdf_dir)
  
  # 4. 过滤低质量细胞
  filtered_data <- filter_low_quality_cells(
    seurat_obj,
    nFeature_min = nFeature_RNA_min,
    nFeature_max = nFeature_RNA_max,
    mt_max = mt_max
  )
  
  # 5. 质控后可视化
  plot_post_qc(filtered_data, group_var, output_pdf_dir)
  
  # 6. 数据标准化
  normalized_data <- normalize_data(filtered_data)
  
  # 7. 保存结果
  save_results(normalized_data, output_rdata)
  
  return(normalized_data)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 从snakemake获取配置
  seurat_obj = as.character(snakemake@input$seurat_obj)
  species = as.character(snakemake@params$species)
  nFeature_RNA_min = as.numeric(snakemake@params$nFeature_RNA_min)
  nFeature_RNA_max = as.numeric(snakemake@params$nFeature_RNA_max)
  mt_max = as.numeric(snakemake@params$mt_max)
  group_var = as.character(snakemake@params$group_var)
  output_pdf_dir = as.character(snakemake@params$output_pdf_dir)
  output_rdata = as.character(snakemake@output$output_rdata)
  
  
  # 确保输出目录存在
  dir.create(dirname(output_pdf_dir), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(output_rdata), recursive = TRUE, showWarnings = FALSE)
  
  # 执行主流程
  final_seurat <- run_qc_pipeline(seurat_obj,
                            species,
                            group_var,
                            nFeature_RNA_min,
                            nFeature_RNA_max,
                            mt_max,
                            output_rdata,
                            output_pdf_dir)
}
