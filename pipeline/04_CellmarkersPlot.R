#' 03_0_cellmarkerplot.R
#' Cell Marker Visualization Pipeline
#' 
#' @description 该脚本用于绘制细胞标记基因表达点图，包括整体点图和分细胞类型点图。
#' @param config 包含以下参数的列表：
#'               - seurat_rdata: 输入Seurat对象路径
#'               - species: 物种类型("human"或"mouse")
#'               - anno_blood_level1: 细胞标记文件路径
#' @param outputs 输出文件列表：
#'               - figure1: 所有标记基因点图路径
#'               - figure2: 分细胞类型点图路径
#' 
#' @return 生成细胞标记基因点图

# 加载依赖包
.libPaths("~/R/library/")
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(readxl)
  library(dplyr)
  library(future)
  library(gdata)
})

#------------------- 核心函数 -------------------#
load_seurat_data <- function(seurat_path) {
  message("[1/5] Loading Seurat object from: ", seurat_path)
  load(seurat_path)
  if (!exists("seurat_obj")) stop("Seurat object 'seurat_obj' not found in RData file")
  message("✓ Data loaded. Cells: ", ncol(seurat_obj), 
          ", Genes: ", nrow(seurat_obj))
  return(seurat_obj)
}

read_cell_markers <- function(marker_file) {
  message("\n[2/5] Reading cell marker file: ", marker_file)
  cellsignal <- read_xlsx(marker_file, sheet = 1, col_names = FALSE)
  
  # 处理标记基因列表
  marker_list <- apply(cellsignal, 1, function(x) {
    genes <- unlist(x[-1])
    genes <- genes[!is.na(genes) & genes != ""]
    unique(genes)
  })
  names(marker_list) <- cellsignal$...1
  
  message("✓ Cell markers loaded. Cell types: ", length(marker_list))
  return(marker_list)
}

plot_all_markers <- function(seurat_obj, marker_list, output_file) {
  message("\n[3/5] Plotting all marker genes")
  
  # 准备所有基因列表
  all_genes <- unlist(marker_list)
  all_genes <- all_genes[!duplicated(all_genes)]
  
  # 绘制点图
  p <- DotPlot(subset(seurat_obj, downsample = min(8000, ncol(seurat_obj))),
               features = all_genes,
               assay = "SCT") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 保存图形
  ggsave(output_file, 
         plot = p, 
         width = 32, 
         height = 15,
         limitsize = FALSE)
  
  message("✓ All markers plot saved: ", output_file)
}

plot_per_celltype <- function(seurat_obj, marker_list, output_file) {
  message("\n[4/5] Plotting markers per cell type")
  
  # 获取缩放数据中的基因
  scaled_genes <- rownames(seurat_obj@assays$SCT@scale.data)
  
  pdf(output_file, width = 8, height = 8)
  for (celltype in names(marker_list)) {
    genes <- marker_list[[celltype]]
    valid_genes <- intersect(scaled_genes, genes)
    
    if (length(valid_genes) > 0) {
      message("  Plotting ", length(valid_genes), " markers for: ", celltype)
      p <- DotPlot(subset(seurat_obj, downsample = min(1000, ncol(seurat_obj))),
                   features = valid_genes,
                   assay = "SCT") +
        labs(title = celltype) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      print(p)
    } else {
      message("  No valid markers found for: ", celltype)
    }
  }
  dev.off()
  
  message("✓ Per-celltype plots saved: ", output_file)
}

#------------------- 主流程函数 -------------------#
run_cellmarker_plot <- function(seurat_path, 
                                marker_file, 
                                all_markers_plot, 
                                per_celltype_plot) {
  # 创建输出目录
  dir.create(dirname(all_markers_plot), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(per_celltype_plot), showWarnings = FALSE, recursive = TRUE)
  
  # 1. 加载数据
  seurat_obj <- load_seurat_data(seurat_path)
  
  # 2. 读取细胞标记
  marker_list <- read_cell_markers(marker_file)
  
  # 3. 绘制所有标记基因点图
  plot_all_markers(seurat_obj, marker_list, all_markers_plot)
  
  # 4. 绘制分细胞类型点图
  plot_per_celltype(seurat_obj, marker_list, per_celltype_plot)
  
  message("\n[5/5] Cell marker visualization completed!")
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 配置并行处理
  options(future.globals.maxSize = 500000 * 1024^2)
  plan("multisession", workers = 5)
  
  # 获取输入参数
  seurat_path <- as.character(snakemake@input$seurat_path)
  marker_file <- as.character(snakemake@params$marker_file)
  all_markers_plot <- as.character(snakemake@output$all_markers_plot)
  per_celltype_plot <- as.character(snakemake@output$per_celltype_plot)
  
  # 执行主流程
  run_cellmarker_plot(
    seurat_path=seurat_path, 
    marker_file=marker_file, 
    all_markers_plot=all_markers_plot, 
    per_celltype_plot=per_celltype_plot
  )
}
