#' 04_1_cellanno_check.R
#' Manual Cell Annotation Verification and Visualization
#' 
#' @description 该脚本加载手动注释的细胞类型信息，更新Seurat对象并进行可视化验证。
#' @param config 包含以下参数的列表：
#'               - seurat_rdata: 输入Seurat对象路径
#'               - manual_annotation: 手动注释文件路径
#' @param outputs 输出文件列表：
#'               - annotation_plot: 注释验证UMAP图路径
#'               - annotated_seurat: 带注释的Seurat对象路径
#' 
#' @return 带手动注释的Seurat对象并保存结果

# 加载依赖包
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(future)
})

#------------------- 核心函数 -------------------#
load_seurat_data <- function(seurat_path) {
  message("[1/4] Loading Seurat object from: ", seurat_path)
  env <- new.env()
  load(seurat_path, envir = env)
  obj_name <- ls(env)[1]
  seurat_obj <- get(obj_name, envir = env)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Loaded object is not a Seurat object")
  }
  
  message("✓ Data loaded. Cells: ", ncol(seurat_obj), 
          ", Clusters: ", length(levels(seurat_obj$seurat_clusters)))
  return(seurat_obj)
}

load_manual_annotations <- function(annotation_path) {
  message("\n[2/4] Loading manual annotations from: ", annotation_path)
  
  if (!file.exists(annotation_path)) {
    stop("Annotation file not found: ", annotation_path)
  }
  
  annotations <- read.csv(annotation_path, header = TRUE)
  
  # 检查必需的列
  required_cols <- c("cluster", "manual_level1")
  missing_cols <- setdiff(required_cols, colnames(annotations))
  
  if (length(missing_cols) > 0) {
    stop("Annotation file missing required columns: ", 
         paste(missing_cols, collapse = ", "))
  }
  
  # 创建命名向量用于重命名
  annotation_vec <- setNames(
    as.character(annotations$manual_level1),
    as.character(annotations$cluster)
  )
  
  message("✓ Annotations loaded. Levels: ", 
          length(unique(annotations$manual_level1)))
  return(annotation_vec)
}

apply_manual_annotations <- function(seurat_obj, annotation_vec) {
  message("\n[3/4] Applying manual annotation_vec to Seurat object")
  
  # 验证聚类ID匹配
  missing_clusters <- setdiff(names(annotation_vec), levels(seurat_obj$seurat_clusters))
  if (length(missing_clusters) > 0) {
    warning("Some clusters in annotation not found in Seurat: ",
            paste(missing_clusters, collapse = ", "))
  }
  
  # 设置当前标识为seurat_clusters
  Idents(seurat_obj) <- "seurat_clusters"
  
  # 重命名聚类
  seurat_obj <- RenameIdents(seurat_obj, annotation_vec)
  
  # 添加新注释到元数据
  seurat_obj$manual_level1 <- Idents(seurat_obj)
  
  message("✓ Annotations applied. New levels: ", 
          length(levels(seurat_obj$manual_level1)))
  return(seurat_obj)
}

visualize_annotations <- function(seurat_obj, output_file) {
  message("\n[4/4] Visualizing manual annotations")
  
  p <- DimPlot(seurat_obj,
               reduction = "umap",
               group.by = "manual_level1",
               label = TRUE,
               label.size = 4,
               repel = TRUE,
               pt.size = 0.5) +
    ggtitle("Manual Cell Type Annotations") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # 保存图形
  ggsave(output_file, 
         plot = p, 
         width = 10, 
         height = 8,
         dpi = 300)
  
  message("✓ Annotation UMAP saved: ", output_file)
}

#------------------- 主流程函数 -------------------#
run_manual_annotation_pipeline <- function(seurat_path, 
                                           annotation_path, 
                                           plot_output, 
                                           seurat_output) {
  # 创建输出目录
  dir.create(dirname(plot_output), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(seurat_output), showWarnings = FALSE, recursive = TRUE)
  
  # 1. 加载Seurat数据
  seurat_obj <- load_seurat_data(seurat_path)
  
  # 2. 加载手动注释
  annotations <- load_manual_annotations(annotation_path)
  
  # 3. 应用注释
  annotated_seurat <- apply_manual_annotations(seurat_obj, annotations)
  
  # 4. 可视化
  visualize_annotations(annotated_seurat, plot_output)
  
  # 5. 保存结果
  save(annotated_seurat, file = seurat_output)
  message("✓ Annotated Seurat object saved: ", seurat_output)
  
  message("\n[SUCCESS] Manual annotation verification completed!")
  return(annotated_seurat)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 配置并行处理
  options(future.globals.maxSize = 500000 * 1024^2)
  plan("multisession", workers = 20)
  
  seurat_path <- as.character(snakemake@input$seurat_path)
  annotation_path <- as.character(snakemake@params$annotation_path)
  plot_output <- as.character(snakemake@output$plot_output)
  seurat_output <- as.character(snakemake@output$seurat_output)
  
  # 执行主流程
  run_manual_annotation_pipeline(
    seurat_path = seurat_path,
    annotation_path = annotation_path,
    plot_output = plot_output,
    seurat_output = seurat_output
  )
}
