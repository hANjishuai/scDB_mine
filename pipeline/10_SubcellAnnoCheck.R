#' 07_0_subcellanno.R
#' Subcluster Annotation and Loom File Preparation
#'
#' @description 该脚本对亚群进行细胞类型注释，并生成pyscenic所需的输入文件
#' @param inputs 输入：
#'              - seurat_rdata: Seurat对象路径
#'              - newcluster: 细胞类型注释文件路径
#'              - res: 聚类分辨率参数
#'              - Condition: 实验条件分组
#' @param outputs 输出：
#'              - figures: UMAP可视化图
#'              - results: 注释后的Seurat对象
#'              - loom_files: pyscenic输入文件（表达矩阵和元数据）
#' 
#' @return 处理后的Seurat对象并保存结果

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

# 增加内存限制
options(future.globals.maxSize = 50000 * 1024^2)  # 50 GB

#------------------- 核心函数 -------------------#
validate_inputs <- function(seurat_path, cluster_path, res_param) {
  message("[1/7] Validating inputs")
  
  # 检查文件存在性
  if (!file.exists(seurat_path)) {
    stop("Seurat object not found: ", seurat_path)
  }
  
  if (!file.exists(cluster_path)) {
    stop("Cluster annotation file not found: ", cluster_path)
  }
  
  # 检查参数有效性
  if (is.null(res_param) || res_param == "") {
    stop("Resolution parameter 'res' must be specified")
  }
  
  
  message("✓ Input validation passed")
  message("  Resolution parameter: ", res_param)
}

load_seurat_data <- function(seurat_path) {
  message("\n[2/7] Loading Seurat object from: ", seurat_path)
  
  tryCatch({
    env <- new.env()
    load(seurat_path, envir = env)
    
    obj_name <- ls(env)[1]
    seurat_obj <- get(obj_name, envir = env)
    
    if (!inherits(seurat_obj, "Seurat")) {
      stop("Loaded object is not a Seurat object")
    }
    
    message("✓ Seurat object loaded")
    message("  Cells: ", ncol(seurat_obj))
    message("  Assays: ", paste(names(seurat_obj@assays), collapse = ", "))
    
    return(seurat_obj)
  }, error = function(e) {
    stop("Failed to load Seurat object: ", e$message)
  })
}

load_cluster_annotations <- function(cluster_path) {
  message("\n[3/7] Loading cluster annotations from: ", cluster_path)
  
  tryCatch({
    cellcluster <- read.csv(cluster_path, header = TRUE)
    
    # 检查必要的列是否存在
    if (!"manual_level2" %in% colnames(cellcluster)) {
      stop("Annotation file missing 'manual_level2' column")
    }
    
    if (!"cluster" %in% colnames(cellcluster)) {
      stop("Annotation file missing 'cluster' column")
    }
    
    # 创建命名向量用于重命名
    newcluster <- cellcluster$manual_level2
    names(newcluster) <- as.character(cellcluster$cluster)
    
    message("✓ Cluster annotations loaded")
    message("  Annotation levels: ", paste(unique(newcluster), collapse = ", "))
    
    return(newcluster)
  }, error = function(e) {
    stop("Failed to load cluster annotations: ", e$message)
  })
}

apply_cluster_annotations <- function(seurat_obj, cluster_mapping, res_param=1) {
  message("\n[4/7] Applying cluster annotations")
  
  tryCatch({
    # 检查分辨率参数是否存在
    res_param = paste0("SCT_snn_res.",res_param)
    if (!res_param %in% colnames(seurat_obj@meta.data)) {
      stop("Resolution parameter '", res_param, "' not found in metadata")
    }
    
    # 设置活动标识
    seurat_obj$seurat_clusters <- seurat_obj[[res_param]][[1]]
    Idents(seurat_obj) <- "seurat_clusters"
    
    # 检查所有聚类是否都有注释
    missing_clusters <- setdiff(levels(Idents(seurat_obj)), names(cluster_mapping))
    
    if (length(missing_clusters) > 0) {
      stop("Missing annotations for clusters: ", paste(missing_clusters, collapse = ", "))
    }
    
    # 应用新注释
    seurat_obj <- RenameIdents(seurat_obj, cluster_mapping)
    
    # 添加到元数据并设置为因子
    seurat_obj$manual_L2 <- Idents(seurat_obj)
    seurat_obj$manual_L2 <- factor(
      seurat_obj$manual_L2,
      levels = unique(cluster_mapping)
    )
    
    message("✓ Annotations applied")
    message("  New cell types: ", paste(levels(seurat_obj$manual_L2), collapse = ", "))
    
    return(seurat_obj)
  }, error = function(e) {
    stop("Failed to apply cluster annotations: ", e$message)
  })
}

visualize_annotations <- function(seurat_obj, output_path) {
  message("\n[5/7] Generating UMAP visualization")
  
  tryCatch({
    cellplot <- DimPlot(
      seurat_obj,
      reduction = "umap",
      label = TRUE,
      pt.size = 0.5,
      label.size = 4
    ) + 
      theme(
        legend.position = "right",
        plot.title = element_text(size = 16, face = "bold")
      ) +
      ggtitle("Cell Type Annotation")
    
    ggsave(
      output_path,
      plot = cellplot,
      width = 12,
      height = 8,
      dpi = 300,create.dir=TRUE
    )
    
    message("✓ UMAP plot saved: ", output_path)
    return(cellplot)
  }, error = function(e) {
    warning("Failed to generate UMAP plot: ", e$message)
    return(NULL)
  })
}

prepare_loom_files <- function(seurat_obj, output_prefix) {
  message("\n[6/7] Preparing loom files for: ", output_prefix)
  
  tryCatch({
    # 1. 表达矩阵
    counts_matrix <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
    
    # 转置为细胞x基因
    transposed_counts <- as.data.frame(as.matrix(t(counts_matrix)))
    
    # 保存表达矩阵
    expression_data <- file.path(output_prefix,"loom_expression.csv")
    write.csv(
      transposed_counts,
      file = expression_data,
      row.names = TRUE
    )
    message("✓ Expression matrix saved: ",expression_data)
    
    # 2. 元数据
    cell_info <- seurat_obj@meta.data[, c("manual_L2","nCount_RNA", "nFeature_RNA")]
    colnames(cell_info) <- c('CellType', 'nGene', 'nUMI')
    
    cell_info_path <- file.path(output_prefix, "loom_metadata.csv")
    write.csv(
      cell_info,
      file = cell_info_path,
      row.names = TRUE
    )
    message("✓ Metadata saved: ", cell_info_path)
    
    return(TRUE)
  }, error = function(e) {
    warning("Failed to prepare loom files for ", output_prefix, ": ", e$message)
    return(FALSE)
  })
}

#------------------- 主流程函数 -------------------#
run_annotation_pipeline <- function(
    seurat_path,
    cluster_path,
    res_param,
    output_path,
    seurat_rdata,
    output_prefix
) {
  # 1. 验证输入
  validate_inputs(seurat_path, cluster_path, res_param)
  
  # 2. 加载Seurat对象
  seurat_obj <- load_seurat_data(seurat_path)
  
  # 3. 加载注释文件
  cluster_mapping <- load_cluster_annotations(cluster_path)
  
  # 4. 应用注释
  seurat_obj <- apply_cluster_annotations(seurat_obj, cluster_mapping, res_param)
  
  # 5. 可视化结果
  visualize_annotations(seurat_obj, output_path)
  
  # 6. 保存注释后的对象
  save(seurat_obj, file = seurat_rdata)
  message("✓ Annotated Seurat object saved: ", seurat_rdata)
  
  # 7. 准备pyscenic文件 - 全数据集
  prepare_loom_files(seurat_obj, output_prefix)

  message("\n[7/7] [SUCCESS] Annotation pipeline completed!")
  return(seurat_obj)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  seurat_path <- as.character(snakemake@input$seurat_path)
  cluster_path <- as.character(snakemake@params$cluster_path)
  res_param <- as.character(snakemake@params$res_param)
  output_prefix <- as.character(snakemake@params$output_prefix)
  output_path <- as.character(snakemake@output$output_path)
  seurat_rdata <- as.character(snakemake@output$seurat_rdata)
 
  # 执行主流程
  run_annotation_pipeline(
    seurat_path = seurat_path,
    cluster_path = cluster_path,
    res_param = res_param,
    output_path = output_path,
    seurat_rdata = seurat_rdata,
    output_prefix = output_prefix
  )
}
