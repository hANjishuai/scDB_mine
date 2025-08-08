#' Subset Specific Cell Types from Seurat Object
#'
#' @description 该脚本根据指定的细胞类型从Seurat对象中提取子集
#' @param config 包含以下参数的列表：
#'              - seurat_rdata: 输入Seurat对象路径
#'              - celltype: 要提取的细胞类型
#' @param outputs 输出文件列表：
#'              - result1: 子集化后的Seurat对象路径
#' 
#' @return 子集化后的Seurat对象并保存结果

# 加载依赖包
suppressPackageStartupMessages({
  library(Seurat)
  library(glue)
  library(dplyr)
})

#------------------- 核心函数 -------------------#
validate_inputs <- function(seurat_path, celltype) {
  message("[1/5] Validating inputs")
  
  if (!file.exists(seurat_path)) {
    stop("Input Seurat file not found: ", seurat_path)
  }
  env <- new.env()
  load(seurat_path, envir = env) 
  obj_name <- ls(env)[1]
  seurat_obj <- get(obj_name, envir = env)
  cell_types <- seurat_obj %>% levels()
  info <- paste0(cell_types,collapse=' ')
  message(glue("Please input celltypes name among: \n {info}"))

  if (is.null(celltype) || celltype == "" || !(celltype %in% cell_types)) {
    stop("Cell type must be specified")
  }
 
  message("✓ Input validation passed")
  message("  Cell type: ", celltype)
}

load_seurat_object <- function(seurat_path) {
  message("\n[2/5] Loading Seurat object from: ", seurat_path)
  
  env <- new.env()
  load(seurat_path, envir = env)
  
  obj_name <- ls(env)[1]
  seurat_obj <- get(obj_name, envir = env)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Loaded object is not a Seurat object")
  }
  
  message("✓ Seurat object loaded")
  message("  Cells: ", ncol(seurat_obj))
  message("  Metadata columns: ", paste(names(seurat_obj@meta.data), collapse = ", "))
  
  return(seurat_obj)
}

subset_by_celltype <- function(seurat_obj, celltype) {
  message("\n[3/5] Subsetting by cell type: ", celltype)
  
  # 检查manual_level1列是否存在
  if (!"manual_level1" %in% colnames(seurat_obj@meta.data)) {
    stop("Metadata column 'manual_L1' not found in Seurat object")
  }
  
  # 设置细胞标识
  Idents(seurat_obj) <- "manual_level1"
  
  # 检查指定的细胞类型是否存在
  available_types <- unique(seurat_obj$manual_level1)
  missing_types <- setdiff(celltype, available_types)
  
  if (length(missing_types) > 0) {
    stop("Specified cell type(s) not found: ", paste(missing_types, collapse = ", "))
  }
  
  # 提取子集
  subcells <- subset(seurat_obj, idents = celltype)
  
  # 更新元数据
  subcells@meta.data$seurat_clusters <- subcells@meta.data$manual_level1
  subcells@meta.data$seurat_clusters <- droplevels(
    subcells@meta.data$seurat_clusters,
    exclude = setdiff(
      levels(subcells@meta.data$seurat_clusters),
      unique(subcells@meta.data$seurat_clusters)
    )
  )
  
  message("✓ Subsetting completed")
  message("  Subset cells: ", ncol(subcells))
  message("  Remaining cell types: ", paste(unique(subcells$manual_level1), collapse = ", "))
  
  return(subcells)
}

save_subset_object <- function(subcells, output_path) {
  message("\n[4/5] Saving subsetted Seurat object")
  
  # 创建输出目录
  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) {
    message("  Creating output directory: ", output_dir)
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  save(subcells, file = output_path)
  message("✓ Subset saved: ", output_path)
}

#------------------- 主流程函数 -------------------#
run_cell_subsetting <- function(seurat_path, celltype, output_path) {
  # 1. 验证输入
  validate_inputs(seurat_path, celltype)
  
  # 2. 加载Seurat对象
  seurat_obj <- load_seurat_object(seurat_path)
  
  # 3. 提取子集
  subcells <- subset_by_celltype(seurat_obj, celltype)
  
  # 4. 保存结果
  save_subset_object(subcells, output_path)
  
  message("\n[5/5] [SUCCESS] Cell subsetting completed!")
  return(subcells)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  seurat_path <- as.character(snakemake@input$seurat_path)
  celltype <- as.character(snakemake@params$celltype)
  output_path <- as.character(snakemake@output$output_path)
  
  # 执行主流程
  run_cell_subsetting(
    seurat_path = seurat_path,
    celltype = celltype,
    output_path = output_path
  )
}
