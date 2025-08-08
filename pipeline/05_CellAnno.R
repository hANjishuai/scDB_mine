#' 05_CellAnno.R
#' SingleR Cell Annotation Pipeline
#' 
#' @description 该脚本使用SingleR进行细胞类型注释，包括自动注释和结果可视化。
#' @param config 包含以下参数的列表：
#'               - species: 物种类型("human"或"mouse")
#'               - seurat_rdata: 输入Seurat对象路径
#' @param outputs 输出文件列表：
#'               - figure1: 注释结果UMAP图路径
#'               - result1: 初始注释结果CSV路径
#'               - result2: 带注释的Seurat对象路径
#' 
#' @return 带注释的Seurat对象并保存结果

# 加载依赖包
suppressPackageStartupMessages({
  library(SingleR)
  library(celldex)
  library(Seurat)
  library(ggplot2)
  library(future)
})

#------------------- 核心函数 -------------------#
load_seurat_data <- function(seurat_path) {
  message("[1/6] Loading Seurat object from: ", seurat_path)
  load(seurat_path)
  if (!exists("seurat_obj")) stop("Seurat object 'seurat_obj' not found in RData file")
  message("✓ Data loaded. Cells: ", ncol(seurat_obj), 
          ", Clusters: ", length(levels(Idents(seurat_obj))))
  return(seurat_obj)
}

load_reference_data <- function(species) {
  message("\n[2/6] Loading reference data for species: ", species)
  
  if (species == "human") {
    ref_data <- celldex::HumanPrimaryCellAtlasData()
    message("✓ Using HumanPrimaryCellAtlas reference")
  } else if (species == "mouse") {
    ref_data <- celldex::ImmGenData()
    message("✓ Using ImmGenData reference")
  } else {
    stop("Unsupported species: ", species)
  }
  
  message("  Reference cell types: ", length(unique(ref_data$label.main)))
  return(ref_data)
}

run_singleR_annotation <- function(seurat_obj, ref_data) {
  message("\n[3/6] Running SingleR annotation")
  
  # 提取表达矩阵
  expr_matrix <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")
  
  # 运行SingleR
  singler_results <- SingleR(
    test = expr_matrix,
    clusters = Idents(seurat_obj),
    ref = ref_data,
    assay.type.test = 1,
    labels = ref_data$label.main
  )
  
  message("✓ Annotation completed. Cluster labels assigned.")
  return(singler_results)
}

add_annotations_to_seurat <- function(seurat_obj, singler_results) {
  message("\n[4/6] Adding annotations to Seurat object")
  
  # 创建注释元数据
  cluster_annotations <- data.frame(
    Cluster = rownames(singler_results),
    SingleR_Label = singler_results$labels
  )
  
  # 将注释添加到Seurat对象
  seurat_obj[["SingleR.cluster.labels"]] <- 
    cluster_annotations$SingleR_Label[match(
      seurat_obj$seurat_clusters, 
      cluster_annotations$Cluster
    )]
  
  message("✓ Annotations added to Seurat metadata")
  return(seurat_obj)
}

visualize_annotations <- function(seurat_obj, output_file) {
  message("\n[5/6] Visualizing annotation results")
  
  # 绘制带注释的UMAP
  p <- DimPlot(seurat_obj,
               reduction = "umap", 
               group.by = "SingleR.cluster.labels",
               label = TRUE,
               label.size = 3,
               repel = TRUE,
               pt.size = 0.5) +
    ggtitle("SingleR Cell Type Annotations") +
    theme(legend.position = "right")
  
  # 保存图形
  ggsave(output_file, 
         plot = p, 
         width = 12, 
         height = 10,
         dpi = 300)
  
  message("✓ Annotation UMAP saved: ", output_file)
}

save_annotation_results <- function(seurat_obj, singler_results, 
                                    annotation_csv, seurat_output) {
  message("\n[6/6] Saving annotation results")
  
  # 保存注释结果
  cluster_labels <- data.frame(
    Cluster = rownames(singler_results),
    SingleR_Label = singler_results$labels
  )
  write.csv(cluster_labels, annotation_csv, row.names = FALSE)
  message("✓ Cluster annotations saved: ", annotation_csv)
  
  # 保存带注释的Seurat对象
  save(seurat_obj, file = seurat_output)
  message("✓ Annotated Seurat object saved: ", seurat_output)
  
  # 用户提示
  message("\nNOTE: For manual annotation refinement, please:")
  message("  1. Open the annotation file: ", annotation_csv)
  message("  2. Rename columns: 'Cluster' to 'cluster', 'SingleR_Label' to 'manual_level1'")
  message("  3. Modify 'manual_level1' column as needed for downstream analysis")
}

#------------------- 主流程函数 -------------------#
run_cell_annotation_pipeline <- function(seurat_path, 
                                         species, 
                                         annotation_plot, 
                                         annotation_csv, 
                                         seurat_output) {
  # 创建输出目录
  dir.create(dirname(annotation_plot), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(annotation_csv), showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(seurat_output), showWarnings = FALSE, recursive = TRUE)
  
  # 1. 加载数据
  seurat_obj <- load_seurat_data(seurat_path)
  
  # 2. 加载参考数据集
  ref_data <- load_reference_data(species)
  
  # 3. 运行SingleR注释
  singler_results <- run_singleR_annotation(seurat_obj, ref_data)
  
  # 4. 添加注释到Seurat对象
  annotated_seurat <- add_annotations_to_seurat(seurat_obj, singler_results)
  
  # 5. 可视化注释结果
  visualize_annotations(annotated_seurat, annotation_plot)
  
  # 6. 保存结果
  save_annotation_results(annotated_seurat, singler_results, 
                          annotation_csv, seurat_output)
  
  message("\n[SUCCESS] Cell annotation pipeline completed!")
  return(annotated_seurat)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 配置并行处理
  options(future.globals.maxSize = 500000 * 1024^2)
  plan("multisession", workers = 20)

  seurat_path <- as.character(snakemake@input$seurat_path)
  species <- as.character(snakemake@params$species)
  annotation_plot <- as.character(snakemake@output$annotation_plot)
  annotation_csv <- as.character(snakemake@output$annotation_csv)
  seurat_output <- as.character(snakemake@output$seurat_output)
  
  # 执行主流程
  run_cell_annotation_pipeline(seurat_path = seurat_path, 
                               species = species, 
                               annotation_plot = annotation_plot, 
                               annotation_csv = annotation_csv, 
                               seurat_output = seurat_output
  )
}
