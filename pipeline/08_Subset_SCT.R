#' 06_0_subset.R
#' Subclustering Pipeline for Specific Cell Types
#'
#' @description 该脚本对特定细胞类型进行亚群分析，包括质量控制、标准化、整合和聚类
#' @param inputs 输入：
#'              - seurat_rdata: Seurat对象路径
#'              - species: 物种类型（human/mouse）
#'              - dims: PCA降维数
#'              - qc_params: 质量控制参数列表
#' @param outputs 输出：
#'              - figures: 质量控制图
#'              - result1: 处理后的Seurat对象
#' 
#' @return 处理后的Seurat对象并保存结果

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(glmGamPoi)
  library(harmony)
})

#------------------- 核心函数 -------------------#
validate_inputs <- function(seurat_path, species, dims, qc_params) {
  message("[1/10] Validating inputs")
  
  if (!file.exists(seurat_path)) {
    stop("Input Seurat file not found: ", seurat_path)
  }
  
  if (!species %in% c("human", "mouse")) {
    stop("Invalid species. Must be 'human' or 'mouse'")
  }
  
  if (!is.numeric(dims) || dims < 1 || dims > 50) {
    stop("Invalid dims value. Must be between 1-50")
  }
  
  required_qc_params <- c("nFeature_RNA_min", "nFeature_RNA_max", "percent_mt_max")
  missing_params <- setdiff(required_qc_params, names(qc_params))
  
  if (length(missing_params) > 0) {
    stop("Missing QC parameters: ", paste(missing_params, collapse = ", "))
  }
  
  message("✓ Input validation passed")
  message("  Species: ", species)
  message("  PCA dims: ", dims)
  message("  QC params: ", paste(names(qc_params), "=", unlist(qc_params), collapse = ", "))
}

load_seurat_data <- function(seurat_path) {
  message("\n[2/10] Loading Seurat object from: ", seurat_path)
  
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
}

create_subset_object <- function(subcells) {
  message("\n[3/10] Creating subsetted Seurat object")
  
  metaData <- subcells@meta.data
  counts <- LayerData(subcells, assay = "RNA", layer = "counts")
  
  seuratdata <- tryCatch(
    {
      CreateSeuratObject(
        counts = counts,
        meta.data = metaData[, c("orig.ident", "manual_level1" )],  # 保留关键元数据列
        project = "Subclustered"
      )
    },
    error = function(e) {
      stop("Failed to create subset object: ", e$message)
    }
  )
  
  DefaultAssay(seuratdata) <- "RNA"
  message("✓ Subset object created")
  
  return(seuratdata)
}

calculate_qc_metrics <- function(seurat_obj, species) {
  message("\n[4/10] Calculating QC metrics for species: ", species)
  
  tryCatch(
    {
      # 计算线粒体基因百分比
      mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
      seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
      
      # 计算核糖体基因百分比
      ribo_pattern <- ifelse(species == "human", "^RP[SL]", "^Rp[sl]")
      rb_genes <- grep(ribo_pattern, rownames(seurat_obj), value = TRUE)
      
      if (length(rb_genes) == 0) {
        warning("No ribosomal genes found with pattern: ", ribo_pattern)
        seurat_obj[["percent.ribo"]] <- 0
      } else {
        counts <- GetAssayData(seurat_obj, layer = "counts")
        ribo_sums <- Matrix::colSums(counts[rb_genes, ])
        total_counts <- Matrix::colSums(counts)
        percent.ribo <- (ribo_sums / total_counts) * 100
        seurat_obj <- AddMetaData(seurat_obj, percent.ribo, col.name = "percent.ribo")
      }
      
      message("✓ QC metrics calculated")
      return(seurat_obj)
    },
    error = function(e) {
      stop("Failed to calculate QC metrics: ", e$message)
    }
  )
}

generate_qc_plots <- function(seurat_obj, figure_paths) {
  message("\n[5/10] Generating QC plots")
  
  tryCatch({
    # 小提琴图
    p1 <- VlnPlot(seurat_obj, alpha = 0.5,
                  features = c("percent.ribo",
                               "percent.mt",
                               "nFeature_RNA"
                               ),
                  ncol = 3,
                  pt.size = 0.01)
    figure1 <- file.path(figure_paths,"qc_VlnPlot.pdf")
    ggsave(
           figure1, 
           plot = p1, 
           width = 12, 
           height = 6, 
           dpi = 300,
           create.dir=TRUE
           )
    message("✓ QC violin plot saved: ", figure1)
    
    # 散点图
    plot1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p <- plot1 + plot2 + plot_layout(ncol = 2)
    
    figure2 <- file.path(figure_paths,"qc_FeatureScatter.pdf")
    ggsave(figure2, plot = p, width = 10, height = 5, dpi = 300, create.dir=TRUE)
    message("✓ Feature scatter plot saved: ", figure2)
  }, error = function(e) {
    warning("Failed to generate QC plots: ", e$message)
  })
}

filter_low_quality_cells <- function(seurat_obj, qc_params,figure_paths) {
  message("\n[6/10] Filtering low-quality cells")
  
  tryCatch({
    seurat_sub <- subset(seurat_obj,
                         subset = nFeature_RNA > qc_params$nFeature_RNA_min &
                           nFeature_RNA < qc_params$nFeature_RNA_max &
                           percent.mt < qc_params$percent_mt_max)
    
    message("✓ Cells after filtering: ", ncol(seurat_sub), 
            " (removed ", ncol(seurat_obj) - ncol(seurat_sub), " cells)")
    
    # 过滤后QC图
    p <- VlnPlot(seurat_sub,
                 features = c("percent.ribo", "percent.mt", "nFeature_RNA", "nCount_RNA"),
                 ncol = 4, pt.size = 0.01)

    figure3 <- file.path(figure_paths,"post_qc_VlnPlot.pdf")
    ggsave(figure3, plot = p, width = 12, height = 6, dpi = 300,create.dir=TRUE)
    message("✓ Post-filter QC violin plot saved: ", figure3)
    
    plot1 <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(seurat_sub, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    p <- plot1 + plot2 + plot_layout(ncol = 2)

    figure4 <- file.path(figure_paths,"post_qc_FeatureScatter.pdf")
    ggsave(figure4, plot = p, width = 10, height = 5, dpi = 300,create.dir=TRUE)
    message("✓ Post-filter feature scatter plot saved: ", figure4)
    
    return(seurat_sub)
  }, error = function(e) {
    stop("Failed to filter cells: ", e$message)
  })
}

run_sctransform <- function(seurat_obj, species) {
  message("\n[7/10] Running SCTransform normalization")
  
  tryCatch({
    # 重新计算线粒体百分比（过滤后）
    mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern, col.name = "percent.mt")
    
    # 运行SCTransform
    seurat_obj <- SCTransform(
      seurat_obj,
      vars.to.regress = "percent.mt",
      method = "glmGamPoi",
      verbose = FALSE
    )
    
    message("✓ SCTransform completed")
    return(seurat_obj)
  }, error = function(e) {
    stop("SCTransform failed: ", e$message)
  })
}

run_pca_and_elbow <- function(seurat_obj, figure_path) {
  message("\n[8/10] Running PCA and elbow plot")
  
  tryCatch({
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    p0 <- ElbowPlot(seurat_obj, ndims = 50)
    figure <- file.path(figure_paths,"elbowPlot.pdf")
    ggsave(figure, plot = p0, width = 8, height = 6, dpi = 300, create.dir=TRUE)
    message("✓ Elbow plot saved: ", figure)
    return(seurat_obj)
  }, error = function(e) {
    stop("PCA failed: ", e$message)
  })
}

run_harmony_integration <- function(seurat_obj, dims) {
  message("\n[9/10] Running Clustering")
  
  tryCatch({
        if(length(unique(seurat_obj$orig.ident))>1){
        seurat_obj <- RunHarmony(seurat_obj, "orig.ident", assay.use = "SCT")
        seurat_obj <- RunUMAP(
          seurat_obj,
          reduction = "harmony",
          dims = 1:dims,
          seed.use = 10,
          verbose = FALSE
        )
        seurat_obj <- FindNeighbors(
          seurat_obj,
          reduction = "harmony",
          dims = 1:dims,
          verbose = FALSE
        )
        seurat_obj <- FindClusters(
          seurat_obj,
          resolution = seq(0, 1, 0.1),
          random.seed = 10,
          verbose = FALSE
        )
        message("✓ Harmony integration and clustering completed")
    }  else {
        seurat_obj <- RunUMAP(
          seurat_obj,
          dims = 1:dims,
          seed.use = 10,
          verbose = FALSE
        )
        seurat_obj <- FindNeighbors(
          seurat_obj,
          dims = 1:dims,
          verbose = FALSE
        )
        seurat_obj <- FindClusters(
          seurat_obj,
          resolution = seq(0, 1, 0.1),
          random.seed = 10,
          verbose = FALSE
        )
        message("✓ clustering completed")
    }    
    message("✓ Harmony integration and clustering completed")
    return(seurat_obj)
  }, error = function(e) {
    stop("Harmony integration failed: ", e$message)
  })
}

save_results <- function(seurat_obj, output_path) {
  message("\n[10/10] Saving results")
  
  tryCatch({
    output_dir <- dirname(output_path)
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
      message("  Created output directory: ", output_dir)
    }
    
    save(seurat_obj, file = output_path)
    message("✓ Final Seurat object saved: ", output_path)
  }, error = function(e) {
    stop("Failed to save results: ", e$message)
  })
}

#------------------- 主流程函数 -------------------#
run_subclustering_pipeline <- function(nFeature_RNA_min,
                                       nFeature_RNA_max,
                                       percent_mt_max,
                                       seurat_rdata,
                                       species,
                                       dims,
                                       figure_paths,
                                       output_path
                                      ) {
  
  qc_params <- list(
    nFeature_RNA_min = as.numeric(nFeature_RNA_min),
    nFeature_RNA_max = as.numeric(nFeature_RNA_max),
    percent_mt_max = as.numeric(percent_mt_max)
  )
  
  # 增加内存限制
  options(future.globals.maxSize = 500000 * 1024^2)  # 500 GB
  
  # 1. 验证输入
  validate_inputs(seurat_rdata, species, dims, qc_params)
  
  # 2. 加载数据
  subcells <- load_seurat_data(seurat_rdata)
  
  # 3. 创建子集对象
  seuratdata <- create_subset_object(subcells)
  
  # 4. 计算QC指标
  seuratdata <- calculate_qc_metrics(seuratdata, species)
  
  # 5. 生成QC图
  generate_qc_plots(seuratdata, figure_paths)
  
  # 6. 过滤低质量细胞
  seuratdata_filtered <- filter_low_quality_cells(seuratdata, qc_params,figure_paths)
  
  # 7. 运行SCTransform
  seuratdata_norm <- run_sctransform(seuratdata_filtered, species)
  
  # 8. 运行PCA和肘部图
  seuratdata_pca <- run_pca_and_elbow(seuratdata_norm, figure_paths)
  
  # 9. 运行Harmony整合
  seuratdata_final <- run_harmony_integration(seuratdata_pca, dims)
  
  # 10. 保存结果
  save_results(seuratdata_final, output_path)
  
  message("\n[SUCCESS] Subclustering pipeline completed!")
  return(seuratdata_final)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  seurat_rdata = as.character(snakemake@input$seurat_rdata)
  species = as.character(snakemake@params$species)
  dims = as.numeric(snakemake@params$dims)
  nFeature_RNA_min = as.numeric(snakemake@params$nFeature_RNA_min)
  nFeature_RNA_max = as.numeric(snakemake@params$nFeature_RNA_max)
  percent_mt_max = as.numeric(snakemake@params$percent_mt_max)
  figure_paths = as.character(snakemake@params$figure_paths)
  output_path = as.character(snakemake@output$output_path)
  
  # 执行主流程
  run_subclustering_pipeline(
    nFeature_RNA_min=nFeature_RNA_min,
    nFeature_RNA_max=nFeature_RNA_max,
    percent_mt_max=percent_mt_max,
    seurat_rdata=seurat_rdata,
    species=species,
    dims=dims,
    figure_paths=figure_paths,
    output_path=output_path
  )
}
