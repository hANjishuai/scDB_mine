
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
    input = list('data/Adult_PeripheralBlood', "rawdatadir" = 'data/Adult_PeripheralBlood'),
    output = list('figure_out/01_rawdata_merge/Adult_PeripheralBlood/integrated.pdf', 'result_out/01_rawdata_merge/Adult_PeripheralBlood/seurat_raw.Rdata', "integrated_pdf" = 'figure_out/01_rawdata_merge/Adult_PeripheralBlood/integrated.pdf', "seurat_rdata" = 'result_out/01_rawdata_merge/Adult_PeripheralBlood/seurat_raw.Rdata'),
    params = list('B1_Adult_PeripheralBlood', 'mouse', "projectname" = 'B1_Adult_PeripheralBlood', "species" = 'mouse'),
    wildcards = list('Adult_PeripheralBlood', "tissue_non10x" = 'Adult_PeripheralBlood'),
    threads = 1,
    log = list(),
    resources = list('tmpdir', "tmpdir" = '/tmp'),
    config = list("tissue" = 'Adult_PeritonealCavity', "input01" = list("rawdatadir" = 'data/{tissue}'), "params01" = list("projectname" = 'B1_{tissue}', "samples" = c('HC_1', 'HC_2'), "species" = 'human'), "output01" = list("integrated_pdf" = 'figure_out/01_rawdata_merge/{tissue}/integrated.pdf', "seurat_rdata" = 'result_out/01_rawdata_merge/{tissue}/seurat_raw.Rdata'), "tissue_non10x" = c('Adult_BoneMarrow', 'Adult_PeripheralBlood', 'Adult_Spleen'), "input01_non" = list("rawdatadir" = 'data/{tissue_non10x}'), "params01_non" = list("projectname_prefix" = 'B1_', "species" = 'mouse'), "output01_non" = list("integrated_pdf" = 'figure_out/01_rawdata_merge/{tissue_non10x}/integrated.pdf', "seurat_rdata" = 'result_out/01_rawdata_merge/{tissue_non10x}/seurat_raw.Rdata'), "tissue_all" = c('Adult_BoneMarrow', 'Adult_PeripheralBlood', 'Adult_Spleen', 'Adult_PeritonealCavity'), "input02" = list("seurat_rdata" = 'result_out/01_rawdata_merge/{tissue_all}/seurat_raw.Rdata'), "params02" = list("species" = 'mouse', "nFeature_RNA_min" = 100L, "nFeature_RNA_max" = 6000L, "mt_max" = 15L, "group_var" = 'orig.ident', "output_pdf_dir" = 'figure_out/02_scRNA_Qc2Sct/{tissue_all}'), "output02" = list("seurat_rdata" = 'result_out/02_scRNA_Qc2Sct/{tissue_all}/seurat_SCT.Rdata')),
    rule = 'rawdata_merge_non10x',
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
#' 01_rawdata_merge_non10x.R
#' Single-cell RNA-seq Raw Data Merging Pipeline
#' 
#' @description 该脚本用于合并多个样本的原始数据，进行整合分析并去除双细胞。
#' @param config 包含以下参数的列表：
#'               - rawdatadir: 原始数据目录
#'               - projectname: 项目名称
#'               - species: 物种 ("human"/"mouse")
#' @param outputs 输出文件列表：
#'               - integrated_pdf: 整合结果PDF路径
#'               - seurat_rdata: Seurat对象保存路径
#' @return 处理后的Seurat对象并保存结果

# 加载依赖包
.libPaths("/home/test/R/library")
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(harmony)
  library(dplyr)
})

#------------------- 核心函数 -------------------#
load_and_merge_data <- function(rawdatadir,projectname) {
  message("\n[1/6] Loading and merging samples from: ", rawdatadir)

  expression_file <- list.files(rawdatadir,full.names=T)[1]
  gene_names_file <- list.files(rawdatadir,full.names=T)[2]

  # 读取基因名
  gene_names <- fread(gene_names_file, header = FALSE)$V1
  
  # 读取表达矩阵（快速高效的大文件读取方法）
  expression_data <- fread(
    expression_file,
    header = TRUE,          # 第一行包含细胞名
    colClasses = "integer", # 指定为整数矩阵
    data.table = FALSE      # 返回数据框而非data.table
  )
  
  # 提取细胞名（第一列）
  expression_matrix <- as.matrix(expression_data)
  # 设置行名
  rownames(expression_matrix) <- gene_names
  
  # 转换为稀疏矩阵（节省内存）
  sparse_matrix <- Matrix(expression_matrix, sparse = TRUE)

  # 创建Seurat对象并添加组织来源
  seurat_obj <- CreateSeuratObject(
    counts = sparse_matrix,
    project = projectname,
    min.cells = 3,      # 至少在3个细胞中表达的基因
    min.features = 200   # 至少检测到200个基因的细胞
  )
  
  # 添加元数据
  seurat_obj@meta.data$tissue_origin <- projectname
  
  # 检查结果
  cat("\n成功创建Seurat对象：\n")
  print(seurat_obj)
  cat("\n组织来源：",  projectname, "\n")
  cat("细胞数量：", ncol(seurat_obj), "\n")
  cat("基因数量：", nrow(seurat_obj), "\n")

  message("✓ Data merged. Total cells: ", ncol(seurat_obj), 
          ", Features: ", nrow(seurat_obj))
  return(seurat_obj)
}

add_qc_metrics <- function(seurat_obj, species) {
  message("\n[2/6] Calculating QC metrics for species: ", species)
  
  # 计算线粒体基因百分比
  mt_pattern <- ifelse(species == "human", "^MT-", "^mt-")
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
  
  message("✓ QC metrics added: percent.mt")
  return(seurat_obj)
}

run_unintegrated_analysis <- function(seurat_obj) {
  message("\n[3/6] Running unintegrated analysis")
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- RunUMAP(
    seurat_obj, 
    dims = 1:30, 
  )
  
  message("✓ Unintegrated analysis completed")
  return(seurat_obj)
}

generate_integration_plots <- function(seurat_obj, output_pdf) {
  message("\n[4/6] Generating integration plots")
  
  umap_i <- DimPlot(seurat_obj, 
                    reduction = "umap", 
                    group.by = "orig.ident") +
    ggtitle("Integrated")
  
  p <- umap_i
  ggsave(output_pdf, plot = p, width = 12, height = 5)
  message("✓ Integration plots saved to: ", output_pdf)
}

remove_doublets <- function(seurat_obj) {
  message("\n[5/6] Removing doublets using scDblFinder")
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

#------------------- 主流程函数 -------------------#
process_raw_merge_pipeline <- function(rawdatadir,
                                       projectname,
                                       species,
                                       integrated_pdf,
                                       seurat_rdata) {

  
  # 1. 加载和合并数据
  seurat_obj <- load_and_merge_data(
    rawdatadir, 
    projectname
  )
  
  # 2. 添加QC指标
  seurat_obj <- add_qc_metrics(seurat_obj, species)
  
  # 3. 未整合分析
  seurat_obj <- run_unintegrated_analysis(seurat_obj)

  # 4. 生成umap图
  generate_integration_plots(seurat_obj, integrated_pdf)
  
  # 5. 去除双细胞
  seurat_obj <- remove_doublets(seurat_obj)
  
  # 6. 保存结果
  message("\n[6/6] Saving results to: ", seurat_rdata)
  save(seurat_obj, file = seurat_rdata)
  message("\n✓ Pipeline completed successfully!")
  
  return(seurat_obj)
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {
  # 从snakemake获取配置
  rawdatadir = snakemake@input$rawdatadir
  projectname = snakemake@params$projectname
  species = snakemake@params$species
  integrated_pdf = snakemake@output$integrated_pdf
  seurat_rdata = snakemake@output$seurat_rdata
  
  # 确保输出目录存在
  dir.create(dirname(integrated_pdf), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(seurat_rdata), recursive = TRUE, showWarnings = FALSE)
  
  # 执行主流程
  process_raw_merge_pipeline(rawdatadir,
                             projectname,
                             species,
                             integrated_pdf,
                             seurat_rdata)
}

