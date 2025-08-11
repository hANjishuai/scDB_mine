# 15_1_monocle2.R
# Monocle2 Pseudotime Analysis Pipeline
#
#' @description 该脚本实现单细胞拟时序分析流程，包括：
#'              - 加载Seurat对象
#'              - 提取所需细胞亚群
#'              - 运行Monocle2拟时序分析
#'              - 保存结果和可视化
#'              
#' @param inputs 输入：
#'              - iSeurat_sub: 输入Seurat对象(.RData)
#'              - pselected_celltype: 关注的细胞类型列表(.xls)
#' 
#' @param outputs 输出：
#'              - oseuratdata: 子集Seurat对象(.RDS)
#'              - otrain_monocle_DEG_df: 时序差异基因结果
#'              - oplot_ordering_genes: 时序基因可视化
#'              - ocds: Monocle2结果对象(.RDS)
#' 
#' @return 生成拟时序分析结果

# 加载必要的库
suppressPackageStartupMessages({
  library(argparse)
  library(Seurat)
  library(monocle)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(ggplotify)
  library(tidydr)
  library(ggrastr)
  library(ggpubr)
})

# 安装缺失包
install_missing_packages <- function() {
  message("[1/5] 检查并安装缺失的R包")
  
  required_packages <- c(
    "AnnotationDbi", "monocle", "ggsci", "tidydr", "ggforce",
    "ggrastr", "ggpubr", "ggplotify", "scales", "ggridges",
    "RColorBrewer", "viridis", "parallel", "readxl"
  )
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message("安装缺失包: ", pkg)
      BiocManager::install(pkg)
    }
  }
  
  # 加载本地函数
  if (!file.exists("src/monocle2_plot.R")) {
    stop("缺少本地函数文件: src/monocle2_plot.R")
  }
  source("src/monocle2_plot.R")
  
  message("✓ 所有依赖包已安装")
}

# 验证输入参数
validate_inputs <- function(args) {
  message("[2/5] 验证输入参数")
  
  required_args <- c(
    "iSeurat_sub", "pmetadat_col", "pcelltype",
    "pselected_celltype", "oseuratdata", "otrain_monocle_DEG_df",
    "oplot_ordering_genes", "ocds"
  )
  
  missing_args <- sapply(required_args, function(arg) {
    if (is.null(args[[arg]])) {
      return(arg)
    }
    return(NULL)
  })
  
  missing_args <- unlist(Filter(Negate(is.null), missing_args))
  
  if (length(missing_args) > 0) {
    stop("缺少必要参数: ", paste(missing_args, collapse = ", "))
  }
  
  # 检查文件是否存在
  if (!file.exists(args$iSeurat_sub)) {
    stop("输入文件不存在: ", args$iSeurat_sub)
  }
  
  if (!file.exists(args$pselected_celltype)) {
    stop("细胞类型文件不存在: ", args$pselected_celltype)
  }
  
  message("✓ 输入验证通过")
}

# 加载和处理Seurat对象
load_and_process_seurat <- function(iSeurat_sub, pmetadat_col, pcelltype, pselected_celltype, oseuratdata) {
  message("[3/5] 加载和处理Seurat对象")
  if(!dir.exists(dirname(oseuratdata))){dir.create(dirname(oseuratdata))}
  tryCatch({
    # 加载Seurat对象
    env <- new.env()
    load(iSeurat_sub, envir = env)
    
    obj_name <- ls(env)[1]
    seurat_obj <- get(obj_name, envir = env)
    if (!exists("seurat_obj")) {
      stop("加载的RData中未找到'seurat_obj'对象")
    }
    
    # 设置细胞标识
    if (!pcelltype %in% colnames(seurat_obj@meta.data)) {
      stop("元数据中不存在指定的细胞类型列: ", pcelltype)
    }
    Idents(seurat_obj) <- pcelltype
    
    # 读取选定的细胞类型
    selected_cells <- read.csv(pselected_celltype) %>% 
      dplyr::pull(selected_celltype) %>% 
      as.character()
    
    # 提取细胞子集
    if (!pmetadat_col %in% colnames(seurat_obj@meta.data)) {
      stop("元数据中不存在指定的子集列: ", pmetadat_col)
    }
    
    seurat_sub <- subset(seurat_obj, subset = !!sym(pmetadat_col) %in% selected_cells)
    
    message("  原始细胞数: ", ncol(seurat_obj))
    message("  子集细胞数: ", ncol(seurat_sub))
    message("  保留细胞类型: ", paste(selected_cells, collapse = ", "))
    
    # 保存子集对象
    saveRDS(seurat_sub, oseuratdata)
    message("✓ Seurat子集已保存: ", oseuratdata)
    
    return(seurat_sub)
  }, error = function(e) {
    stop("Seurat处理失败: ", e$message)
  })
}

# 运行Monocle2分析
run_monocle_analysis <- function(seurat_sub, formats) {
  message("[4/5] 运行Monocle2分析")
  
  tryCatch({
    # 拟时序需要用rnaslot
    DefaultAssay(seurat_sub) <- "RNA"
    SeuratObject <- seurat_sub
    monocle_list <- list()
    
    # 提取表达矩阵
    expr_matrix <- as(LayerData(SeuratObject, assay = "RNA", layer = "counts"), "sparseMatrix")
    
    # 提取表型信息
    p_data <- SeuratObject@meta.data
    
    # 提取基因信息
    f_data <- data.frame(
      gene_short_name = row.names(SeuratObject),
      row.names = row.names(SeuratObject)
    ) 
    
    # 构建CDS对象
    pd <- new('AnnotatedDataFrame', data = p_data)
    fd <- new('AnnotatedDataFrame', data = f_data)
    cds <- newCellDataSet(
      expr_matrix,
      phenoData = pd,
      featureData = fd,
      lowerDetectionLimit = 0.5,
      expressionFamily = negbinomial.size()
    ) 
    
    # 估计size factor和离散度
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    
    # 过滤低表达基因
    cds <- detectGenes(cds, min_expr = 0.1)
    expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
    
    # 差异表达分析
    if(!is.null(formats)){
    residualMFS <- paste0("~", formats)
    diff <- differentialGeneTest(
      cds[expressed_genes,],
      fullModelFormulaStr = residualMFS
    )
    
    deg <- subset(diff, qval < 0.01)
    deg <- deg[order(deg$qval, decreasing = FALSE),]
    monocle_list[["monocle_deg"]] <- deg
    
    # 设置排序基因
    ordergene <- rownames(deg)
    cds <- setOrderingFilter(cds, ordergene)
    }
    # 可视化轨迹构建基因
    plot_ordering_genes <- plot_ordering_genes(cds)
    monocle_list[["plot_ordering_genes"]] <- plot_ordering_genes
    
    # 降维
    message("  正在降维...")
    cds <- reduceDimension(
      cds,
      max_components = 2,
      method = "DDRTree"
    )
    
    # 拟时序排列细胞
    message("  正在排列细胞...")
    cds <- orderCells(cds)
    
    monocle_list[["monocle_result"]] <- cds
    message("✓ Monocle2分析完成")
    
    return(monocle_list)
  }, error = function(e) {
    stop("Monocle2分析失败: ", e$message)
  })
}

# 保存结果
save_monocle_results <- function(results, outputs) {
  message("[5/5] 保存分析结果")
  
  tryCatch({
    # 保存差异基因结果
    deg <- results$monocle_deg
    write.table(
      deg,
      file = outputs$otrain_monocle_DEG_df,
      col.names = TRUE,
      row.names = FALSE,
      sep = "\t",
      quote = FALSE
    )
    message("✓ 时序差异基因已保存: ", outputs$otrain_monocle_DEG_df)
    
    # 保存基因排序图
    plot <- results$plot_ordering_genes
    plot <- as.ggplot(plot)
    ggsave(
      outputs$oplot_ordering_genes,
      plot,
      width = 8,
      height = 6,
      dpi = 300,
      create.dir = T  
    )
    message("✓ 基因排序图已保存: ", outputs$oplot_ordering_genes)
    
    # 保存CDS对象
    cds <- results$monocle_result
    saveRDS(cds, file = outputs$ocds)
    message("✓ CDS对象已保存: ", outputs$ocds)
    
  }, error = function(e) {
    stop("结果保存失败: ", e$message)
  })
}

# 主分析流程
run_monocle_pipeline <- function(args) {
  # 安装依赖
  install_missing_packages()
  
  # 验证输入
  validate_inputs(args)
  
  # 提取输入参数
  outputs <- list(
    oseuratdata = args$oseuratdata,
    otrain_monocle_DEG_df = args$otrain_monocle_DEG_df,
    oplot_ordering_genes = args$oplot_ordering_genes,
    ocds = args$ocds
  )
    
  # 加载和处理Seurat数据
  seurat_sub <- load_and_process_seurat(
    args$iSeurat_sub,
    args$pmetadat_col,
    args$pcelltype,
    args$pselected_celltype,
    args$oseuratdata
  )
  
  # 运行Monocle2分析
  results <- run_monocle_analysis(seurat_sub, formats = args$pmetadat_col)
  
  
  # 保存结果
  save_monocle_results(results, outputs)
  
  message("\n[SUCCESS] Monocle2分析流程完成!")
}

# 命令行参数解析
if (exists("snakemake")) {
  # Snakemake模式
  args <- list(
    iSeurat_sub = snakemake@input$iSeurat_sub,
    pmetadat_col = snakemake@params$pmetadat_col,
    pcelltype = snakemake@params$pcelltype,
    pselected_celltype = snakemake@params$pselected_celltype,
    oseuratdata = snakemake@output$oseuratdata,
    otrain_monocle_DEG_df = snakemake@output$otrain_monocle_DEG_df,
    oplot_ordering_genes = snakemake@output$oplot_ordering_genes,
    ocds = snakemake@output$ocds
  )
  run_monocle_pipeline(args)
}

#args <- list(
#    iSeurat_sub = "result_out/17_Enrichment/Adult_Spleen/seurat_rdata.Rdata",
#    pmetadat_col = "manual_L2",
#    pcelltype = "manual_L2",
#    pselected_celltype = "result_out/17_Enrichment/Adult_Spleen/monocle_selected_cell.csv",
#    oseuratdata = "result_out/18_monocle2/Adult_Spleen/seurat_rdata.Rdata",
#    otrain_monocle_DEG_df = "result_out/18_monocle2/Adult_Spleen/otrain_monocle_DEG_df.xls",
#    oplot_ordering_genes = "result_out/18_monocle2/Adult_Spleen/oplot_ordering_genes.pdf",
#    ocds = "result_out/18_monocle2/Adult_PeritonealCavAdult_Spleenity/ocds.rds"
#)
