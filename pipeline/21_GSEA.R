# 14_0_GSEA_plot.R
# GSEA Analysis with Seurat Integration
#
#' @description 该脚本实现：
#'              - 加载Seurat对象
#'              - 对指定的细胞类型(manual_L2)执行差异基因分析
#'              - 使用选定通路进行GSEA分析
#'              - 可视化并保存GSEA结果
#'              
#' @param inputs 输入：
#'              - iseurat: Seurat对象(.rds)
#'              - ipathway_ID: 通路ID文件(.xls)
#' 
#' @param outputs 输出：
#'              - pgsea_list: GSEA可视化结果目录
#'              - ogsea_egmt: GSEA结果对象(.RData)
#' 
#' @return 生成每个细胞类型的GSEA分析结果和可视化图表
suppressPackageStartupMessages({
  library(Seurat)
  library(msigdbr)
  library(clusterProfiler)
  library(GseaVis)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(readxl)
  library(ggplot2)
})

# 验证输入参数
validate_inputs <- function(args) {
  message("[1/5] 验证输入参数")
  
  required_args <- c("iseurat", "ipathway_ID", "pgsea_list", "ogsea_egmt")
  missing_args <- setdiff(required_args, names(args))
  
  if (length(missing_args) > 0) {
    stop("缺少必要参数: ", paste(missing_args, collapse = ", "))
  }
  
  # 检查文件是否存在gsea_results
  file_args <- c("iseurat", "ipathway_ID")
  for (arg in file_args) {
    if (!file.exists(args[[arg]])) {
      stop("输入文件不存在: ", args[[arg]])
    }
  }
  
  # 创建输出目录
  output_dirs <- c(args$pgsea_list, dirname(args$ogsea_egmt))
  for (odir in output_dirs) {
    if (!dir.exists(odir)) {
      dir.create(odir, recursive = TRUE)
      message("✓ 创建输出目录: ", odir)
    }
  }
  
  message("✓ 输入验证通过")
}

# 加载数据和计算差异基因
prepare_data <- function(args) {
  message("\n[2/5] 加载数据和计算差异基因")
  
  tryCatch({
    # 加载Seurat对象
    seurat_obj <- readRDS(args$iseurat)
    message("✓ 加载Seurat对象 | 细胞数量: ", ncol(seurat_obj))
    
    # 设置细胞类型标识
    if (!"manual_L2" %in% names(seurat_obj@meta.data)) {
      stop("Seurat对象中不存在 'manual_L2' 元数据列")
    }
    
    Idents(seurat_obj) <- "manual_L2"
    cell_types <- levels(Idents(seurat_obj))
    message("✓ 细胞类型: ", paste(cell_types, collapse = ", "))
    
    # 加载通路ID
    pathway_df <- readxl::read_xls(args$ipathway_ID, sheet = "Sheet3")
    pathway_ids <- pathway_df %>% 
      dplyr::select(ID) %>% 
      unique() %>%
      pull() %>%
      as.character()
    
    if (length(pathway_ids) == 0) {
      stop("未找到有效通路ID")
    }
    
    message("✓ 通路ID数量: ", length(pathway_ids))
    
    return(list(
      seurat_obj = seurat_obj,
      cell_types = cell_types,
      pathway_ids = pathway_ids
    ))
  }, error = function(e) {
    stop("数据准备失败: ", e$message)
  })
}

# 执行差异基因分析和GSEA
run_celltype_gsea <- function(data, output_dir) {
  message("\n[3/5] 执行细胞类型GSEA分析")
  
  tryCatch({
    all_gsea_results <- list()
    
    for (cell_type in data$cell_types) {
      message("\n>>> 分析细胞类型: ", cell_type)
      
      # 计算差异基因
      deg_df <- tryCatch({
        FindMarkers(
          data$seurat_obj,
          ident.1 = cell_type,
          ident.2 = NULL,  # 对比所有其他细胞
          logfc.threshold = 0.1,
          min.pct = 0.1,
          test.use = "wilcox",
          only.pos = FALSE,
          verbose = FALSE
        ) %>%
          rownames_to_column(var = "gene") %>%
          mutate(avg_log2FC = avg_log2FC)  # 确保列名一致
      }, error = function(e) {
        warning("差异基因分析失败 (", cell_type, "): ", e$message)
        return(NULL)
      })
      
      if (is.null(deg_df) || nrow(deg_df) == 0) {
        message("  无显著差异基因，跳过")
        next
      }
      
      # 创建基因列表
      genelist <- deg_df$avg_log2FC
      names(genelist) <- deg_df$gene
      genelist <- sort(genelist, decreasing = TRUE)
      genelist <- genelist[genelist != 0]  # 移除零值
      
      message("  差异基因数量: ", length(genelist))
      
      # 获取基因集
      genesets <- get_genesets(data$pathway_ids, species = "Mus musculus")
      
      if (length(genesets) == 0) {
        warning("  未获取到任何有效基因集")
        next
      }
      
      # 为每个通路执行GSEA
      cell_gsea <- run_gsea_for_cell(genelist, genesets, cell_type, output_dir)
      
      if (!is.null(cell_gsea)) {
        all_gsea_results[[cell_type]] <- cell_gsea
      }
    }
    
    message("✓ 完成 ", length(all_gsea_results), "/", length(data$cell_types), 
            " 个细胞类型的GSEA分析")
    
    return(all_gsea_results)
  }, error = function(e) {
    stop("细胞类型GSEA分析失败: ", e$message)
  })
}

# 获取基因集
get_genesets <- function(pathway_ids, species = "Mus musculus") {
  tryCatch({
    GO_df_ALL <- msigdbr(species = species,
                         category = "C5",
                         subcategory = "GO:BP")
    
    if (nrow(GO_df_ALL) == 0) {
      warning("未获取到任何基因集数据")
      return(NULL)
    }
    
    GO_df <- dplyr::select(GO_df_ALL, gene_symbol, gs_exact_source, gs_name)
    
    # 创建基因集列表
    geneset_list <- lapply(pathway_ids, function(id) {
      GO_SUB <- GO_df %>% 
        dplyr::filter(gs_exact_source == id) %>%
        dplyr::distinct(gene_symbol, .keep_all = TRUE)
      
      if (nrow(GO_SUB) == 0) {
        warning("未找到通路ID对应的基因: ", id)
        return(NULL)
      }
      
      # 返回命名基因向量
      setNames(list(unique(GO_SUB$gene_symbol)), 
               unique(GO_SUB$gs_name))
    })
    
    # 移除空元素并扁平化
    geneset_list <- Filter(Negate(is.null), geneset_list)
    flattened_list <- unlist(geneset_list, recursive = FALSE)
    
    message("✓ 获取 ", length(flattened_list), " 个基因集")
    
    return(flattened_list)
  }, error = function(e) {
    warning("基因集获取失败: ", e$message)
    return(NULL)
  })
}

# 为单个细胞类型执行GSEA
run_gsea_for_cell <- function(genelist, genesets, cell_type, output_dir) {
  tryCatch({
    cell_results <- list()
    success_count <- 0
    
    for (i in seq_along(genesets)) {
      pathway_name <- names(genesets)[i]
      pathway_genes <- genesets[[i]]
      
      # 准备TERM2GENE数据框
      term2gene <- data.frame(
        term = rep(pathway_name, length(pathway_genes)),
        gene = pathway_genes
      )
      
      # 检查基因标识符匹配
      matched_genes <- intersect(names(genelist), term2gene$gene)
      
      if (length(matched_genes) < 5) {
        next  # 跳过匹配不足的通路
      }
      
      # 执行GSEA
      egmt <- tryCatch({
        GSEA(
          geneList = genelist,
          TERM2GENE = term2gene,
          verbose = FALSE,
          nPerm = 1000,
          minGSSize = 5,
          maxGSSize = 500,
          pvalueCutoff = 1,
          eps = 1e-50
        )
      }, error = function(e) {
        return(NULL)
      })
      
      if (is.null(egmt) || nrow(egmt@result) == 0) {
        next
      }
      
      cell_results[[pathway_name]] <- egmt
      success_count <- success_count + 1
      
      # 可视化
      plot_gsea(egmt, pathway_name, cell_type, output_dir)
    }
    
    message("  ✓ ", cell_type, ": 完成 ", success_count, "/", length(genesets), " 个通路的GSEA分析")
    
    return(cell_results)
  }, error = function(e) {
    warning("细胞类型GSEA失败 (", cell_type, "): ", e$message)
    return(NULL)
  })
}

# 可视化GSEA结果
plot_gsea <- function(egmt, pathway_name, cell_type, output_dir) {
  tryCatch({
    # 创建细胞类型目录
    cell_dir <- file.path(output_dir, cell_type)
    if (!dir.exists(cell_dir)) {
      dir.create(cell_dir, recursive = TRUE)
    }
    
    # 生成输出路径
    output_path <- file.path(cell_dir, paste0(pathway_name, "_gsearesult.pdf"))
    
    # 获取富集分数
    es <- egmt@result$enrichmentScore
    
    # 创建绘图
    gsea_plot <- gseaNb(
      object = egmt,
      geneSetID = egmt@result$ID,
      subPlot = 2,
      addPval = TRUE,
      pvalX = 0.75,
      pvalY = ifelse(es > 0, 0.75, 0.25),
      pCol = 'black',
      pvalSize = 4,
      pHjust = 0
    )
    
    # 保存PDF
    ggsave(
      filename = output_path,
      plot = gsea_plot,
      width = 8.5,
      height = 4.7
    )
  }, error = function(e) {
    warning("GSEA可视化失败 (", cell_type, " - ", pathway_name, "): ", e$message)
  })
}

# 保存结果
save_results <- function(gsea_results, output_path) {
  message("\n[4/5] 保存GSEA结果")
   data 
  tryCatch({
    if (length(gsea_results) == 0) {
      warning("无有效GSEA结果可保存")
      return(NULL)
    }
    
    save(gsea_results, file = output_path)
    message("✓ GSEA结果已保存: ", output_path)
    message("  包含的细胞类型数量: ", length(gsea_results))
  }, error = function(e) {
    stop("结果保存失败: ", e$message)
  })
}

# 主分析流程
run_gsea_pipeline <- function(args) {
  # 验证输入
  validate_inputs(args)
  
  # 准备数据
  data <- prepare_data(args)
  
  # 执行GSEA分析
  gsea_results <- run_celltype_gsea(data, args$pgsea_list)
  
  # 保存结果
  save_results(gsea_results, args$ogsea_egmt)
  
  message("\n[SUCCESS] GSEA分析完成!")
}

# 命令行参数解析
if (exists("snakemake")) {
  # Snakemake模式
  args <- list(
    iseurat = snakemake@input$iseurat,
    ipathway_ID = snakemake@input$ipathway_ID,
    pgsea_list = snakemake@params$pgsea_list,
    ogsea_egmt = snakemake@output$ogsea_egmt
  )
  run_gsea_pipeline(args)
} 
