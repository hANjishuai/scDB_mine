# 14_1_singlegene_showexpress.R
# 单基因表达热图可视化流程
#
#' @description 该脚本实现单基因在细胞亚群中的表达热图可视化，包括：
#'              - 加载Seurat对象和通路ID列表
#'              - 按细胞亚群计算基因平均表达量
#'              - 绘制基因相对表达热图
#'              
#' @param inputs 输入：
#'              - iSeurat: Seurat对象(.RData)
#'              - ipathway_ID: 通路ID文件(.xls)
#'              - psingle_gene_list: 热图输出目录
#' 
#' @param outputs 输出：
#'              - osingle_gene_expression: 热图对象(.RData)
#' 
#' @return 生成单基因表达热图
suppressPackageStartupMessages({
  library(Seurat)
  library(msigdbr)
  library(ComplexHeatmap)
  library(circlize)
  library(readr)
  library(ggplotify)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(argparse)
  library(readxl)
})

# 验证输入参数
validate_inputs <- function(args) {
  message("\n[1/5] 验证输入参数")
  
  required_args <- c(
    "iSeurat", "ipathway_ID", "osingle_gene_expression",
    "pmin", "pmax", "psize_w", "psize_h", "psingle_gene_list"
  )
  
  missing_args <- setdiff(required_args, names(args))
  if (length(missing_args) > 0) {
    stop("缺少必要参数: ", paste(missing_args, collapse = ", "))
  }
  
  # 检查文件是否存在2
  if (!file.exists(args$iSeurat)) {
    stop("Seurat文件不存在: ", args$iSeurat)
  }
  if (!file.exists(args$ipathway_ID)) {
    stop("通路ID文件不存在: ", args$ipathway_ID)
  }
  
  # 检查数值参数
  num_args <- c("pmin", "pmax", "psize_w", "psize_h")
  for (arg in num_args) {
    if (!is.numeric(args[[arg]]) || is.na(args[[arg]])) {
      stop("参数必须为数值: ", arg)
    }
  }
  
  # 创建输出目录
  output_dirs <- unique(dirname(args$osingle_gene_expression))
  if (!dir.exists(output_dirs)) {
    dir.create(output_dirs, recursive = TRUE)
    message("✓ 创建输出目录: ", output_dirs)
  }
  
  if (!dir.exists(args$psingle_gene_list)) {
    dir.create(args$psingle_gene_list, recursive = TRUE)
    message("✓ 创建热图输出目录: ", args$psingle_gene_list)
  }
  
  message("✓ 输入验证通过")
}

# 加载Seurat对象
load_seurat_data <- function(seurat_path) {
  message("\n[2/5] 加载Seurat对象")
  
  tryCatch({
    seuratdata <- readRDS(seurat_path)
    
    message("✓ Seurat对象加载成功")
    message("  细胞数量: ", ncol(seuratdata))
    message("  基因数量: ", nrow(seuratdata))
    
    return(seuratdata)
  }, error = function(e) {
    stop("Seurat对象加载失败: ", e$message)
  })
}

# 加载通路ID
load_pathway_ids <- function(pathway_file) {
  message("\n[3/5] 加载通路ID")
  
  tryCatch({
    pathway_df <- read_xls(pathway_file, sheet = "Sheet4")
    
    if (!"Special_ID" %in% colnames(pathway_df)) {
      stop("通路文件缺少 'Special_ID' 列")
    }
    
    pathway_ids <- pathway_df %>% 
      dplyr::select(Special_ID) %>% 
      unique() %>% 
      unlist() %>% 
      as.character()
    
    if (length(pathway_ids) == 0) {
      stop("未找到有效通路ID")
    }
    
    message("✓ 加载通路数量: ", length(pathway_ids))
    return(pathway_ids)
  }, error = function(e) {
    stop("通路ID加载失败: ", e$message)
  })
}

# 获取基因集
get_genesets <- function(pathway_list,
                         species = "Mus musculus",
                         category = "C5",
                         subcategory = "GO:BP") {
  tryCatch({
    GO_df_ALL <- msigdbr(species = species,
                         category = category,
                         subcategory = subcategory)
    
    GO_df <- dplyr::select(GO_df_ALL, gene_symbol, gs_exact_source, gs_name)
    
    GO_match <- lapply(pathway_list, function(x) {
      GO_results <- lapply(x, function(i) {
        GO_SUB <- GO_df %>% filter(gs_exact_source == i)
        GO_SUB <- GO_SUB[!duplicated(GO_SUB$gene_symbol), ]
        split(GO_SUB$gene_symbol, GO_SUB$gs_name)
      })
      unlist(GO_results, recursive = FALSE)
    })
    
    unique_genes <- unique(unlist(GO_match))
    message("✓ 获取唯一基因数量: ", length(unique_genes))
    return(unique_genes)
  }, error = function(e) {
    stop("基因集获取失败: ", e$message)
  })
}

# 处理单个细胞类型
process_celltype <- function(subcell, celltype, genesets, args) {
  message("\n处理细胞类型: ", celltype)
  
  tryCatch({
#    # 设置细胞标识和元数据
#    Idents(subcell) <- "condition"
#    subcell$condition <- factor(subcell$condition, 
#                               levels = c("HC", "DLE"),
#                               labels = c("HC", "DLE"))
    # 计算平均表达量
    dat <- AverageExpression(subcell, 
                            assays = "SCT", 
                            slot = "data")[[1]]
    
    # 过滤零表达基因
    dat <- dat[rowSums(dat) > 0, , drop = FALSE]
    dat <- as.matrix(dat)
    
    # 筛选目标基因
    target_genes <- rownames(dat)[rownames(dat) %in% genesets]
    if (length(target_genes) == 0) {
      message("⚠️ 无匹配基因，跳过该细胞类型")
      return(NULL)
    }
    
    heatmap_data <- dat[target_genes, , drop = FALSE]
    
    if (nrow(heatmap_data) == 0) {
      message("⚠️ 无显著变化基因，跳过该细胞类型")
      return(NULL)
    }
    
    message("✓ 有效基因数量: ", nrow(heatmap_data))
    return(heatmap_data)
  }, error = function(e) {
    stop("细胞类型处理失败: ", celltype, " - ", e$message)
  })
}

# 创建热图
create_heatmap <- function(heatmap_data, celltype, args) {
  message("  创建热图...")
  
  tryCatch({
    # 设置颜色映射
    col_fun <- colorRamp2(c(args$pmin, 0, args$pmax), 
                         c("navy", "#e8f6e5", "#d43627"))
    
    # 自动调整行名字体大小
    font_size <- ifelse(nrow(heatmap_data) > 50, 6, 8)
    
    # 创建热图
    heatmap_obj <- Heatmap(
      scale(heatmap_data),
      row_names_max_width = unit(12, "cm"),
      cluster_columns = FALSE,
      cluster_rows = TRUE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontface = "bold", fontsize = font_size),
      col = col_fun,
      show_heatmap_legend = TRUE,
      border = FALSE,
      gap = unit(0.5, "mm"),
      rect_gp = gpar(col = "black", lwd = 0.1),
      heatmap_legend_param = list(
        title = "",
        title_position = "lefttop-rot",
        grid_height = unit(30, "mm"),
        grid_width = unit(6, "mm")
      ),
      column_names_gp = gpar(fontface = "bold", fontsize = 11),
      column_names_centered = FALSE,
      column_names_rot = 90,
      row_title_rot = 90
    )
    
    # 转换为ggplot对象
    plot <- as.ggplot(heatmap_obj)
    
    # 创建细胞类型专属目录
    celltype_dir <- file.path(args$psingle_gene_list, celltype)
    if (!dir.exists(celltype_dir)) {
      dir.create(celltype_dir, recursive = TRUE)
    }
    
    # 保存热图
    output_file <- file.path(celltype_dir, "single_gene_show_express.pdf")
    ggsave(output_file,
           plot = plot,create.dir=T,
           width = args$psize_w,
           height = args$psize_h)
    
    message("✓ 热图已保存: ", output_file)
    return(plot)
  }, error = function(e) {
    stop("热图创建失败: ", e$message)
  })
}

# 主分析流程
run_analysis_pipeline <- function(args) {
  # 验证输入
  validate_inputs(args)
  
  # 加载数据
  seuratdata <- load_seurat_data(args$iSeurat)
  pathway_ids <- load_pathway_ids(args$ipathway_ID)
  
  # 设置细胞标识
  Idents(seuratdata) <- args$pident
  # 获取基因集
  genesets <- get_genesets(pathway_ids)
  
  data <- process_celltype(seuratdata, celltype="Total", genesets, args)
  plot_T <- create_heatmap(data, celltype="Total", args)


  # 处理每个细胞类型
  plot_list <- list()
  celltypes <- levels(seuratdata)
  message("\n[4/5] 分析细胞类型 (总数: ", length(celltypes), ")")
  
  for (celltype in celltypes) {
    tryCatch({
      message("\n--------------------------------------------------")
      message("开始处理: ", celltype)
      
      # 提取细胞子集
      subcell <- subset(seuratdata, manual_L2 == celltype)
      
      # 处理数据
      heatmap_data <- process_celltype(subcell, celltype, genesets, args)
      if (is.null(heatmap_data)) next
      
      # 创建热图
      plot <- create_heatmap(heatmap_data, celltype, args)
      plot_list[[celltype]] <- plot
      
      message("✓ 完成: ", celltype)
    }, error = function(e) {
      warning("处理 ", celltype, " 时出错: ", e$message)
    })
  }
  
  # 保存结果
  message("\n[5/5] 保存结果")
  save(plot_list, file = args$osingle_gene_expression)
  message("✓ 热图对象已保存: ", args$osingle_gene_expression)
  message("\n[SUCCESS] 单基因表达分析完成!")
}

# 命令行参数解析
if (exists("snakemake")) {
  # Snakemake模式
  args <- list(
    iSeurat = snakemake@input$iSeurat,
    ipathway_ID = snakemake@input$ipathway_ID,
    pident = snakemake@params$pident,
    pmin = snakemake@params$pmin,
    pmax = snakemake@params$pmax,
    psize_w = snakemake@params$psize_w,
    psize_h = snakemake@params$psize_h,
    psingle_gene_list = snakemake@params$psingle_gene_list,
    osingle_gene_expression = snakemake@output$osingle_gene_expression
  )
  run_analysis_pipeline(args)
}
