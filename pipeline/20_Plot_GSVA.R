# 13_1_GSVA_heatmap.R
# GSVA Heatmap Visualization Pipeline
#
#' @description 该脚本实现GSVA结果的热图可视化流程，包括：
#'              - 加载GSVA结果和基因集数据
#'              - 筛选通路中差异表达的基因
#'              - 绘制热图展示通路富集情况
#'              
#' @param inputs 输入：
#'              - ibulk_gsva: GSVA富集结果(.csv)
#'              - igenesets: 基因集文件(.rds)
#'              - iDEG_path: 差异基因分析结果(.csv)
#' 
#' @param outputs 输出：
#'              - obulk_gsva_heatmap_data: 热图数据(.csv)
#'              - obulk_gsva_heatmap: 热图结果(PDF)
#'              - ideg_list: 通路中差异表达基因的目录
#' 
#' @return 生成GSVA热图可视化结果
suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(argparse)
})

# 验证输入参数
validate_inputs <- function(args) {
  message("[1/6] 验证输入参数")
  
  required_args <- c(
    "ibulk_gsva", "igenesets", "iDEG_path",
    "pmin", "pmax", "psize_w", "psize_h",
    "obulk_gsva_heatmap_data", "obulk_gsva_heatmap", "ideg_list"
  )
  
  missing_args <- setdiff(required_args, names(args))
  
  if (length(missing_args) > 0) {
    stop("缺少必要参数: ", paste(missing_args, collapse = ", "))
  }
  
  # 检查文件是否存在
  file_args <- c("ibulk_gsva", "igenesets", "iDEG_path")
  for (arg in file_args) {
    if (!file.exists(args[[arg]])) {
      stop("输入文件不存在: ", args[[arg]])
    }
  }
  
  # 检查数值参数
  num_args <- c("pmin", "pmax", "psize_w", "psize_h")
  for (arg in num_args) {
    if (!is.numeric(args[[arg]]) || is.na(args[[arg]])) {
      stop("参数必须为数值: ", arg)
    }
  }
  
  # 创建输出目录
  output_dirs <- unique(dirname(unlist(args[c("obulk_gsva_heatmap_data", "obulk_gsva_heatmap", "ideg_list")])))
  for (odir in output_dirs) {
    if (!dir.exists(odir)) {
      dir.create(odir, recursive = TRUE)
      message("✓ 创建输出目录: ", odir)
    }
  }
  
  message("✓ 输入验证通过")
}

# 加载GSVA结果
load_gsva_data <- function(ibulk_gsva) {
  message("\n[2/6] 加载GSVA结果")
  
  tryCatch({
    bulk_gsva <- read_csv(ibulk_gsva, show_col_types = FALSE)
    
    if (nrow(bulk_gsva) == 0 || ncol(bulk_gsva) < 2) {
      stop("GSVA结果为空或格式不正确")
    }
    
    # 处理重复的通路名称
    if (any(duplicated(bulk_gsva[[1]]))) {
      message("检测到重复的通路名称，保留第一个出现项")
      bulk_gsva <- bulk_gsva %>% 
        distinct(!!sym(names(bulk_gsva)[1]), .keep_all = TRUE)
    }
    
    # 转换为矩阵格式
    gsva_matrix <- bulk_gsva %>%
      column_to_rownames(var = names(bulk_gsva)[1]) %>%
      as.matrix()
    
    message("✓ GSVA矩阵维度: ", nrow(gsva_matrix), " 个通路 x ", ncol(gsva_matrix), " 个样本")
    
    return(list(raw = bulk_gsva, matrix = gsva_matrix))
  }, error = function(e) {
    stop("GSVA结果加载失败: ", e$message)
  })
}

# 保存热图数据
save_heatmap_data <- function(gsva_data, output_path) {
  message("\n[3/6] 保存热图数据")
  
  tryCatch({
    # 转换回数据框格式
    heatmap_df <- gsva_data$matrix %>%
      as.data.frame() %>%
      rownames_to_column(var = "Pathway")
    
    write_csv(heatmap_df, output_path)
    message("✓ 热图数据已保存: ", output_path)
    
    return(heatmap_df)
  }, error = function(e) {
    stop("热图数据保存失败: ", e$message)
  })
}

# 筛选通路中的差异表达基因
filter_pathway_deg <- function(genesets_path, deg_path, output_dir) {
  message("\n[4/6] 筛选通路中的差异表达基因")
  
  tryCatch({
    # 加载基因集
    genesets <- readRDS(genesets_path)
    
    if (length(genesets) == 0) {
      stop("基因集为空")
    }
    
    # 创建基因-通路映射
    gene_pathway <- lapply(names(genesets), function(path) {
      data.frame(
        pathway = path,
        gene = genesets[[path]]
      )
    }) %>% bind_rows()
    
    # 加载差异表达基因
    degs <- read_csv(deg_path, show_col_types = FALSE)
    
    if (nrow(degs) == 0) {
      warning("差异表达基因文件为空")
      return(NULL)
    }
    
    # 提取文件名
    base_name <- tools::file_path_sans_ext(basename(deg_path))
    output_file <- file.path(output_dir, paste0(base_name, "_DEG_inpathway.csv"))
    
    # 合并数据
    degs_merge <- degs %>%
      inner_join(gene_pathway, by = c("gene" = "gene"))
    
    # 保存结果
    write_csv(degs_merge, output_file)
    message("✓ 通路差异基因已保存: ", output_file)
    message("  筛选到 ", nrow(degs_merge), " 个通路相关差异基因")
    
    return(degs_merge)
  }, error = function(e) {
    stop("通路差异基因筛选失败: ", e$message)
  })
}

# 创建热图
create_heatmap <- function(gsva_matrix, pmin, pmax, width, height, output_path) {
  message("\n[5/6] 创建热图")
  
  tryCatch({
    # 设置颜色映射
    col_fun <- colorRamp2(
      c(pmin, 0, pmax), 
      c("navy", "#e8f6e5", "#d43627")
    )
    
    # 自动调整行名字体大小
    n_rows <- nrow(gsva_matrix)
    font_size <- ifelse(n_rows > 50, 6, 
                       ifelse(n_rows > 30, 8, 10))
    
    # 创建热图对象
    heatmap_obj <- Heatmap(
      gsva_matrix,
      name = "GSVA Score",
      row_names_max_width = unit(12, "cm"),
      cluster_columns = FALSE,
      cluster_rows = FALSE,
      show_row_names = TRUE,
      row_names_gp = gpar(fontface = "bold", fontsize = font_size),
      col = col_fun,
      show_heatmap_legend = TRUE,
      border = TRUE,
      border_gp = gpar(col = "black", lwd = 0.5),
      gap = unit(1, "mm"),
      rect_gp = gpar(col = "gray90", lwd = 0.1),
      heatmap_legend_param = list(
        title = "GSVA Score",
        title_position = "lefttop",
        grid_height = unit(20, "mm"),
        grid_width = unit(5, "mm"),
        legend_direction = "horizontal"
      ),
      column_names_gp = gpar(fontface = "bold", fontsize = 10),
      column_names_rot = 45,
      column_names_centered = FALSE
    )
    
    # 保存PDF
    pdf(output_path, width = width, height = height)
    draw(heatmap_obj, heatmap_legend_side = "bottom")
    dev.off()
    
    message("✓ 热图已保存: ", output_path)
    message("  热图尺寸: ", width, " x ", height, " 英寸")
    
    return(heatmap_obj)
  }, error = function(e) {
    stop("热图创建失败: ", e$message)
  })
}

# 主分析流程
run_heatmap_pipeline <- function(args) {
  # 验证输入
  validate_inputs(args)
  
  # 加载GSVA结果
  gsva_data <- load_gsva_data(args$ibulk_gsva)
  
  # 保存热图数据
  heatmap_df <- save_heatmap_data(
    gsva_data, 
    args$obulk_gsva_heatmap_data
  )
  
  # 筛选通路中的差异表达基因
  pathway_deg <- filter_pathway_deg(
    args$igenesets,
    args$iDEG_path,
    args$ideg_list
  )
  
  # 创建热图
  heatmap_obj <- create_heatmap(
    gsva_data$matrix,
    args$pmin,
    args$pmax,
    args$psize_w,
    args$psize_h,
    args$obulk_gsva_heatmap
  )
  
  message("\n[SUCCESS] GSVA热图分析完成!")
}

# 命令行参数解析
if (exists("snakemake")) {
  # Snakemake模式
  args <- list(
    ibulk_gsva = snakemake@input$ibulk_gsva,
    igenesets = snakemake@input$igenesets,
    iDEG_path = snakemake@input$iDEG_path,
    pmin = snakemake@params$pmin,
    pmax = snakemake@params$pmax,
    psize_w = snakemake@params$psize_w,
    psize_h = snakemake@params$psize_h,
    obulk_gsva_heatmap_data = snakemake@output$obulk_gsva_heatmap_data,
    obulk_gsva_heatmap = snakemake@output$obulk_gsva_heatmap,
    ideg_list = snakemake@params$ideg_list
  )
  run_heatmap_pipeline(args)
}
