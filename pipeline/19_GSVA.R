# 13_0_GSVA.R
# GSVA Analysis Pipeline
#
#' @description 该脚本实现GSVA（Gene Set Variation Analysis）分析流程，包括：
#'              - 加载基因集和表达矩阵
#'              - 执行GSVA分析
#'              - 保存富集结果
#'              
#' @param inputs 输入：
#'              - igenesets: 基因集文件(.rds)
#'              - idats: 表达矩阵文件(.rds)
#' 
#' @param outputs 输出：
#'              - obulk_gsva: GSVA富集结果(.csv)
#' 
#' @return 生成GSVA富集分析结果
suppressPackageStartupMessages({
  library(GSVA)
  library(BiocParallel)
  library(dplyr)
  library(purrr)
  library(parallel)
})

# 验证输入参数
validate_inputs <- function(args) {
  message("[1/4] 验证输入参数")
  
  required_args <- c("igenesets", "idats", "obulk_gsva")
  missing_args <- setdiff(required_args, names(args))
  
  if (length(missing_args) > 0) {
    stop("缺少必要参数: ", paste(missing_args, collapse = ", "))
  }
  
  # 检查文件是否存在
  file_args <- c("igenesets", "idats")
  for (arg in file_args) {
    if (!file.exists(args[[arg]])) {
      stop("输入文件不存在: ", args[[arg]])
    }
  }
  
  # 创建输出目录
  output_dir <- dirname(args$obulk_gsva)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("✓ 创建输出目录: ", output_dir)
  }
  
  message("✓ 输入验证通过")
}

# 加载数据
load_data <- function(args) {
  message("\n[2/4] 加载数据")
  
  tryCatch({
    # 加载基因集
    genesets <- readRDS(args$igenesets)
    if (length(genesets) == 0) {
      stop("基因集为空，请检查输入文件")
    }
    message("✓ 加载基因集: ",length(genesets),"个基因集")
    
    # 加载表达矩阵
    expr_matrix <- readRDS(args$idats)
    if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
      stop("表达矩阵为空，请检查输入文件")
    }
    message("✓ 加载表达矩阵: ", nrow(expr_matrix), " 个基因 x ", ncol(expr_matrix), " 个样本")
    
    return(list(genesets = genesets, expr_matrix = expr_matrix))
  }, error = function(e) {
    stop("数据加载失败: ", e$message)
  })
}

# 执行GSVA分析
run_gsva <- function(data) {
  message("\n[3/4] 执行GSVA分析")
  
  tryCatch({
    # 确定核心数（留一个核心给系统）
    cores <- max(1, detectCores() - 1)
    message("✓ 使用核心数: ", cores)
    
    # 执行GSVA
    start_time <- Sys.time()
    param <- gsvaParam(
        exprData = as.matrix(data$expr_matrix),
        geneSets = data$genesets,
        minSize = 5,
        maxSize = 500,
        kcdf = "Gaussian"
    )
    parallel.sz = 96
    if(parallel.sz > 1){
        if (.Platform$OS.type == "windows"){
            BiocParallel::register(SnowParam(workers = parallel.sz))
        } else {
            BiocParallel::register(MulticoreParam(workers = parallel.sz))
        }
    }
    gsva_res <- gsva(param, verbose= TRUE)
    end_time <- Sys.time()
    
    duration <- round(difftime(end_time, start_time, units = "mins"), 1)
    message("✓ GSVA分析完成! 耗时: ", duration, " 分钟")
    
    return(gsva_res)
  }, error = function(e) {
    stop("GSVA分析失败: ", e$message)
  })
}

# 保存结果
save_results <- function(gsva_res, output_path) {
  message("\n[4/4] 保存结果")
  
  tryCatch({
    # 转换为数据框
    gsva_df <- data.frame(
      Genesets = rownames(gsva_res),
      gsva_res,
      check.names = FALSE,
      row.names = NULL
    )
    
    # 保存结果
    write.csv(gsva_df, file = output_path, row.names = FALSE)
    message("✓ GSVA结果已保存: ", output_path)
    message("  结果维度: ", nrow(gsva_df), " 行 x ", ncol(gsva_df), " 列")
  }, error = function(e) {
    stop("结果保存失败: ", e$message)
  })
}

# 主分析流程
run_gsva_pipeline <- function(args) {
  # 验证输入
  validate_inputs(args)
  
  # 加载数据
  data <- load_data(args)
  
  # 执行GSVA分析
  gsva_result <- run_gsva(data)
  
  # 保存结果
  save_results(gsva_result, args$obulk_gsva)
  
  message("\n[SUCCESS] GSVA分析完成!")
}

# 命令行参数解析
if (exists("snakemake")) {
  # Snakemake模式
  args <- list(
    igenesets = snakemake@input$igenesets,
    idats = snakemake@input$idats,
    obulk_gsva = snakemake@output$obulk_gsva
  )
  run_gsva_pipeline(args)
} else {
  # 命令行参数模式
  # 定义命令行参数规范
  spec <- matrix(c(
    "igenesets", "i", 1, "character",
    "idats", "d", 1, "character",
    "obulk_gsva", "o", 1, "character"
  ), byrow = TRUE, ncol = 4)
  
  # 解析命令行参数
  parsed_args <- getopt(spec)
  
  # 检查必要参数
  if (any(sapply(parsed_args[c("igenesets", "idats", "obulk_gsva")], is.null))) {
    stop("igenesets, idats, and obulk_gsva files must be specified.")
  }
  
  # 运行流程
  run_gsva_pipeline(parsed_args)
}
