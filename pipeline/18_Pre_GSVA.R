# 10_SubcellAnnoCheck_GSVAprep.R
# GSVA Preparation Pipeline
#
#' @description 该脚本实现GSVA分析前的数据准备，包括：
#'              - 加载Seurat对象
#'              - 计算细胞类型的平均表达矩阵
#'              - 获取选定通路的基因集
#'              
#' @param inputs 输入：
#'              - iSeurat: Seurat对象文件(.RData)
#'              - pident: 用于设置细胞身份的列名
#'              - ipathway_ID: 选定通路ID文件(.xls)
#' 
#' @param outputs 输出：
#'              - odat: 平均表达矩阵(.rds)
#'              - ogenesets: 基因集对象(.rds)
#' 
#' @return 生成GSVA分析所需数据
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(msigdbr)
  library(readxl)
})

# 加载本地函数
if (!file.exists("src/gsva_prep_functions.R")) {
  stop("缺少本地函数文件: src/gsva_prep_functions.R")
}
source("src/gsva_prep_functions.R")

# 验证输入参数
validate_inputs <- function(args) {
  message("[1/4] 验证输入参数")
  
  required_args <- c(
    "iSeurat", "pident", "ipathway_ID", 
    "odat", "ogenesets"
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
  file_args <- c("iSeurat", "ipathway_ID")
  for (arg in file_args) {
    if (!file.exists(args[[arg]])) {
      stop("输入文件不存在: ", args[[arg]])
    }
  }
  
  # 创建输出目录
  output_dirs <- unique(dirname(unlist(args[c("odat", "ogenesets")])))
  for (odir in output_dirs) {
    if (!dir.exists(odir)) {
      dir.create(odir, recursive = TRUE)
      message("✓ 创建输出目录: ", odir)
    }
  }
  
  message("✓ 输入验证通过")
}

# 加载Seurat对象
load_seurat <- function(iSeurat) {
  message("\n[2/4] 加载Seurat对象")
  
  tryCatch({
    env <- new.env()
    load(iSeurat, envir = env)
    obj_name <- ls(env)[1]
    seurat_obj <- get(obj_name, envir = env)
    
    if (!inherits(seurat_obj, "Seurat")) {
      stop("加载的对象不是Seurat类型")
    }
    
    message("✓ Seurat对象已加载: ", iSeurat)
    message("  对象包含细胞数: ", ncol(seurat_obj))
    message("  对象包含基因数: ", nrow(seurat_obj))
    
    return(seurat_obj)
  }, error = function(e) {
    stop("Seurat对象加载失败: ", e$message)
  })
}

# 计算平均表达矩阵
calculate_average_expression <- function(seurat_obj, pident) {
  message("\n[3/4] 计算平均表达矩阵")
  
  tryCatch({
    # 设置细胞身份
    if (!pident %in% colnames(seurat_obj@meta.data)) {
      stop("指定的身份列 '", pident, "' 不在元数据中")
    }
    
    Idents(seurat_obj) <- pident
    message("✓ 设置细胞身份为: ", pident)
    message("  细胞类型数量: ", length(unique(seurat_obj@active.ident)))
    
    # 计算平均表达
    dat <- AverageExpression(
      seurat_obj, 
      assays = "SCT", 
      slot = "data"
    )[[1]]
    
    # 过滤零表达基因
    dat <- dat[rowSums(dat) > 0, ]
    dat <- as.matrix(dat)
    
    message("✓ 平均表达矩阵维度: ", nrow(dat), " x ", ncol(dat))
    
    if (ncol(dat) < 2) {
      warning("细胞类型数量小于2，可能无法进行后续分析")
    }
    
    return(dat)
  }, error = function(e) {
    stop("平均表达计算失败: ", e$message)
  })
}

# 获取基因集
prepare_genesets <- function(ipathway_ID) {
  message("\n[4/4] 准备基因集")
  
  tryCatch({
    # 读取选定的通路ID
    ID_selected <- readxl::read_xls(ipathway_ID, sheet = "Sheet3") %>% 
      dplyr::select(ID) %>% 
      unique() %>%
      pull()
    
    if (length(ID_selected) == 0) {
      stop("未找到有效通路ID，请检查输入文件")
    }
    
    message("✓ 读取通路ID数量: ", length(ID_selected))
    
    # 获取基因集
    genesets <- get_genesets(ID_selected)
    
    if (length(genesets) == 0) {
      stop("未获取到基因集，请检查通路ID")
    }
    
    message("✓ 获取基因集数量: ", length(genesets))
    message("  示例基因集: ", names(genesets)[1])
    
    return(genesets)
  }, error = function(e) {
    stop("基因集准备失败: ", e$message)
  })
}

# 主分析流程
run_gsva_prep_pipeline <- function(args) {
  # 验证输入
  validate_inputs(args)
  
  # 加载Seurat对象
  seurat_obj <- load_seurat(args$iSeurat)
  
  # 计算平均表达
  avg_expr <- calculate_average_expression(seurat_obj, args$pident)
  
  # 准备基因集
  genesets <- prepare_genesets(args$ipathway_ID)
  
  # 保存结果
  message("\n保存结果...")
  saveRDS(avg_expr, file = args$odat)
  message("✓ 平均表达矩阵已保存: ", args$odat)
  
  saveRDS(genesets, file = args$ogenesets)
  message("✓ 基因集已保存: ", args$ogenesets)
  
  message("\n[SUCCESS] GSVA数据准备完成!")
}

# 命令行参数解析
if (exists("snakemake")) {
  # Snakemake模式
  args <- list(
    iSeurat = snakemake@input$iSeurat,
    pident = snakemake@params$pident,
    ipathway_ID = snakemake@input$ipathway_ID,
    odat = snakemake@output$odat,
    ogenesets = snakemake@output$ogenesets
  )
  run_gsva_prep_pipeline(args)
} else {
  # 直接运行模式（用于测试）
  args <- list(
    iSeurat = "result_out/10_SubcellAnnoCheck/Adult_Spleen/seurat_rdata.Rdata",
    pident = "manual_L2",
    ipathway_ID = "db/pyscenic/mouse/Regulon_TFs.xls",
    odat = "result_out/11_GSVAprep/gsva_dats.rds",
    ogenesets = "result_out/11_GSVAprep/gsva_genesets.rds"
  )
  run_gsva_prep_pipeline(args)
}

