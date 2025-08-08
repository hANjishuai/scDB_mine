#' 09_0_pyscenic_plot.R
#' SCENIC Analysis: Regulon Specificity Score Calculation
#'
#' @description 该脚本计算SCENIC分析中的Regulon特异性评分（RSS）
#' @param inputs 输入：
#'              - iScenic_loom: SCENIC分析生成的loom文件路径
#'              - iSeurat: Seurat对象文件路径
#' @param params 参数：
#'              - pcelltype: 细胞类型元数据列名
#'              - pgroup: 分组元数据列名
#'              - pnF_RNA: 基因数量元数据列名
#'              - pnC_RNA: UMI数量元数据列名
#' @param outputs 输出：
#'              - orss: 完整RSS结果文件
#'              - orss_activate: 激活态Regulon的RSS结果
#'              - orss_inactivate: 抑制态Regulon的RSS结果
#' 
#' @return 保存RSS分析结果到指定文件

suppressPackageStartupMessages({
  library(SCopeLoomR)
  library(AUCell)
  library(SCENIC)
  library(dplyr)
  library(stringr)
})

#------------------- 核心函数 -------------------#
validate_inputs <- function(loom_path, seurat_path) {
  message("[1/6] Validating inputs")
  
  # 检查文件存在性
  if (!file.exists(loom_path)) {
    stop("SCENIC loom file not found: ", loom_path)
  }
  
  if (!file.exists(seurat_path)) {
    stop("Seurat object file not found: ", seurat_path)
  }
  
  message("✓ Input validation passed")
}

load_scenic_data <- function(loom_path) {
  message("\n[2/6] Loading SCENIC loom data from: ", loom_path)
  
  tryCatch({
    # 打开loom文件并提取必要数据
    scenic_loom <- open_loom(loom_path)
    
    regulons_incidMat <- get_regulons(scenic_loom, column.attr.name = "Regulons")
    regulons <- regulonsToGeneLists(regulons_incidMat)
    regulonAUC <- get_regulons_AUC(scenic_loom, column.attr.name = 'RegulonsAUC')
    
    message("✓ SCENIC data loaded successfully")
    message("  Regulons count: ", length(regulons))
    message("  Cells in AUC matrix: ", ncol(regulonAUC))
    
    return(list(
      regulons = regulons,
      regulonAUC = regulonAUC
    ))
  }, error = function(e) {
    stop("Failed to load SCENIC loom data: ", e$message)
  })
}

load_and_subset_seurat <- function(seurat_path) {
  message("\n[3/6] Loading Seurat object from: ", seurat_path)
  
  tryCatch({
    env <- new.env()
    load(seurat_path, envir = env)  # 加载Seurat对象

    obj_name <- ls(env)[1]
    seurat_obj <- get(obj_name, envir = env)

    # 检查对象是否存在
    if (!exists("seurat_obj")) {
      stop("Seurat object 'seurat_obj' not found in the loaded file")
    }
    
    message("✓ Seurat object processed")
    message("  Cells retained: ", ncol(seurat_obj))
    
    return(seurat_obj)
  }, error = function(e) {
    stop("Failed to load or subset Seurat object: ", e$message)
  })
}

prepare_cell_metadata <- function(seurat_obj, celltype_col, group_col, nF_col, nC_col) {
  message("\n[4/6] Preparing cell metadata")
  
  tryCatch({
    # 提取元数据并重命名列
    cellinfo <- seurat_obj@meta.data[, c(celltype_col, group_col, nF_col, nC_col)]
    colnames(cellinfo) <- c('celltype', 'group', 'nGene', 'nUMI')
    # 这里的celltype很重要哦！！！！
    
    # 检查是否有缺失值
    if (any(is.na(cellinfo$celltype))) {
      warning("NA values found in celltype column, they will be removed in RSS calculation")
    }
    
    message("✓ Cell metadata prepared")
    message("  Celltypes: ", paste(unique(cellinfo$celltype), collapse = ", "))
    
    return(cellinfo)
  }, error = function(e) {
    stop("Failed to prepare cell metadata: ", e$message)
  })
}

.H <- function(pVect){
  pVect <- pVect[pVect>0] # /sum(pVect) ??
  - sum(pVect * log2(pVect))
}

# Jensen-Shannon Divergence (JSD)
calcJSD <- function(pRegulon, pCellType){
  (.H((pRegulon+pCellType)/2)) - ((.H(pRegulon)+.H(pCellType))/2)
}
 
# Regulon specificity score (RSS)
.calcRSS.oneRegulon <- function(pRegulon, pCellType){
  jsd <- calcJSD(pRegulon, pCellType)
  1 - sqrt(jsd)
}

calcRSS <- function(AUC, cellAnnotation, cellTypes=NULL){
  if(any(is.na(cellAnnotation))) stop("NAs in annotation")
  if(any(class(AUC)=="aucellResults")) AUC <- getAUC(AUC)
  normAUC <- AUC/rowSums(AUC)
  if(is.null(cellTypes)) cellTypes <- unique(cellAnnotation)
  # 
  ctapply <- lapply
  if(require('BiocParallel')) ctapply <- bplapply
  
  rss <- ctapply(cellTypes, function(thisType)
    sapply(rownames(normAUC), function(thisRegulon)
    {
      pRegulon <- normAUC[thisRegulon,]
      pCellType <- as.numeric(cellAnnotation==thisType)
      pCellType <- pCellType/sum(pCellType)
      .calcRSS.oneRegulon(pRegulon, pCellType)
    })
  )
  rss <- do.call(cbind, rss)
  colnames(rss) <- cellTypes
  return(rss)
}

calculate_rss <- function(regulonAUC, cellinfo) {
  message("\n[5/6] Calculating Regulon Specificity Scores (RSS)")
  celltype_col = "celltype"
  cellTypes <- as.data.frame(subset(cellinfo,select = celltype_col))
  selectedResolution <- celltype_col

  tryCatch({
    # 计算RSS
    rss <- calcRSS(
      AUC = getAUC(regulonAUC),
      cellAnnotation =cellTypes[colnames(regulonAUC), celltype_col]
    )
    
    # 移除含有NA的行
    rss <- na.omit(rss)
    
    message("✓ RSS calculation completed")
    message("  Regulons in RSS matrix: ", nrow(rss))
    
    return(rss)
  }, error = function(e) {
    stop("RSS calculation failed: ", e$message)
  })
}

save_rss_results <- function(rss, results_path) {
  message("\n[6/6] Saving RSS results")
  if(!dir.exists(results_path)){dir.create(results_path)}
  tryCatch({
    # 分离激活态和抑制态regulons
    rss_activate <- rss[grepl("\\+", rownames(rss)), ]
    rss_inactivate <- rss[grepl("\\-", rownames(rss)), ]
    
    # 保存结果
    full_path <- file.path(results_path,"rss_all.Rdata")
    activate_path <- file.path(results_path,"rss_activate.Rdata")
    inactivate_path <- file.path(results_path,"rss_inactivate.Rdata")

    save(rss_activate, file = activate_path)
    save(rss_inactivate, file = inactivate_path)
    save(rss, file = full_path)
    
    message("✓ Results saved successfully")
    message("  Activated regulons: ", nrow(rss_activate))
    message("  Inactivated regulons: ", nrow(rss_inactivate))
    message("  Full RSS saved to: ", full_path)
  }, error = function(e) {
    stop("Failed to save RSS results: ", e$message)
  })
}

#------------------- 主流程函数 -------------------#
run_rss_pipeline <- function(
    loom_path,
    seurat_path,
    celltype_col,
    group_col,
    nF_col,
    nC_col,
    results_path,
    seurat_res_path
) {

  # 1. 验证输入
  validate_inputs(loom_path, seurat_path)
  
  # 3. 加载SCENIC数据
  scenic_data <- load_scenic_data(loom_path)
  
  # 4. 加载并处理Seurat对象
  seurat_obj <- load_and_subset_seurat(seurat_path)
  
  # 5. 准备细胞元数据
  cell_metadata <- prepare_cell_metadata(
    seurat_obj, 
    celltype_col, 
    group_col, 
    nF_col, 
    nC_col
  )
  
  # 6. 计算RSS
  rss_results <- calculate_rss(
    scenic_data$regulonAUC,
    cell_metadata
  )
  
  # 7. 保存结果
  save_rss_results(
    rss_results,
    results_path
  )
  save(seurat_obj,file = seurat_res_path)
  message("\n[SUCCESS] RSS pipeline completed!")
}

#------------------- Snakemake集成入口 -------------------#
if (exists("snakemake")) {

  loom_path = as.character(snakemake@input$loom_path)
  seurat_path = as.character(snakemake@input$seurat_path)
  celltype_col = as.character(snakemake@params$celltype_col)
  group_col = as.character(snakemake@params$group_col)
  nF_col = as.character(snakemake@params$nF_col)
  nC_col = as.character(snakemake@params$nC_col)
  results_path = as.character(snakemake@params$results_path)
  seurat_res_path = as.character(snakemake@output$seurat_res_path)

  run_rss_pipeline(
    loom_path = loom_path,
    seurat_path = seurat_path,
    celltype_col = celltype_col,
    group_col = group_col,
    nF_col = nF_col,
    nC_col = nC_col,
    results_path = results_path,
    seurat_res_path = seurat_res_path
  )
}
