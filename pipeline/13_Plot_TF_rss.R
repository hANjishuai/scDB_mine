#' 10_pyscenic_TF_rss_plot.R
#' SCENIC Analysis: Regulon Visualization Pipeline
#'
#' @description 该脚本实现SCENIC分析结果的可视化流程，包括：
#'              - 细胞特异性TF气泡图
#'              - TF表达量热图
#'              - Regulon活性热图
#'              - Seurat点图和特征图
#'              - RSS排名图
#'              
#' @param inputs 输入：
#'              - iScenic_loom: SCENIC分析生成的loom文件路径
#'              - irss_result: RSS计算结果文件路径
#'              - iSeurat: Seurat对象文件路径
#'              - icellinfo: 细胞元数据文件路径
#'              - iRegulon: 关注的TF列表文件路径
#' @param outputs 输出：
#'              - oresult1: TF RSS气泡图
#'              - oresult2: TF表达量热图
#'              - oresult3: Regulon活性热图
#'              - oresult4: Seurat点图
#'              - oresult5: Seurat特征图
#'              - oresult6: RSS排名图
#'              - oregulon_activity: 分组Scaled的Regulon活性
#' 
#' @return 生成7种可视化结果并保存到指定路径

suppressPackageStartupMessages({
  library(dplyr)
  library(readxl)
  library(ggrepel)
  library(Seurat)
  library(SCopeLoomR)
  library(SCENIC)
  library(AUCell)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  library(reshape2)
})

#------------------- 核心函数 -------------------#
validate_inputs <- function(inputs) {
  message("[1/8] Validating inputs")
  
  required_inputs <- c(
    "iScenic_loom", "irss_result", "iSeurat", 
    "icellinfo", "iRegulon"
  )
  
  missing <- sapply(required_inputs, function(x) {
    if (is.null(inputs[[x]]) || !file.exists(inputs[[x]])) {
      return(x)
    }
    return(NULL)
  })
  
  missing <- unlist(Filter(Negate(is.null), missing))
  
  if (length(missing) > 0) {
    stop("Missing or invalid input files: ", paste(missing, collapse = ", "))
  }
  
  message("✓ Input validation passed")
}

load_rss_data <- function(irss_result) {
  message("\n[2/8] Loading RSS results from: ", irss_result)
  
  tryCatch({
    load(irss_result)  # 加载后变量名为rss
    if (!exists("rss")) {
      stop("RSS data 'rss' not found in: ", irss_result)
    }
    
    message("✓ RSS data loaded")
    message("  Regulons: ", nrow(rss), ", Cell types: ", ncol(rss))
    return(rss)
  }, error = function(e) {
    stop("Failed to load RSS data: ", e$message)
  })
}

load_scenic_data <- function(loom_path) {
  message("\n[3/8] Loading SCENIC loom data from: ", loom_path)
  
  tryCatch({
    scenic_loom <- open_loom(loom_path)
    regulonAUC <- get_regulons_AUC(scenic_loom, column.attr.name = 'RegulonsAUC')
    
    message("✓ SCENIC data loaded")
    message("  Cells in AUC matrix: ", ncol(regulonAUC))
    
    return(regulonAUC)
  }, error = function(e) {
    stop("Failed to load SCENIC loom data: ", e$message)
  })
}

load_seurat_data <- function(seurat_path) {
  message("\n[4/8] Loading Seurat object from: ", seurat_path)
  
  tryCatch({
    env <- new.env()
    load(seurat_path, envir = env)
    
    obj_name <- ls(env)[1]
    seurat_obj <- get(obj_name, envir = env)
        
    message("✓ Seurat object loaded")
    message("  Cells: ", ncol(seurat_obj), ", Features: ", nrow(seurat_obj))
    
    return(seurat_obj)
  }, error = function(e) {
    stop("Failed to load Seurat object: ", e$message)
  })
}

load_cell_metadata <- function(cellinfo_path) {
  message("\n[5/8] Loading cell metadata from: ", cellinfo_path)
  
  tryCatch({
    cellinfo <- read.csv(cellinfo_path)
    required_cols <- c("CellType")
    
    if (!all(required_cols %in% colnames(cellinfo))) {
      stop("Missing required columns in cellinfo: ", 
           paste(setdiff(required_cols, colnames(cellinfo)), collapse = ", "))
    }
    
    message("✓ Cell metadata loaded")
    message("  Cells: ", nrow(cellinfo), ", Cell types: ", 
            paste(unique(cellinfo$CellType), collapse = ", "))
    
    return(cellinfo)
  }, error = function(e) {
    stop("Failed to load cell metadata: ", e$message)
  })
}

load_regulon_lists <- function(regulon_path) {
  message("\n[6/8] Loading regulon lists from: ", regulon_path)
  
  tryCatch({
    regulon_sheets <- excel_sheets(regulon_path)
    
    if (length(regulon_sheets) < 2) {
      stop("Regulon file must contain at least 2 sheets")
    }
    
    tf_regulon <- read_xls(regulon_path, sheet = 1)
    gene_regulon <- read_xls(regulon_path, sheet = 2)
    
    message("✓ Regulon lists loaded")
    message("  TF regulons: ", nrow(tf_regulon), ", Gene regulons: ", nrow(gene_regulon))
    
    return(list(
      tf = tf_regulon$regulon,
      gene = gene_regulon$Tfs
    ))
  }, error = function(e) {
    stop("Failed to load regulon lists: ", e$message)
  })
}

calculate_regulon_activity <- function(regulonAUC, cellinfo) {
  message("\n[7/8] Calculating regulon activity by group")
  
  tryCatch({
    cellTypes <- cellinfo %>% select('CellType')
    selectedResolution <- "CellType"
    cellsPerGroup <- split(rownames(cellTypes), cellTypes[, selectedResolution])
    
    regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells) {
      rowMeans(getAUC(regulonAUC)[, as.numeric(cells)])
    })
    
    regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup), center = TRUE, scale = TRUE))
    regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled) %>% 
      as.data.frame()
    
    message("✓ Regulon activity calculated")
    message("  Cell groups: ", ncol(regulonActivity_byGroup_Scaled),
            ", Regulons: ", nrow(regulonActivity_byGroup_Scaled))
    
    return(regulonActivity_byGroup_Scaled)
  }, error = function(e) {
    stop("Regulon activity calculation failed: ", e$message)
  })
}

#------------------- 可视化函数 -------------------#
generate_plots <- function(rss, regulonAUC, seurat_obj, 
                           regulon_lists, regulon_activity, outputs_file_dir,output_rank_RSS) {
  message("\n[8/8] Generating visualizations")
  
  tryCatch({
    # 加载自定义绘图函数
    source("src/Pyscenic_TF_rss_plot.R")
    if(!dir.exists(outputs_file_dir)){dir.create(path=outputs_file_dir,recursive = T)}
    
    # 01: TF RSS气泡图
    message("  Generating TF RSS bubble plot...")
    oresult1 <- file.path(outputs_file_dir,"Rss_bubble.pdf")
    pdf(file = oresult1, width = 3.5, height = 5.5)
    p1 <- TF_rss_droplot(rss)
    print(p1)
    dev.off()
    
    # 02: TF表达量热图
    message("  Generating TF expression heatmap...")
    rss_data <- p1$plot$data
    rss_data <- reshape2::dcast(rss_data, Topic ~ cellType, value.var = 'Z')
    rownames(rss_data) <- rss_data[, 1]
    rss_data <- rss_data[, -1]
    
    group_colors <- c("#e04f61", "#f7f398", "#e07c95", "#322dfd",
                     "#23c7d6", "#7491b7", "#8fcc4a", "#91CDC8", "#f2311233",
                     "#cfa560", "#c2606d", "#F0988C", "#ed1223", "#9E9E9E")
    names(group_colors) <- colnames(rss_data)
    
    oresult2 <- file.path(outputs_file_dir,"TF_exprewss.pdf")
    pdf(oresult2, width = 4.6, height = 6.5)
    p2 <- TF_exprewss_heatmap(
      Data = rss_data,
#      annotation_cols = data.frame(Celltypes = colnames(rss_data)),
      annotation_color = list(Celltypes = group_colors)
    )
    print(p2)
    dev.off()
    
    # 03: Regulon活性热图 (仅当提供TF列表时)
    if (!is.null(regulon_lists$tf)) {
      message("  Generating regulon activity heatmap...")
      tf_c <- regulon_lists$tf
      anno_heatmap <- regulon_activity %>%
        filter(rownames(.) %in% tf_c) %>%
        mutate(index = which(rownames(regulon_activity) %in% tf_c)) %>%
        tibble::rownames_to_column("mark") %>%
        mutate(mark = gsub('\\(.*\\)', '', mark)) %>%
        dplyr::select(index, mark)
      
      lab <- rowAnnotation(
        ano = anno_mark(
          at = anno_heatmap$index,
          labels = anno_heatmap$mark,
          labels_gp = gpar(fontsize = 8)
        )
      )
      oresult3 <- file.path(outputs_file_dir,"TF_regulonActivity.pdf")      
      pdf(oresult3, width = 3, height = 7)
      p3 <- TF_regulonActivity_heatmap(
        data = regulon_activity,
        right_annotatio = lab
      )
      print(p3)
      dev.off()
    } else {
      message("  Skipping regulon activity heatmap (no TF list provided)")
    }
    
    # 04-05: Seurat点图和特征图 (仅当提供TF和基因列表时)
    if (!is.null(regulon_lists$tf) && !is.null(regulon_lists$gene)) {
      message("  Generating Seurat dotplot and featureplot...")
      tf_c <- regulon_lists$tf
      gene_c <- regulon_lists$gene
      
      # 将AUC矩阵添加到Seurat元数据
      next_regulonAUC <- regulonAUC[, match(colnames(seurat_obj), colnames(regulonAUC))]
      seurat_obj@meta.data <- cbind(
        seurat_obj@meta.data,
        t(SummarizedExperiment::assay(next_regulonAUC[regulonAUC@NAMES, ]))
      )
      
      p_list <- TF_seurat_dotplot(
        seuratdata = seurat_obj,
#        group = 'orig.ident',
        TF_lists = tf_c,
        gene_c = gene_c
      )
      oresult4 <- file.path(outputs_file_dir,"TF_bubble.pdf")
      oresult5 <- file.path(outputs_file_dir,"TF_dimplot.pdf")
      ggsave(oresult4, plot = p_list[[1]], width = 7, height = 3)
      ggsave(oresult5, plot = p_list[[2]], width = 17, height = 12)
    } else {
      message("  Skipping Seurat plots (missing TF or gene list)")
    }
    
    # 06: RSS排名图
    message("  Generating RSS rank plot...")
    rss_rankplot <- TF_rssRank_plot(rss)
    oresult6 <- file.path(outputs_file_dir,"TF_rank_RSS.pdf")
    plot <- rss_rankplot$rankplot
    ggsave(oresult6,plot=plot,width=10,height=10)
    saveRDS(rss_rankplot, output_rank_RSS)
    
    message("✓ All visualizations generated successfully")
    
  }, error = function(e) {
    stop("Visualization generation failed: ", e$message)
  })
}

#------------------- 主流程函数 -------------------#
run_visualization_pipeline <- function(inputs,
                                       outputs,
                                       outputs_file_dir
                                      ) {
  # 验证输入
  validate_inputs(inputs)
  if(!dir.exists(outputs_file_dir)){dir.create(outputs_file_dir, recursive= T)}
  if(!dir.exists(dirname(outputs$oregulon_activity))){dir.create(dirname(outputs$oregulon_activity), recursive= T)}
  # 加载数据
  rss_data <- load_rss_data(inputs$irss_result)
  scenic_data <- load_scenic_data(inputs$iScenic_loom)
  seurat_data <- load_seurat_data(inputs$iSeurat)
  cell_metadata <- load_cell_metadata(inputs$icellinfo)
  regulon_lists <- load_regulon_lists(inputs$iRegulon)
  
  # 计算regulon活性
  regulon_activity <- calculate_regulon_activity(scenic_data, cell_metadata)
  write.csv(regulon_activity, outputs$oregulon_activity)
  
  # 生成可视化
  generate_plots(
    rss = rss_data,
    regulonAUC = scenic_data,
    seurat_obj = seurat_data,
    regulon_lists = regulon_lists,
    regulon_activity = regulon_activity,
    outputs_file_dir = outputs_file_dir,
    output_rank_RSS = outputs$output_rank_RSS
  )
  
  message("\n[SUCCESS] Visualization pipeline completed!")
}

#------------------- 命令行/Snakemake集成 -------------------#
if (exists("snakemake")){
  # Snakemake模式
  inputs = list(
    iScenic_loom = snakemake@input$iScenic_loom,
    irss_result = snakemake@input$irss_result,
    iSeurat = snakemake@input$iSeurat,
    icellinfo = snakemake@input$icellinfo,
    iRegulon = snakemake@input$iRegulon
  )

  outputs_file_dir <- snakemake@params$outputs_file_dir

  outputs = list(
    oregulon_activity = snakemake@output$oregulon_activity,
    output_rank_RSS = snakemake@output$output_rank_RSS
  )
  

  run_visualization_pipeline(inputs,
                             outputs,
                             outputs_file_dir)
}


