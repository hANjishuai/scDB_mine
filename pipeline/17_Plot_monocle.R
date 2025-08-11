# 15_2_monocle2_plot.R
# Monocle2 Visualization Pipeline
#
#' @description 该脚本实现Monocle2可视化分析流程，包括：
#'              - 加载CDS和Seurat对象
#'              - 生成伪时间热图
#'              - 绘制基因表达随时间变化图
#'              - 绘制通路富集变化图
#'              
#' @param inputs 输入：
#'              - icds: Monocle CDS对象文件(.RDS)
#'              - iseu: Seurat对象文件(.RDS)
#'              - igenes: 关注的基因列表文件(.xls)
#'              - itopn: 差异表达基因文件(.RDS)
#'              - istates_de: 状态相关基因文件(.RDS)
#'              - iBEAM_res_list: BEAM分析结果(.RDS)
#'              - pmetadat_col: 用于分组的元数据列名
#'              - poutdir: 输出目录
#'              - pgemt_list_path: 基因模块列表路径
#' 
#' @param outputs 输出：
#'              - oGenexpress_changes_over_pseudotime: 基因表达随时间变化图(PDF)
#'              - owrap_plots: 通路富集图(PDF)
#'              - oplot_pseudotimeHeatmap: 伪时间热图(PDF)
#'              - oplots_listn_RDS: 通路富集结果(.RDS)
#' 
#' @return 生成Monocle2可视化结果
suppressPackageStartupMessages({
  library(ggplotify)
  library(readxl)
  library(stringr)
  library(dplyr)
  library(ggplot2)
})

# 加载本地函数
if (!file.exists("src/monocle2_plot.R")) {
  stop("缺少本地函数文件: src/monocle2_plot.R")
}
source("src/monocle2_plot.R")

# 验证输入参数
validate_inputs <- function(args) {
  message("[1/6] 验证输入参数")
  
  required_args <- c(
    "icds", "iseu", "igenes", "itopn", "istates_de", "iBEAM_res_list",
    "pmetadat_col", "poutdir", "pgemt_list_path",
    "oGenexpress_changes_over_pseudotime", 
    "oplot_pseudotimeHeatmap", "oplots_listn_RDS"
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
  file_args <- c("icds", "iseu", "igenes", "itopn", "istates_de", "iBEAM_res_list")
  for (arg in file_args) {
    if (!file.exists(args[[arg]])) {
      stop("输入文件不存在: ", args[[arg]])
    }
  }
  
  # 创建输出目录
  if (!dir.exists(args$poutdir)) {
    dir.create(args$poutdir, recursive = TRUE)
    message("✓ 创建输出目录: ", args$poutdir)
  }
  
  message("✓ 输入验证通过")
}

# 加载数据
load_data <- function(args) {
  message("\n[2/6] 加载数据")
  
  tryCatch({    
    # 加载CDS对象
    icds <- readRDS(args$icds)
    message("✓ CDS对象已加载: ", args$icds)
    
    # 加载Seurat对象
    iseu <- readRDS(args$iseu)
    message("✓ Seurat对象已加载: ", args$iseu)
    
    # 加载基因列表
    tryCatch({
      # 尝试读取Sheet2
      igenes <- readxl::read_xls(args$igenes, sheet = "Sheet2") %>% 
        dplyr::pull(Tfs) %>% 
        unique() %>% 
        as.character()
    }, error = function(e) {
      message("尝试读取工作表'Sheet2'失败，尝试默认工作表")
      # 回退到第一个工作表
      igenes <- readxl::read_xls(args$igenes, sheet = 1) %>% 
        dplyr::pull(Tfs) %>% 
        unique() %>% 
        as.character()
    })
    
    if (length(igenes) == 0) {
      warning("基因列表为空，请检查输入文件")
    } else {
      message("✓ 基因列表已加载: ", length(igenes), " 个基因")
    }
    
    # 加载差异表达基因
    itopn <- readRDS(args$itopn)
    message("✓ 差异表达基因已加载: ", args$itopn)
    
    # 加载状态相关基因
    istates_de <- readRDS(args$istates_de)
    message("✓ 状态相关基因已加载: ", args$istates_de)
    
    # 加载BEAM分析结果
    iBEAM_res_list <- readRDS(args$iBEAM_res_list)
    message("✓ BEAM分析结果已加载: ", args$iBEAM_res_list)
    
    return(list(
      icds = icds,
      iseu = iseu,
      igenes = igenes,
      itopn = itopn,
      istates_de = istates_de,
      iBEAM_res_list = iBEAM_res_list    
      ))
  }, error = function(e) {
    stop("数据加载失败: ", e$message)
  })
}

# 生成伪时间热图
generate_pseudotime_heatmap <- function(data, outputs) {
  message("\n[3/6] 生成伪时间热图")
  
  tryCatch({
    plot <- tryCatch({
      plot_pseudotimeHeatmap(
        cds = data$icds,
        topn = data$itopn,
        states_de = data$istates_de,
        special_gene = data$igenes,
        ngene = 20,
        n_cluster = 3
      )
    }, error = function(e) {
      message("伪时间热图生成失败: ", e$message)
      return(NULL)
    }, warning = function(w) {
      message("伪时间热图生成警告: ", w$message)
      return(NULL)
    })
    
    if (!is.null(plot)) {
      ggsave(
        outputs$oplot_pseudotimeHeatmap,
        plot = plot,
        width = 8,
        height = 10,create.dir=TRUE
      )
      message("✓ 伪时间热图已保存: ", outputs$oplot_pseudotimeHeatmap)
    } else {
      warning("伪时间热图生成失败，跳过保存")
    }
    
    return(plot)
  }, error = function(e) {
    stop("伪时间热图生成失败: ", e$message)
  })
}

# 生成分支热图 (保留但注释掉)
generate_branched_heatmaps <- function(data, outputs) {
  message("\n[跳过] 生成分支热图 (当前未启用)")
  # 原脚本中此功能被注释掉，保留结构以备后续扩展
   for(i in seq(1, length(data$iBEAM_res_list), 1)) {
     BEAM_res <- data$iBEAM_res_list[[i]]
     plot <- plot_genes_branchedHeatmap(
       data$icds,
       BEAM_res,
       special_gene = data$igenes,
       branch_p = i
     )
     plot <- as.ggplot(plot)
     ggsave(
       filename = paste0(outputs$poutdir, "/","node_", i, ".pdf"),
       plot = plot,
       width = 6,
       height = 11,create.dir=T
     )
   }
}

# 生成基因表达变化图
generate_expression_changes <- function(data, pmetadat_col, outputs) {
  message("\n[4/6] 生成基因表达随时间变化图")
  
  tryCatch({
    plot <- Genexpress_changes_over_pseudotime(
      data$icds,
      genes = data$igenes,
      meta_col = pmetadat_col
    )
    
    ggsave(
      outputs$oGenexpress_changes_over_pseudotime,
      plot = plot,
      width = 10,
      height = 8,create.dir=TRUE
    )
    message("✓ 基因表达变化图已保存: ", outputs$oGenexpress_changes_over_pseudotime)
    
    return(plot)
  }, error = function(e) {
    stop("基因表达变化图生成失败: ", e$message)
  })
}

# 生成通路富集图
generate_pathway_enrichment <- function(data, pgemt_list_path, pmetadat_col, outputs) {
  message("\n[5/6] 生成通路富集图")
  
  tryCatch({
    # 检查基因模块列表路径
    if (!file.exists(pgemt_list_path)) {
      warning("基因模块列表文件不存在: ", pgemt_list_path)
      return(NULL)
    }
    
    plots_listn <- pathwayenrich_changes_over_pseudotime(
      data$icds,
      data$iseu,
      pgemt_list_path,
      meta_col = pmetadat_col
    )
    for(i in seq_along(plots_listn)){
        ggsave(
            file.path(outputs$poutdir,paste0(names(plots_listn)[i],".pdf")),
            plot = plots_listn[[i]],
            width = 11.5,
            height = 10,create.dir=TRUE
          )
    message("✓ 通路富集图已保存: ", outputs$poutdir)

    }
    
    saveRDS(plots_listn, file = outputs$oplots_listn_RDS)
    message("✓ 通路富集结果已保存: ", outputs$oplots_listn_RDS)
    
    return(plots_listn)
  }, error = function(e) {
    stop("通路富集图生成失败: ", e$message)
  })
}

# 主分析流程
run_visualization_pipeline <- function(args) {
  # 验证输入
  validate_inputs(args)
  
  # 加载数据
  data <- load_data(args)
  
  # 准备输出列表
  outputs <- list(
    oplot_pseudotimeHeatmap = args$oplot_pseudotimeHeatmap,
    oGenexpress_changes_over_pseudotime = args$oGenexpress_changes_over_pseudotime,
#    owrap_plots = args$owrap_plots,
    oplots_listn_RDS = args$oplots_listn_RDS,
    poutdir = args$poutdir
  )
  
  # 生成伪时间热图
  heatmap_plot <- generate_pseudotime_heatmap(data, outputs)
  
  # 生成基因表达变化图
  expression_plot <- generate_expression_changes(data, args$pmetadat_col, outputs)
  
  generate_branched_heatmaps(data, outputs)

  # 生成通路富集图
  enrichment_results <- generate_pathway_enrichment(
    data, 
    args$pgemt_list_path, 
    args$pmetadat_col, 
    outputs
  )
  
  message("\n[SUCCESS] Monocle2可视化分析完成!")
}

# 命令行参数解析
if (exists("snakemake")){
  # Snakemake模式
  args <- list(
    icds = snakemake@input$icds,
    iseu = snakemake@input$iseu,
    igenes = snakemake@input$igenes,
    itopn = snakemake@input$itopn,
    istates_de = snakemake@input$istates_de,
    iBEAM_res_list = snakemake@input$iBEAM_res_list,
    pmetadat_col = snakemake@params$pmetadat_col,
    poutdir = snakemake@params$poutdir,
    pgemt_list_path = snakemake@params$pgemt_list_path,
    oGenexpress_changes_over_pseudotime = snakemake@output$oGenexpress_changes_over_pseudotime,
#    owrap_plots = snakemake@output$owrap_plots,
    oplot_pseudotimeHeatmap = snakemake@output$oplot_pseudotimeHeatmap,
    oplots_listn_RDS = snakemake@output$oplots_listn_RDS
  )
  run_visualization_pipeline(args)
}

#args <- list(
#    icds = "result_out/18_monocle2/Adult_Spleen/ocds.rds",
#    iseu = "result_out/18_monocle2/Adult_Spleen/seurat_rdata.rds",
#    igenes = "db/pyscenic/mouse/Regulon_TFs.xls",
#    itopn = "result_out/19_pre_Plot_monocle/Adult_Spleen/otopn.rds",
#    istates_de = "result_out/19_pre_Plot_monocle/Adult_Spleen/ostates_de.rds",
#    iBEAM_res_list = "result_out/19_pre_Plot_monocle/Adult_Spleen/oBEAM_res_list.rds",
#    pmetadat_col = "manual_L2",
#    poutdir = "result_out/20_Plot_monocle",
#    pgemt_list_path = "db/pyscenic/mouse/Regulon_TFs.xls",
#    oGenexpress_changes_over_pseudotime = "figure_out/20_Plot_monocle/oGenexpress_changes_over_pseudotime.pdf",
#    owrap_plots = "figure_out/20_Plot_monocle/owrap_plots.pdf",
#    oplot_pseudotimeHeatmap = "figure_out/20_Plot_monocle/oplot_pseudotimeHeatmap.pdf",
#    oplots_listn_RDS = "result_out/20_Plot_monocle/oplots_listn_RDS.rds"
#)
