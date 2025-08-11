# 15_2_monocle2_plot_pre.R
# Monocle2 Visualization Preparation Pipeline
#
#' @description 该脚本准备Monocle2可视化所需的数据和图表，包括：
#'              - 加载CDS和Seurat对象
#'              - 提取随拟时序和状态变化的基因
#'              - 生成轨迹图、密度图等可视化
#'              - 保存中间结果
#'              
#' @param inputs 输入：
#'              - icds: Monocle CDS对象文件(.RDS)
#'              - iseu: Seurat对象文件(.RDS)
#'              - igenes: 关注的基因列表文件(.xls)
#'              - pmetadat_col: 用于分组的元数据列名
#' 
#' @param outputs 输出：
#'              - ocds_DGT_pseudotimegenes: 随拟时序变化的基因(.RDS)
#'              - ostates_de: 随状态变化的基因(.RDS)
#'              - oBEAM_res_list: BEAM分析结果(.RDS)
#'              - otopn: 差异表达基因(.RDS)
#'              - oplot_monocle2_trajectory_Group: 分组轨迹图(PDF)
#'              - oplot_monocle2_trajectory_State: 状态轨迹图(PDF)
#'              - oplot_monocle2_density: 伪时间密度图(PDF)
#' 
#' @return 生成可视化所需的数据和图表

suppressPackageStartupMessages({
  library(ggplotify)
  library(readxl)
  library(stringr)
  library(dplyr)
  library(ggplot2)
})

# 加载本地函数
if (!file.exists("src/monocle2_plot.R")) {
  stop("缺少本地函数文件:src/monocle2_plot.R")
}
source("src/monocle2_plot.R")

# 验证输入参数
validate_inputs <- function(args) {
  message("[1/7] 验证输入参数")
  
  required_args <- c(
    "icds", "iseu", "igenes", "pmetadat_col",
    "ocds_DGT_pseudotimegenes", "ostates_de",
    "oBEAM_res_list", "otopn",
     "oplot_monocle2_trajectory_State", "oplot_monocle2_density"
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
  file_args <- c("icds", "iseu", "igenes")
  for (arg in file_args) {
    if (!file.exists(args[[arg]])) {
      stop("输入文件不存在: ", args[[arg]])
    }
  }

  sapply(args,function(x) dir.create(dirname(x),recursive=T))

  message("✓ 输入验证通过")
}

# 加载数据
load_data <- function(args) {
  message("\n[2/7] 加载数据")
  
  tryCatch({
    # 加载CDS对象
    icds <- readRDS(args$icds)
    message("✓ CDS对象已加载: ", args$icds)
    
    # 加载Seurat对象
    seuratdata <- readRDS(args$iseu)
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
    
    return(list(
      icds = icds,
      seuratdata = seuratdata,
      igenes = igenes
    ))
  }, error = function(e) {
    stop("数据加载失败: ", e$message)
  })
}

# 准备核心分析数据
prepare_core_data <- function(icds, outputs) {
  message("\n[3/7] 准备核心分析数据")
  
  tryCatch({
    # 获取随拟时序和状态变化的基因
    preparation <- plot_gene_express_mode(icds)
    
    # 保存结果
    saveRDS(preparation$cds_DGT_pseudotimegenes, outputs$ocds_DGT_pseudotimegenes)
    message("✓ 时序相关基因已保存: ", outputs$ocds_DGT_pseudotimegenes)
    
    saveRDS(preparation$states_de, outputs$ostates_de)
    message("✓ 状态相关基因已保存: ", outputs$ostates_de)
    
    return(preparation)
  }, error = function(e) {
    stop("核心数据准备失败: ", e$message)
  })
}

# 生成轨迹图
generate_trajectory_plots <- function(icds, pmetadat_col, outputs) {
  message("\n[4/7] 生成轨迹图")
  
  tryCatch({
    # 按状态分组
    plot_state <- plot_monocle2_trajectory(icds, bycolor = "State")
    ggsave(
      outputs$oplot_monocle2_trajectory_State,
      plot = plot_state,
      width = 8,
      height = 10,create.dir=T
    )
    message("✓ 状态轨迹图已保存: ", outputs$oplot_monocle2_trajectory_State)
    
    return(list(
      plot_state = plot_state
    ))
  }, error = function(e) {
    stop("轨迹图生成失败: ", e$message)
  })
}

# 运行BEAM分析
run_beam_analysis <- function(icds, plot_state, outputs) {
  message("\n[5/7] 运行BEAM分析")
  
  tryCatch({
    BEAM_res_list <- plot_BEAM_anaylsis(icds, plot_state, all_node = NULL)
    saveRDS(BEAM_res_list, outputs$oBEAM_res_list)
    message("✓ BEAM分析结果已保存: ", outputs$oBEAM_res_list)
    
    return(BEAM_res_list)
  }, error = function(e) {
    stop("BEAM分析失败: ", e$message)
  })
}

# 准备差异表达基因
prepare_deg_data <- function(seuratdata, pmetadat_col, outputs) {
  message("\n[6/7] 准备差异表达基因")
  
  tryCatch({
    topn <- seurat_deg_topn(seu = seuratdata, ident = pmetadat_col)
    saveRDS(topn, outputs$otopn)
    message("✓ 差异表达基因已保存: ", outputs$otopn)
    
    return(topn)
  }, error = function(e) {
    stop("差异基因准备失败: ", e$message)
  })
}

# 生成密度图
generate_density_plot <- function(icds, pmetadat_col, outputs) {
  message("\n[7/7] 生成伪时间密度图")
  
  tryCatch({
    density_plot <- plot_monocle2_density(icds, bycolor = pmetadat_col)
    ggsave(
      outputs$oplot_monocle2_density,
      plot = density_plot,
      width = 8,
      height = 10,create.dir=T
    )
    message("✓ 伪时间密度图已保存: ", outputs$oplot_monocle2_density)
    
    return(density_plot)
  }, error = function(e) {
    stop("密度图生成失败: ", e$message)
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
    ocds_DGT_pseudotimegenes = args$ocds_DGT_pseudotimegenes,
    ostates_de = args$ostates_de,
    oBEAM_res_list = args$oBEAM_res_list,
    otopn = args$otopn,
    oplot_monocle2_trajectory_State = args$oplot_monocle2_trajectory_State,
    oplot_monocle2_density = args$oplot_monocle2_density
  )
  
  # 准备核心数据
  preparation <- prepare_core_data(data$icds, outputs)
  
  # 生成轨迹图
  trajectory_plots <- generate_trajectory_plots(data$icds, args$pmetadat_col, outputs)
  
  # 运行BEAM分析
  beam_results <- run_beam_analysis(data$icds, trajectory_plots$plot_state, outputs)
  
  # 准备差异基因
  deg_data <- prepare_deg_data(data$seuratdata, args$pmetadat_col, outputs)
  
  # 生成密度图
  density_plot <- generate_density_plot(data$icds, args$pmetadat_col, outputs)
  
  message("\n[SUCCESS] Monocle2可视化准备流程完成!")
}

# 命令行参数解析
if (exists("snakemake")){
  # Snakemake模式
  args <- list(
    icds = snakemake@input$icds,
    iseu = snakemake@input$iseu,
    igenes = snakemake@input$igenes,
    pmetadat_col = snakemake@params$pmetadat_col,
    ocds_DGT_pseudotimegenes = snakemake@output$ocds_DGT_pseudotimegenes,
    ostates_de = snakemake@output$ostates_de,
    oBEAM_res_list = snakemake@output$oBEAM_res_list,
    otopn = snakemake@output$otopn,
    oplot_monocle2_trajectory_State = snakemake@output$oplot_monocle2_trajectory_State,
    oplot_monocle2_density = snakemake@output$oplot_monocle2_density
  )
  run_visualization_pipeline(args)
}

#args <- list(
#    icds = "result_out/18_monocle2/Adult_Spleen/ocds.rds", 
#    iseu = "result_out/18_monocle2/Adult_Spleen/seurat_rdata.rds", 
#    igenes = "db/pyscenic/mouse/Regulon_TFs.xls", 
#    pmetadat_col = "manual_L2",    
#    ocds_DGT_pseudotimegenes = "result_out/19_pre_Plot_monocle/Adult_Spleen/ocds_DGT_pseudotimegenes.rds",
#    ostates_de = "result_out/19_pre_Plot_monocle/Adult_Spleen/ostates_de.rds",
#    oBEAM_res_list = "result_out/19_pre_Plot_monocle/Adult_Spleen/oBEAM_res_list.rds",
#    otopn = "result_out/19_pre_Plot_monocle/Adult_Spleen/otopn.rds",
#    #oplot_monocle2_trajectory_Group = "result_out/19_pre_Plot_monocle/Adult_Spleen/oplot_monocle2_trajectory_Group.pdf",
#    oplot_monocle2_trajectory_State = "result_out/19_pre_Plot_monocle/Adult_Spleen/oplot_monocle2_trajectory_State.pdf",
#    oplot_monocle2_density = "result_out/19_pre_Plot_monocle/Adult_Spleen/oplot_monocle2_density.pdf"
#)
