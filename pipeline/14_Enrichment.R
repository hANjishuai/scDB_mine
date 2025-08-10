#' 11_enrich_analysis.R
#' Gene Enrichment Analysis Pipeline
#'
#' @description 该脚本实现基因富集分析流程，包括：
#'              - 细胞亚群差异表达分析
#'              - 富集数据准备
#'              - GO富集分析
#'              - 基因表达-富集通路联合可视化
#'              - 条件特异性分析（当Condition为"ALL"时）
#'              
#' @param inputs 输入：
#'              - iSeurat: Seurat对象文件路径(.RData)
#'              - iRegulon: 关注的TF列表文件路径(.xls)
#'              - pidents: 用作细胞标识的元数据列名
#' 
#' @param outputs 输出：
#'              - oenrich_plot: 富集通路可视化PDF
#'              - ofindallmarkers_df: 全差异基因CSV
#'              - ofindallmarkers_filter_df: 过滤后差异基因CSV
#'              - oSTdata_rds: 富集分析数据RDS
#'              - oenrich_df: 富集结果CSV
#' 
#' @return 生成富集分析结果和可视化图形

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ClusterGVis)
#  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(stringr)
  library(readxl)
})

#------------------- 核心函数 -------------------#
validate_inputs <- function(inputs) {
  message("[1/7] Validating inputs")
  
  required_inputs <- c("iSeurat", "iRegulon")
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

load_seurat_data <- function(seurat_path, pidents) {
  message("\n[2/7] Loading Seurat object from: ", seurat_path)
  
  tryCatch({
    env <- new.env()
    load(seurat_path, envir = env)
    
    obj_name <- ls(env)[1]
    seurat_obj <- get(obj_name, envir = env)
    
    if (!pidents %in% colnames(seurat_obj@meta.data)) {
      stop("Invalid ident column: '", pidents, "' not found in metadata")
    }
    
    Idents(seurat_obj) <- pidents
    message("✓ Seurat object loaded")
    message("  Cells: ", ncol(seurat_obj), 
            ", Idents: ", paste(levels(seurat_obj), collapse = ", "))
    
    return(seurat_obj)
  }, error = function(e) {
    stop("Failed to load Seurat object: ", e$message)
  })
}

run_differential_expression <- function(seurat_obj, outputs_file_dir) {
  message("\n[3/7] Running differential expression analysis")
  if(!dir.exists(outputs_file_dir)){dir.create(outputs_file_dir,recursive=T)}
  tryCatch({
    # 全差异基因分析
    pbmc.markers.all <- Seurat::FindAllMarkers(
      seurat_obj,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 0.25
    )
    ofindallmarkers_df = file.path(outputs_file_dir,"ofindallmarkers_df.csv")
    write.csv(pbmc.markers.all, ofindallmarkers_df)
    
    # 过滤差异基因
    pbmc.markers <- pbmc.markers.all %>%
      dplyr::group_by(cluster) %>%
      dplyr::top_n(n = 800, wt = avg_log2FC)

    ofindallmarkers_filter_df <- file.path(outputs_file_dir,'ofindallmarkers_filter_df.csv')
    write.csv(pbmc.markers, ofindallmarkers_filter_df)
    
    message("✓ Differential expression completed")
    message("  Total markers: ", nrow(pbmc.markers.all),
            ", Filtered markers: ", nrow(pbmc.markers))
    
    return(pbmc.markers)
  }, error = function(e) {
    stop("Differential expression failed: ", e$message)
  })
}

prepare_enrichment_data <- function(seurat_obj, 
                                    pbmc.markers,  
                                    outputs_file_dir) {
  message("\n[4/7] Preparing enrichment data")
  if(!dir.exists(outputs_file_dir)){dir.create(outputs_file_dir)}

  tryCatch({
    st.data <- ClusterGVis::prepareDataFromscRNA(
      object = seurat_obj,
      assays = "SCT",
      diffData = pbmc.markers,
      keep.uniqGene = FALSE,
      showAverage = TRUE
    )
    oSTdata_rds <- file.path(outputs_file_dir,"oSTdata.RDS")
    saveRDS(st.data, oSTdata_rds)
    
    message("✓ Enrichment data prepared")
    message("  Genes: ", nrow(st.data$wide.res),
            ", Clusters: ", length(unique(st.data$long.res$cluster)))
    
    return(st.data)
  }, error = function(e) {
    stop("Enrichment data preparation failed: ", e$message)
  })
}

run_enrichment_analysis <- function(st.data, outputs_file_dir) {
  message("\n[5/7] Running enrichment analysis")
  
  tryCatch({
    enrich <- ClusterGVis::enrichCluster(
      object = st.data,
      OrgDb = org.Mm.eg.db,
      type = "BP",
      organism = "hsa",
      pvalueCutoff = 0.5,
      topn = 5,
      seed = 5201314
    )
    oenrich_df <- file.path(outputs_file_dir,"oenrich_df.csv")
    write.csv(enrich, oenrich_df)
    
    message("✓ Enrichment analysis completed")
    message("  Enriched terms: ", nrow(enrich),
            ", Gene clusters: ", length(unique(enrich$group)))
    
    return(enrich)
  }, error = function(e) {
    stop("Enrichment analysis failed: ", e$message)
  })
}

visualize_enrichment <- function(st.data, enrich, 
                                regulon_genes, 
                                outputs_file_dir,
                                figure_dir,
                                ocolor_palette=NULL) {
  message("\n[6/7] Visualizing enrichment results")
  if(!dir.exists(figure_dir)){dir.create(figure_dir,recursive=T)}
  tryCatch({
    # 提取关注的基因
    label_gene_c <- NULL
    for (i in regulon_genes) {
      a <- st.data$wide.res %>% filter(.,grepl(i,st.data$wide.res[,'gene'])) %>% dplyr::select("gene") 
      label_gene_c <- c(label_gene_c,a)
    }
    label_gene_c <- unlist(label_gene_c)
    
    # 创建可视化
    ncluster_gene <- length(unique(enrich$group))

    if(is.null(ocolor_palette)){
        color_palette <- rep(jjAnno::useMyCol("calm", n = ncluster_gene), 
                            each = 5)#(这里的5是指top_n)
        color_palette <- setNames(color_palette,enrich$group)
        ocolor_palette <- file.path(outputs_file_dir,"ocolor_palette.csv") 
        color_table <- data.frame(Clusters=names(color_palette) %>% unique,
                                  colors=color_palette %>% unique)
        write.csv(color_table,ocolor_palette)
    } else {
        color_palette <- read.csv(ocolor_palette)
        color <- color_palette$colors %>% setNames(.,color_palette$Clusters)         
        color_palette <- rep(color,each = 5)#(这里的5是指top_n)
    }

    oenrich_plot <- file.path(figure_dir,"oenrich_plot.pdf")
    pdf(oenrich_plot, height = 15, width = 14)
    ClusterGVis::visCluster(
      object = st.data,
      plot.type = "both",
      add.sampleanno = FALSE,
      column_names_rot = 90,
      show_row_dend = FALSE,
      markGenes = label_gene_c,
      annoTerm.data = enrich,
      textbar.pos = c(0.89, 0.15),
      markGenes.side = "left",
      line.side = "left",
      cluster.order = 1:ncluster_gene,
      go.col = color_palette,
      add.bar = TRUE
    )
    dev.off()
    
    message("✓ Visualization saved to: ",oenrich_plot)
  }, error = function(e) {
    stop("Visualization failed: ", e$message)
  })
}

#------------------- 主流程函数 -------------------#
run_enrich_analysis_pipeline <- function(inputs, 
                                        outputs_file_dir,
                                        figure_dir,
                                        ocolor_palette = NULL,
                                        rds_output) {
  # 输入验证
  validate_inputs(inputs)
  
  # 加载数据
  seurat_obj <- load_seurat_data(inputs$iSeurat, inputs$pidents)

  # 加载调控基因列表
  regulon <- read_xls(inputs$iRegulon, sheet = 2)
  regulon_genes <- do.call("$",list(regulon,"Tfs"))# 想要关注的基因
  
  # 主分析流程
  pbmc.markers <- run_differential_expression(seurat_obj, outputs_file_dir)
  st.data <- prepare_enrichment_data(seurat_obj, pbmc.markers, outputs_file_dir)
  enrich <- run_enrichment_analysis(st.data, outputs_file_dir)
  visualize_enrichment(st.data, enrich, regulon_genes,outputs_file_dir,figure_dir,ocolor_palette)
  
  save(seurat_obj,file=rds_output)
  message("\n[SUCCESS] Enrichment analysis pipeline completed!")
}

#------------------- 命令行/Snakemake集成 -------------------#
if (exists("snakemake")){
  # Snakemake模式
  inputs <- list(
    iSeurat = snakemake@input$iSeurat,
    iRegulon = snakemake@params$iRegulon,
    pidents = snakemake@params$pidents
  )
  outputs_file_dir = snakemake@params$outputs_file_dir
  figure_dir =  snakemake@params$figure_dir
#  ocolor_palette = snakemake@params$ocolor_palette
  rds_output = snakemake@output$rds_output

  
  # 运行主流程
  run_enrich_analysis_pipeline(inputs = inputs, 
                               outputs_file_dir = outputs_file_dir,
                               figure_dir = figure_dir,
                               ocolor_palette = NULL,
                               rds_output = rds_output
                               )
}


