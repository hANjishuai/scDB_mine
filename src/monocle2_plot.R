# monocle2_plot.R
# Monocle2 Visualization Functions with Enhanced Robustness

# 加载必要的库
suppressPackageStartupMessages({
  library(monocle)
  library(ggsci)
  library(patchwork)
  library(Seurat)
  library(ggplot2)
  library(tidydr)
  library(stringr)
  library(ggforce)
  library(ggrastr)
  library(ggpubr)
  library(fgsea)
  library(scales)
  library(dplyr)
  library(ggridges)
  library(RColorBrewer)
  library(viridis)
})

plot_monocle2_trajectory <- function(cds, bycolor = "State") {
  #' 绘制发育轨迹图
  #' 
  #' @param cds Monocle CDS对象
  #' @param bycolor 分组颜色依据的元数据列名
  #' @return ggplot对象
  
  tryCatch({
    plot <- plot_cell_trajectory(
      cds, 
      color_by = bycolor,
      size = 1,
      show_backbone = TRUE
    ) +
      scale_color_npg() +
      guides(colour = guide_legend("")) +
      theme(
        legend.text = element_text(size = 12, face = "bold"),
        legend.key.size = unit(20, "pt"),
        legend.position = "right"
      )
    
    message("如需分面展示，可添加: + facet_wrap(facets = 'condition')")
    return(plot)
  }, error = function(e) {
    stop("轨迹图绘制失败: ", e$message)
  })
}

plot_monocle2_density <- function(cds, bycolor = "State") {
  #' 绘制发育轨迹密度图
  #' 
  #' @param cds Monocle CDS对象
  #' @param bycolor 分组颜色依据的元数据列名
  #' @return ggplot对象
  
  tryCatch({
    df <- pData(cds)
    plot <- ggplot(df, aes(Pseudotime, colour = !!sym(bycolor), fill = !!sym(bycolor))) + 
      geom_density(bw = 0.5, linewidth = 1, alpha = 0.5) + 
      theme_classic2() +
      guides(fill = guide_legend(""), colour = guide_legend(""))
    
    message("如需分面展示，可添加: + facet_wrap(facets = 'condition')")
    return(plot)
  }, error = function(e) {
    stop("密度图绘制失败: ", e$message)
  })
}

plot_gene_express_mode <- function(cds) {
  #' 识别随拟时序和状态变化的基因
  #' 
  #' @param cds Monocle CDS对象
  #' @return 包含时序相关基因和状态相关基因的列表
  
  tryCatch({
    plot_result <- list()
    
    # 过滤低表达基因
    expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
    if (length(expressed_genes) == 0) {
      stop("没有找到在任何细胞中表达的基因")
    }
    
    # 识别随伪时间变化的基因
    cds_DGT_pseudotimegenes <- tryCatch({
      differentialGeneTest(
        cds[expressed_genes, ],
        fullModelFormulaStr = "~sm.ns(Pseudotime)"
      )
    }, error = function(e) {
      warning("伪时间差异基因检测失败: ", e$message)
      return(NULL)
    })
    
    if (!is.null(cds_DGT_pseudotimegenes)) {
      cds_DGT_pseudotimegenes <- cds_DGT_pseudotimegenes[order(cds_DGT_pseudotimegenes$qval), ]
    }
    
    # 识别随状态变化的基因
    states_de <- tryCatch({
      differentialGeneTest(
        cds[expressed_genes, ],
        fullModelFormulaStr = "~State"
      )
    }, error = function(e) {
      warning("状态差异基因检测失败: ", e$message)
      return(NULL)
    })
    
    if (!is.null(states_de)) {
      states_de <- states_de[order(states_de$qval), ]
    }
    
    plot_result[["cds_DGT_pseudotimegenes"]] <- cds_DGT_pseudotimegenes
    plot_result[["states_de"]] <- states_de
    
    return(plot_result)
  }, error = function(e) {
    stop("基因表达模式分析失败: ", e$message)
  })
}

plot_BEAM_anaylsis <- function(cds, plot_trajectory, all_node = NULL) {
  #' 执行分支点分析
  #' 
  #' @param cds Monocle CDS对象
  #' @param plot_trajectory 轨迹图对象
  #' @param all_node 指定分析的分支点
  #' @return BEAM分析结果列表
  
  tryCatch({
    BEAM_res_list <- list()
    
    if (is.null(all_node)) {
      nodes <- tryCatch({
        plot_trajectory$plot_env$branch_point_df$branch_point_idx
      }, error = function(e) {
        warning("无法从轨迹图中提取分支点: ", e$message)
        return(integer(0))
      })
      
      if (length(nodes) == 0) {
        warning("未找到分支点，尝试默认分支点1")
        nodes <- 1
      }
    } else {
      nodes <- as.integer(all_node)
    }
    
    for (i in nodes) {
      message("分析分支点 ", i, " ...")
      BEAM_res <- tryCatch({
        BEAM(
          cds,
          branch_point = i,
          cores = 1,
          progenitor_method = 'duplicate'
        )
      }, error = function(e) {
        warning("分支点 ", i, " 分析失败: ", e$message)
        return(NULL)
      })
      
      if (!is.null(BEAM_res)) {
        BEAM_res_list[[as.character(i)]] <- BEAM_res
      }
    }
    
    if (length(BEAM_res_list) == 0) {
      warning("所有分支点分析均失败")
    }
    
    return(BEAM_res_list)
  }, error = function(e) {
    stop("分支点分析失败: ", e$message)
  })
}

seurat_deg_topn <- function(seu, ident, ntop = 20) {
  #' 识别差异表达基因
  #' 
  #' @param seu Seurat对象
  #' @param ident 分组依据的元数据列名
  #' @param ntop 每组取前N个差异基因
  #' @return 差异表达基因数据框
  
  tryCatch({
    message("确认Idents已设置为: ", ident)
    Idents(seu) <- ident
    
    marker <- tryCatch({
      FindAllMarkers(
        seu,
        only.pos = TRUE,
        logfc.threshold = 0.5
      )
    }, error = function(e) {
      stop("差异表达分析失败: ", e$message)
    })
    
    marker <- marker[which(marker$p_val_adj < 0.05), ]
    
    if (nrow(marker) == 0) {
      warning("未找到显著差异表达基因")
      return(data.frame())
    }
    
    topn <- marker %>% 
      group_by(cluster) %>% 
      top_n(n = ntop, wt = avg_log2FC)
    
    return(topn)
  }, error = function(e) {
    stop("差异基因分析失败: ", e$message)
  })
}

plot_pseudotimeHeatmap <- function(cds, topn, states_de, special_gene = NULL, ngene = NULL, n_cluster = 3) {
  #' 修复颜色映射问题的伪时间热图函数
  #' 
  #' @param cds Monocle CDS对象
  #' @param topn 差异表达基因数据框
  #' @param states_de 状态相关基因数据框
  #' @param special_gene 需要特别关注的基因列表
  #' @param ngene 限制显示的基因数量
  #' @param n_cluster 基因聚类数
  #' @return 热图对象或替代热图
  
  tryCatch({
    # 1. 获取基因列表
    if (!is.null(states_de) && nrow(states_de) > 0 && "gene_short_name" %in% colnames(states_de)) {
      state_genes <- states_de[rownames(states_de) %in% topn$gene, "gene_short_name"] %>% 
        as.character() %>% 
        unique() %>% 
        na.omit()
    } else {
      state_genes <- unique(topn$gene)
    }
    
    # 添加特别关注的基因
    if (!is.null(special_gene)) {
      valid_special_genes <- special_gene[special_gene %in% rownames(fData(cds))]
      message("添加关注的基因: ", paste(valid_special_genes, collapse = ", "))
      state_genes <- unique(c(valid_special_genes, state_genes))
    }
    
    # 过滤存在的基因
    valid_genes <- state_genes[state_genes %in% rownames(fData(cds))]
    if (length(valid_genes) == 0) {
      stop("没有有效的基因用于热图绘制")
    }
    
    # 限制基因数量
    if (!is.null(ngene) && ngene > 0) {
      valid_genes <- head(valid_genes, min(ngene, length(valid_genes)))
    }
    message("最终用于热图的基因数量: ", length(valid_genes))
    
    # 2. 提取表达数据
    expr_data <- exprs(cds)[valid_genes, , drop = FALSE]
    
    # 3. 处理伪时间
    pseudotime <- pData(cds)$Pseudotime
    if (any(is.na(pseudotime))) {
      pseudotime[is.na(pseudotime)] <- mean(pseudotime, na.rm = TRUE)
    }
    
    # 4. 按伪时间排序
    time_order <- order(pseudotime)
    expr_data <- expr_data[, time_order]
    pseudotime <- pseudotime[time_order]
    
    # 5. 关键修复：正确标准化和颜色映射
    # 计算每行（基因）的z-score
    expr_scale <- t(apply(expr_data, 1, function(x) {
      if (sd(x) == 0) {
        rep(0, length(x)) # 避免除以0
      } else {
        (x - mean(x)) / sd(x)
      }
    }))
    
    # 设置颜色映射的断点
    breaks <- seq(-2, 2, length.out = 100)
    
    # 6. 使用pheatmap创建热图
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      message("使用pheatmap创建热图")
      
      # 创建列注释（伪时间）
      annotation_col <- data.frame(Pseudotime = pseudotime)
      rownames(annotation_col) <- colnames(expr_scale)
      
      # 创建热图
      heatmap_plot <- pheatmap::pheatmap(
        expr_scale,
        cluster_cols = FALSE,
        show_colnames = FALSE,
        show_rownames = ifelse(length(valid_genes) <= 50, TRUE, FALSE),
        annotation_col = annotation_col,
        annotation_colors = list(
          Pseudotime = grDevices::colorRampPalette(c("navy", "yellow", "red"))(50)
        ),
        color = grDevices::colorRampPalette(c("navy", "white", "firebrick3"))(100),
        breaks = breaks,  # 关键：设置明确的断点
        cutree_rows = n_cluster,
        silent = TRUE,
        main = "Pseudotime Expression Heatmap"
      )
      
      return(heatmap_plot)
    }
    
    # 7. 使用基础的ggplot2作为备选
    message("使用ggplot2创建热图")
    
    # 准备数据
    plot_data <- reshape2::melt(
      expr_scale,  # 使用标准化后的数据
      varnames = c("Gene", "Cell"),
      value.name = "Expression"
    )
    
    plot_data$Pseudotime <- pseudotime[plot_data$Cell]
    
    # 按伪时间排序基因
    gene_means <- plot_data %>% 
      group_by(Gene) %>% 
      summarise(mean_expr = mean(Expression)) %>% 
      arrange(mean_expr)
    
    plot_data$Gene <- factor(plot_data$Gene, levels = gene_means$Gene)
    
    # 创建热图
    heatmap_plot <- ggplot(plot_data, aes(x = Pseudotime, y = Gene, fill = Expression)) +
      geom_tile() +
      scale_fill_gradientn(
        colors = c("navy", "blue", "white", "red", "firebrick3"),
        values = scales::rescale(c(-2, -1, 0, 1, 2)),
        limits = c(-2, 2),  # 关键：设置明确的颜色范围
        na.value = "grey80"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(angle = 90, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "right"
      ) +
      labs(
        x = "Pseudotime", 
        y = "Gene",
        title = "Pseudotime Expression Heatmap",
        fill = "Z-score"
      )
    
    return(heatmap_plot)
    
  }, error = function(e) {
    warning("热图方法失败: ", e$message)
    return(NULL)
  })
}



#plot_pseudotimeHeatmap <- function(cds, topn, states_de, special_gene = NULL, ngene = NULL, n_cluster = 3) {
#  #' 绘制伪时间热图 - 修复版
#  #' 
#  #' @param cds Monocle CDS对象
#  #' @param topn 差异表达基因数据框
#  #' @param states_de 状态相关基因数据框
#  #' @param special_gene 需要特别关注的基因列表
#  #' @param ngene 限制显示的基因数量
#  #' @param n_cluster 基因聚类数
#  #' @return 热图对象
#  
#  tryCatch({
#    # 验证输入数据
#    if (is.null(states_de) || nrow(states_de) == 0) {
#      stop("状态相关基因数据为空")
#    }
#    
#    if (is.null(topn) || nrow(topn) == 0) {
#      stop("差异表达基因数据为空")
#    }
#    
#    # 获取状态相关基因
#    if (!"gene_short_name" %in% colnames(states_de)) {
#      warning("states_de缺少gene_short_name列，使用行名代替")
#      state_genes <- rownames(states_de)
#    } else {
#      state_genes <- states_de[rownames(states_de) %in% topn$gene, "gene_short_name"] %>% 
#        as.character() %>% 
#        unique() %>% 
#        na.omit()
#    }
#    
#    if (length(state_genes) == 0) {
#      warning("未找到共同的基因，使用前100个状态基因")
#      state_genes <- rownames(states_de)[1:min(100, nrow(states_de))]
#    }
#    
#    # 添加特别关注的基因
#    if (!is.null(special_gene)) {
#      valid_special_genes <- special_gene[special_gene %in% rownames(fData(cds))]
#      message("添加关注的基因: ", paste(valid_special_genes, collapse = ", "))
#      state_genes <- unique(c(valid_special_genes, state_genes))
#    }
#    
#    # 过滤不存在的基因
#    valid_genes <- state_genes[state_genes %in% rownames(fData(cds))]
#    if (length(valid_genes) == 0) {
#      stop("没有有效的基因用于热图绘制")
#    }
#    
#    # 限制基因数量
#    if (!is.null(ngene) && ngene > 0) {
#      valid_genes <- head(valid_genes, min(ngene, length(valid_genes)))
#    }
#    
#    message("最终用于热图的基因数量: ", length(valid_genes))
#    
#    # 检查并清理表达数据
#    expr_mat <- exprs(cds)[valid_genes, , drop = FALSE]
#    
#    # 移除有NA或Inf的基因
#    na_genes <- apply(expr_mat, 1, function(x) any(is.na(x) | is.infinite(x)))
#    if (any(na_genes)) {
#      warning("移除包含NA/Inf的基因: ", paste(valid_genes[na_genes], collapse = ", "))
#      valid_genes <- valid_genes[!na_genes]
#      
#      if (length(valid_genes) == 0) {
#        stop("所有基因都包含无效值")
#      }
#    }
#    
#    # 移除零方差基因
#    zero_var <- apply(expr_mat, 1, var) == 0
#    if (any(zero_var)) {
#      warning("移除零方差基因: ", paste(valid_genes[zero_var], collapse = ", "))
#      valid_genes <- valid_genes[!zero_var]
#      
#      if (length(valid_genes) == 0) {
#        stop("所有基因都是零方差")
#      }
#    }
#    
#    # 确保伪时间没有NA
#    pseudotime <- pData(cds)$Pseudotime
#    if (any(is.na(pseudotime))) {
#      warning("伪时间包含NA值，使用均值填充")
#      pseudotime[is.na(pseudotime)] <- mean(pseudotime, na.rm = TRUE)
#      pData(cds)$Pseudotime <- pseudotime
#    }
#    
#    # 尝试绘制热图 - 添加稳定性选项
#    heatmap_plot <- tryCatch({
#      monocle::plot_pseudotime_heatmap(
#        cds[valid_genes, ],
#        num_cluster = n_cluster,
#        show_rownames = TRUE,
#        return_heatmap = TRUE,
#        hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100),
#        cores = 1,  # 禁用并行计算
#        use_gene_short_name = TRUE
#      )
#    }, error = function(e) {
#      # 尝试简化方法作为备选
#      warning("标准热图方法失败，尝试替代方法: ", e$message)
#      
#      # 使用更稳定的平滑方法
#      monocle::plot_pseudotime_heatmap(
#        cds[valid_genes, ],
#        num_cluster = n_cluster,
#        show_rownames = TRUE,
#        return_heatmap = TRUE,
#        hmcols = colorRampPalette(c("navy", "white", "firebrick3"))(100),
#        cores = 1,
#        use_gene_short_name = TRUE,
#        trend_formula = "~ splines::ns(Pseudotime, df=3)"  # 使用更稳定的平滑函数
#      )
#    })
#    
#    return(heatmap_plot)
#  }, error = function(e) {
#    stop("热图绘制失败: ", e$message)
#  })
#}

plot_genes_branchedHeatmap <- function(cds, BEAM_res, special_gene = NULL, ngene = 50, n_cluster = 3, branch_p = 1) {
  #' 绘制分支热图
  #' 
  #' @param cds Monocle CDS对象
  #' @param BEAM_res BEAM分析结果
  #' @param special_gene 特别关注的基因
  #' @param ngene 显示的基因数量
  #' @param n_cluster 基因聚类数
  #' @param branch_p 分支点编号
  #' @return 热图对象
  
  tryCatch({
    if (is.null(BEAM_res) || nrow(BEAM_res) == 0) {
      stop("BEAM分析结果为空")
    }
    
    # 获取显著基因
    sig_genes <- rownames(subset(BEAM_res, qval < 1e-4))
    if (length(sig_genes) == 0) {
      warning("没有qval < 1e-4的显著基因")
      sig_genes <- rownames(BEAM_res)[1:min(ngene, nrow(BEAM_res))]
    }
    
    # 限制基因数量
    Show_gene <- head(sig_genes, min(ngene, length(sig_genes)))
    
    # 添加特别关注的基因
    if (!is.null(special_gene)) {
      valid_special_genes <- special_gene[special_gene %in% rownames(fData(cds))]
      message("添加关注的基因: ", paste(valid_special_genes, collapse = ", "))
      Show_gene <- unique(c(valid_special_genes, Show_gene))
    }
    
    # 过滤不存在的基因
    valid_genes <- Show_gene[Show_gene %in% rownames(fData(cds))]
    if (length(valid_genes) == 0) {
      stop("没有有效的基因用于分支热图")
    }
    
    # 绘制分支热图
    plot <- monocle::plot_genes_branched_heatmap(
      cds[valid_genes, ],
      branch_point = branch_p,
      num_clusters = n_cluster,
      cores = 1,
      hmcols = colorRampPalette(rev(brewer.pal(9, "PRGn")))(65),
      branch_colors = c("#979797", "#F05662", "#7990C8"),
      use_gene_short_name = TRUE,
      show_rownames = TRUE,
      return_heatmap = TRUE
    )
    
    return(plot$ph_res)
  }, error = function(e) {
    stop("分支热图绘制失败: ", e$message)
  })
}

Genexpress_changes_over_pseudotime <- function(cds, genes, meta_col = 'orig.ident') {
  #' 绘制基因表达随时间变化图
  #' 
  #' @param cds Monocle CDS对象
  #' @param genes 关注的基因列表
  #' @param meta_col 分组元数据列名
  #' @return ggplot对象
  
  tryCatch({
    # 验证基因列表
    valid_genes <- genes[genes %in% rownames(fData(cds))]
    if (length(valid_genes) == 0) {
      stop("没有有效的基因可用于绘图")
    }
    
    message("拟合基因: ", paste(valid_genes, collapse = ", "))
    
    # 创建基因表达矩阵
    gene_exp <- as.data.frame(t(log2(exprs(cds)[valid_genes, , drop = FALSE] + 1)))
    colnames(gene_exp) <- valid_genes
    
    # 合并元数据
    meta_data <- pData(cds)[, c(meta_col, "Pseudotime"), drop = FALSE]
    combined_data <- cbind(meta_data, gene_exp)
    
    # 转换为长格式
    data_long <- reshape2::melt(
      combined_data,
      id.vars = c(meta_col, "Pseudotime"),
      variable.name = 'gene',
      value.name = 'value'
    )
    
    # 分组数量
    num_groups <- length(unique(data_long[[meta_col]]))
    
    # 创建颜色调色板
    if (num_groups <= 8) {
      colors <- RColorBrewer::brewer.pal(max(3, num_groups), "Set1")[1:num_groups]
    } else {
      colors <- viridis::viridis(num_groups)
    }
    
    # 绘制图形
    plot <- ggplot(data_long, aes(x = Pseudotime, y = value, color = .data[[meta_col]])) +
      geom_smooth(aes(fill = .data[[meta_col]]), se = FALSE) +
      xlab('Pseudotime') + 
      ylab('Relative expression') +
      facet_wrap(~ gene, scales = "free_y") +
      theme_bw() +
      theme(
        axis.text = element_text(color = 'black', size = 10),
        axis.title = element_text(color = 'black', size = 12),
        strip.text = element_text(color = 'black', size = 10),
        legend.position = "bottom"
      ) +
      scale_color_manual(values = colors) +
      scale_fill_manual(values = colors)
    
    return(plot)
  }, error = function(e) {
    stop("基因表达变化图绘制失败: ", e$message)
  })
}

pathwayenrich_changes_over_pseudotime <- function(cds, seuratdata, gmt_list_path, meta_col = 'orig.ident') {
  #' 绘制通路富集随时间变化图
  #' 
  #' @param cds Monocle CDS对象
  #' @param seuratdata Seurat对象
  #' @param gmt_list_path GMT文件路径
  #' @param meta_col 分组元数据列名
  #' @return 包含绘图结果的列表
  
  tryCatch({
    # 检查GMT文件路径
    if (!dir.exists(gmt_list_path)) {
      stop("GMT文件目录不存在: ", gmt_list_path)
    }
    
    pathway_files <- list.files(gmt_list_path, full.names = TRUE, pattern = "\\.gmt$")
    if (length(pathway_files) == 0) {
      stop("未找到GMT文件: ", gmt_list_path)
    }
    
    # 添加伪时间到Seurat对象
    common_cells <- intersect(colnames(cds), colnames(seuratdata))
    if (length(common_cells) == 0) {
      stop("CDS和Seurat对象中没有共同的细胞")
    }
    
    seuratdata <- seuratdata[, common_cells]
    seuratdata$pseudotime <- cds$Pseudotime[match(common_cells, colnames(cds))]
    
    plots_list <- list()
    for (pathway_file in pathway_files) {
      pathway_name <- tools::file_path_sans_ext(basename(pathway_file))
      
      tryCatch({
        pathways <- fgsea::gmtPathways(pathway_file)
        if (length(pathways) == 0) {
          warning("空通路文件: ", pathway_file)
          next
        }
        
        # 计算通路评分
        seuratdata <- tryCatch({
          Seurat::AddModuleScore(
            object = seuratdata,
            features = pathways,
            name = pathway_name,
            pool = rownames(seuratdata),
            nbin = 24,
            ctrl = min(vapply(pathways, length, integer(1)))
        )}, error = function(e) {
          warning("通路评分计算失败: ", pathway_file, " - ", e$message)
          return(NULL)
        })
        
        if (is.null(seuratdata)) next
        
        # 提取数据
        score_col <- grep(paste0("^", pathway_name, "\\d+$"), colnames(seuratdata@meta.data), value = TRUE)
        if (length(score_col) == 0) {
          warning("未找到通路评分列: ", pathway_name)
          next
        }
        
        plot_data <- seuratdata@meta.data[, c(meta_col, "pseudotime", score_col)]
        colnames(plot_data)[3] <- "pathway_score"
        
        # 绘制通路图
        plot <- ggplot(plot_data, aes(x = pseudotime, y = pathway_score, color = .data[[meta_col]])) +
          geom_smooth(se = FALSE) +
          labs(title = pathway_name, x = "Pseudotime", y = "Pathway Score") +
          theme_minimal()
        
        plots_list[[pathway_name]] <- plot
      }, error = function(e) {
        warning("通路处理失败: ", pathway_file, " - ", e$message)
      })
    }
    
    if (length(plots_list) == 0) {
      warning("未生成任何通路图")
      return(list())
    }
    
    return(plots_list)
  }, error = function(e) {
    stop("通路富集分析失败: ", e$message)
  })
}
