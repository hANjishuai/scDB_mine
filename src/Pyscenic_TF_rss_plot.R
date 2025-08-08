# Pyscenic_TF_rss_plot.R
## 以下均是画图打包函数，使用时，请加载到运行脚本中，不要轻易改动！
## source("/pathway/to/Pyscenic_TF_rss_plot.R")
## 创建时间：2025.02.17 创建人：蓟方涵
library(SCopeLoomR)
library(cowplot)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(cowplot)
library(stringr)
library(ComplexHeatmap)
library(data.table)
library(ggplot2)
library(pheatmap)
library(ggheatmap)
library(reshape2)

#可视化细胞特异性TF气泡图,反映regulon的特异性
TF_rss_droplot <- function(
        Matrix, #行为转录因子，列为细胞种类
#        zThreshold = 2,
        cluster_columns = FALSE,
        order_rows = TRUE,
        thr=0.1,
        varName = "cellType",
        col.low = '#330066',
        col.mid = '#66CC66',
        col.high = '#FFCC33'
){
    p <-plotRSS(
            Matrix[,order(colnames(Matrix))],
#            zThreshold = zThreshold,
            cluster_columns = cluster_columns,
            order_rows = order_rows,
            thr=thr,
            varName = varName,
            col.low = col.low,
            col.mid = col.mid,
            col.high = col.high)
    return(p)
}

#可视化细胞特异性TF热图，反映TF的表达量
TF_exprewss_heatmap <- function(
    Data, #TF_rss_droplot函数返回的对象：object$plot$data，
    color=colorRampPalette(c('#91CDC8','white',"#F0988C"))(100),
    cluster_rows = F,
    cluster_cols = F,
    scale = "row",
    annotation_cols = NULL,
    annotation_color = NULL,
    legendName = "Relative value",
    #shape = 'circle',
    border = T
){
# 作图
    p <- ggheatmap(Data,
               color=color,
               cluster_rows = cluster_rows,
               cluster_cols = cluster_cols,
               scale = scale,
               annotation_cols =  annotation_cols,
               annotation_color = annotation_color,
               legendName = legendName,
             #  shape = 'circle',
               border = border) %>% 
    ggheatmap_theme(1,theme =list(
                   theme(axis.text.x = element_text(angle = 90,face = "bold",
                                                    vjust =1,hjust = 0.8))
                                                    ))
    return(p)
}

#TF_AUC与seurat结合
TF_seurat_dotplot <- function(
    seuratdata, #seurat对象
    regulonAUC, #pyscenic第三步分析中AUC结果
    TF_lists, #感兴趣的或者比较重要的转录因子向量
    group = NULL, #Seurat对象的metadata中的分类变量，一般是细胞种类或者疾病状态
    direction = "horizontal",
    position = "bottom",
    gene_c = NULL #关注的基因
){
p1 <- DotPlot(seuratdata, features = TF_lists, group.by = group)+
        theme_bw()+
        theme(panel.grid = element_blank(),
                axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
        theme(legend.direction = direction ,
                legend.position = position)+
        labs(x=NULL,y=NULL)
p2 <- NULL
if(!is.null(gene_c)){
    p2 <- FeaturePlot(seuratdata,gene_c)
    print("featureplot已完成！")
    }
return(list(p1,p2))
}

#可视化所有细胞特异性TF的表达热图
TF_regulonActivity_heatmap <- function(
    data,#经过标准化的AUC矩阵
    name="Regulon activity",
    show_heatmap_legend = F,
    cluster_rows = T,
    show_row_dend = F,
    cluster_columns = F,
    fontsize=5,
    show_row_names = F,
    right_annotation = NULL #rowAnnotation注释
){
    hm <- draw(ComplexHeatmap::Heatmap(data, 
                                   name=name,
                                   show_heatmap_legend = show_heatmap_legend,
                                   cluster_rows = cluster_rows,
                                   show_row_dend = show_row_dend,
                                   cluster_columns = cluster_columns,
                                   column_names_gp = grid::gpar(fontsize=fontsize),
                                   show_row_names = show_row_names,
                                   right_annotation = right_annotation)) 
    return(hm)
}


#rank可视化rss
TF_rssRank_plot <- function(
    rssobj,
    spec_celltype=NULL,#关注的细胞类型
    seqnum=6
){
#rss特异性TF结果
B_rss <- as.data.frame(rssobj)
if(!is.null(spec_celltype)){
    celltype <- spec_celltype
} else {
    celltype <- colnames(B_rss)
}

rssRanklist <- list()
for(i in 1:length(celltype)) {
  #提取数据
  data_rank_plot <- cbind(as.data.frame(rownames(B_rss)),
                          as.data.frame(B_rss[,celltype[i]]))
  
  colnames(data_rank_plot) <- c("TF", "celltype")
  data_rank_plot=na.omit(data_rank_plot)#去除NA
  data_rank_plot <- data_rank_plot[order(data_rank_plot$celltype,decreasing=T),]#降序排列
  data_rank_plot$rank <- seq(1, nrow(data_rank_plot))#添加排序
  
  p <- ggplot(data_rank_plot, aes(x=rank, y=celltype)) + 
    geom_point(size=1, shape=16, color="#1F77B4",alpha =0.4)+
    geom_point(data = data_rank_plot[1:seqnum,],
               size=1, color='#DC050C')+ #选择前6个标记，自行按照需求选择
    theme_bw()+
    theme(axis.title = element_text(colour = 'black', size = 12),
          axis.text = element_text(colour = 'black', size = 10),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    labs(x='Regulons Rank', y='Specificity Score',title =celltype[i])+
    geom_text_repel(data= data_rank_plot[1:6,],max.overlaps =15 ,
                    aes(label=TF), color="black", size=3.5, fontface="italic", 
                    arrow = arrow(ends="first", length = unit(0.01, "npc")), box.padding = 0.2,
                    point.padding = 0.3, segment.color = 'black', 
                    segment.size = 0.3, force = 1, max.iter = 3e3)
  rssRanklist[[celltype[i]]] <- p
}
# 获取rssRanklist的元素数量
n <- length(rssRanklist)
# 计算每行每列的图片数量（尽量接近正方形）
side <- ceiling(sqrt(n))
# 使用plot_grid函数排列图片
p_list <- plot_grid(plotlist = rssRanklist, ncol = side, nrow = side)
result_list <- list(rankplot_Rdata=rssRanklist,rankplot=p_list)
return(result_list)
}