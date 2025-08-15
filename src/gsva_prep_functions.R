get_genesets <- function(pathway_list,
                         species = "Mus musculus",
                         category = "C5",
                         subcategory = "GO:BP") {
  
  # 获取完整的基因集数据库
  GO_df_ALL <- msigdbr(
    species = species,
    category = category,
    subcategory = subcategory
  )  
  
  # 选择需要的列
  GO_df <- dplyr::select(
    GO_df_ALL,
    gene_symbol,
    gs_exact_source,
    gs_name
  )
  
  # 匹配通路ID - 修复结构问题
  geneset_list <- lapply(pathway_list, function(i) {
    GO_SUB <- GO_df %>% 
      filter(gs_exact_source == i) %>%
      distinct(gene_symbol, .keep_all = TRUE)
    
    # 返回基因向量，使用通路名称作为列表元素名
    setNames(list(unique(GO_SUB$gene_symbol)), 
             unique(GO_SUB$gs_name))
  })
  
  # 扁平化列表结构
  flattened_list <- unlist(geneset_list, recursive = FALSE)
  
  # 检查结果
  if (length(flattened_list) == 0) {
    warning("未匹配到任何基因集，请检查通路ID")
  }
  
  return(flattened_list)
}
