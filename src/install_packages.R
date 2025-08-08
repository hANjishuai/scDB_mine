#!/usr/bin/env Rscript

# 设置超时时间（秒），解决大包下载问题
# 如果用的是系统R，要下载到系统的R库中请运行以下命令：
# .libPaths("~/R/library/") # 如果这是你系统R制定的包位置
# pak::pkg_install("xxx")
# CRAN 包：直接写包名
# pak::pkg_install("ggplot2")
# 
# Bioconductor 包：加 Bioconductor/ 前缀
# pak::pkg_install("Bioconductor/SingleR")
# 
# GitHub 包：加 用户名/仓库
# pak::pkg_install("satijalab/seurat")

options(timeout = 600)

# 安装 CRAN 包
cran_packages <- c(
  "Seurat", "patchwork", "ggplot2", "DoubletFinder", "ggrepel", "stringr", "dplyr",
  "gdata", "readxl", "argparse", "cowplot", "plotly", "data.table", "pheatmap",
  "ggheatmap", "reshape2", "tibble", "getopt", "parallel", "circlize", "ggplotify",
  "tidydr", "ggforce", "ggrastr", "ggpubr", "scales", "ggridges", "viridis", "future"
)

# 安装 Bioconductor 包
bioc_packages <- c(
  "glmGamPoi", "SingleR", "celldex", "SCopeLoomR", "AUCell","RColorBrewer",
  "BiocParallel", "ComplexHeatmap", "org.Hs.eg.db", "org.Mm.eg.db", "clusterProfiler",
  "enrichplot", "DOSE", "GO.db","pathview", "msigdbr", "GSVA", "monocle"
)

# 安装 GitHub 包
github_packages <- list(
  "SeuratData" = "satijalab/seurat-data",
  "harmony" = "immunogenomics/harmony",
  "ClusterGVis" = "junjunlab/ClusterGVis",
  "GseaVis" = "junjunlab/GseaVis"
)

# 1. 安装 CRAN 包
to_install_cran <- cran_packages[!cran_packages %in% installed.packages()[, "Package"]]
if (length(to_install_cran) > 0) {
  message("安装 CRAN 包: ", paste(to_install_cran, collapse = ", "))
  install.packages(
    to_install_cran,
#    lib = conda_lib_path,
#    repos = "https://cloud.r-project.org",
    dependencies = TRUE
  )
} else {
  message("所有 CRAN 包已安装")
}

# 2. 安装 Bioconductor 包
if (length(bioc_packages) > 0) {
  # 检查并安装 BiocManager
  if (!require("BiocManager", character.only = TRUE, quietly = TRUE)) {
    message("安装 BiocManager...")
    install.packages("BiocManager")
  }
  
  to_install_bioc <- bioc_packages[!bioc_packages %in% installed.packages()[, "Package"]]
  if (length(to_install_bioc) > 0) {
    message("安装 Bioconductor 包: ", paste(to_install_bioc, collapse = ", "))
    BiocManager::install(
      to_install_bioc,
#     lib = conda_lib_path,
      site_repository = "https://bioconductor.org/packages/release/bioc",
      update = FALSE,
      ask = FALSE
    )
  } else {
    message("所有 Bioconductor 包已安装")
  }
}

# 3. 安装 GitHub 包
if (!require("remotes", character.only = TRUE)) {
  message("安装 remotes 包...")
  install.packages("remotes")
}

to_install_github <- names(github_packages)[!names(github_packages) %in% installed.packages()[, "Package"]]
if (length(to_install_github) > 0) {
  message("安装 GitHub 包: ", paste(to_install_github, collapse = ", "))
  for (pkg in to_install_github) {
    message("安装 ", pkg, " (", github_packages[[pkg]], ")")
    remotes::install_github(
      github_packages[[pkg]],
#      lib = conda_lib_path,
#      dependencies = TRUE,
#      upgrade = "never"
    )
  }
} else {
  message("所有 GitHub 包已安装")
}

# 4. 验证安装
all_packages <- c(cran_packages, bioc_packages, names(github_packages))
success <- sapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(require(pkg, character.only = TRUE, quietly = TRUE))
})

if (all(success)) {
  message("✅ 所有包成功安装并加载！")
} else {
  failed_pkgs <- names(success)[!success]
  message("❌ 以下包安装失败: ", paste(failed_pkgs, collapse = ", "))
  message("请尝试手动安装: ")
#  message("install.packages(c('", paste(failed_pkgs, collapse = "','"), "'), lib = '", conda_lib_path, "')")
}

message("安装路径: ", .libPaths())
