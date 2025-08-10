# Snakefile

# 加载配置文件
configfile: "config.yaml"

# 定义通配符
# 获取所有组织类型（Non 10X）
tissue_non10x = config["tissue_non10x"]
tissue_all = config["tissue_all"]

# 定义规则
# 读取表达矩阵（10X format）
#rule all_rawdata_merge:
#    input:
#        touch("result_out/01_rawdata_merge/Adult_PeritonealCavity/completed.txt")
#
#rule rawdata_merge:
#    input:
#        r_script = "pipeline/01_rawdata_merge.R"
#    output:
#        touch("result_out/01_rawdata_merge/Adult_PeritonealCavity/completed.txt")
#    params:
#        rawdatadir = config['input01']["rawdatadir"].format(tissue=config["tissue"]),
#        projectname = config['params01']['projectname'].format(tissue=config["tissue"]),
#        samples = config['params01']['samples'],
#        species = config['params01']['species'],
#        integrated_pdf = config['output01']['integrated_pdf'].format(tissue=config["tissue"]),
#        seurat_rdata = config['output01']['seurat_rdata'].format(tissue=config["tissue"])
#    script:
#        "pipeline/01_rawdata_merge.R"
#
## 读取表达矩阵（Non 10X format）
#rule all_rawdata_merge_non10x:
#    input:
#        expand(config["output01_non"]["integrated_pdf"], tissue_non10x=tissue_non10x),
#        expand(config["output01_non"]["seurat_rdata"], tissue_non10x=tissue_non10x)
#
#rule rawdata_merge_non10x:
#    input:
#        rawdatadir = lambda wildcards: config["input01_non"]["rawdatadir"].format(tissue_non10x=wildcards.tissue_non10x)
#    output:
#        integrated_pdf = config["output01_non"]["integrated_pdf"],
#        seurat_rdata = config["output01_non"]["seurat_rdata"]
#    params:
#        projectname = lambda wildcards: config["params01_non"]["projectname_prefix"] + wildcards.tissue_non10x,
#        species = config['params01_non']['species']
#    script:
#        "pipeline/01_rawdata_merge_non10x.R"
#
## 质控+标准化表达矩阵
#rule all_scRNA_Qc2Sct:
#    input:
#        expand(config["output02"]["seurat_rdata"], tissue_all=tissue_all)
#
#rule scRNA_Qc2Sct:
#    input:
#        seurat_obj = lambda wildcards: config["input02"]["seurat_rdata"].format(tissue_all=wildcards.tissue_all)
#    output:
#        output_rdata = config["output02"]["seurat_rdata"]
#    params:
#        species = config["params02"]["species"],
#        nFeature_RNA_min = config["params02"]["nFeature_RNA_min"],
#        nFeature_RNA_max = config["params02"]["nFeature_RNA_max"],
#        mt_max = config["params02"]["mt_max"],
#        group_var = config["params02"]["group_var"],
#        output_pdf_dir = lambda wildcards: config["params02"]["output_pdf_dir"].format(tissue_all=wildcards.tissue_all)
#    script:
#        "pipeline/02_scRNA_Qc2Sct.R"
#
## 聚类分群
#rule all_cluster:
#    input:
#        expand(config["output03"]["seurat_rdata"], tissue_all=tissue_all)
#
#rule cluster:
#    input:
#        seurat_path = lambda wildcards: config["input03"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#    output:
#        seurat_rdata = config["output03"]["seurat_rdata"]
#    params:
#        orig_ident = config["params03"]["orig_ident"],
#        result_out_dir = config["params03"]["result_out_dir"],
#        figure_out_dir = config["params03"]["figure_out_dir"],
#        cluster_resolution = config["params03"]["cluster_resolution"]
#    script:
#        "pipeline/03_cluster.R"
#
## 亚群标识基因展示
#rule all_CellmarkersPlot:
#    input:
#        expand(config["output04"]["all_markers_plot"], tissue_all=tissue_all),
#        expand(config["output04"]["per_celltype_plot"], tissue_all=tissue_all)
#
#rule CellmarkersPlot:
#    input:
#        seurat_path = lambda wildcards: config["input04"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#    output:
#        all_markers_plot = config["output04"]["all_markers_plot"],
#        per_celltype_plot = config["output04"]["per_celltype_plot"]
#    params:
#        marker_file = lambda wildcards: config["params04"]["marker_file"].format(tissue_all=wildcards.tissue_all)
#    script:
#        "pipeline/04_CellmarkersPlot.R"
#
## 亚群注释level1
#rule all_CellAnno:
#    input:
#        expand(config["output05"]["annotation_plot"],tissue_all=tissue_all),
#        expand(config["output05"]["annotation_csv"],tissue_all=tissue_all),
#        expand(config["output05"]["seurat_output"],tissue_all=tissue_all)
#
#rule CellAnno:
#    input:
#        seurat_path = lambda wildcards: config["input05"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#    output:
#        annotation_plot = config["output05"]["annotation_plot"],
#        annotation_csv = config["output05"]["annotation_csv"],
#        seurat_output = config["output05"]["seurat_output"]
#    params:
#        species =  config["params05"]["species"]
#    script:
#        "pipeline/05_CellAnno.R"
#
## 亚群注释level1_check
#rule all_CellAnnoCheck:
#    input:
#        expand(config["output06"]["plot_output"],tissue_all=tissue_all),
#        expand(config["output06"]["seurat_output"],tissue_all=tissue_all)
#
#rule CellAnnoCheck:
#    input:
#        seurat_path = lambda wildcards: config["input06"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#    params:
#        annotation_path = lambda wildcards: config["params06"]["annotation_path"].format(tissue_all=wildcards.tissue_all)
#    output:
#        plot_output = config["output06"]["plot_output"],
#        seurat_output = config["output06"]["seurat_output"]
#    script:
#        "pipeline/06_CellAnnoCheck.R"
#
## 提取感兴趣的亚群细胞
#rule all_SubsetCells:
#    input:
#        expand(config["output07"]["output_path"],tissue_all=tissue_all)
#
#rule SubsetCells:
#    input:
#        seurat_path = lambda wildcards: config["input07"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#
#    params:
#        celltype = config["params07"]["celltype"]
#
#    output:
#        output_path = config["output07"]["output_path"]
#
#    script:
#        "pipeline/07_SubsetCells.R"
#
## 亚群细胞进行SCT流程
#rule all_Subset_SCT:
#    input:
#        expand(config["output08"]["output_path"],tissue_all=tissue_all)
#
#rule Subset_SCT:
#    input:
#        seurat_rdata = lambda wildcards: config["input08"]["seurat_rdata"].format(tissue_all=wildcards.tissue_all)
#    output:
#        output_path = config["output08"]["output_path"]
#    params:
#        nFeature_RNA_min = config["params08"]["nFeature_RNA_min"],
#        nFeature_RNA_max = config["params08"]["nFeature_RNA_max"],
#        percent_mt_max = config["params08"]["percent_mt_max"],
#        species = config["params08"]["species"],
#        dims = config["params08"]["dims"],
#        figure_paths = config["params08"]["figure_paths"]
#    script:
#        "pipeline/08_Subset_SCT.R"
#
## 亚群细胞注释
#rule all_Subset_CellAnno:
#    input:
#        expand(config["output09"]["annotation_plot"],tissue_all=tissue_all),
#        expand(config["output09"]["annotation_csv"],tissue_all=tissue_all),
#        expand(config["output09"]["seurat_output"],tissue_all=tissue_all)
#
#rule Subset_CellAnno:
#    input:
#        seurat_path = lambda wildcards: config["input09"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#    params:
#        species = config["params09"]["species"]
#    output:
#        annotation_plot = config["output09"]["annotation_plot"],
#        annotation_csv = config["output09"]["annotation_csv"],
#        seurat_output = config["output09"]["seurat_output"]
#    script:
#        "pipeline/05_CellAnno.R"
#
## 生成loom文件的前置文件
#rule all_SubcellAnnoCheck:
#    input:
#        expand(config["output10"]["seurat_rdata"],tissue_all=tissue_all),
#        expand(config["output10"]["output_path"],tissue_all=tissue_all)
#
#rule SubcellAnnoCheck:
#    input:
#        seurat_path= lambda wildcards: config["input10"]["seurat_path"].format(tissue_all=wildcards.tissue_all)
#
#    params:
#        cluster_path = lambda wildcards: config["params10"]["cluster_path"].format(tissue_all=wildcards.tissue_all),
#        res_param = config["params10"]["res_param"],
#        output_prefix = lambda wildcards: config["params10"]["output_prefix"].format(tissue_all=wildcards.tissue_all)
#
#    output:
#        seurat_rdata = config["output10"]["seurat_rdata"],
#        output_path = config["output10"]["output_path"]
#
#    script:
#        "pipeline/10_SubcellAnnoCheck.R"
#
## 生成loom文件
#rule all_trans2loom:
#    input:
#        expand(config['output11']['loom'],tissue_all=tissue_all)
#
#rule trans2loom:
#    input:
#        expression = config['input11']['expression']
#    output:
#        loom = config['output11']['loom']
#    shell:
#        "python pipeline/11_trans2loom.py {input.expression} {output.loom}" 
#
## pyscenic_grn
#rule all_pyscenic_grn:
#    input:
#        expand(config["output12"]["adj"], tissue_all=tissue_all) 
#
#rule pyscenic_grn:
#    input:
#        loomfile = lambda wildcards: config["input12"]["loomfile"].format(tissue_all=wildcards.tissue_all)# 输入loom文件
#
#    output:
#        adj = config["output12"]["adj"] # 输出adjacency表
#
#    params:
#        num_workers = config["params12"]["num_workers"], # 并行工作线程数
#        tf_list = config["params12"]["tf_list"]
#
#    shell:
#        """
#        pyscenic grn \
#            --num_workers {params.num_workers} \
#            --sparse \
#            --method grnboost2 \
#            --output {output.adj} \
#            {input.loomfile} \
#            {params.tf_list} 
#        """
#
## pyscenic_cistarget
#rule all_RcisTarget:
#    input:
#        expand(config["output13"]["regulons"],tissue_all=tissue_all)
#
#rule RcisTarget:
#    input:
#        loomfile = lambda wildcards: config["input13"]["loomfile"].format(tissue_all=wildcards.tissue_all),
#        grn_result = lambda wildcards: config["input13"]["grn_result"].format(tissue_all=wildcards.tissue_all)
#
#    params:
#        motifs = config["params13"]["motifs"],
#        feather = config["params13"]["feather"],
#        num_workers = config["params13"]["num_workers"]
#
#    output:
#        regulons = config["output13"]["regulons"],
#
#    shell:
#        """
#        pyscenic ctx \
#            --num_workers {params.num_workers} \
#            --output {output.regulons} \
#            --expression_mtx_fname {input.loomfile} \
#            --all_modules \
#            --mask_dropouts \
#            --mode "dask_multiprocessing" \
#            --min_genes 10 \
#            --annotations_fname {params.motifs} \
#            {input.grn_result} \
#            {params.feather}
#        """
#
## pyscenic_aucell
#rule all_aucell:
#    input:
#        expand(config['output14']['Scenic_loom'],tissue_all=tissue_all)
#
#rule aucell:
#    input:
#        loomfile = lambda wildcards: config['input14']["loomfile"].format(tissue_all=wildcards.tissue_all),
#        regulons = lambda wildcards: config['input14']["regulons"].format(tissue_all=wildcards.tissue_all)
#    
#    output:
#        Scenic_loom = config['output14']['Scenic_loom']
#
#    params:
#        num_workers = config['params14']['num_workers']
#
#    shell:
#        """
#        pyscenic aucell \
#            --num_workers {params.num_workers} \
#            --output {output.Scenic_loom} \
#            {input.loomfile} \
#            {input.regulons}
#        """
#
## pyscenic绘图前准备数据
#rule all_Plot_pyscenic:
#    input:
#        expand(config['output15']['seurat_res_path'],tissue_all=tissue_all)
#
#rule Plot_pyscenic:
#    input:
#        loom_path = lambda wildcards: config['input15']['loom_path'].format(tissue_all=wildcards.tissue_all),
#        seurat_path = lambda wildcards: config['input15']['seurat_path'].format(tissue_all=wildcards.tissue_all),
#    
#    output:
#        seurat_res_path = config['output15']['seurat_res_path']
#
#    params:
#        celltype_col = config['params15']['celltype_col'],
#        group_col = config['params15']['group_col'],
#        nF_col = config['params15']['nF_col'],
#        nC_col = config['params15']['nC_col'],
#        results_path = lambda wildcards: config['params15']['results_path'].format(tissue_all=wildcards.tissue_all)
#
#    script:
#        "pipeline/12_Plot_pyscenic.R"
#
## pyscenic绘图
#rule all_Plot_pyscenic_R:
#    input:
#        expand(config["output16"]["oregulon_activity"],tissue_all=tissue_all),
#        expand(config["output16"]["output_rank_RSS"],tissue_all=tissue_all)
#    
#rule Plot_pyscenic_R:
#    input:
#        iScenic_loom = lambda wildcards: config["input16"]["iScenic_loom"].format(tissue_all=wildcards.tissue_all),
#        irss_result = lambda wildcards: config["input16"]["irss_result"].format(tissue_all=wildcards.tissue_all),
#        iSeurat = lambda wildcards: config["input16"]["iSeurat"].format(tissue_all=wildcards.tissue_all),
#        icellinfo = lambda wildcards: config["input16"]["icellinfo"].format(tissue_all=wildcards.tissue_all),
#        iRegulon = config["input16"]["iRegulon"]
#
#    output:
#        oregulon_activity = config["output16"]["oregulon_activity"],
#        output_rank_RSS = config["output16"]["output_rank_RSS"]
#
#    params:
#        outputs_file_dir = lambda wildcards: config["params16"]["outputs_file_dir"].format(tissue_all=wildcards.tissue_all),
#
#    script:
#        "pipeline/13_Plot_TF_rss.R"

# 功能注释
rule all_Enrichment:
    input:
        expand(config['output17']["rds_output"],tissue_all=tissue_all)

rule Enrichment:
    input:
        iSeurat = lambda wildcards: config["input17"]["iSeurat"].format(tissue_all=wildcards.tissue_all)

    output:
        rds_output = config['output17']["rds_output"]
        
    params:
        iRegulon = config['params17']['iRegulon'],
        pidents = config['params17']['pidents'],
        outputs_file_dir = lambda wildcards: config['params17']['outputs_file_dir'].format(tissue_all=wildcards.tissue_all),
#        ocolor_palette = lambda wildcards: config['params17']['ocolor_palette'].format(tissue_all=wildcards.tissue_all),
        figure_dir = lambda wildcards: config['params17']['figure_dir'].format(tissue_all=wildcards.tissue_all)

    script:
        "pipeline/14_Enrichment.R"



