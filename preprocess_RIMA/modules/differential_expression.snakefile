#!/usr/bin/env python

#-------------------------------DESeq2 for differential gene expression-----------------------#
import itertools
import pandas as pd


def diffexpr_targets(wildcards):
    ls=[]
    ls.append("files/ssgsea/ssgsea_all_terms_score.txt")
    #ls.append("files/multiqc/DESeq2/DESeq2_Volcano_Plot_mqc.png")
    for design in config["designs"]:
        meta = config["metasheet"]
        design_file = "./" + meta
        design_meta = pd.read_csv(design_file, index_col=0, sep=',')
        design_meta = design_meta.rename(columns = {design: 'Condition'})
        comps = [i for i in list(set(design_meta["Condition"].tolist()))]
        combinations = list(itertools.combinations(comps,2))
        if len(comps) >1:
            if config["comparison"] == "between":
                compare_list = combinations
            if config["comparison"] == "loop":
                compare_list = comps
            for cp in compare_list:
                if config["comparison"] == "between":
                    cp = sorted(list(cp))
                    cp_list = cp[0]+"_VS_"+cp[1]
                else:
                    cp_list = cp+"_VS_others"
                ls.append("images/deseq2/%s/%s_%s_diff_volcano_plot.png" % (design,design,cp_list))
                ls.append("files/deseq2/%s/%s_%s_DESeq2_ConvertID.txt" % (design,design,cp_list))
                ls.append("analysis/deseq2/%s/%s_%s_DESeq2_raw.txt" % (design,design,cp_list))
                ls.append("files/multiqc/DESeq2/%s_%s_DESeq2_sub.txt" % (design,cp_list))
                #ls.append("files/multiqc/differentialexpr/%s_%s/GO_BP_terms.txt" % (design,cp_list)
                #ls.append("files/multiqc/differentialexpr/%s_%s/GO_MF_terms.txt" % (design,cp_list)
                #ls.append("files/multiqc/differentialexpr/%s_%s/GO_CC_terms.txt" % (design,cp_list)
                #ls.append("files/multiqc/differentialexpr/%s_%s/KEGG_terms.txt" % (design,cp_list))
    return ls


rule differentialexpr_all:
    input:
       diffexpr_targets


rule deseq2_differential_genes:
    input:
        files = expand("analysis/salmon/{sample}/{sample}.quant.sf",sample=config["samples"])
    output:
        deseq2_res = "files/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt",
        deseq2_raw = "analysis/deseq2/{design}/{design}_{compare}_DESeq2_raw.txt",
        multiqc = "files/multiqc/DESeq2/{design}_{compare}_DESeq2_sub.txt"
    log:
        "logs/deseq2/{design}/{design}_{compare}.deseq2.log"
    params:
        comps = config["comparison"],
        filelist = lambda wildcards, input: ','.join(str(i) for i in list({input.files})[0]),
        batch = config["batch"],
        out_path = "analysis/deseq2/{design}/{design}",
        file_path = "files/deseq2/{design}/",
        tx_annot = "static/deseq2/tx2gene.csv",
        condition = config["designs"],
        meta = config["metasheet"],
        multiqc = "files/multiqc/DESeq2/{design}",
        source = "salmon",
        path="set +eu;source activate %s" % config['stat_root']
    message:
        "Running DESeq2 on the samples"
    benchmark:
        "benchmarks/deseq2/{design}/{design}_{compare}.deseq2.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/DESeq2.R --input {params.filelist} --batch {params.batch} --typ\
e {params.source} --comparison {params.comps} --meta {params.meta} --tx2gene {params.tx_annot} -\
-condition {params.condition} --outpath {params.out_path} --multiqc {params.multiqc}"
        " && mv {params.out_path}*DESeq2_ConvertID.txt {params.file_path}"

rule volcano_plot:
    input:
        "files/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt"
    output:
       plot =  "images/deseq2/{design}/{design}_{compare}_diff_volcano_plot.png",
       #multiqc_plot = "files/multiqc/DESeq2/{design}_{cDESeq2_Volcano_Plot_mqc.png"
    params:
        path="set +eu;source activate %s" % config['stat_root'],
        multiqc_folder = "files/multiqc/DESeq2/"
    log:
        "logs/deseq2/{design}/{design}_{compare}.deseq2.volcano.log"
    message:
        "Running Volcano plot on the DESeq2 result"
    benchmark:
        "benchmarks/deseq2/{design}/{design}_{compare}.deseq2.volcano.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/volcano_plot.R --deseq2_mat {input} --outdir {output.plot}"
        " && cp {output.plot} {params.multiqc_folder}"


rule gsea_plot:
    input:
        "files/deseq2/{design}/{design}_{compare}_DESeq2_ConvertID.txt"
    output:
        go_bp = "files/multiqc/differentialexpr/{design}_{compare}/GO_BP_terms.txt",
        go_mf = "files/multiqc/differentialexpr/{design}_{compare}/GO_MF_terms.txt",
        go_cc = "files/multiqc/differentialexpr/{design}_{compare}/GO_CC_terms.txt",
        kegg = "files/multiqc/differentialexpr/{design}_{compare}/KEGG_terms.txt"
    log:
        "logs/gsea/.gsea.log"
    params:
        files_path = "files/multiqc/differentialexpr/{design}_{compare}/",
        gsea_pcut = "0.1",
        gsea_minsize = "5",
        gsea_permutation = "1000",
        path="set +eu;source activate %s" % config['stat_root']
    message:
        "Running GSEA on the samples"
    benchmark:
        "benchmarks/gsea/gsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/gsea_plot.R --deseq2_mat {input} --pcut {params.gsea_pcut} --minsize {params.gsea_minsize} --npermutation {params.gsea_permutation} --outdir {params.files_path}"
        " && mv {params.out_path}*txt {params.files_path}"

rule ssgsea_score:
    input:
        "analysis/batchremoval/tpm_convertID.batch"
    output:
        score = "files/ssgsea/ssgsea_all_terms_score.txt",
    log:
        "logs/ssgsea/ssgsea.log"
    params:
        gmt = config['gmt_path'],
        comparison = lambda wildcards: ','.join(str(i) for i in list(config['ssgsea_comparisons'])),
        meta = config['metasheet'],
        top_n = config['ssgsea_top_n_terms'],
        outpath = "files/ssgsea/",
        path= "set +eu;source activate %s" % config['stat_root'],
        order= config['sam_ord']
    message:
        "Running single sample gene set enrichment analysis"
    benchmark:
        "benchmarks/ssgsea/ssgsea.benchmark"
    conda: "../envs/stat_perl_r.yml"
    shell:
        "{params.path}; Rscript src/ssgsea_plot.R -e {input} -g {params.gmt} -c {params.comparison} -m {params.meta} -n {params.top_n} -o {params.outpath} --order {params.order}"
