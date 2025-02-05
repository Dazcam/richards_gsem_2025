rule ldsr_make_annot:
    # Input can be bed file with gene boundaries or gene set with separate gene coord file
    input:   gene_set = "../results/01GENESETS/LDSR/{CELL_TYPE}.{GENE_WINDOW}.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_Phase3_plinkfiles/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/02LDSR/annotation_files/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz"
    resources: threads = 1, mem_mb = 10000
#    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for Stiletti_superclust: {wildcards.CELL_TYPE}, {wildcards.GENE_WINDOW}, Chr {wildcards.CHR}"
    log:     "../results/00LOG/02LDSR/ldsr_make_annot.Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}.log"
    shell:
             """
             eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
             conda activate ldsc
             python ../resources/ldsr/make_annot.py \
             --bed-file {input.gene_set} \
             --windowsize 0 \
             --bimfile {input.bim_file} \
             --annot-file {output} 2> {log}
             
             """
        
rule ldsr_ld_scores:
    input:   annot = "../results/02LDSR/annotation_files/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz",
             bfile_folder = "../resources/ldsr/reference_files/1000G_Phase3_plinkfiles",
             snps = "../resources/ldsr/reference_files/hm3_no_MHC.list.txt"
    output:  "../results/02LDSR/annotation_files/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.l2.ldscore.gz"
    resources: threads = 1, mem_mb = 10000       
#    conda:   "../envs/ldsr.yml"
    params:  bfile = "../resources/ldsr/reference_files/1000G_Phase3_plinkfiles/1000G.EUR.QC.{CHR}",
             ldscores = "../results/02LDSR/annotation_files/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{CHR}",
    message: "Running LDSR Phase 3 for {wildcards.CELL_TYPE}, {wildcards.GENE_WINDOW}, CHR {wildcards.CHR}" 
    log:     "../results/00LOG/02LDSR/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}_ldsc.log"
    shell:
             """
             eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
             conda activate ldsc
             python ../resources/ldsr/ldsc.py --thin-annot --l2 --bfile {params.bfile} --ld-wind-cm 1 \
             --annot {input.annot} --out {params.ldscores} --print-snps {input.snps} 2> {log}
             """

rule ldsr_stratified_baseline_v12:
    input:   GWAS = "../resources/gwas/ldsr_ready/{GWAS}_hg19_ldsr_ready_sumstats.gz",
             LDSR = expand("../results/02LDSR/annotation_files/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.l2.ldscore.gz", CELL_TYPE = config["CELL_TYPES"], GENE_WINDOW = config["GENE_WINDOW"], CHR = range(1,23))
    output:  "../results/02LDSR/SLDSR_baseline_v1.2/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2.results"
    resources: threads = 1, mem_mb = 20000   
#    conda:   "../envs/ldsr.yml"
    params:  weights = "../resources/ldsr/reference_files/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.",
             baseline = "../resources/ldsr/reference_files/1000G_Phase3_baseline_v1.2_ldscores/baseline.",
             frqfile = "../resources/ldsr/reference_files/1000G_Phase3_frq/1000G.EUR.QC.",
             LD_anns = "../results/02LDSR/annotation_files/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.",
             out_file = "../results/02LDSR/SLDSR_baseline_v1.2/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2"
    message: "Running SLDSR with {wildcards.CELL_TYPE} {wildcards.GENE_WINDOW} and {wildcards.GWAS} GWAS"
    log:     "../results/00LOG/02LDSR/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}.baseline.v1.2_partHerit.log"
    shell:
             """
             eval "$(/apps/languages/miniforge3/24.3.0-0/bin/conda shell.bash hook)"
             conda activate ldsc
             python ../resources/ldsr/ldsc.py --h2 {input.GWAS} --w-ld-chr {params.weights} \
             --ref-ld-chr {params.baseline},{params.LD_anns} --overlap-annot \
             --frqfile-chr {params.frqfile} --out {params.out_file} --print-coefficients 2> {log}
             """

rule ldsr_stratified_summary:
    input:   expand("../results/02LDSR/SLDSR_baseline_v1.2/Stiletti_superclust.{CELL_TYPE}.{GENE_WINDOW}.{GWAS}_baseline.v1.2.results", CELL_TYPE = config["CELL_TYPES"], GENE_WINDOW = config["GENE_WINDOW"], GWAS = config["GWAS"])
    output:  "../results/02LDSR/SLDSR_baseline_v1.2/SLDSR_{GWAS}_baseline.v1.2_summary.tsv"
    message: "Creating summary file for {wildcards.GWAS} GWAS"
    params:  dir = "../results/02LDSR/SLDSR_baseline_v1.2/",
             cell_types = "../resources/sheets/superclust_celltypes.tsv"
    log:     "../results/00LOG/02LDSR/Stiletti_superclust.{GWAS}_baseline.v1.2_partHerit.summary.log"
    shell:
             """

             
             head -1 {params.dir}Stiletti_superclust.Amygdala_excitatory.100UP100DOWN.SZ_baseline.v1.2.results > {output}
             File={params.cell_types}
             Lines=$(cat $File)
             for Line in $Lines
             do
             grep L2_1 ../results/02LDSR/SLDSR_baseline_v1.2/Stiletti_superclust."$Line".{wildcards.GWAS}_baseline.v1.2.results | sed "s/L2_1/$Line/g" >> {output} 2> {log}
             done

             """
