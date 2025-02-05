rule ldsr_setup:
    # Input can be bed file with gene boundaries or gene se                               
    input:   gene_set = "../results/02GENE_LISTS/shi_bc/LDSR/{CELL_TYPE}.{GENE_WINDOW}.bed",
             bim_file = "../resources/ldsr/reference_files/1000G_EUR_Phase3_plink/1000G.EUR.QC.{CHR}.bim"
    output:  "../results/05LDSR/annotation_files/snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.{CHR}.annot.gz"
    conda:   "../envs/ldsr.yml"
    message: "Creating annotation files for snRNAseq: {wildcards.CELL_TYPE}, {wildcards.GENE_WINDOW}, Chr {wildcards.CHR}"
    log:     "../results/00LOG/05LDSR/ldsr_make_annot.snRNAseq.{CELL_TYPE}.{GENE_WINDOW}.Chr{CHR}.log"
    shell:
             """

             git clone https://github.com/bulik/ldsc.git {params.resources}

             """
