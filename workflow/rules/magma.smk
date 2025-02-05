rule magma_gene_set_analysis:
    input:   genes = "../resources/gwas/magma_ready/{GWAS}_hg19_magma_ready_35UP10DOWN.genes.raw",
             data  = "../results/01GENESETS/MAGMA/stiletti_superclust_top10pc.txt"
    output:  "../results/03MAGMA/{GWAS}_hg19_magma_35UP10DOWN.gsa.out"
    resources: threads = 1, mem_mb = 10000  
    params:  out = "../results/03MAGMA/{GWAS}_hg19_magma_35UP10DOWN"
    message: "Running MAGMA gene set analysis step for {wildcards.GWAS}, Gene window: 35UP10DOWN"
    log:     "../results/00LOG/03MAGMA/snRNAseq.GE.gene_set_analysis.{GWAS}.35UP10DOWN.log"
    shell:
             """

             module load magma/1.10
             magma --gene-results {input.genes} --set-annot {input.data} --out {params.out} &> {log}

             """

