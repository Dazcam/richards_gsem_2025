#--------------------------------------------------------------------------------------
#
#    Richards 2025 - Generate Stiletti Supercluster Top 10% specificity scores
#
#--------------------------------------------------------------------------------------

## Info  ------------------------------------------------------------------------------
# Top 10% defined in TDEP col as 1
# Note that hg19_100kb_before_start end position didn't match

# Load pkgs & funcs
library(biomaRt)
library(tidyverse)
source(paste0(script_dir, 'functions.R'))

# Set variables
root_dir <- '~/Desktop/richards_gsem_2025/'
resources_dir <- paste0(root_dir, 'resources/')
script_dir <- paste0(root_dir, 'workflow/scripts/')
specificity_dir <- paste0(resources_dir, 'Source_Data_v241210/Specificity_per_cell_type/')
geneset_dir <- paste0(root_dir, 'results/01GENESETS/')


# Load data
download_supercluster_spec_scores(resources_dir)
superclust_spec_scores <- read_tsv(paste0(specificity_dir, 'Siletti_Supercluster_expression_specificity_TDEP_label.tsv.gz'))
supercluster_top10pc <- superclust_spec_scores |> filter(TDEP == 1)

# Get lookup: 88,802 genes in hg19
hg19_lookup <- get_biomart_gene_lookup(genome_build = 'hg19')
hg38_lookup <- get_biomart_gene_lookup(genome_build = 'hg38')

hg19_lookup_subset <- hg19_lookup |>
  dplyr::mutate(hg19_chr = paste0('chr', chromosome_name)) |>
  dplyr::select(hg19_chr, start_position, end_position, ensembl_gene_id, hgnc_symbol)

# Join
supercluster_results <- supercluster_top10pc |> 
  group_by(Supercluster) |> 
  nest() |> 
  mutate(
    matched_data = map(data, ~ .x |> 
                         inner_join(hg19_lookup_subset,  
                                    by = join_by(hg19_chr, ENSGID == ensembl_gene_id)) |> 
                         distinct()),  
    unmatched_data = map(data, ~ .x |> 
                           anti_join(hg19_lookup_subset,  
                                     by = join_by(hg19_chr, ENSGID == ensembl_gene_id)))
  ) |> 
  dplyr::select(Supercluster, matched_data, unmatched_data) |> 
  split(~Supercluster) 

# Check genes in top10pc not in lookup
results_list <- list()
for (cluster in names(supercluster_results)) {
  message('Checking for genes not in hg19 for: ', cluster)
    
    # Count overlapping and non-overlapping genes
    overlapping_genes <- nrow(supercluster_results[[cluster]]$matched_data[[1]])
    non_overlapping_genes <- nrow(supercluster_results[[cluster]]$unmatched_data[[1]])
    
    # Find non-overlapping genes that are in hg19
    unmatched_genes_in_hg19 <- supercluster_results[[cluster]]$unmatched_data[[1]] |>
      inner_join(hg19_lookup, join_by(ENSGID == ensembl_gene_id))
    
    # Find non-overlapping genes that are in hg38
    unmatched_genes_in_hg38 <- supercluster_results[[cluster]]$unmatched_data[[1]] |>
      inner_join(hg38_lookup, join_by(ENSGID == ensembl_gene_id))
    
    # Store in list
    results_list[[cluster]] <- tibble(
      Supercluster = cluster,
      Overlapping_genes_hg19 = overlapping_genes,
      Non_overlapping_genes_hg19 = non_overlapping_genes,
      Non_overlapping_genes_hg19_pc =  non_overlapping_genes / overlapping_genes * 100,
      Non_overlapping_genes_in_hg19_lookup = nrow(unmatched_genes_in_hg19),
      Non_overlapping_genes_in_hg38_lookup = nrow(unmatched_genes_in_hg38)
    )
  }
  
# Convert the list to a single dataframe
summary_df <- bind_rows(results_list)
print(summary_df, n = Inf)
  

# Check for patches, sex chromosomes
# Define valid chromosome list
valid_chromosomes <- paste0("chr", 1:22)

# Check each Supercluster for invalid chromosomes
invalid_chromosomes <- map(supercluster_results, ~ .x$matched_data[[1]] |> 
                             distinct(hg19_chr) |> 
                             filter(!hg19_chr %in% valid_chromosomes))

# Filter results to only Superclusters that have invalid chromosomes
invalid_chromosomes <- keep(invalid_chromosomes, ~ nrow(.x) > 0)

# Print results
if (length(invalid_chromosomes) > 0) {
  for (cluster in names(invalid_chromosomes)) {
    message("\nSupercluster: ", cluster)
    print(invalid_chromosomes[[cluster]])
  }
} else {
  message("No superclusters contain invalid chromosomes.")
}

matched_list <- map(supercluster_results, ~ .x$matched_data[[1]])


# Remove MHC genes
# Test 2 methods, by gene name, by hard filter: Identical output!
mhc_coords <- c("chr6", "25000000", "35000000") # hg19 coords - extended region (PGC3 SCZ ref)
mart_hg19 <- useMart('ENSEMBL_MART_ENSEMBL', host = 'https://grch37.ensembl.org')
mart_hg19 <- useDataset('hsapiens_gene_ensembl', mart_hg19)
mhc_genes <- getBM(attributes = c('ensembl_gene_id', "external_gene_name", "chromosome_name", "start_position", "end_position"),
                   filters = c("chromosome_name","start","end"),
                   values = list(chromosome = mhc_coords[1], start = mhc_coords[2], end = mhc_coords[3]), 
                   mart = mart_hg19)


matched_list_mhc_genes <- map(matched_list, ~ .x |> 
                                anti_join(mhc_genes, join_by(ENSGID == ensembl_gene_id))) # Remove genes using the MHC gene list

matched_list_mhc_coords <- map(matched_list, ~ .x |> 
                                 filter(!(hg19_chr == mhc_coords[1] & 
                                            ((start_position > as.numeric(mhc_coords[2]) & start_position < as.numeric(mhc_coords[3])) | 
                                               (end_position > as.numeric(mhc_coords[2]) & end_position < as.numeric(mhc_coords[3]))))))


comparison_results <- map_df(names(matched_list), function(cluster) {
  original_count <- nrow(matched_list[[cluster]])
  gene_list_count <- nrow(matched_list_mhc_genes[[cluster]])
  coord_filter_count <- nrow(matched_list_mhc_coords[[cluster]])
  
  tibble(
    Supercluster = cluster,
    Original = original_count,
    After_Gene_List = gene_list_count,
    After_Coord_Filter = coord_filter_count,
    Removed_Gene_List = original_count - gene_list_count,
    Removed_Coord_Filter = original_count - coord_filter_count
  )
})
print(comparison_results, n = Inf)

# Write gene set files for MAGMA and LDSR ---------------------------------------------
# MAGMA
for(i in names(matched_list_mhc_genes)) {
  
  ensg_ids <- matched_list_mhc_genes[[i]] |> 
    dplyr::pull(ENSGID) 

  cat(i, " ", paste(ensg_ids, collapse = " "), "\n", 
      file = paste0(geneset_dir, 'MAGMA/stiletti_superclust_top10pc.txt'), sep = '', append = TRUE)
  
}

# LDSR
for (cell_type in names(matched_list_mhc_genes)) {
  
  window <- '100UP100DOWN'
  matched_list_mhc_genes[[cell_type]] |> 
    mutate(start = pmax(0, start_position - 100000), # Avoid negative start
           end = end_position + 100000) |>  
    dplyr::select(chr = hg19_chr, start, end, ensembl_id = ENSGID) |>  
    write_tsv(paste0(geneset_dir, 'LDSR/', cell_type, '.', window, '.bed'), 
              col_names = FALSE)
}

# Check extensions
# matched_list_mhc_genes[[cell_type]] |> 
#   mutate(start = pmax(0, start_position - 100000), # Avoid negative start
#          end = end_position + 100000) |>
#   dplyr::select(start, end, start_position, end_position) |>
#   mutate(start_diff = start - start_position,
#          end_diff = end - end_position)

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
