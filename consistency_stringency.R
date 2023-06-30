library(biomaRt)
library(ggplot2)
library(cowplot)
library(ggvenn)
library(rtracklayer)


##############
# Parameters #
##############
# dataset ID
dataset <- "GSE136913"
# dataset <- "GSE140420"
# dataset <- "GSE179101"
work_dir <- "./"
raw_dir <- sprintf("%s/00-Raw_data", work_dir)
input_dir <- sprintf("%s/20-Expression_analysis/input/%s", work_dir, dataset)
output_dir <- sprintf("%s/20-Expression_analysis/output/10-Consistency_stringency/%s", work_dir, dataset)
src_dir <- sprintf("%s/20-Expression_analysis/src", work_dir)

# reference genomes
ensembl_104 <- list(ref_genome="NCBIRefSeq106_Ensembl104GTF", biomaRt_attribute="ensembl_gene_id", AnnotationDb_keytype="ENSEMBL", gtf_file=sprintf("%s/reference_genomes/Rattus_norvegicus.Rnor_6.0.104.NCBIseqname.gtf.gz", raw_dir))
ncbi_refseq_106 <- list(ref_genome="NCBIRefSeq106_NCBIRefSeq106GTF", biomaRt_attribute="entrezgene_accession", AnnotationDb_keytype="SYMBOL", gtf_file=sprintf("%s/reference_genomes/GCF_000001895.5_Rnor_6.0_genomic.gtf.gz", raw_dir))
ncbi_refseq_108 <- list(ref_genome="NCBIRefSeq108_NCBIRefSeq108GTF", biomaRt_attribute="entrezgene_accession", AnnotationDb_keytype="SYMBOL", gtf_file=sprintf("%s/reference_genomes/GCF_015227675.2_mRatBN7.2_genomic.gtf.gz", raw_dir))
ensembl_105 <- list(ref_genome="NCBIRefSeq108_Ensembl105GTF", biomaRt_attribute="ensembl_gene_id", AnnotationDb_keytype="ENSEMBL", gtf_file=sprintf("%s/reference_genomes/Rattus_norvegicus.mRatBN7.2.105.NCBIseqname.gtf.gz", raw_dir))
genome_list <- list(ensembl_104, ncbi_refseq_106, ncbi_refseq_108, ensembl_105)
names(genome_list) <- c("Ensembl 104", "NCBI RefSeq 106", "NCBI RefSeq 108", "Ensembl 105")

# n value in consistency/stringency formula
consistency_stringency_n_value_to_add <- 1

# consistency/stringency threshold: 5% difference of expression, i.e. consistency/stringency threshold: -log10(log2(1.1))=1.152493
threshold <- -log(log(1.05, 2), 10)


#############
# Functions #
#############
source(sprintf("%s/ref_genome_eval_functions.R", src_dir))

############
# Analysis #
############

# GSE136913, GSE140420 and GSE179101 all samples comparisons
for(i in 1:(length(names(genome_list))-1)) {
  exp_matrix_1_name <- names(genome_list)[i]
  exp_matrix_1_ref_genome <- genome_list[[exp_matrix_1_name]]$ref_genome
  exp_matrix_1_attribute <- genome_list[[exp_matrix_1_name]]$biomaRt_attribute
  exp_matrix_1_gtf_file <- genome_list[[exp_matrix_1_name]]$gtf_file
  # load expression matrices
  load(sprintf("%s/%s_%s_exp_matrix.RData", input_dir, dataset, exp_matrix_1_ref_genome))
  exp_matrix_1 <- exp_matrix
  exp_matrix_1_metadata <- metadata_df
  for(j in (i+1):length(names(genome_list))) {
    exp_matrix_2_name <- names(genome_list)[j]
    exp_matrix_2_ref_genome <- genome_list[[exp_matrix_2_name]]$ref_genome
    exp_matrix_2_attribute <- genome_list[[exp_matrix_2_name]]$biomaRt_attribute
    exp_matrix_2_gtf_file <- genome_list[[exp_matrix_2_name]]$gtf_file
    # load expression matrices
    load(sprintf("%s/%s_%s_exp_matrix.RData", input_dir, dataset, exp_matrix_2_ref_genome))
    exp_matrix_2 <- exp_matrix
    exp_matrix_2_metadata <- metadata_df
    
    comparison <- sprintf("%s_%s", gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
    print(sprintf("comparison: %s", comparison))
    comparison_output_dir <- sprintf("%s/%s", output_dir, comparison)
    if (dataset == "GSE140420") {
      comparison_output_dir <- sprintf("%s/all_samples", comparison_output_dir)
    }
    comparison_output_dir <- sprintf("%s/n_%d", comparison_output_dir, consistency_stringency_n_value_to_add)
    # comparison_output_dir <- sprintf("%s/n_%d/abs_ratios", comparison_output_dir, consistency_stringency_n_value_to_add)
    if (! dir.exists(comparison_output_dir)) {
      dir.create(comparison_output_dir, recursive=TRUE, mode="0775")
    }
    
    # get common samples
    exp_matrix_1_samples <- rownames(exp_matrix_1_metadata)
    exp_matrix_2_samples <- rownames(exp_matrix_2_metadata)
    samples <- exp_matrix_1_samples[exp_matrix_1_samples %in% exp_matrix_2_samples]
    # expression matrices with common samples
    exp_matrix_1_expressed <- expression_matrix_common_samples_expressed_features(exp_matrix_1, samples)
    exp_matrix_2_expressed <- expression_matrix_common_samples_expressed_features(exp_matrix_2, samples)
    
    # mean expression density plot
    exp_matrix_1_mean_exp_log <- apply(log(exp_matrix_1_expressed+1, 2), 1, mean)
    exp_matrix_2_mean_exp_log <- apply(log(exp_matrix_2_expressed+1, 2), 1, mean)
    df2ggplot <- data.frame(exp=c(exp_matrix_1_mean_exp_log, exp_matrix_2_mean_exp_log), genome=c(rep(exp_matrix_1_name, length(exp_matrix_1_mean_exp_log)), rep(exp_matrix_2_name, length(exp_matrix_2_mean_exp_log))))
    mean_expression_output_basename <- sprintf("%s_%s_mean_expression", dataset, comparison)
    pdf(sprintf("%s/%s.pdf", comparison_output_dir, mean_expression_output_basename))
    percentiles_df <- density_percentiles_plot(df2ggplot, "exp", "genome", "Mean log2 expression")
    dev.off()
    write.csv(percentiles_df, file=sprintf("%s/%s_percentiles.csv", comparison_output_dir, mean_expression_output_basename), quote=FALSE, row.names=FALSE)
    
    if (exp_matrix_1_attribute != exp_matrix_2_attribute) {
      # Ensembl vs NCBI Refseq expression matrices comparison
      if (exp_matrix_1_attribute == "ensembl_gene_id") {
        # exp_matrix_1 is the Ensembl expression matrix
        consistency_stringency_pipeline_ensembl_ncbi_gene_ids(exp_matrix_1_expressed, exp_matrix_2_expressed, exp_matrix_1_gtf_file, exp_matrix_2_gtf_file, exp_matrix_1_name, exp_matrix_2_name, exp_matrix_1_attribute, exp_matrix_2_attribute, consistency_stringency_n_value_to_add, comparison, dataset, comparison_output_dir)
      } else {
        if (exp_matrix_2_attribute == "ensembl_gene_id") {
          # exp_matrix_2 is the Ensembl expression matrix
          consistency_stringency_pipeline_ensembl_ncbi_gene_ids(exp_matrix_2_expressed, exp_matrix_1_expressed, exp_matrix_2_gtf_file, exp_matrix_1_gtf_file, exp_matrix_2_name, exp_matrix_1_name, exp_matrix_2_attribute, exp_matrix_1_attribute, consistency_stringency_n_value_to_add, comparison, dataset, comparison_output_dir)
        }
      }
    } else {
      # comparison of expression matrices  with the same gene identifiers
      consistency_stringency_pipeline_same_gene_id_type(exp_matrix_1_expressed, exp_matrix_2_expressed, exp_matrix_1_gtf_file, exp_matrix_2_gtf_file, exp_matrix_1_name, exp_matrix_2_name, consistency_stringency_n_value_to_add, comparison, dataset, comparison_output_dir)
    }
  }
}

# GSE140420 per region comparisons
for(i in 1:(length(names(genome_list))-1)) {
  exp_matrix_1_name <- names(genome_list)[i]
  exp_matrix_1_ref_genome <- genome_list[[exp_matrix_1_name]]$ref_genome
  exp_matrix_1_attribute <- genome_list[[exp_matrix_1_name]]$biomaRt_attribute
  exp_matrix_1_gtf_file <- genome_list[[exp_matrix_1_name]]$gtf_file
  # load expression matrices
  load(sprintf("%s/%s_%s_exp_matrix.RData", input_dir, dataset, exp_matrix_1_ref_genome))
  exp_matrix_1 <- exp_matrix
  exp_matrix_1_metadata <- metadata_df
  for(j in (i+1):length(names(genome_list))) {
    exp_matrix_2_name <- names(genome_list)[j]
    exp_matrix_2_ref_genome <- genome_list[[exp_matrix_2_name]]$ref_genome
    exp_matrix_2_attribute <- genome_list[[exp_matrix_2_name]]$biomaRt_attribute
    exp_matrix_2_gtf_file <- genome_list[[exp_matrix_2_name]]$gtf_file
    # load expression matrices
    load(sprintf("%s/%s_%s_exp_matrix.RData", input_dir, dataset, exp_matrix_2_ref_genome))
    exp_matrix_2 <- exp_matrix
    exp_matrix_2_metadata <- metadata_df
    
    comparison <- sprintf("%s_%s", gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
    print(sprintf("comparison: %s", comparison))
    comparison_output_dir <- sprintf("%s/%s", output_dir, comparison)
    
    for (one_region in levels(exp_matrix_1_metadata$Region)) {
      print(one_region)
      region_dataset <- sprintf("%s_%s", dataset, one_region)
      region_output_dir <- sprintf("%s/%s/n_%d", comparison_output_dir, one_region, consistency_stringency_n_value_to_add)
      if (! dir.exists(region_output_dir)) {
        dir.create(region_output_dir, recursive=TRUE, mode="0775")
      }
      
      # get common samples
      exp_matrix_1_region_samples <- rownames(exp_matrix_1_metadata[which(exp_matrix_1_metadata$Region==one_region),])
      exp_matrix_2_region_samples <- rownames(exp_matrix_2_metadata[which(exp_matrix_2_metadata$Region==one_region),])
      samples <- exp_matrix_1_region_samples[exp_matrix_1_region_samples %in% exp_matrix_2_region_samples]
      # expression matrices with common samples
      exp_matrix_1_region_exp_matrix_expressed <- expression_matrix_common_samples_expressed_features(exp_matrix_1, samples)
      exp_matrix_2_region_exp_matrix_expressed <- expression_matrix_common_samples_expressed_features(exp_matrix_2, samples)
      
      # mean expression density plot
      exp_matrix_1_mean_exp_log <- apply(log(exp_matrix_1_region_exp_matrix_expressed+1, 2), 1, mean)
      exp_matrix_2_mean_exp_log <- apply(log(exp_matrix_2_region_exp_matrix_expressed+1, 2), 1, mean)
      df2ggplot <- data.frame(exp=c(exp_matrix_1_mean_exp_log, exp_matrix_2_mean_exp_log), genome=c(rep(exp_matrix_1_name, length(exp_matrix_1_mean_exp_log)), rep(exp_matrix_2_name, length(exp_matrix_2_mean_exp_log))))
      mean_expression_output_basename <- sprintf("%s_%s_mean_expression", region_dataset, comparison)
      pdf(sprintf("%s/%s.pdf", region_output_dir, mean_expression_output_basename))
      percentiles_df <- density_percentiles_plot(df2ggplot, "exp", "genome", "Mean log2 expression")
      dev.off()
      write.csv(percentiles_df, file=sprintf("%s/%s_percentiles.csv", region_output_dir, mean_expression_output_basename), quote=FALSE, row.names=FALSE)
      
      if (exp_matrix_1_attribute == "ensembl_gene_id" | exp_matrix_2_attribute == "ensembl_gene_id") {
        # Ensembl vs NCBI Refseq expression matrices comparison
        if (exp_matrix_1_attribute == "ensembl_gene_id") {
          # exp_matrix_1 is the Ensembl expression matrix
          consistency_stringency_pipeline_ensembl_ncbi_gene_ids(exp_matrix_1_region_exp_matrix_expressed, exp_matrix_2_region_exp_matrix_expressed, exp_matrix_1_gtf_file, exp_matrix_2_gtf_file, exp_matrix_1_name, exp_matrix_2_name, exp_matrix_1_attribute, exp_matrix_2_attribute, consistency_stringency_n_value_to_add, comparison, region_dataset, region_output_dir)
        } else {
          if (exp_matrix_2_attribute == "ensembl_gene_id") {
            # exp_matrix_2 is the Ensembl expression matrix
            consistency_stringency_pipeline_ensembl_ncbi_gene_ids(exp_matrix_2_region_exp_matrix_expressed, exp_matrix_1_region_exp_matrix_expressed, exp_matrix_2_gtf_file, exp_matrix_1_gtf_file, exp_matrix_2_name, exp_matrix_1_name, exp_matrix_2_attribute, exp_matrix_1_attribute, consistency_stringency_n_value_to_add, comparison, region_dataset, region_output_dir)
          }
        }
      } else {
        # comparison of expression matrices  with the same gene identifiers
        consistency_stringency_pipeline_same_gene_id_type(exp_matrix_1_region_exp_matrix_expressed, exp_matrix_2_region_exp_matrix_expressed, exp_matrix_1_gtf_file, exp_matrix_2_gtf_file, exp_matrix_1_name, exp_matrix_2_name, consistency_stringency_n_value_to_add, comparison, region_dataset, region_output_dir)
      }
    }
  }
}

