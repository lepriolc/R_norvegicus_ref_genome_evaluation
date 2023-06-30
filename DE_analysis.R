library(edgeR)
require(org.Rn.eg.db)
library(biomaRt)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(rtracklayer)


##############
# Parameters #
##############
# dataset <- "GSE136913"
# dataset <- "GSE140420"
dataset <- "GSE179101"
work_dir <- "./"
raw_dir <- sprintf("%s/00-Raw_data", work_dir)
input_dir <- sprintf("%s/20-Expression_analysis/input/%s", work_dir, dataset)
output_dir <- sprintf("%s/20-Expression_analysis/output/20-DE_analysis/%s", work_dir, dataset)
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

# designs
## GSE136913
# design_formula <- "~Instrument+Litter+Condition"
# condition_name <- "Condition"
# condition_ref <- NA
# second_condition_name <- NA
# design_output_name_append <- NULL

## GSE140420
### Region, blocking: Age
# design_formula <- "~Age+Region"
# condition_name <- "Region"
# condition_ref <- NA
# second_condition_name <- NA
# design_output_name_append <- NULL

## GSE179101
design_formula <- "~Subfield"
condition_name <- "Subfield"
condition_ref <- NA
second_condition_name <- NA
design_output_name_append <- NULL

# design output name
design_output_name <- sub("~|~0\\+", "", design_formula)
design_output_name <- gsub("\\+", "_", design_output_name)
design_output_name <- gsub("\\.", "", design_output_name)

# intercept boolean
intercept <- ifelse(grepl("^~0+", design_formula), FALSE, TRUE)

# fold-change thresholds
fc_thresholds <- c(1, 1.1, 1.25, 1.5, 2)

# Gene Ontology evidence categories: http://geneontology.org/docs/guide-go-evidence-codes/
evidence_categories <- data.frame(evidence=c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "HTP", "HDA", "HMP", "HGI", "HEP", "IBA", "IBD", "IKR", "IRD", "ISS", "ISO", "ISA", "ISM", "IGC", "RCA", "TAS", "NAS", "IC", "ND", "IEA"), category=c(rep("experimental", 11), rep("phylogeny", 4), rep("computational", 6), rep("author statement", 2), rep("curator statement", 2), "electronic"))

# plot theme
plot_theme <- theme_bw() +
  theme(panel.border=element_rect(color="grey50")) +
  theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=8)) +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=7), strip.text=element_text(size=8))


#############
# Functions #
#############
source(sprintf("%s/ref_genome_eval_functions.R", src_dir))


############
# Analysis #
############

#####################
# DE and GO analysis
for (one_genome in genome_list) {
  one_ref_genome <- one_genome$ref_genome
  print(one_ref_genome)
  one_ref_genome_in_path <- gsub(" ", "_", one_ref_genome)
  one_ref_genome_output_dir <- sprintf("%s/%s", output_dir, one_ref_genome_in_path)
  if (! dir.exists(one_ref_genome_output_dir)) {
    dir.create(one_ref_genome_output_dir, recursive=TRUE, mode="0775")
  }
  output_basename <- sprintf("%s_%s", dataset, one_ref_genome_in_path)
  
  # load expression matrix and metadata table
  load(sprintf("%s/%s_%s_exp_matrix.RData", input_dir, dataset, one_ref_genome))
  # samples <- rownames(metadata_df)
  ## remove '__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned' and '__too_low_aQual'
  exp_matrix <- exp_matrix[rownames(exp_matrix)[! grepl("^_", rownames(exp_matrix))],]
  ## remove non expressed features
  expressed_features <- rownames(exp_matrix)[unname(apply(exp_matrix, 1, mean) > 0)]
  exp_matrix_expressed <- exp_matrix[expressed_features,]
  ## write metadata table
  write.csv(metadata_df, file=sprintf("%s/%s_metadata.csv", one_ref_genome_output_dir, output_basename), quote=FALSE)
  
  # Fit a GLM and test for DE genes
  ## condition
  design_formula_variables <- sub("~|~0\\+", "", design_formula)
  design_variables <- unlist(strsplit(design_formula_variables, "\\+"))
  design_last_variable <- design_variables[length(design_variables)]
  if (grepl("\\.", design_last_variable)) {
    ### all combinations of multiple factors
    multiple_factor_combination <- TRUE
    condition_factor <- factor(paste(metadata_df[[unlist(strsplit(design_last_variable, "\\."))[1]]], make.names(metadata_df[[unlist(strsplit(design_last_variable, "\\."))[2]]]), sep="."))
  } else {
    multiple_factor_combination <- FALSE
    condition_factor <- metadata_df[[design_last_variable]]
    if (! is.na(condition_ref)) {
      condition_factor <- relevel(condition_factor, ref=condition_ref)
    }
  }
  
  ## design matrix
  print(sprintf("design: %s", design_formula))
  if (multiple_factor_combination) {
    design <- model.matrix(~0+condition_factor, data=metadata_df)
    ### rename design matrix columns: no intercept
    colnames(design) <- levels(condition_factor)
  } else {
    design <- model.matrix(as.formula(design_formula), data=metadata_df)
    ### rename design matrix columns: keep Intercept and remove reference level, i.e. first level, for both Region and Age
    design_colnames <- c()
    if (intercept) {
      design_colnames <- c(design_colnames, "Intercept")
    }
    for (one_variable in design_variables[1:length(design_variables)-1]) {
      design_colnames <- c(design_colnames, make.names(levels(metadata_df[[one_variable]])[-1]))
    }
    design_colnames <- c(design_colnames, make.names(levels(condition_factor)[-1]))
    colnames(design) <- design_colnames
  }
  
  ## contrasts
  if (is.na(second_condition_name)) {
    ### no second condition factor
    if (length(levels(condition_factor)) == 2) {
      #### only 2 levels for condition factor
      contrasts <- build_1vs1_1vsall_contrasts(design, condition_factor)
    } else {
      #### more than 2 levels for condition factor
      constrast_list <- all_combinations_mutliple_factors_contrasts(condition_factor, intercept, design_output_name_append, FALSE)
      contrasts <- makeContrasts(contrasts=constrast_list$contrasts, levels=design)
      colnames(contrasts) <- constrast_list$names
    }
    contrasts_output_name <- design_output_name
  } else {
    ### second condition factor
    condition <- metadata_df[[condition_name]]
    if (! is.na(condition_ref)) {
      condition <- relevel(condition, ref=condition_ref)
    }
    second_condition <- metadata_df[[second_condition_name]]
    first <- ifelse(second_condition_name == unlist(strsplit(design_last_variable, "\\."))[1], TRUE, FALSE)
    all_contrasts_vector <- c()
    all_contrasts_names_vector <- c()
    for (one_level in levels(second_condition)) {
      constrast_list <- all_combinations_mutliple_factors_contrasts(condition, FALSE, one_level, first)
      all_contrasts_vector <- c(all_contrasts_vector, constrast_list$contrasts)
      all_contrasts_names_vector <- c(all_contrasts_names_vector, constrast_list$names)
    }
    contrasts <- makeContrasts(contrasts=all_contrasts_vector, levels=design)
    colnames(contrasts) <- all_contrasts_names_vector
    contrasts_output_name <- sprintf("%s%s", design_output_name, design_output_name_append)
  }
  
  design_output_dir <- sprintf("%s/%s", one_ref_genome_output_dir, design_output_name)
  if (! dir.exists(design_output_dir)) {
    dir.create(design_output_dir, recursive=TRUE, mode="0775")
  }
  pipeline_basename <- sprintf("%s_%s", output_basename, contrasts_output_name)
  DE_analysis_pipeline(exp_matrix_expressed, condition_factor, design, contrasts, fc_thresholds, one_genome$AnnotationDb_keytype, one_genome$biomaRt_attribute, design_output_dir, pipeline_basename)
}
#####################


###############################################################
# DE analysis result comparison according to reference genomes
library(biomaRt)
library(dplyr)

if (is.null(design_output_name_append)) {
  design_output_name_append <- ""
}

# contrast names
## GSE136913
# contrast_output_names <- c("kainatevscontrol")
## GSE140420
### Region, blocking: Age
# contrast_output_names <- c("CA3vsCA1", "DGvsCA1", "CA1vsCA3DG", "DGvsCA3", "CA3vsCA1DG", "DGvsCA1CA3")
## GSE179101
contrast_output_names <- c("CA2vsCA1", "CA3vsCA1", "DGvsCA1", "CA1vsCA2CA3DG", "CA3vsCA2", "DGvsCA2", "CA2vsCA1CA3DG", "DGvsCA3", "CA3vsCA1CA2DG", "DGvsCA1CA2CA3")

# ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)
ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", ensemblRedirect = FALSE)

for(i in 1:(length(names(genome_list))-1)) {
  exp_matrix_1_name <- names(genome_list)[i]
  exp_matrix_1_ref_genome <- genome_list[[exp_matrix_1_name]]$ref_genome
  exp_matrix_1_attribute <- genome_list[[exp_matrix_1_name]]$biomaRt_attribute
  for(j in (i+1):length(names(genome_list))) {
    exp_matrix_2_name <- names(genome_list)[j]
    exp_matrix_2_ref_genome <- genome_list[[exp_matrix_2_name]]$ref_genome
    exp_matrix_2_attribute <- genome_list[[exp_matrix_2_name]]$biomaRt_attribute
    comparison <- sprintf("%s_%s", gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
    print(sprintf("comparison: %s", comparison))
    for (one_contrast_output_name in contrast_output_names) {
      print(sprintf("contrast: %s", one_contrast_output_name))
      for (one_fc_threshold in fc_thresholds) {
        print(sprintf("FC: %s", one_fc_threshold))
        comparison_FC_output_dir <- sprintf("%s/comparisons/%s/%s/%s/FC_%s/DE_analysis", output_dir, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
        if (! dir.exists(comparison_FC_output_dir)) {
          dir.create(comparison_FC_output_dir, recursive=TRUE, mode="0775")
        }
        comparison_FC_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, comparison, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
        
        # get DE results
        ## annotation 1
        exp_matrix_1_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s", output_dir, exp_matrix_1_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
        exp_matrix_1_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, exp_matrix_1_ref_genome, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
        exp_matrix_1_de_results_file <- sprintf("%s/%s_pval.csv", exp_matrix_1_ref_genome_output_dir, exp_matrix_1_output_basename)
        exp_matrix_1_de_results_df <- read.csv(exp_matrix_1_de_results_file, quote="", row.names=1)
        exp_matrix_1_DE_genes <- rownames(exp_matrix_1_de_results_df[which(exp_matrix_1_de_results_df$FDR < 0.05),])
        ## annotation 2
        exp_matrix_2_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s", output_dir, exp_matrix_2_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
        exp_matrix_2_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, exp_matrix_2_ref_genome, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
        exp_matrix_2_de_results_file <- sprintf("%s/%s_pval.csv", exp_matrix_2_ref_genome_output_dir, exp_matrix_2_output_basename)
        exp_matrix_2_de_results_df <- read.csv(exp_matrix_2_de_results_file, quote="", row.names=1)
        exp_matrix_2_DE_genes <- rownames(exp_matrix_2_de_results_df[which(exp_matrix_2_de_results_df$FDR < 0.05),])
        
        if (length(exp_matrix_1_DE_genes) != 0 & length(exp_matrix_2_DE_genes) != 0) {
          # gene ID matching
          gene_ID_matching_output_dir <- sprintf("%s/gene_ID_matching", comparison_FC_output_dir)
          if (! dir.exists(gene_ID_matching_output_dir)) {
            dir.create(gene_ID_matching_output_dir, recursive=TRUE, mode="0775")
          }
          ## different types of gene IDs
          if (exp_matrix_1_attribute != exp_matrix_2_attribute) {
            if (exp_matrix_1_attribute == "ensembl_gene_id") {
              ensembl_name <- exp_matrix_1_name
              ensembl_attribute <- exp_matrix_1_attribute
              ensembl_DE_genes <- exp_matrix_1_DE_genes
              ensembl_output_basename <- exp_matrix_1_output_basename
              ncbi_name <- exp_matrix_2_name
              ncbi_attribute <- exp_matrix_2_attribute
              ncbi_DE_genes <- exp_matrix_2_DE_genes
              ncbi_output_basename <- exp_matrix_2_output_basename
            } else {
              if (exp_matrix_2_attribute == "ensembl_gene_id") {
                ensembl_name <- exp_matrix_2_name
                ensembl_attribute <- exp_matrix_2_attribute
                ensembl_DE_genes <- exp_matrix_2_DE_genes
                ensembl_output_basename <- exp_matrix_2_output_basename
                ncbi_name <- exp_matrix_1_name
                ncbi_attribute <- exp_matrix_1_attribute
                ncbi_DE_genes <- exp_matrix_1_DE_genes
                ncbi_output_basename <- exp_matrix_1_output_basename
              }
            }
            
            print("Ensembl IDs and NCBI gene names")
            ### get Entrez accessions corresponding to Ensembl gene IDs
            gene.data.ensembl <- getBM(attributes=c(ensembl_attribute, ncbi_attribute), filters = ensembl_attribute, values = ensembl_DE_genes, mart = ensembl)
            output_basename <- sprintf("%s_useMart_getBM", ensembl_output_basename)
            #### getBM() function may include not requested IDs: remove them if any
            gene.data.ensembl.keep <- remove_not_requested_getBM_ids(gene.data.ensembl, ensembl_attribute, ensembl_DE_genes, gene_ID_matching_output_dir, output_basename)
            write.csv(gene.data.ensembl, file=sprintf("%s/%s.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
            ### get Ensembl gene IDs corresponding to Entrez accessions
            gene.data.ncbi <- getBM(attributes=c(ensembl_attribute, ncbi_attribute), filters = ncbi_attribute, values = ncbi_DE_genes, mart = ensembl)
            output_basename <- sprintf("%s_useMart_getBM", ncbi_output_basename)
            #### getBM() function may include not requested IDs: remove them if any
            gene.data.ncbi.keep <- remove_not_requested_getBM_ids(gene.data.ncbi, ncbi_attribute, ncbi_DE_genes, gene_ID_matching_output_dir, output_basename)
            write.csv(gene.data.ncbi, file=sprintf("%s/%s.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
            
            ## gene ID matching
            ### Ensembl gene IDs to Entrez gene accessions matching
            print("Ensembl to Entrez gene ID matching")
            output_basename <- sprintf("%s_useMart_getBM", ensembl_output_basename)
            ensembl2entrez_df <- Ensembl_NCBI_ID_matching(gene.data.ensembl.keep, ensembl_attribute, ncbi_attribute, ensembl_DE_genes, ncbi_DE_genes, gene_ID_matching_output_dir, output_basename)
            ### Entrez gene accessions to Ensembl gene IDs matching
            print("Entrez to Ensembl gene ID matching")
            output_basename <- sprintf("%s_useMart_getBM", exp_matrix_2_output_basename)
            entrez2ensembl_df <- Ensembl_NCBI_ID_matching(gene.data.ncbi.keep, ncbi_attribute, ensembl_attribute, ncbi_DE_genes, ensembl_DE_genes, gene_ID_matching_output_dir, output_basename)
            
            if (exp_matrix_1_attribute == "ensembl_gene_id") {
              exp_matrix_1_to_2_match_df <- ensembl2entrez_df
              exp_matrix_2_to_1_match_df <- entrez2ensembl_df
            } else {
              if (exp_matrix_2_attribute == "ensembl_gene_id") {
                exp_matrix_1_to_2_match_df <- entrez2ensembl_df
                exp_matrix_2_to_1_match_df <- ensembl2entrez_df
              }
            }
          } else {
            ## same type of gene IDs
            print(sprintf("%s to %s gene ID matching", exp_matrix_1_name, exp_matrix_2_name))
            # exp_matrix_1_to_2_match_df <- same_gene_id_type_matching(exp_matrix_1_DE_genes, exp_matrix_2_DE_genes, rownames(exp_matrix_2_de_results_df))
            exp_matrix_1_to_2_match_df <- same_gene_id_type_matching(exp_matrix_1_DE_genes, exp_matrix_2_DE_genes)
            print(sprintf("%s to %s gene ID matching", exp_matrix_2_name, exp_matrix_1_name))
            # exp_matrix_2_to_1_match_df <- same_gene_id_type_matching(exp_matrix_2_DE_genes, exp_matrix_1_DE_genes, rownames(exp_matrix_1_de_results_df))
            exp_matrix_2_to_1_match_df <- same_gene_id_type_matching(exp_matrix_2_DE_genes, exp_matrix_1_DE_genes)
          }
          
          ## add annotation column
          annotation <- rep(exp_matrix_1_name, dim(exp_matrix_1_to_2_match_df)[1])
          exp_matrix_1_to_2_match_df <- cbind(exp_matrix_1_to_2_match_df, annotation)
          write.csv(exp_matrix_1_to_2_match_df, file=sprintf("%s/%s_useMart_getBM_matching_categories.csv", gene_ID_matching_output_dir, exp_matrix_1_output_basename), quote=FALSE, row.names=FALSE)
          annotation <- rep(exp_matrix_2_name, dim(exp_matrix_2_to_1_match_df)[1])
          exp_matrix_2_to_1_match_df <- cbind(exp_matrix_2_to_1_match_df, annotation)
          write.csv(exp_matrix_2_to_1_match_df, file=sprintf("%s/%s_useMart_getBM_matching_categories.csv", gene_ID_matching_output_dir, exp_matrix_2_output_basename), quote=FALSE, row.names=FALSE)
          ## stats data frame
          nb_all_exp_matrix_1_ids <- length(exp_matrix_1_DE_genes)
          one_df2barplot <- gene_counts_per_ID_matching_category(exp_matrix_1_to_2_match_df, nb_all_exp_matrix_1_ids)
          annotation_IDs_df2barplot <- one_df2barplot
          nb_all_exp_matrix_2_ids <- length(exp_matrix_2_DE_genes)
          one_df2barplot <- gene_counts_per_ID_matching_category(exp_matrix_2_to_1_match_df, nb_all_exp_matrix_2_ids)
          annotation_IDs_df2barplot <- rbind(annotation_IDs_df2barplot, one_df2barplot)
          output_name <- sprintf("%s_gene_ID_matching", comparison_FC_output_basename)
          write.csv(annotation_IDs_df2barplot, file=sprintf("%s/%s_barplot_data.csv", comparison_FC_output_dir, output_name), quote=FALSE, row.names=FALSE)
          
          ## sanity check
          print("sanity check: verify that there are the same number of ID matchings for both matchings")
          ### different types of gene IDs
          if (exp_matrix_1_attribute != exp_matrix_2_attribute) {
            ensembl_matching_number <- sum(annotation_IDs_df2barplot[which(annotation_IDs_df2barplot$annotation==ensembl_name & annotation_IDs_df2barplot$category_2 %in% c("1:1_1", "1:N_1", "N:1_1", "N:P_1")), "value"]) + sum(gene.data.ensembl.keep[which(gene.data.ensembl.keep[[ensembl_attribute]] %in% as.character(ensembl2entrez_df[which(ensembl2entrez_df$category_2=="1:N_P"),"feature"])), ncbi_attribute] %in% ncbi_DE_genes) + sum(gene.data.ensembl.keep[which(gene.data.ensembl.keep[[ensembl_attribute]] %in% as.character(ensembl2entrez_df[which(ensembl2entrez_df$category_2=="N:P_Q"),"feature"])), ncbi_attribute] %in% ncbi_DE_genes)
            ncbi_matching_number <- sum(annotation_IDs_df2barplot[which(annotation_IDs_df2barplot$annotation==ncbi_name & annotation_IDs_df2barplot$category_2 %in% c("1:1_1", "1:N_1", "N:1_1", "N:P_1")), "value"]) + sum(gene.data.ncbi.keep[which(gene.data.ncbi.keep[[ncbi_attribute]] %in% as.character(entrez2ensembl_df[which(entrez2ensembl_df$category_2=="1:N_P"),"feature"])), ensembl_attribute] %in% ensembl_DE_genes) + sum(gene.data.ncbi.keep[which(gene.data.ncbi.keep[[ncbi_attribute]] %in% as.character(entrez2ensembl_df[which(entrez2ensembl_df$category_2=="N:P_Q"),"feature"])), ensembl_attribute] %in% ensembl_DE_genes)
            print(sprintf("number of %s to %s ID matchings: %d", ensembl_name, ncbi_name, ensembl_matching_number))
            print(sprintf("number of %s to %s ID matchings: %d", ncbi_name, ensembl_name, ncbi_matching_number))
          } else {
            ### same type of gene IDs
            exp_matrix_1_to_2_matching_number <- sum(annotation_IDs_df2barplot[which(annotation_IDs_df2barplot$annotation==exp_matrix_1_name & annotation_IDs_df2barplot$category_2 == "1:1_1"), "value"])
            exp_matrix_2_to_1_matching_number <- sum(annotation_IDs_df2barplot[which(annotation_IDs_df2barplot$annotation==exp_matrix_2_name & annotation_IDs_df2barplot$category_2 == "1:1_1"), "value"])
            print(sprintf("number of %s to %s ID matchings: %d", exp_matrix_1_name, exp_matrix_2_name, exp_matrix_1_to_2_matching_number))
            print(sprintf("number of %s to %s ID matchings: %d", exp_matrix_2_name, exp_matrix_1_name, exp_matrix_2_to_1_matching_number))
          }
          
          ## plots
          print("Gene ID matching plots")
          gene_ID_matching_barplots(annotation_IDs_df2barplot, comparison_FC_output_dir, output_name)
          
          ## matches with expressed but non-DE matching genes: only deal with 1:1_0, 1:N_0, 1:N_1, 1:N_P, N:1_0, N:P_0, N:P_1 categories
          print("Matches with non DE genes")
          matching_output_name <- sprintf("%s_gene_ID_matching_with_non_DE", comparison_FC_output_basename)
          categories_look_for_nonDE <- c("1:1_0", "1:N_0", "1:N_1", "1:N_P", "N:1_0", "N:P_0", "N:P_1", "N:P_Q")
          ### exp matrix 1 to exp matrix 2 ID matching
          print(sprintf("%s to %s gene ID matching", exp_matrix_1_name, exp_matrix_2_name))
          exp_matrix_2_non_DE_genes <- rownames(exp_matrix_2_de_results_df[which(exp_matrix_2_de_results_df$FDR > 0.05),])
          output_basename <- sprintf("%s_useMart_getBM", exp_matrix_1_output_basename)
          exp_matrix_1_to_2_match_with_nonDE_df <- gene_ID_matching_add_second_match(exp_matrix_1_to_2_match_df, categories_look_for_nonDE, exp_matrix_2_DE_genes, exp_matrix_2_non_DE_genes, gene_ID_matching_output_dir, output_basename)
          write.csv(exp_matrix_1_to_2_match_with_nonDE_df, file=sprintf("%s/%s_matching_categories_with_non_DE.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
          
          ### exp matrix 2 to exp matrix 1 ID matching
          print(sprintf("%s to %s gene ID matching", exp_matrix_2_name, exp_matrix_1_name))
          exp_matrix_1_non_DE_genes <- rownames(exp_matrix_1_de_results_df[which(exp_matrix_1_de_results_df$FDR > 0.05),])
          output_basename <- sprintf("%s_useMart_getBM", exp_matrix_2_output_basename)
          exp_matrix_2_to_1_match_with_nonDE_df <- gene_ID_matching_add_second_match(exp_matrix_2_to_1_match_df, categories_look_for_nonDE, exp_matrix_1_DE_genes, exp_matrix_1_non_DE_genes, gene_ID_matching_output_dir, output_basename)
          write.csv(exp_matrix_2_to_1_match_with_nonDE_df, file=sprintf("%s/%s_matching_categories_with_non_DE.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
          
          ### plots
          print("Gene ID matching plots")
          nb_all_exp_matrix_1_ids <- length(exp_matrix_1_DE_genes)
          one_df2barplot <- gene_counts_per_ID_matching_category(exp_matrix_1_to_2_match_with_nonDE_df, nb_all_exp_matrix_1_ids)
          annotation_IDs_with_nonDE_df2barplot <- one_df2barplot
          nb_all_exp_matrix_2_ids <- length(exp_matrix_2_DE_genes)
          one_df2barplot <- gene_counts_per_ID_matching_category(exp_matrix_2_to_1_match_with_nonDE_df, nb_all_exp_matrix_2_ids)
          annotation_IDs_with_nonDE_df2barplot <- rbind(annotation_IDs_with_nonDE_df2barplot, one_df2barplot)
          write.csv(annotation_IDs_with_nonDE_df2barplot, file=sprintf("%s/%s_barplot_data.csv", comparison_FC_output_dir, matching_output_name), quote=FALSE, row.names=FALSE)
          gene_ID_matching_barplots(annotation_IDs_with_nonDE_df2barplot, comparison_FC_output_dir, matching_output_name)
          
          ## sanity check
          print("sanity check: verify that there are the same number of ID matchings for both matchings with matches with non DE genes")
          ### different types of gene IDs
          if (exp_matrix_1_attribute != exp_matrix_2_attribute) {
            ensembl_matching_number <- sum(annotation_IDs_with_nonDE_df2barplot[which(annotation_IDs_with_nonDE_df2barplot$annotation==ensembl_name & annotation_IDs_with_nonDE_df2barplot$category_2 %in% c("1:1_0_1", "1:1_1", "1:N_0_1", "1:N_1_0", "N:1_0_1", "N:1_1", "N:P_0_1", "N:P_1_0")), "value"]) + sum(gene.data.ensembl.keep[which(gene.data.ensembl.keep[[ensembl_attribute]] %in% as.character(exp_matrix_1_to_2_match_with_nonDE_df[which(exp_matrix_1_to_2_match_with_nonDE_df$category_2 %in% c("1:N_0_P", "1:N_1_1", "1:N_1_P", "1:N_P_0", "1:N_P_1", "1:N_P_Q")),"feature"])), ncbi_attribute] %in% rownames(exp_matrix_2_de_results_df)) + sum(gene.data.ensembl.keep[which(gene.data.ensembl.keep[[ensembl_attribute]] %in% as.character(exp_matrix_1_to_2_match_with_nonDE_df[which(exp_matrix_1_to_2_match_with_nonDE_df$category_2 %in% c("N:P_0_Q", "N:P_1_1", "N:P_1_Q", "N:P_Q_0", "N:P_Q_1", "N:P_Q_R")),"feature"])), ncbi_attribute] %in% rownames(rownames(exp_matrix_2_de_results_df)))
            ncbi_matching_number <- sum(annotation_IDs_with_nonDE_df2barplot[which(annotation_IDs_with_nonDE_df2barplot$annotation==ncbi_name & annotation_IDs_with_nonDE_df2barplot$category_2 %in% c("1:1_0_1", "1:1_1", "1:N_0_1", "1:N_1_0", "N:1_0_1", "N:1_1", "N:P_0_1", "N:P_1_0")), "value"]) + sum(gene.data.ncbi.keep[which(gene.data.ncbi.keep[[ncbi_attribute]] %in% as.character(exp_matrix_2_to_1_match_with_nonDE_df[which(exp_matrix_2_to_1_match_with_nonDE_df$category_2 %in% c("1:N_0_P", "1:N_1_1", "1:N_1_P", "1:N_P_0", "1:N_P_1", "1:N_P_Q")),"feature"])), ensembl_attribute] %in% rownames(exp_matrix_1_de_results_df)) + sum(gene.data.ncbi.keep[which(gene.data.ncbi.keep[[ncbi_attribute]] %in% as.character(exp_matrix_2_to_1_match_with_nonDE_df[which(exp_matrix_2_to_1_match_with_nonDE_df$category_2 %in% c("N:P_0_Q", "N:P_1_1", "N:P_1_Q", "N:P_Q_0", "N:P_Q_1", "N:P_Q_R")),"feature"])), ensembl_attribute] %in% rownames(exp_matrix_1_de_results_df))
            print(sprintf("number of %s to %s ID matchings: %d", ensembl_name, ncbi_name, ensembl_matching_number))
            print(sprintf("number of %s to %s ID matchings: %d", ncbi_name, ensembl_name, ncbi_matching_number))
          } else {
            ### same type of gene IDs
            exp_matrix_1_to_2_matching_number <- sum(annotation_IDs_with_nonDE_df2barplot[which(annotation_IDs_with_nonDE_df2barplot$annotation==exp_matrix_1_name & annotation_IDs_with_nonDE_df2barplot$category_2 %in% c("1:1_0_1", "1:1_1")), "value"])
            exp_matrix_2_to_1_matching_number <- sum(annotation_IDs_with_nonDE_df2barplot[which(annotation_IDs_with_nonDE_df2barplot$annotation==exp_matrix_2_name & annotation_IDs_with_nonDE_df2barplot$category_2 %in% c("1:1_0_1", "1:1_1")), "value"])
            print(sprintf("number of %s to %s ID matchings: %d", exp_matrix_1_name, exp_matrix_2_name, exp_matrix_1_to_2_matching_number))
            print(sprintf("number of %s to %s ID matchings: %d", exp_matrix_2_name, exp_matrix_1_name, exp_matrix_2_to_1_matching_number))
          }
          
          ## get common and specific DE genes
          ### common DE genes
          common_matching_df <- common_features_in_matching_df(exp_matrix_1_to_2_match_with_nonDE_df, exp_matrix_2_to_1_match_with_nonDE_df, gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
          output_name <- sprintf("%s_common", matching_output_name)
          write.csv(common_matching_df, file=sprintf("%s/%s.csv", comparison_FC_output_dir, output_name), quote=FALSE, row.names=FALSE)
          ### specific DE genes
          #### exp matrix 1
          exp_matrix_1_specific_matching_df <- specific_features_in_matching_df(exp_matrix_1_to_2_match_with_nonDE_df, gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
          output_name <- sprintf("%s_%s_specific", matching_output_name, gsub(" ", "", exp_matrix_1_name))
          write.csv(exp_matrix_1_specific_matching_df, file=sprintf("%s/%s.csv", comparison_FC_output_dir, output_name), quote=FALSE, row.names=FALSE)
          #### exp matrix 2
          exp_matrix_2_specific_matching_df <- specific_features_in_matching_df(exp_matrix_2_to_1_match_with_nonDE_df, gsub(" ", "", exp_matrix_2_name), gsub(" ", "", exp_matrix_1_name))
          output_name <- sprintf("%s_%s_specific", matching_output_name, gsub(" ", "", exp_matrix_2_name))
          write.csv(exp_matrix_2_specific_matching_df, file=sprintf("%s/%s.csv", comparison_FC_output_dir, output_name), quote=FALSE, row.names=FALSE)
          
          if (dim(common_matching_df)[1] > 0) {
            ## DE gene p-value scatter plot: common and specific genes
            output_name <- sprintf("%s_DE_pval_per_specificity", comparison_FC_output_basename)
            DE_pval_scatter_plot(exp_matrix_1_de_results_df, exp_matrix_2_de_results_df, common_matching_df, exp_matrix_1_specific_matching_df, exp_matrix_2_specific_matching_df, exp_matrix_1_name, exp_matrix_2_name, comparison_FC_output_dir, output_name)
            
            ## consistency/stringency for DE analysis results
            ### read consistency/stringency results file
            print("Consistency/stringency for DE analysis results")
            #### read files
            consistency_stringency_dir <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/20-Expression_analysis/output/10-Consistency_stringency/%s/%s", work_dir, dataset, comparison)
            if (dataset == "GSE140420") {
              consistency_stringency_dir <- sprintf("%s/all_samples", consistency_stringency_dir)
            }
            consistency_stringency_file <- sprintf("%s/n_%d/%s_%s_consistency_stringency.csv", consistency_stringency_dir, consistency_stringency_n_value_to_add, dataset, comparison)
            consistency_stringency_df <- read.csv(consistency_stringency_file, quote="")
            # consistency_stringency_file <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/20-Expression_analysis/output/10-Consistency_stringency/%s/%s/n_%d/%s_%s_consistency_stringency.csv", work_dir, dataset, comparison, consistency_stringency_n_value_to_add, dataset, comparison)
            # consistency_stringency_df <- read.csv(consistency_stringency_file, quote="")
            #### computations
            # consistency_stringency_df <- consistency_stringency_computation(ncbi_matrix_1, ncbi_matrix_1_name, ncbi_matrix_2, ncbi_matrix_2_name, ncbi_matrix_1_to_2_match_df)
            # write.csv(consistency_stringency_df, file=sprintf("%s/%s_consistency_stringency.csv", out_dir, output_basename), quote=FALSE, row.names=FALSE)
            output_name <- sprintf("%s_consistency_stringency", comparison_FC_output_basename)
            consistency_stringency_plots_for_DE_analysis(consistency_stringency_df, common_matching_df, exp_matrix_1_specific_matching_df, exp_matrix_2_specific_matching_df, exp_matrix_1_name, exp_matrix_2_name, comparison_FC_output_dir, output_name)
            ### gene statistics
            # gtf_file_1 <- genome_list[[exp_matrix_1_name]]$gtf_file
            # gtf_file_2 <- genome_list[[exp_matrix_2_name]]$gtf_file
            gene_stats_file <- sprintf("%s/n_%d/%s_%s_gene_stats.csv", consistency_stringency_dir, consistency_stringency_n_value_to_add, dataset, comparison)
            gene_stats_df <- read.csv(gene_stats_file, quote="")
            # consistency_gene_stats_df <- consistency_stringency_gene_stats_for_DE_analysis(consistency_stringency_df, common_matching_df, exp_matrix_1_specific_matching_df, exp_matrix_2_specific_matching_df, gtf_file_1, gtf_file_2, exp_matrix_1_name, exp_matrix_2_name)
            consistency_gene_stats_df <- consistency_stringency_gene_stats_for_DE_analysis_2(consistency_stringency_df, gene_stats_df, common_matching_df, exp_matrix_1_specific_matching_df, exp_matrix_2_specific_matching_df, exp_matrix_1_name, exp_matrix_2_name)
            consistency_stringency_gene_stats_output_name <- sprintf("%s_gene_stats", output_name)
            write.csv(consistency_gene_stats_df, file=sprintf("%s/%s.csv", comparison_FC_output_dir, consistency_stringency_gene_stats_output_name), quote=FALSE, row.names=FALSE)
            for (threshold in c(1, -log10(log2(1.05)), -log10(log2(1.1)))) {
              #### consistency/stringency categories
              if (threshold == 1) {
                consistency_stringency_category_name <- "consistency_stringency_category"
                threshold_in_path <- sprintf("%d", threshold)
                threshold_to_print <- sprintf("%d", threshold)
              } else {
                consistency_stringency_category_name <- ifelse(threshold == -log10(log2(1.05)), "consistency_stringency_category_2", "consistency_stringency_category_3")
                threshold_in_path <- sub("[.]", "_", sprintf("%.2f", threshold))
                threshold_to_print <- sprintf("%.2f", threshold)
              }
              consistency_stringency_gene_stats_plots_output_basename <- sprintf("%s_thres_%s_plots", consistency_stringency_gene_stats_output_name, threshold_in_path)
              plot_title_base <- sprintf("%s and %s - consistency/stringency threshold: %s", exp_matrix_1_name, exp_matrix_2_name, threshold_to_print)
              consistency_stringency_gene_stats_plots_for_DE_analysis(consistency_gene_stats_df, exp_matrix_1_name, exp_matrix_2_name, consistency_stringency_category_name, plot_title_base, comparison_FC_output_dir, consistency_stringency_gene_stats_plots_output_basename)
            }
          }
          
          # GO term enrichment analysis result comparison according to reference genomes
          comparison_FC_output_dir <- sprintf("%s/comparisons/%s/%s/%s/FC_%s/GO_analysis", output_dir, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
          if (! dir.exists(comparison_FC_output_dir)) {
            dir.create(comparison_FC_output_dir, recursive=TRUE, mode="0775")
          }
          comparison_FC_output_name <- sprintf("%s_%s_%s%s_%s_FC_%s_GO_analysis", dataset, comparison, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold)) # NEW for GSE136913
          
          # Get GO analysis results
          ## annotation 1
          exp_matrix_1_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s/GO_analysis", output_dir, exp_matrix_1_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
          exp_matrix_1_go_output_basename <- sprintf("%s_GO", exp_matrix_1_output_basename) # NEW for GSE136913
          exp_matrix_1_go_results_file <- sprintf("%s/%s_analysis.tsv", exp_matrix_1_ref_genome_output_dir, exp_matrix_1_go_output_basename)
          exp_matrix_1_go_results_df <- read.table(exp_matrix_1_go_results_file, header=TRUE, sep="\t", quote="", row.names=1)
          exp_matrix_1_GO_terms_up <- rownames(exp_matrix_1_go_results_df[which(exp_matrix_1_go_results_df$Ont=="BP" & exp_matrix_1_go_results_df$P.Up < 0.05),])
          exp_matrix_1_GO_terms_down <- rownames(exp_matrix_1_go_results_df[which(exp_matrix_1_go_results_df$Ont=="BP" & exp_matrix_1_go_results_df$P.Down < 0.05),])
          
          ## annotation 2
          exp_matrix_2_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s/GO_analysis", output_dir, exp_matrix_2_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
          exp_matrix_2_go_output_basename <- sprintf("%s_GO", exp_matrix_2_output_basename) # NEW for GSE136913
          exp_matrix_2_go_results_file <- sprintf("%s/%s_analysis.tsv", exp_matrix_2_ref_genome_output_dir, exp_matrix_2_go_output_basename)
          exp_matrix_2_go_results_df <- read.table(exp_matrix_2_go_results_file, header=TRUE, sep="\t", quote="", row.names=1)
          exp_matrix_2_GO_terms_up <- rownames(exp_matrix_2_go_results_df[which(exp_matrix_2_go_results_df$Ont=="BP" & exp_matrix_2_go_results_df$P.Up < 0.05),])
          exp_matrix_2_GO_terms_down <- rownames(exp_matrix_2_go_results_df[which(exp_matrix_2_go_results_df$Ont=="BP" & exp_matrix_2_go_results_df$P.Down < 0.05),])
          
          # Entrez IDs used for GO analysis
          ## annotation 1
          exp_matrix_1_go_entrez_id_file <- sprintf("%s/ID_conversion/%s_Entrez_ID_conversion_retained.csv", exp_matrix_1_ref_genome_output_dir, exp_matrix_1_go_output_basename)
          exp_matrix_1_go_entrez_id_df <- read.csv(exp_matrix_1_go_entrez_id_file, quote="", row.names=1)
          exp_matrix_1_entrez_ID_name <- ifelse(genome_list[[exp_matrix_1_name]]$AnnotationDb_keytype == "SYMBOL", "ENTREZID", "entrezgene_id")
          exp_matrix_1_go_entrez_ids <- exp_matrix_1_go_entrez_id_df[[exp_matrix_1_entrez_ID_name]]
          ## annotation 2
          exp_matrix_2_go_entrez_id_file <- sprintf("%s/ID_conversion/%s_Entrez_ID_conversion_retained.csv", exp_matrix_2_ref_genome_output_dir, exp_matrix_2_go_output_basename)
          exp_matrix_2_go_entrez_id_df <- read.csv(exp_matrix_2_go_entrez_id_file, quote="", row.names=1)
          exp_matrix_2_entrez_ID_name <- ifelse(genome_list[[exp_matrix_2_name]]$AnnotationDb_keytype == "SYMBOL", "ENTREZID", "entrezgene_id")
          exp_matrix_2_go_entrez_ids <- exp_matrix_2_go_entrez_id_df[[exp_matrix_2_entrez_ID_name]]
          ## matchings
          exp_matrix_1_to_2_match_df <- feature_matching(exp_matrix_1_go_entrez_ids, exp_matrix_2_go_entrez_ids)
          write.csv(exp_matrix_1_to_2_match_df, file=sprintf("%s/%s_%s_matching.csv", comparison_FC_output_dir, comparison_FC_output_name, gsub(" ", "", exp_matrix_1_name)), quote=FALSE, row.names=FALSE)
          exp_matrix_2_to_1_match_df <- feature_matching(exp_matrix_2_go_entrez_ids, exp_matrix_1_go_entrez_ids)
          write.csv(exp_matrix_2_to_1_match_df, file=sprintf("%s/%s_%s_matching.csv", comparison_FC_output_dir, comparison_FC_output_name, gsub(" ", "", exp_matrix_2_name)), quote=FALSE, row.names=FALSE)
          ## matching plot
          plot_title <- "Matching for Entrez IDs used for GO analysis"
          feature_matching_plot(exp_matrix_1_to_2_match_df, exp_matrix_2_to_1_match_df, exp_matrix_1_name, exp_matrix_2_name, plot_title, comparison_FC_output_dir, sprintf("%s_Entrez_ID_matching", comparison_FC_output_name))
          ## evidences for specific Entrez IDs
          ### annotation 1
          exp_matrix_1_specific_entrez2GO_df <- GO_terms_per_Entrez_IDs(as.character(exp_matrix_1_to_2_match_df[which(exp_matrix_1_to_2_match_df$category=="specific"), "feature"]), evidence_categories)
          write.csv(exp_matrix_1_specific_entrez2GO_df, file=sprintf("%s/%s_Entrez_ID_%s_specific_GO_evidence.csv", comparison_FC_output_dir, comparison_FC_output_name, gsub(" ", "", exp_matrix_1_name)), quote=FALSE, row.names=FALSE)
          nb_GO_terms_per_Entrez_ID <- unlist(lapply(levels(exp_matrix_1_specific_entrez2GO_df$Entrez_ID), function(x, df=exp_matrix_1_specific_entrez2GO_df){ return(length(df[which(df$Entrez_ID==x), "GO_ID"])) }))
          df2ggplot <- data.frame(count=nb_GO_terms_per_Entrez_ID, annotation=rep(exp_matrix_1_name, length(nb_GO_terms_per_Entrez_ID)))
          counts <- as.numeric(unname(table(exp_matrix_1_specific_entrez2GO_df$evidence)))
          evidence_df <- data.frame(evidence=names(table(exp_matrix_1_specific_entrez2GO_df$evidence)), count=counts, pct=counts/sum(counts)*100)
          evidence_df <- merge(evidence_df, evidence_categories, by="evidence")
          evidence_df$category <- factor(evidence_df$category, levels=c("experimental", "phylogeny", "computational", "author statement", "curator statement", "electronic"))
          evidence_df <- evidence_df[order(evidence_df$category),]
          evidence_df$evidence <- factor(evidence_df$evidence, levels=evidence_df$evidence)
          evidence_df <- cbind(evidence_df, annotation=rep(exp_matrix_1_name, dim(evidence_df)[1]))
          evidence_df2ggplot <- evidence_df
          ### annotation 2
          exp_matrix_2_specific_entrez2GO_df <- GO_terms_per_Entrez_IDs(as.character(exp_matrix_2_to_1_match_df[which(exp_matrix_2_to_1_match_df$category=="specific"), "feature"]), evidence_categories)
          write.csv(exp_matrix_2_specific_entrez2GO_df, file=sprintf("%s/%s_Entrez_ID_%s_specific_GO_evidence.csv", comparison_FC_output_dir, comparison_FC_output_name, gsub(" ", "", exp_matrix_2_name)), quote=FALSE, row.names=FALSE)
          nb_GO_terms_per_Entrez_ID <- unlist(lapply(levels(exp_matrix_2_specific_entrez2GO_df$Entrez_ID), function(x, df=exp_matrix_2_specific_entrez2GO_df){ return(length(df[which(df$Entrez_ID==x), "GO_ID"])) }))
          df2ggplot <- rbind(df2ggplot, data.frame(count=nb_GO_terms_per_Entrez_ID, annotation=rep(exp_matrix_2_name, length(nb_GO_terms_per_Entrez_ID))))
          counts <- as.numeric(unname(table(exp_matrix_2_specific_entrez2GO_df$evidence)))
          evidence_df <- data.frame(evidence=names(table(exp_matrix_2_specific_entrez2GO_df$evidence)), count=counts, pct=counts/sum(counts)*100)
          evidence_df <- merge(evidence_df, evidence_categories, by="evidence")
          evidence_df$category <- factor(evidence_df$category, levels=c("experimental", "phylogeny", "computational", "author statement", "curator statement", "electronic"))
          evidence_df <- evidence_df[order(evidence_df$category),]
          evidence_df$evidence <- factor(evidence_df$evidence, levels=evidence_df$evidence)
          evidence_df <- cbind(evidence_df, annotation=rep(exp_matrix_2_name, dim(evidence_df)[1]))
          evidence_df2ggplot <- rbind(evidence_df2ggplot, evidence_df)
          ### plots
          pdf(sprintf("%s/%s_Entrez_ID_specific.pdf", comparison_FC_output_dir, comparison_FC_output_name))
          histogram <- ggplot(df2ggplot, aes(x=count)) +
            geom_histogram() +
            facet_wrap(~annotation, scales="free") +
            labs(x="Number of GO terms", y="Count") +
            plot_theme
          barplot <- ggplot(evidence_df2ggplot, aes(x=evidence, y=count, fill=category)) +
            geom_bar(stat="identity") +
            geom_text(aes(label=sprintf("%.1f%%", pct)), size=2.5, position=position_stack(vjust=0.5)) +
            facet_wrap(~annotation, scales="free") +
            labs(x="Evidence", y="Count") +
            plot_theme +
            theme(legend.position="bottom")
          multiplot <- ggdraw() +
            draw_plot(histogram, 0, 0.5, 1, 0.45) +
            draw_plot(barplot, 0, 0, 1, 0.45) +
            draw_label("Specific Entrez IDs used for GO analysis", x=0.02, y=0.99, hjust=0, vjust=1, size=12)
          print(multiplot)
          dev.off()
          
          # GO ID matching
          for (dysregulation in c("Up", "Down")) {
            for (pval_correction in c("none", "BH", "BY")) {
              dysregulation_output_name <- sprintf("%s_%s", comparison_FC_output_name, dysregulation)
              pval_colname <- sprintf("P.%s", dysregulation)
              if (pval_correction != "none") {
                dysregulation_output_name <- sprintf("%s_%s", dysregulation_output_name, pval_correction)
                pval_colname <- sprintf("%s.%s", pval_colname, pval_correction)
              }
              # get enriched GO terms for both annotations
              exp_matrix_1_GO_terms <- rownames(exp_matrix_1_go_results_df[which(exp_matrix_1_go_results_df$Ont=="BP" & exp_matrix_1_go_results_df[[pval_colname]] < 0.05),])
              exp_matrix_2_GO_terms <- rownames(exp_matrix_2_go_results_df[which(exp_matrix_2_go_results_df$Ont=="BP" & exp_matrix_2_go_results_df[[pval_colname]] < 0.05),])
              if (length(exp_matrix_1_GO_terms) != 0 & length(exp_matrix_2_GO_terms) != 0) {
                # GO ID matching
                ## exp matrix 1
                ### common and specific GO terms
                print(sprintf("%s GO ID matching: %s and %s", dysregulation, exp_matrix_1_name, exp_matrix_2_name))
                exp_matrix_1_to_2_match_df <- feature_matching(exp_matrix_1_GO_terms, exp_matrix_2_GO_terms)
                write.csv(exp_matrix_1_to_2_match_df, file=sprintf("%s/%s_%s_matching.csv", comparison_FC_output_dir, dysregulation_output_name, gsub(" ", "", exp_matrix_1_name)), quote=FALSE, row.names=FALSE)
                ### GO analysis results for specific enriched GO terms
                go_specific_terms <- get_GO_analysis_specific_results(exp_matrix_1_to_2_match_df, exp_matrix_1_go_results_df, sprintf("P.%s", dysregulation))
                write.table(go_specific_terms, file=sprintf("%s/%s_%s_GO_analysis_specific.tsv", comparison_FC_output_dir, dysregulation_output_name, gsub(" ", "", exp_matrix_1_name)), sep="\t", quote=FALSE, row.names=TRUE)
                
                ## exp matrix 2
                print(sprintf("%s GO ID matching: %s and %s", dysregulation, exp_matrix_2_name, exp_matrix_1_name))
                exp_matrix_2_to_1_match_df <- feature_matching(exp_matrix_2_GO_terms, exp_matrix_1_GO_terms)
                write.csv(exp_matrix_2_to_1_match_df, file=sprintf("%s/%s_%s_matching.csv", comparison_FC_output_dir, dysregulation_output_name, gsub(" ", "", exp_matrix_2_name)), quote=FALSE, row.names=FALSE)
                ### GO analysis results for specific enriched GO terms
                go_specific_terms <- get_GO_analysis_specific_results(exp_matrix_2_to_1_match_df, exp_matrix_2_go_results_df, sprintf("P.%s", dysregulation))
                write.table(go_specific_terms, file=sprintf("%s/%s_%s_GO_analysis_specific.tsv", comparison_FC_output_dir, dysregulation_output_name, gsub(" ", "", exp_matrix_2_name)), sep="\t", quote=FALSE, row.names=TRUE)
                
                # GO ID matching plot
                print(sprintf("%s GO ID matching plot, p-value correction: %s", dysregulation, pval_correction))
                plot_title <- sprintf("%s GO IDs per matching category, p-value correction: %s", dysregulation, pval_correction)
                feature_matching_plot(exp_matrix_1_to_2_match_df, exp_matrix_2_to_1_match_df, exp_matrix_1_name, exp_matrix_2_name, plot_title, comparison_FC_output_dir, sprintf("%s_GO_term_matching", dysregulation_output_name))
              } else {
                if (length(exp_matrix_1_GO_terms_up) == 0) {
                  print(sprintf("No %s enriched GO terms for %s", dysregulation, exp_matrix_1_name))
                }
                if (length(exp_matrix_2_GO_terms_up) == 0) {
                  print(sprintf("No %s enriched GO term for %s", dysregulation, exp_matrix_2_name))
                }
              }
            }
          }
        } else {
          if (length(exp_matrix_1_DE_genes) == 0) {
            print(sprintf("No DE genes for %s", exp_matrix_1_name))
          }
          if (length(exp_matrix_2_DE_genes) == 0) {
            print(sprintf("No DE genes for %s", exp_matrix_2_name))
          }
        }
      }
    }
  }
}
###############################################################


