library(ggplot2)
library(cowplot)
library(ggExtra)


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

# consistency/stringency categories
threshold <- -log(log(1.05, 2), 10)
threshold_in_path <- sub("[.]", "_", sprintf("%.2f", threshold))

# plot theme
plot_theme <- theme_bw() +
  theme(panel.border=element_rect(color="grey50")) +
  theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=8)) +
  theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=7), strip.text=element_text(size=8))
plot_theme_figure <- theme_bw() +
  theme(panel.border=element_rect(color="grey50")) +
  theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=11)) +
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10), strip.text=element_text(size=10))


#############
# Functions #
#############
source(sprintf("%s/ref_genome_eval_functions.R", src_dir))

specific_GO_term_matching_features <- function(specific_genes_1, match_df, entrez_col, match_df_col_1, id_col_1, de_df_1, de_df_2, exon_df, exp_df, name_1, name_2) {
  # annotation_1_DE_specific, match_df, entrez_colname, annotation_1_matching_colname_1, annotation_1_id_colname, annotation_1_de_df, annotation_2_de_df, exon_df, exp_df, annotation_1_name, annotation_2_name
  name_1_col <- gsub(" ", "", name_1)
  name_2_col <- gsub(" ", "", name_2)
  exonlog2ratio_col <- ifelse(sprintf("exon_mergelength_%s_%s_log2ratio", name_1_col, name_2_col) %in% colnames(exon_df), sprintf("exon_mergelength_%s_%s_log2ratio", name_1_col, name_2_col), sprintf("exon_mergelength_%s_%s_log2ratio", name_2_col, name_1_col))
  ## exon length ratio to plot: NCBI RefSeq 108 vs Ensembl 105
  exon_df$exon_mergelength_ratio <- -exon_df[[exonlog2ratio_col]]
  
  # get Entrez IDs with matching IDs in annotation 2
  specific_genes_1_matching <- match_df[which(match_df[[entrez_col]] %in% specific_genes_1 & ! is.na(match_df[[match_df_col_1]])), entrez_col]
  if (length(specific_genes_1_matching) > 0) {
    ## feature IDs
    features <- as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_matching), id_col_1])
    ## matching IDs
    matching_features <- as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_matching), match_df_col_1])
    ## DE results for both annotation
    matching_features_df2ggplot <- data.frame(Entrez_ID=as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_matching), entrez_col]), feature=features, matching_feature=matching_features, DE_pval=de_df_1[as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_matching), id_col_1]), "FDR"], DE_pval_matching=de_df_2[matching_features, "FDR"])
    ## exon length statistics
    matching_features_exon_length_df <- exon_df[which(as.character(exon_df[[sprintf("gene_ID_%s", name_2_col)]]) %in% matching_features),]
    matching_features_df2ggplot <- merge(matching_features_df2ggplot, matching_features_exon_length_df[, c(sprintf("gene_ID_%s", name_2_col), sprintf("mean_exp_%s", name_1_col), sprintf("mean_exp_%s", name_2_col), sprintf("exon_mergelength_%s", name_1_col), sprintf("exon_mergelength_%s", name_2_col), "exon_mergelength_ratio", "consistency_stringency_category_2")], by.x="matching_feature", by.y=sprintf("gene_ID_%s", name_2_col))
    #colnames(matching_features_df2ggplot)[colnames(matching_features_df2ggplot)==exonlog2ratio_col] <- "exon_mergelength_ratio"
    colnames(matching_features_df2ggplot)[colnames(matching_features_df2ggplot)=="consistency_stringency_category_2"] <- "consistency_stringency_category"
    matching_features_df2ggplot$consistency_stringency_category <- factor(matching_features_df2ggplot$consistency_stringency_category)
    matching_features_df2ggplot <- cbind(matching_features_df2ggplot, matching_ID_type=rep("Entrez ID", dim(matching_features_df2ggplot)[1]), annotation=rep(name_1, dim(matching_features_df2ggplot)[1]))
  } else {
    matching_features_df2ggplot <- data.frame()
  }
  
  # get Entrez IDs without matching IDs in annotation 2
  specific_genes_1_no_matching_annot2 <- match_df[which(match_df[[entrez_col]] %in% specific_genes_1 & is.na(match_df[[match_df_col_1]])), entrez_col]
  if (length(specific_genes_1_no_matching_annot2) > 0) {
    # ## Entrez IDs without matching IDs in annotation 2 but with matching IDs in annotation 1
    # specific_genes_1_no_matching_annot2_matching_annot1 <- match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2 & ! is.na(match_df[[match_df_col_1]])), entrez_col]
    # if (length(specific_genes_1_no_matching_annot2_matching_annot1) > 0) {
    #   ### feature IDs
    #   features <- as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2_matching_annot1 & ! is.na(match_df[[match_df_col_1]])), id_col_1])
    #   ### matching features
    #   matching_features <- as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2_matching_annot1 & ! is.na(match_df[[match_df_col_1]])), match_df_col_1])
    #   ### DE results for both annotation
    #   matching_features_df2ggplot_no_matching_annot2_matching_annot1 <- data.frame(Entrez_ID=as.factor(specific_genes_1_no_matching_annot2_matching_annot1), feature=features, matching_feature=matching_features, DE_pval=de_df_1[as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2_matching_annot1), id_col_1]), "FDR"], DE_pval_matching=de_df_2[matching_features, "FDR"])
    #   ### exon length statistics
    #   matching_features_exon_length_df <- exon_df[which(as.character(exon_df[[sprintf("gene_ID_%s", name_2_col)]]) %in% matching_features),]
    #   matching_features_df2ggplot_no_matching_annot2_matching_annot1 <- merge(matching_features_df2ggplot_no_matching_annot2_matching_annot1, matching_features_exon_length_df[, c(sprintf("gene_ID_%s", name_2_col), sprintf("mean_exp_%s", name_1_col), sprintf("mean_exp_%s", name_2_col), sprintf("exon_mergelength_%s", name_1_col), sprintf("exon_mergelength_%s", name_2_col), exonlog2ratio_col, "consistency_stringency_category_2")], by.x="matching_feature", by.y=sprintf("gene_ID_%s", name_2_col))
    #   colnames(matching_features_df2ggplot_no_matching_annot2_matching_annot1)[colnames(matching_features_df2ggplot_no_matching_annot2_matching_annot1)==exonlog2ratio_col] <- "exon_mergelenth_ratio"
    #   matching_features_df2ggplot_no_matching_annot2_matching_annot1 <- cbind(matching_features_df2ggplot_no_matching_annot2_matching_annot1, matching_ID_type=rep("different Entrez IDs", dim(matching_features_df2ggplot_no_matching_annot2_matching_annot1)[1]), annotation=rep(name_1, dim(matching_features_df2ggplot_no_matching_annot2_matching_annot1)[1]))
    #   colnames(matching_features_df2ggplot_no_matching_annot2_matching_annot1)[colnames(matching_features_df2ggplot_no_matching_annot2_matching_annot1)=="consistency_stringency_category_2"] <- "consistency_stringency_category"
    #   matching_features_df2ggplot_no_matching_annot2_matching_annot1$consistency_stringency_category <- factor(matching_features_df2ggplot_no_matching_annot2_matching_annot1$consistency_stringency_category)
    #   if (dim(matching_features_df2ggplot)[1] == 0) {
    #     matching_features_df2ggplot <- matching_features_df2ggplot_no_matching_annot2_matching_annot1
    #   } else {
    #     matching_features_df2ggplot <- rbind(matching_features_df2ggplot, matching_features_df2ggplot_no_matching_annot2_matching_annot1[, colnames(matching_features_df2ggplot)])
    #   }
    # }
    
    # ### Entrez IDs without matching IDs in annotation 2 and in annotation 1: Entrez IDs are specific genes to annotation 1
    # specific_genes_1_no_matching_annot2_no_matching_annot1 <- match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2 & is.na(match_df[[match_df_col_1]])), entrez_col]
    # nb_specific_genes_1_no_matching_annot2_no_matching_annot1 <- length(specific_genes_1_no_matching_annot2_no_matching_annot1)
    # if (nb_specific_genes_1_no_matching_annot2_no_matching_annot1 > 0) {
      nb_specific_genes_1_no_matching_annot2 <- length(specific_genes_1_no_matching_annot2)
      na_values <- rep(NA, nb_specific_genes_1_no_matching_annot2)
      ### feature IDs
      features <- as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2), id_col_1])
      matching_features_df2ggplot_no_matching_annot2_no_matching_annot1 <- data.frame(Entrez_ID=as.factor(specific_genes_1_no_matching_annot2), feature=features, matching_feature=na_values, DE_pval=de_df_1[as.character(match_df[which(match_df[[entrez_col]] %in% specific_genes_1_no_matching_annot2), id_col_1]), "FDR"], DE_pval_matching=na_values)
      mean_exp_df <- exp_df[which(exp_df$annotation==name_1 & exp_df$feature %in% features),]
      matching_features_df2ggplot_no_matching_annot2_no_matching_annot1 <- cbind(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1, mean_exp_1=mean_exp_df[which(!duplicated(mean_exp_df$feature)), "mean_exp"], mean_exp_2=na_values, exon_mergelength_1=na_values, exon_mergelength_2=na_values, exon_mergelength_ratio=na_values, consistency_stringency_category=na_values, matching_ID_type=rep("specific", nb_specific_genes_1_no_matching_annot2), annotation=rep(name_1, nb_specific_genes_1_no_matching_annot2))
      colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)[colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)=="mean_exp_1"] <- sprintf("mean_exp_%s", gsub(" ", "", name_1))
      colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)[colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)=="mean_exp_2"] <- sprintf("mean_exp_%s", gsub(" ", "", name_2))
      colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)[colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)=="exon_mergelength_1"] <- sprintf("exon_mergelength_%s", gsub(" ", "", name_1))
      colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)[colnames(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1)=="exon_mergelength_2"] <- sprintf("exon_mergelength_%s", gsub(" ", "", name_2))
      matching_features_df2ggplot_no_matching_annot2_no_matching_annot1$consistency_stringency_category <- factor(matching_features_df2ggplot_no_matching_annot2_no_matching_annot1$consistency_stringency_category)
      if (dim(matching_features_df2ggplot)[1] == 0) {
        matching_features_df2ggplot <- matching_features_df2ggplot_no_matching_annot2_no_matching_annot1
      } else {
        matching_features_df2ggplot <- rbind(matching_features_df2ggplot, matching_features_df2ggplot_no_matching_annot2_no_matching_annot1[, colnames(matching_features_df2ggplot)])
      }
    # }
  }
  
  return(matching_features_df2ggplot)
}

specific_GO_term_matching_features_outputs <- function(specific_go_term_df, matching_go_term_df, annotation_1_DE_genes_Entrez, annotation_2_DE_genes_Entrez, annotation_1_N_genes_Entrez, annotation_2_N_genes_Entrez, match_df, annotation_1_id_colname, annotation_2_id_colname, annotation_1_matching_colname_1, annotation_2_matching_colname_1, annotation_1_de_df, annotation_2_de_df, exon_df, exp_df, annotation_1_name, annotation_2_name, dys, out_dir, out_name) {
  # annotation_1_go_specific_results_df, annotation_2_go_results_df, annotation_1_DE_genes_Entrez, annotation_2_DE_genes_Entrez, annotation_1_expressed_genes_Entrez, annotation_2_expressed_genes_Entrez, annotation_entrez_ID_df, annotation_1_id_colname, annotation_2_id_colname, annotation_1_matching_colname_1, annotation_2_matching_colname_1, annotation_1_de_results_df, annotation_2_de_results_df, exon_length_df, mean_exp_df, annotation_1_name, annotation_2_name, dysregulation, specific_GO_term_output_dir, specific_GO_output_name
  
  # annotation 1: annotation with specific GO terms
  # annotation 2: the other annotation
  go_term_vector <- rownames(specific_go_term_df)
  
  all_GO_terms_matching_features_df <- data.frame()
  pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=9)
  for (GO_term in go_term_vector) {
    print(sprintf("GO term: %s", GO_term))
    annotation_1_matching_features_df <- data.frame()
    annotation_2_matching_features_df <- data.frame()
    
    # get Entrez IDs related to GO term
    GO_Entrez_IDs <- as.integer(GO2EntrezID(GO_term))
    
    # Entrez IDs related related to GO term in DE genes
    ## annotation 1
    annotation_1_DE_genes_Entrez_IDs_GO_term <- GO_Entrez_IDs[GO_Entrez_IDs %in% annotation_1_DE_genes_Entrez]
    ## annotation 2
    annotation_2_DE_genes_Entrez_IDs_GO_term <- GO_Entrez_IDs[GO_Entrez_IDs %in% annotation_2_DE_genes_Entrez]
    ## specific Entrez IDs
    ### annotation 1 specific
    annotation_1_DE_specific <- annotation_1_DE_genes_Entrez_IDs_GO_term[which(! annotation_1_DE_genes_Entrez_IDs_GO_term %in% annotation_2_DE_genes_Entrez_IDs_GO_term)]
    nb_annotation_1_DE_specific <- length(annotation_1_DE_specific)
    ### annotation 2 specific
    annotation_2_DE_specific <- annotation_2_DE_genes_Entrez_IDs_GO_term[which(! annotation_2_DE_genes_Entrez_IDs_GO_term %in% annotation_1_DE_genes_Entrez_IDs_GO_term)]
    nb_annotation_2_DE_specific <- length(annotation_2_DE_specific)
    
    # annotation 1 specific DE genes
    if (nb_annotation_1_DE_specific > 0) {
      entrez_colname <- "ENTREZID"
      annotation_1_DE_matching_features_df <- specific_GO_term_matching_features(annotation_1_DE_specific, match_df, entrez_colname, annotation_1_matching_colname_1, annotation_1_id_colname, annotation_1_de_df, annotation_2_de_df, exon_df, exp_df, annotation_1_name, annotation_2_name)
      annotation_1_DE_matching_features_df <- cbind(annotation_1_DE_matching_features_df, category=rep(dys, dim(annotation_1_DE_matching_features_df)[1]))
      annotation_1_matching_features_df <- rbind(annotation_1_matching_features_df, annotation_1_DE_matching_features_df)
    }
    
    # annotation 2 specific DE genes
    if (nb_annotation_2_DE_specific > 0) {
      entrez_colname <- "ENTREZID"
      annotation_2_DE_matching_features_df <- specific_GO_term_matching_features(annotation_2_DE_specific, match_df, entrez_colname, annotation_2_matching_colname_1, annotation_2_id_colname, annotation_2_de_df, annotation_1_de_df, exon_df, exp_df, annotation_2_name, annotation_1_name)
      annotation_2_DE_matching_features_df <- cbind(annotation_2_DE_matching_features_df, category=rep(dys, dim(annotation_2_DE_matching_features_df)[1]))
      annotation_2_matching_features_df <- rbind(annotation_2_matching_features_df, annotation_2_DE_matching_features_df)
    }
    
    # Entrez IDs related to GO term in background genes
    ## annotation 1
    annotation_1_DE_genes_Entrez_IDs_GO_term <- GO_Entrez_IDs[GO_Entrez_IDs %in% annotation_1_N_genes_Entrez]
    ## annotation 2
    annotation_2_DE_genes_Entrez_IDs_GO_term <- GO_Entrez_IDs[GO_Entrez_IDs %in% annotation_2_N_genes_Entrez]
    # specific Entrez IDs
    ## annotation 1 specific
    annotation_1_N_specific <- annotation_1_DE_genes_Entrez_IDs_GO_term[which(! annotation_1_DE_genes_Entrez_IDs_GO_term %in% annotation_2_DE_genes_Entrez_IDs_GO_term)]
    nb_annotation_1_N_specific <- length(annotation_1_N_specific)
    ## annotation 2 specific
    annotation_2_N_specific <- annotation_2_DE_genes_Entrez_IDs_GO_term[which(! annotation_2_DE_genes_Entrez_IDs_GO_term %in% annotation_1_DE_genes_Entrez_IDs_GO_term)]
    nb_annotation_2_N_specific <- length(annotation_2_N_specific)
    
    # annotation 1 specific background genes
    if (nb_annotation_1_N_specific > 0) {
      entrez_colname <- "ENTREZID"
      annotation_1_N_matching_features_df <- specific_GO_term_matching_features(annotation_1_N_specific, match_df, entrez_colname, annotation_1_matching_colname_1, annotation_1_id_colname, annotation_1_de_df, annotation_2_de_df, exon_df, exp_df, annotation_1_name, annotation_2_name)
      annotation_1_N_matching_features_df <- cbind(annotation_1_N_matching_features_df, category=rep("background", dim(annotation_1_N_matching_features_df)[1]))
      annotation_1_matching_features_df <- rbind(annotation_1_matching_features_df, annotation_1_N_matching_features_df)
    }
    # annotation 2 specific background genes
    if (nb_annotation_2_N_specific > 0) {
      entrez_colname <- "ENTREZID"
      annotation_2_N_matching_features_df <- specific_GO_term_matching_features(annotation_2_N_specific, match_df, entrez_colname, annotation_2_matching_colname_1, annotation_2_id_colname, annotation_2_de_df, annotation_1_de_df, exon_df, exp_df, annotation_2_name, annotation_1_name)
      annotation_2_N_matching_features_df <- cbind(annotation_2_N_matching_features_df, category=rep("background", dim(annotation_2_N_matching_features_df)[1]))
      annotation_2_matching_features_df <- rbind(annotation_2_matching_features_df, annotation_2_N_matching_features_df)
    }
    
    # plots
    if (dim(annotation_1_matching_features_df)[1] > 0 | dim(annotation_2_matching_features_df)[1] > 0) {
      if (dim(annotation_1_matching_features_df)[1] > 0) {
        if (dim(annotation_2_matching_features_df)[1] > 0) {
          matching_features_df <- rbind(annotation_1_matching_features_df, annotation_2_matching_features_df[, colnames(annotation_1_matching_features_df)])
        } else {
          matching_features_df <- annotation_1_matching_features_df
        }
      } else {
        if (dim(annotation_2_matching_features_df)[1] > 0) {
          matching_features_df <- annotation_2_matching_features_df
        }
      }
      
      category_vector <- annotation_vector <- matching_type_vector <- count_vector <- c()
      for (one_category in levels(matching_features_df$category)) {
        one_category_df <- matching_features_df[which(matching_features_df$category==one_category),]
        one_category_df <- droplevels(one_category_df)
        for (one_annotation in levels(one_category_df$annotation)) {
          one_annotation_df <- one_category_df[which(one_category_df$annotation==one_annotation),]
          one_annotation_df <- droplevels(one_annotation_df)
          for (one_matching_type in levels(one_annotation_df$matching_ID_type)) {
            nb_matching_type <- dim(one_annotation_df[which(one_annotation_df$matching_ID_type==one_matching_type),])[1]
            matching_type_vector <- c(matching_type_vector, one_matching_type)
            count_vector <- c(count_vector, nb_matching_type)
          }
          annotation_vector <- c(annotation_vector, rep(one_annotation, length(levels(one_annotation_df$matching_ID_type))))
          category_vector <- c(category_vector, rep(one_category, length(levels(one_annotation_df$matching_ID_type))))
        }
      }
      df2ggplot <- data.frame(category=category_vector, annotation=annotation_vector, matching_type=matching_type_vector, count=count_vector)
      
      
      GO_term_name <- as.character(specific_go_term_df[GO_term, "Term"])
      annotation_1_GO_term_N <- specific_go_term_df[GO_term, "N"]
      annotation_1_GO_term_DE <- specific_go_term_df[GO_term, dys]
      annotation_1_GO_term_pval <- specific_go_term_df[GO_term, sprintf("P.%s.BY", dys)]
      annotation_2_GO_term_N <- matching_go_term_df[GO_term, "N"]
      annotation_2_GO_term_DE <- matching_go_term_df[GO_term, dys]
      annotation_2_GO_term_pval <- matching_go_term_df[GO_term, sprintf("P.%s.BY", dys)]
      
      plot_title <- sprintf("%s: %s\n%s: N=%d, %s=%d, p-value=%.2e\n%s: N=%d, %s=%d, p-value=%.2e\nDE genes: %d %s specific DE genes, %d %s specific DE genes", GO_term, GO_term_name, annotation_1_name, annotation_1_GO_term_N, dys, annotation_1_GO_term_DE, annotation_1_GO_term_pval, annotation_2_name, annotation_2_GO_term_N, dys, annotation_2_GO_term_DE, annotation_2_GO_term_pval, nb_annotation_1_DE_specific, annotation_1_name, nb_annotation_2_DE_specific, annotation_2_name)
      
      # counts per annotation and matching ID type
      barplot <- ggplot(df2ggplot, aes(x=annotation, y=count, fill=matching_type)) +
        geom_bar(stat="identity") +
        geom_text(aes(label=sprintf("%d", count)), size=3, position=position_stack(vjust=0.5)) +
        facet_wrap(~category, nrow=2, scales="free_y") +
        labs(x="Annotation", y="# IDs", fill="ID matching type") +
        plot_theme +
        theme(legend.position="bottom") +
        guides(fill=guide_legend(nrow=length(levels(df2ggplot$matching_type)), byrow=TRUE))
      
      # specific features
      ## scatter plot: mean expression
      specific_features_df <- matching_features_df[which(matching_features_df$matching_ID_type=="specific"),]
      if (dim(specific_features_df)[1] > 0) {
        specific_features_plot <- TRUE
        specific_features_df <- droplevels(specific_features_df)
        annotation_1_specific_features_df <- specific_features_df[which(specific_features_df$annotation==annotation_1_name),]
        if (dim(annotation_1_specific_features_df)[1]) {
          annotation_1_specific_features_df <- droplevels(annotation_1_specific_features_df)
          specific_features_df2ggplot <- data.frame(DE_pval=annotation_1_specific_features_df$DE_pval, mean_exp=annotation_1_specific_features_df[, sprintf("mean_exp_%s", gsub(" ", "", annotation_1_name))], annotation=annotation_1_specific_features_df$annotation, category=annotation_1_specific_features_df$category)
        } else {
          specific_features_df2ggplot <- data.frame()
        }
        annotation_2_specific_features_df <- specific_features_df[which(specific_features_df$annotation==annotation_2_name),]
        if (dim(annotation_2_specific_features_df)[1]) {
          annotation_2_specific_features_df <- droplevels(annotation_2_specific_features_df)
          specific_features_df2ggplot <- rbind(specific_features_df2ggplot, data.frame(DE_pval=annotation_2_specific_features_df$DE_pval, mean_exp=annotation_2_specific_features_df[, sprintf("mean_exp_%s", gsub(" ", "", annotation_2_name))], annotation=annotation_2_specific_features_df$annotation, category=annotation_2_specific_features_df$category))
        }
        
        scatterplot_mean_exp_specific <- ggplot(specific_features_df2ggplot, aes(x=DE_pval, y=mean_exp, color=annotation)) +
          geom_point(size=2) +
          xlim(0, 1) +
          geom_vline(xintercept=0, alpha=0.5) + 
          geom_hline(yintercept=0, alpha=0.5) +
          geom_vline(xintercept=0.05, linetype="dashed", alpha=0.5) +
          facet_wrap(~category, nrow=1) +
          labs(x=sprintf("DE p-value", annotation_1_name), y="Mean expression", color="Annotation") +
          plot_theme +
          theme(legend.position="bottom", legend.box="vertical")
      } else {
        specific_features_plot <- FALSE
      }
      
      # matching features
      ## scatter plot: matching ID DE p-value and exon mergelenth ratio
      matching_features_df2ggplot <- matching_features_df[which(! is.na(matching_features_df$DE_pval_matching)),]
      matching_features_df2ggplot <- droplevels(matching_features_df2ggplot)
      levels(matching_features_df2ggplot$consistency_stringency_category) <- unlist(lapply(levels(matching_features_df2ggplot$consistency_stringency_category), function(x) { return(sub("_", " ", x)) }))
      if (dim(matching_features_df2ggplot)[1] > 0) {
        scatterplot <- ggplot(matching_features_df2ggplot, aes(x=DE_pval_matching, y=exon_mergelength_ratio, color=annotation, shape=consistency_stringency_category)) +
          geom_point(size=2) +
          xlim(0, 1) +
          geom_vline(xintercept=0, alpha=0.5) + 
          geom_hline(yintercept=0, alpha=0.5) +
          geom_vline(xintercept=0.05, linetype="dashed", alpha=0.5) +
          facet_wrap(~category, nrow=2) +
          labs(x=sprintf("matching annotation DE p-value", annotation_1_name), y=sprintf("%s/%s exon length log2ratio", annotation_1_name, annotation_2_name), color="Annotation", shape="Consistency/stringency") +
          plot_theme +
          theme(legend.position="bottom", legend.box="vertical") +
          guides(shape=guide_legend(nrow=length(levels(matching_features_df2ggplot$consistency_stringency_category)), byrow=TRUE))
        
        scatterplot_mean_exp_1 <- ggplot(matching_features_df2ggplot, aes_string(x="DE_pval_matching", y=sprintf("mean_exp_%s", gsub(" ", "", annotation_1_name)), color="annotation", shape="consistency_stringency_category")) +
          geom_point(size=2) +
          xlim(0, 1) +
          geom_vline(xintercept=0, alpha=0.5) + 
          geom_hline(yintercept=0, alpha=0.5) +
          geom_vline(xintercept=0.05, linetype="dashed", alpha=0.5) +
          facet_wrap(~category, nrow=2) +
          labs(x=sprintf("matching annotation DE p-value", annotation_1_name), y=sprintf("%s mean expression", annotation_1_name), color="Annotation", shape="Consistency/stringency") +
          plot_theme +
          theme(legend.position="bottom", legend.box="vertical") +
          guides(shape=guide_legend(nrow=length(levels(matching_features_df2ggplot$consistency_stringency_category)), byrow=TRUE))
        
        scatterplot_mean_exp_2 <- ggplot(matching_features_df2ggplot, aes_string(x="DE_pval_matching", y=sprintf("mean_exp_%s", gsub(" ", "", annotation_2_name)), color="annotation", shape="consistency_stringency_category")) +
          geom_point(size=2) +
          xlim(0, 1) +
          geom_vline(xintercept=0, alpha=0.5) + 
          geom_hline(yintercept=0, alpha=0.5) +
          geom_vline(xintercept=0.05, linetype="dashed", alpha=0.5) +
          facet_wrap(~category, nrow=2) +
          labs(x=sprintf("matching annotation DE p-value", annotation_1_name), y=sprintf("%s mean expression", annotation_2_name), color="Annotation", shape="Consistency/stringency") +
          plot_theme +
          theme(legend.position="bottom", legend.box="vertical") +
          guides(shape=guide_legend(nrow=length(levels(matching_features_df2ggplot$consistency_stringency_category)), byrow=TRUE))
        
        multiplot1 <- ggdraw() +
          draw_plot(barplot, 0, 0, 0.5, 0.9) +
          draw_plot(scatterplot, 0.5, 0, 0.5, 0.9) +
          draw_label(plot_title, x=0.02, y=0.99, hjust=0, vjust=1, size=12)
        
        if (specific_features_plot) {
          multiplot2 <- ggdraw() +
            draw_plot(scatterplot_mean_exp_1, 0, 0, 0.33, 0.9) +
            draw_plot(scatterplot_mean_exp_2, 0.33, 0, 0.33, 0.9) +
            draw_plot(scatterplot_mean_exp_specific, 0.66, 0, 0.33, 0.9) +
            draw_label(plot_title, x=0.02, y=0.99, hjust=0, vjust=1, size=12)
        } else {
          multiplot2 <- ggdraw() +
            draw_plot(scatterplot_mean_exp_1, 0, 0, 0.5, 0.9) +
            draw_plot(scatterplot_mean_exp_2, 0.5, 0, 0.5, 0.9) +
            draw_label(plot_title, x=0.02, y=0.99, hjust=0, vjust=1, size=12)
        }
        print(multiplot1)
        print(multiplot2)
      } else {
        if (specific_features_plot) {
          multiplot1 <- ggdraw() +
            draw_plot(barplot, 0, 0, 0.5, 0.9) +
            draw_plot(scatterplot_mean_exp_specific, 0.5, 0, 0.5, 0.9) +
            draw_label(plot_title, x=0.02, y=0.99, hjust=0, vjust=1, size=12)
        } else {
          multiplot1 <- ggdraw() +
            draw_plot(barplot, 0, 0, 1, 0.9) +
            draw_label(plot_title, x=0.02, y=0.99, hjust=0, vjust=1, size=12)
        }
        print(multiplot1)
      }
      all_GO_terms_matching_features_df <- rbind(all_GO_terms_matching_features_df, cbind(GO_term=rep(GO_term, dim(matching_features_df)[1]), matching_features_df))
    }
  }
  dev.off()
  write.csv(all_GO_terms_matching_features_df, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  write.csv(all_GO_terms_matching_features_df[which(all_GO_terms_matching_features_df$category != "background"),], file=sprintf("%s/%s_supp_table.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  return(all_GO_terms_matching_features_df)
}

specific_GO_term_matching_features_outputs_all_plots <- function(go_specific_matching_features_df, out_dir, out_name) {
  
  go_specific_matching_features_df2ggplot <- go_specific_matching_features_df[which(! is.na(go_specific_matching_features_df$DE_pval_matching)),]
  go_specific_matching_features_df2ggplot <- droplevels(go_specific_matching_features_df2ggplot)
  levels(go_specific_matching_features_df2ggplot$consistency_stringency_category) <- unlist(lapply(levels(go_specific_matching_features_df2ggplot$consistency_stringency_category), function(x) { return(sub("_", " ", x)) }))
  levels(go_specific_matching_features_df2ggplot$annotation) <- unlist(lapply(levels(go_specific_matching_features_df2ggplot$annotation), function(x) { return(sub("NCBI ", "", x)) }))
  go_specific_matching_features_df2ggplot$annotation <- factor(go_specific_matching_features_df2ggplot$annotation, levels=sort(levels(go_specific_matching_features_df2ggplot$annotation)))
  
  ## remove duplicated lines: same genes found for multiple different GO terms
  go_specific_matching_features_df2ggplot_nodup <- go_specific_matching_features_df2ggplot[! duplicated(go_specific_matching_features_df2ggplot[,2:dim(go_specific_matching_features_df2ggplot)[2]]),]
  ## exon length ratio to plot: NCBI RefSeq 108 vs Ensembl 105
  annotation_1_name <- "NCBI RefSeq 108"
  annotation_2_name <- "Ensembl 105"
  go_specific_matching_features_df2ggplot_nodup$exon_mergelenth_ratio2plot <- log(go_specific_matching_features_df2ggplot_nodup$exon_mergelength_NCBIRefSeq108/go_specific_matching_features_df2ggplot_nodup$exon_mergelength_Ensembl105, 2)
  write.csv(go_specific_matching_features_df2ggplot_nodup, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=9)
  # exon length ratio and matching annotation DE p-value scatter plot with density margins for specific DE genes of all specific GO terms
  scatterplot <- ggplot(go_specific_matching_features_df2ggplot_nodup, aes(x=DE_pval_matching, y=exon_mergelenth_ratio2plot)) +
    geom_point(aes(color=annotation), alpha=0) + # fake call to geom_point() to make ggMarginal() plot densities only according to color and not according to both color and shape
    geom_point(aes(color=annotation, shape=consistency_stringency_category), size=2) +
    xlim(0, 1) +
    geom_vline(xintercept=0, alpha=0.5) + 
    geom_hline(yintercept=0, alpha=0.5) +
    geom_vline(xintercept=0.05, linetype="dashed", alpha=0.5) +
    labs(x=sprintf("other genome\ndifferential expression p-value", sub("NCBI ", "", annotation_1_name)), y=sprintf("%s/%s\nexon length log2 ratio", sub("NCBI ", "", annotation_1_name), annotation_2_name), color="Specificity", shape="Category") +
    plot_theme_figure +
    theme(legend.position="bottom") +
    guides(shape=guide_legend(nrow=2, byrow=TRUE))
  ## get plot ranges
  scatter_y_scale <- layer_scales(scatterplot)$y$range$range
  scatter_y_min <- scatter_y_scale[1]
  scatter_y_max <- scatter_y_scale[2]
  scatter_y_range <- scatter_y_max - scatter_y_scale[1]
  
  # exon length log ratio distributions
  exon_mergelength_density <- ggplot(go_specific_matching_features_df2ggplot_nodup, aes(x=exon_mergelenth_ratio2plot, color=annotation)) +
    geom_line(stat="density") +
    geom_vline(xintercept=0, linetype="dashed", alpha=0.5) +
    labs(title="Exon length ratio per DE specificity", x=sprintf("%s/%s exon length log2ratio", sub("NCBI ", "", annotation_1_name), annotation_2_name), y="Density", color="Annotation") +
    plot_theme +
    theme(legend.position="bottom", legend.box="vertical")
  ## get plot ranges
  density_x_scale <- layer_scales(exon_mergelength_density)$x$range$range
  density_x_min <- density_x_scale[1]
  density_x_max <- density_x_scale[2]
  density_x_range <- density_x_max - density_x_min
  density_y_scale <- layer_scales(exon_mergelength_density)$y$range$range
  density_y_max <- density_y_scale[2]
  density_y_range <- density_y_max - density_y_scale[1]
  
  # exon length ratio distribution tests: build data frame to write p-values on plots
  annotation_vector <- label_vector <- density_x_vector <- density_y_vector <- scatter_x_vector <- scatter_y_vector <- c()
  i <- 1
  for (one_annotation in levels(go_specific_matching_features_df2ggplot_nodup$annotation)) {
    # annotation_vector <- c(annotation_vector, one_annotation)
    annotation_vector <- c(annotation_vector, rep(one_annotation, 2))
    annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation <- go_specific_matching_features_df2ggplot_nodup[which(go_specific_matching_features_df2ggplot_nodup$annotation == one_annotation),]
    nb_genes <- dim(annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation)[1]
    if (nb_genes >= 3) {
      differences <- annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1_name))]] - annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2_name))]]
      normality_test_pval <- shapiro.test(differences)$p.value
      # test_alternative <- ifelse(mean(annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation$exon_mergelenth_ratio2plot) < 0, "less", "greater")
      # alternative: greater (i.e. longer exon in RefSeq 108)
      wilcoxon_test_result <- wilcox.test(annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1_name))]], annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2_name))]], paired=TRUE, alternative="greater")
      label <- sprintf("p = %.2e", wilcoxon_test_result$p.value)
      label_vector <- c(label_vector, label)
      # t_test_result <- t.test(annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1_name))]], annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2_name))]], paired=TRUE, alternative=test_alternative)
      # alternative: less (i.e. shorter exon in RefSeq 108)
      wilcoxon_test_result <- wilcox.test(annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1_name))]], annotation_1_go_specific_matching_features_df2ggplot_nodup_one_annotation[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2_name))]], paired=TRUE, alternative="less")
      label <- sprintf("p = %.2e", wilcoxon_test_result$p.value)
      label_vector <- c(label_vector, label)
    } else {
      label_vector <- c(label_vector, rep("NA", 2))
    }
    density_x_vector <- c(density_x_vector, c(density_x_min+0.15*density_x_range, density_x_min+0.15*density_x_range))
    density_y_vector <- c(density_y_vector, c(density_y_max-0.05*density_y_range-((i-1)*0.1*density_y_range), density_y_max-0.05*density_y_range-(i*0.1*density_y_range)))
    scatter_x_vector <- c(scatter_x_vector, c(0.87, 0.87))
    scatter_y_vector <- c(scatter_y_vector, c(scatter_y_max-0.01*scatter_y_range-((i-1)*0.03*scatter_y_range), scatter_y_min+0.01*scatter_y_range-((i-1)*0.03*scatter_y_range)))
    i <- i + 2
  }
  density_annotation_df <- data.frame(annotation=annotation_vector, x=density_x_vector, y=density_y_vector, label=label_vector)
  scatter_annotation_df <- data.frame(annotation=annotation_vector, x=scatter_x_vector, y=scatter_y_vector, label=label_vector)
  ## write p-values on plots
  exon_mergelength_density <- exon_mergelength_density +
    geom_text(data=density_annotation_df, aes(x=x, y=y, label=label, color=annotation), size=3.5, show.legend=FALSE)
  print(exon_mergelength_density)
  scatterplot <- scatterplot +
    geom_text(data=scatter_annotation_df, aes(x=x, y=y, label=label, color=annotation), size=3.5, show.legend=FALSE)
  print(scatterplot)
  scatterplot2return <- scatterplot
  ## remove legend
  scatterplot <- scatterplot + theme(legend.position="none")
  
  # add distributions to scatter plot
  scatterplot2 <- ggMarginal(scatterplot, type="density", groupColour=TRUE)
  plot.new()
  print(scatterplot2)
  dev.off()
  
  return(list(scatter_plot=scatterplot2return, scatter_plot_marginal=scatterplot2, nb_categories=length(levels(go_specific_matching_features_df2ggplot_nodup$consistency_stringency_category)), nb_annotations=length(levels(go_specific_matching_features_df2ggplot_nodup$annotation))))
}

############
# Analysis #
############
# NCBI RefSeq 108 vs Ensembl 105 comparison
annotation_1_name <- "NCBI RefSeq 108"
annotation_2_name <- "Ensembl 105"
annotation_1_ref_genome <- genome_list[[annotation_1_name]]$ref_genome
annotation_2_ref_genome <- genome_list[[annotation_2_name]]$ref_genome
comparison <- sprintf("%s_%s", gsub(" ", "", annotation_1_name), gsub(" ", "", annotation_2_name))
print(sprintf("comparison: %s", comparison))

# design and contrast name
## GSE136913
# design_formula <- "~Instrument+Litter+Condition"
# design_output_name_append <- NULL
# contrast_output_names <- c("kainatevscontrol")
## GSE140420
### Region, blocking: Age
# design_formula <- "~Age+Region"
# design_output_name_append <- NULL
# contrast_output_names <- c("CA1vsCA3DG", "CA3vsCA1DG", "DGvsCA1CA3")
## GSE179101
design_formula <- "~Subfield"
design_output_name_append <- NULL
contrast_output_names <- c("CA1vsCA2CA3DG", "CA2vsCA1CA3DG", "CA3vsCA1CA2DG", "DGvsCA1CA2CA3")

fc_thresholds <- c(1, 1.1, 1.25, 1.5, 2)



if (is.null(design_output_name_append)) {
  design_output_name_append <- ""
}
design_output_name <- sub("~|~0\\+", "", design_formula)
design_output_name <- gsub("\\+", "_", design_output_name)
design_output_name <- gsub("\\.", "", design_output_name)


for (one_contrast_output_name in contrast_output_names) {
  print(sprintf("contrast: %s", one_contrast_output_name))
  for (one_fc_threshold in fc_thresholds) {
    print(sprintf("FC: %s", one_fc_threshold))
    
    # GO term enrichment analysis result comparison according to reference genomes
    comparison_FC_output_dir <- sprintf("%s/comparisons/%s/%s/%s/FC_%s", output_dir, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
    comparison_FC_output_name <- sprintf("%s_%s_%s%s_%s_FC_%s_GO_analysis", dataset, comparison, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
    specific_GO_term_output_dir <- sprintf("%s/GO_analysis/specific_GO_terms", comparison_FC_output_dir)
    if (! dir.exists(specific_GO_term_output_dir)) {
      dir.create(specific_GO_term_output_dir, recursive=TRUE, mode="0775")
    }
    
    scatter_plot_list <- list()
    for (dysregulation in c("Up", "Down")) {
      print(sprintf("dysregulation: %s", dysregulation))
      
      annotation_1_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s", output_dir, annotation_1_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
      annotation_1_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, annotation_1_ref_genome, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
      annotation_2_ref_genome_output_dir <- sprintf("%s/%s/%s/%s/FC_%s", output_dir, annotation_2_ref_genome, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
      annotation_2_output_basename <- sprintf("%s_%s_%s%s_%s_FC_%s", dataset, annotation_2_ref_genome, design_output_name, design_output_name_append, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
      
      # Check that GO analysis results files exist
      annotation_1_ref_genome_go_output_dir <- sprintf("%s/GO_analysis", annotation_1_ref_genome_output_dir)
      annotation_1_go_results_file <- sprintf("%s/%s_GO_analysis.tsv", annotation_1_ref_genome_go_output_dir, annotation_1_output_basename)
      annotation_2_ref_genome_go_output_dir <- sprintf("%s/GO_analysis", annotation_2_ref_genome_output_dir)
      annotation_2_go_results_file <- sprintf("%s/%s_GO_analysis.tsv", annotation_2_ref_genome_go_output_dir, annotation_2_output_basename)
      if (file.exists(annotation_1_go_results_file) & file.exists(annotation_2_go_results_file)) {
        # Get GO analysis results and specific enriched GO term analysis results
        ## annotation 1
        ### GO analysis results
        annotation_1_go_results_df <- read.table(annotation_1_go_results_file, header=TRUE, sep="\t", quote="", row.names=1)
        ### specific enriched GO term analysis results
        annotation_1_go_specific_results_file <- sprintf("%s/GO_analysis/%s_%s_BY_%s_GO_analysis_specific.tsv", comparison_FC_output_dir, comparison_FC_output_name, dysregulation, gsub(" ", "", annotation_1_name))
        if (file.exists(annotation_1_go_specific_results_file)) {
          annotation_1_go_specific_results_df <- read.table(annotation_1_go_specific_results_file, header=TRUE, sep="\t", quote="", row.names=1)
        }
        ## annotation 2
        ### GO analysis results
        annotation_2_go_results_df <- read.table(annotation_2_go_results_file, header=TRUE, sep="\t", quote="", row.names=1)
        ### specific enriched GO term analysis results
        annotation_2_go_specific_results_file <- sprintf("%s/GO_analysis/%s_%s_BY_%s_GO_analysis_specific.tsv", comparison_FC_output_dir, comparison_FC_output_name, dysregulation, gsub(" ", "", annotation_2_name))
        if (file.exists(annotation_2_go_specific_results_file)) {
          annotation_2_go_specific_results_df <- read.table(annotation_2_go_specific_results_file, header=TRUE, sep="\t", quote="", row.names=1)
        }
        
        # get DE results
        ## annotation 1
        annotation_1_de_results_file <- sprintf("%s/%s_pval.csv", annotation_1_ref_genome_output_dir, annotation_1_output_basename)
        annotation_1_de_results_df <- read.csv(annotation_1_de_results_file, quote="", row.names=1)
        if (dysregulation == "Up") {
          annotation_1_DE_genes <- rownames(annotation_1_de_results_df[which(annotation_1_de_results_df$FDR < 0.05 & annotation_1_de_results_df$logFC > 0),])
        } else {
          if (dysregulation == "Down") {
            annotation_1_DE_genes <- rownames(annotation_1_de_results_df[which(annotation_1_de_results_df$FDR < 0.05 & annotation_1_de_results_df$logFC < 0),])
          }
        }
        ## annotation 2
        annotation_2_de_results_file <- sprintf("%s/%s_pval.csv", annotation_2_ref_genome_output_dir, annotation_2_output_basename)
        annotation_2_de_results_df <- read.csv(annotation_2_de_results_file, quote="", row.names=1)
        if (dysregulation == "Up") {
          annotation_2_DE_genes <- rownames(annotation_2_de_results_df[which(annotation_2_de_results_df$FDR < 0.05 & annotation_2_de_results_df$logFC > 0),])
        } else {
          if (dysregulation == "Down") {
            annotation_2_DE_genes <- rownames(annotation_2_de_results_df[which(annotation_2_de_results_df$FDR < 0.05 & annotation_2_de_results_df$logFC < 0),])
          }
        }
        
        # gene IDs used for GO analysis
        ## annotation 1
        annotation_1_entrez_ID_file <- sprintf("%s/ID_conversion/%s_GO_Entrez_ID_conversion_retained.csv", annotation_1_ref_genome_go_output_dir, annotation_1_output_basename)
        annotation_1_entrez_ID_df <- read.csv(annotation_1_entrez_ID_file, quote="")
        ## annotation 2
        annotation_2_entrez_ID_file <- sprintf("%s/ID_conversion/%s_GO_Entrez_ID_conversion_retained.csv", annotation_2_ref_genome_go_output_dir, annotation_2_output_basename)
        annotation_2_entrez_ID_df <- read.csv(annotation_2_entrez_ID_file, quote="")
        ## merge
        annotation_entrez_ID_df <- merge(annotation_1_entrez_ID_df, annotation_2_entrez_ID_df, by.x="ENTREZID", by.y="entrezgene_id", all=TRUE)
        
        # gene ID matching
        consistency_dir <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/20-Expression_analysis/output/10-Consistency_stringency/%s/%s", work_dir, dataset, comparison)
        if (dataset == "GSE140420") {
          consistency_dir <- sprintf("%s/all_samples", consistency_dir)
        }
        consistency_dir <- sprintf("%s/n_%d", consistency_dir, consistency_stringency_n_value_to_add)
        annotation_1_to_2_match_file <- sprintf("%s/gene_ID_matching/%s_%s_useMart_getBM_matching_categories.csv", consistency_dir, dataset, gsub(" ", "", annotation_1_name))
        annotation_1_to_2_match_df <- read.csv(annotation_1_to_2_match_file, quote="")
        
        annotation_2_to_1_match_file <- sprintf("%s/gene_ID_matching/%s_%s_useMart_getBM_matching_categories.csv", consistency_dir, dataset, gsub(" ", "", annotation_2_name))
        annotation_2_to_1_match_df <- read.csv(annotation_2_to_1_match_file, quote="")
        
        # mean expression
        mean_exp_file <- sprintf("%s/%s_%s_gene_ID_matching_density_data.csv", consistency_dir, dataset, comparison)
        mean_exp_df <- read.csv(mean_exp_file, quote="")
        
        # exon length
        exon_length_file <- sprintf("%s/consistency_stringency_gene_stats_plots/%s_%s_consistency_stringency_gene_stats_plots_thres_%s.csv", consistency_dir, dataset, comparison, threshold_in_path)
        exon_length_df <- read.csv(exon_length_file, quote="")
        
        # Entrez IDs
        ## annotation 1 Entrez IDs of genes used for GO analysis and DE genes
        annotation_1_GO_analysis_Entrez_IDs <- annotation_1_entrez_ID_df$ENTREZID
        annotation_1_DE_genes_Entrez <- annotation_1_entrez_ID_df[which(annotation_1_entrez_ID_df$SYMBOL %in% annotation_1_DE_genes), "ENTREZID"]
        annotation_1_expressed_genes_Entrez <- annotation_1_entrez_ID_df[which(annotation_1_entrez_ID_df$SYMBOL %in% rownames(annotation_1_de_results_df)), "ENTREZID"]
        ## annotation 2 Entrez IDs of genes used for GO analysis and DE genes
        annotation_2_GO_analysis_Entrez_IDs <- annotation_2_entrez_ID_df$entrezgene_id
        annotation_2_DE_genes_Entrez <- annotation_2_entrez_ID_df[which(annotation_2_entrez_ID_df$ensembl_gene_id %in% annotation_2_DE_genes), "entrezgene_id"]
        annotation_2_expressed_genes_Entrez <- annotation_2_entrez_ID_df[which(annotation_2_entrez_ID_df$ensembl_gene_id %in% rownames(annotation_2_de_results_df)), "entrezgene_id"]
        
        # IDs
        ## IDs used in exp_matrix_1
        annotation_1_id_colname <- "SYMBOL"
        # ## IDs used in annotation_entrez_ID_df to match IDs used in exp_matrix_2 (not used)
        # annotation_1_matching_colname_1 <- "ensembl_gene_id"
        ## IDs used in annotation_entrez_ID_df and used in exp_matrix_2 to match IDs used in exp_matrix_2
        annotation_1_matching_colname_1 <- "ensembl_gene_id"
        ## IDs used in exp_matrix_2
        annotation_2_id_colname <- "ensembl_gene_id"
        # ## IDs used in annotation_entrez_ID_df to match IDs used in exp_matrix_1
        # annotation_2_matching_colname_1 <- "entrezgene_accession"
        ## IDs used in annotation_entrez_ID_df and used in exp_matrix_1 to match IDs used in exp_matrix_1
        annotation_2_matching_colname_1 <- "SYMBOL"
        
        # annotation 1 specific GO terms (RefSeq)
        if (file.exists(annotation_1_go_specific_results_file)) {
          print(sprintf("annotation: %s", annotation_1_name))
          specific_GO_output_name <- sprintf("%s_%s_%s_specific_matching_features", comparison_FC_output_name, dysregulation, gsub(" ", "", annotation_1_name))
          # specific_GO_term_matching_features_outputs(annotation_1_go_specific_results_df, annotation_2_go_results_df, annotation_1_DE_genes_Entrez, annotation_2_DE_genes_Entrez, annotation_1_expressed_genes_Entrez, annotation_2_expressed_genes_Entrez, annotation_entrez_ID_df, annotation_1_id_colname, annotation_2_id_colname, annotation_1_matching_colname_1, annotation_1_matching_colname_2, annotation_2_matching_colname_1, annotation_2_matching_colname_2, annotation_1_de_results_df, annotation_2_de_results_df, exon_length_df, mean_exp_df, annotation_1_name, annotation_2_name, dysregulation, specific_GO_term_output_dir, specific_GO_output_name)
          annotation_1_go_specific_matching_features_df <- specific_GO_term_matching_features_outputs(annotation_1_go_specific_results_df, annotation_2_go_results_df, annotation_1_DE_genes_Entrez, annotation_2_DE_genes_Entrez, annotation_1_expressed_genes_Entrez, annotation_2_expressed_genes_Entrez, annotation_entrez_ID_df, annotation_1_id_colname, annotation_2_id_colname, annotation_1_matching_colname_1, annotation_2_matching_colname_1, annotation_1_de_results_df, annotation_2_de_results_df, exon_length_df, mean_exp_df, annotation_1_name, annotation_2_name, dysregulation, specific_GO_term_output_dir, specific_GO_output_name)
          # all matching features plots
          if (dim(annotation_1_go_specific_matching_features_df[which(! is.na(annotation_1_go_specific_matching_features_df$DE_pval_matching)),])[1] > 0) {
            scatter_plot_list[[sprintf("%s_%s", dysregulation, gsub(" ", "", annotation_1_name))]] <- specific_GO_term_matching_features_outputs_all_plots(annotation_1_go_specific_matching_features_df, specific_GO_term_output_dir, sprintf("%s_all_plots", specific_GO_output_name))
          } else {
            empty_plot <- ggplot() + theme_void() + geom_text(aes(0, 0, label="no specific GO terms"))
            scatter_plot_list[[sprintf("%s_%s", dysregulation, gsub(" ", "", annotation_1_name))]] <- list(scatter_plot=empty_plot, scatter_plot_marginal=empty_plot, nb_categories=0, nb_annotations=0)
          }
        }
        # annotation 2 specific GO terms (Ensembl)
        if (file.exists(annotation_2_go_specific_results_file)) {
          print(sprintf("annotation: %s", annotation_2_name))
          specific_GO_output_name <- sprintf("%s_%s_%s_specific_matching_features", comparison_FC_output_name, dysregulation, gsub(" ", "", annotation_2_name))
          # specific_GO_term_matching_features_outputs(annotation_2_go_specific_results_df, annotation_1_go_results_df, annotation_2_DE_genes_Entrez, annotation_1_DE_genes_Entrez, annotation_2_expressed_genes_Entrez, annotation_1_expressed_genes_Entrez, annotation_entrez_ID_df, annotation_2_id_colname, annotation_1_id_colname, annotation_2_matching_colname_1, annotation_2_matching_colname_2, annotation_1_matching_colname_1, annotation_1_matching_colname_2, annotation_2_de_results_df, annotation_1_de_results_df, exon_length_df, mean_exp_df, annotation_2_name, annotation_1_name, dysregulation, specific_GO_term_output_dir, specific_GO_output_name)
          annotation_2_go_specific_matching_features_df <- specific_GO_term_matching_features_outputs(annotation_2_go_specific_results_df, annotation_1_go_results_df, annotation_2_DE_genes_Entrez, annotation_1_DE_genes_Entrez, annotation_2_expressed_genes_Entrez, annotation_1_expressed_genes_Entrez, annotation_entrez_ID_df, annotation_2_id_colname, annotation_1_id_colname, annotation_2_matching_colname_1, annotation_1_matching_colname_1, annotation_2_de_results_df, annotation_1_de_results_df, exon_length_df, mean_exp_df, annotation_2_name, annotation_1_name, dysregulation, specific_GO_term_output_dir, specific_GO_output_name)
          # all matching features plots
          if (dim(annotation_2_go_specific_matching_features_df[which(! is.na(annotation_2_go_specific_matching_features_df$DE_pval_matching)),])[1] > 0) {
            scatter_plot_list[[sprintf("%s_%s", dysregulation, gsub(" ", "", annotation_2_name))]] <- specific_GO_term_matching_features_outputs_all_plots(annotation_2_go_specific_matching_features_df, specific_GO_term_output_dir, sprintf("%s_all_plots", specific_GO_output_name))
          } else {
            empty_plot <- ggplot() + theme_void() + geom_text(aes(0, 0, label="no specific GO terms"))
            scatter_plot_list[[sprintf("%s_%s", dysregulation, gsub(" ", "", annotation_2_name))]] <- list(scatter_plot=empty_plot, scatter_plot_marginal=empty_plot, nb_categories=0, nb_annotations=0)
          }
        }
      }
    }
    
    # figure
    ## get legend: plots with the more consistency/stringency categories
    more_categories_plot <- NA
    nb_categories <- 0
    for (one_name in names(scatter_plot_list)) {
      if (scatter_plot_list[[one_name]]$nb_categories > nb_categories) {
        more_categories_plot <- scatter_plot_list[[one_name]]$scatter_plot
        nb_categories <- scatter_plot_list[[one_name]]$nb_categories
      }
    }
    figure_legend <- get_legend(more_categories_plot)
    
    ## draw figure
    pdf(sprintf("%s/%s.pdf", specific_GO_term_output_dir, sprintf("%s_specific_figure", comparison_FC_output_name)), width=9)
    multiplot <- ggdraw() +
      draw_plot(scatter_plot_list[[sprintf("Up_%s",  gsub(" ", "", annotation_1_name))]]$scatter_plot_marginal, 0, 0.55, 0.5, 0.45) +
      draw_plot(scatter_plot_list[[sprintf("Up_%s",  gsub(" ", "", annotation_2_name))]]$scatter_plot_marginal, 0.5, 0.55, 0.5, 0.45) +
      draw_plot(scatter_plot_list[[sprintf("Down_%s",  gsub(" ", "", annotation_1_name))]]$scatter_plot_marginal, 0, 0.1, 0.5, 0.45) +
      draw_plot(scatter_plot_list[[sprintf("Down_%s",  gsub(" ", "", annotation_2_name))]]$scatter_plot_marginal, 0.5, 0.1, 0.5, 0.45) +
      draw_plot(figure_legend, 0, 0, 1, 0.1) +
      draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.55, 0.55), size=15)
    print(multiplot)
    dev.off()
  }
}


