# Functions for reference genomes evaluation

expression_matrix_common_samples_expressed_features <- function(exp_mat, samp) {
  exp_mat_common_samples <- exp_mat[, samp]
  ### remove '__alignment_not_unique', '__ambiguous', '__no_feature', '__not_aligned' and '__too_low_aQual'
  exp_mat_common_samples <- exp_mat_common_samples[rownames(exp_mat_common_samples)[! grepl("^_", rownames(exp_mat_common_samples))],]
  ### remove non expressed features
  expressed_features <- rownames(exp_mat_common_samples)[unname(apply(exp_mat_common_samples, 1, sum) > 0)]
  exp_mat_expressed_features <- exp_mat_common_samples[expressed_features,]
  return(exp_mat_expressed_features)
}

density_percentiles_plot <- function(data_df, data_name, category_name, plot_title) {
  # percentiles
  df2return_colnames <- c()
  df2return <- data.frame(percentile=seq(0,1,0.05))
  for (one_level in levels(data_df[, category_name])) {
    percentiles <- quantile(as.data.frame(data_df)[which(data_df[, category_name]==one_level), data_name], probs=seq(0,1,0.05))
    df2return <- cbind(df2return, percentiles)
    df2return_colnames <- c(df2return_colnames, sprintf("%s_%s", data_name, one_level))
  }
  df2return$percentile <- names(percentiles)
  colnames(df2return)[2:length(colnames(df2return))] <- df2return_colnames
  
  ## for plot
  median_df2ggplot <- data_df %>% group_by(!!sym(category_name)) %>% summarize(median=median(!!sym(data_name)))
  lowerquartile_df2ggplot <- data_df %>% group_by(!!sym(category_name)) %>% summarize(quartile=unname(quantile(!!sym(data_name), probs=seq(0,1,0.25))["25%"]))
  upperquartile_df2ggplot <- data_df %>% group_by(!!sym(category_name)) %>% summarize(quartile=unname(quantile(!!sym(data_name), probs=seq(0,1,0.25))["75%"]))
  firstdecile_df2ggplot <- data_df %>% group_by(!!sym(category_name)) %>% summarize(decile=unname(quantile(!!sym(data_name), probs=seq(0,1,0.1))["10%"]))
  lastdecile_df2ggplot <- data_df %>% group_by(!!sym(category_name)) %>% summarize(decile=unname(quantile(!!sym(data_name), probs=seq(0,1,0.1))["90%"]))
  
  # plot theme
  plot_theme <- theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8))
  
  p <- ggplot(data_df, aes_string(x=data_name, color=category_name)) +
    geom_line(stat="density") +
    geom_vline(data=median_df2ggplot, aes_string(xintercept="median", color=category_name), linetype="longdash", alpha=0.7) +
    geom_vline(data=lowerquartile_df2ggplot, aes_string(xintercept="quartile", color=category_name), linetype="dashed", alpha=0.7) +
    geom_vline(data=upperquartile_df2ggplot, aes_string(xintercept="quartile", color=category_name), linetype="dashed", alpha=0.7) +
    geom_vline(data=firstdecile_df2ggplot, aes_string(xintercept="decile", color=category_name), linetype="dotted", alpha=0.7) +
    geom_vline(data=lastdecile_df2ggplot, aes_string(xintercept="decile", color=category_name), linetype="dotted", alpha=0.7) +
    labs(title=plot_title, x="Value", y="Density", color="Category") +
    plot_theme
  print(p)
  
  return(df2return)
}

consistency_category_per_DE_specificity_plot <- function(data_df, annotation_1_to_2_match_with_nonDE_df, annotation_2_to_1_match_with_nonDE_df, annotation_1, annotation_2, category_colname, out_dir, out_name) {
  
  # consistency/stringency category percentage per DE specificity
  ## vectors for consistency/stringency category percentage barplot data frame
  specificity_vector <- consistency_stringency_category_vector <- value_vector <- pct_vector <- c()
  df2write <- data.frame()
  ## common DE genes
  gene_category <- "common DE genes"
  common_DE_genes <- as.character(annotation_1_to_2_match_with_nonDE_df[which(annotation_1_to_2_match_with_nonDE_df$category_2=="1:1_1"), "feature"])
  sub_data_df <- data_df[which(data_df[[sprintf("gene_ID_%s", gsub(" ", "", annotation_1))]] %in% common_DE_genes),]
  df2write <- rbind(df2write, sub_data_df)
  ### nb genes per consistency/stringency categories
  nb_genes <- unname(table(sub_data_df[[category_colname]])[consistency_stringency_categories])
  pct <- nb_genes/sum(nb_genes)*100
  ### vectors for consistency/stringency category percentage per DE specificty barplot data frame
  specificity_vector <- c(specificity_vector, rep(gene_category, length(consistency_stringency_categories)))
  consistency_stringency_category_vector <- c(consistency_stringency_category_vector, consistency_stringency_categories)
  value_vector <- c(value_vector, nb_genes)
  pct_vector <- c(pct_vector, pct)
  ## specific DE genes
  ### exp matrix 1
  gene_category <- sprintf("%s specific DE genes", annotation_1)
  specific_DE_genes_1 <- as.character(annotation_1_to_2_match_with_nonDE_df[which(annotation_1_to_2_match_with_nonDE_df$category_2=="1:1_0_1"), "feature"])
  sub_data_df <- data_df[which(data_df[[sprintf("gene_ID_%s", gsub(" ", "", annotation_1))]] %in% specific_DE_genes_1),]
  df2write <- rbind(df2write, sub_data_df)
  #### nb genes per consistency/stringency categories
  nb_genes <- unname(table(sub_data_df[[category_colname]])[consistency_stringency_categories])
  pct <- nb_genes/sum(nb_genes)*100
  #### vectors for consistency/stringency category percentage per DE specificty barplot data frame
  specificity_vector <- c(specificity_vector, rep(gene_category, length(consistency_stringency_categories)))
  consistency_stringency_category_vector <- c(consistency_stringency_category_vector, consistency_stringency_categories)
  value_vector <- c(value_vector, nb_genes)
  pct_vector <- c(pct_vector, pct)
  ### exp matrix 2
  gene_category <- sprintf("%s specific DE genes", annotation_2)
  specific_DE_genes_2 <- as.character(annotation_2_to_1_match_with_nonDE_df[which(annotation_2_to_1_match_with_nonDE_df$category_2=="1:1_0_1"), "feature"])
  sub_data_df <- data_df[which(data_df[[sprintf("gene_ID_%s", gsub(" ", "", annotation_2))]] %in% specific_DE_genes_2),]
  df2write <- rbind(df2write, sub_data_df)
  #### nb genes per consistency/stringency categories
  nb_genes <- unname(table(sub_data_df[[category_colname]])[consistency_stringency_categories])
  pct <- nb_genes/sum(nb_genes)*100
  #### vectors for consistency/stringency category percentage per DE specificty barplot data frame
  specificity_vector <- c(specificity_vector, rep(gene_category, length(consistency_stringency_categories)))
  consistency_stringency_category_vector <- c(consistency_stringency_category_vector, consistency_stringency_categories)
  value_vector <- c(value_vector, nb_genes)
  pct_vector <- c(pct_vector, pct)
  
  write.csv(df2write, file=sprintf("%s/%s_consistency_stringency.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  ## consistency/stringency category percentage per DE specificty barplot data frame
  category_pct_df2ggplot <- data.frame(specificity=specificity_vector, category=consistency_stringency_category_vector, value=value_vector, pct=pct_vector)
  category_pct_df2ggplot$category <- factor(category_pct_df2ggplot$category, levels=rev(c("consistent_stringent", "consistent_nonstringent", "inconsistent_stringent", "inconsistent_nonstringent")))
  levels(category_pct_df2ggplot$category) <- unlist(lapply(levels(category_pct_df2ggplot$category), function(x) { return(sub("_", " ", x)) }))
  levels(category_pct_df2ggplot$category) <- unlist(lapply(levels(category_pct_df2ggplot$category), function(x) { return(sub("nonstringent", "non stringent", x)) }))
  category_pct_df2ggplot$specificity <- factor(category_pct_df2ggplot$specificity, levels=rev(c("common DE genes", sprintf("%s specific DE genes", exp_matrix_1_name), sprintf("%s specific DE genes", exp_matrix_2_name))))
  levels(category_pct_df2ggplot$specificity) <- unlist(lapply(levels(category_pct_df2ggplot$specificity), function(x) { return(sub(" DE genes", "", x)) }))
  levels(category_pct_df2ggplot$specificity) <- unlist(lapply(levels(category_pct_df2ggplot$specificity), function(x) { return(sub("specific", "", x)) }))
  levels(category_pct_df2ggplot$specificity) <- unlist(lapply(levels(category_pct_df2ggplot$specificity), function(x) { return(sub("NCBI ", "", x)) }))
  # write.csv(category_pct_df2ggplot, file=sprintf("%s/%s_consistency_stringency_category_pct.csv", figures_fc_output_dir, comparison_output_basename), quote=FALSE, row.names=FALSE)
  pairwise_category_pct_plot <- ggplot(data=category_pct_df2ggplot, aes(x=specificity, y=pct, fill=category)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=value), size=2.5, position=position_stack(vjust=0.5)) +
    coord_flip() +
    labs(title="Consistency/stringency categories per DE specificty", x="Specificity", y="Percentage", fill="Category") +
    scale_fill_manual(values=rev(unname(consistency_stringency_plot_colors))) +
    # scale_fill_manual(values=c(brewer.pal(9, "Set1")[2], brewer.pal(9, "Set1")[4], brewer.pal(9, "Set1")[1])) +
    plot_theme +
    # theme(legend.position="bottom", legend.box="horizontal") +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(reverse=TRUE))
  
  write.csv(category_pct_df2ggplot, file=sprintf("%s/%s_consistency_stringency_DE_specificity_df2ggplot.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  return(pairwise_category_pct_plot)
}

exon_mergelength_density_plot <- function(data_df, annotation_1_to_2_match_with_nonDE_df, annotation_2_to_1_match_with_nonDE_df, annotation_1, annotation_2, category_colname, log2ratio_colname, out_dir, out_name) {
  
  ## data frame for inconsistent gene exon mergelenth density plots
  exon_mergelength_df2ggplot <- data.frame()
  ### common DE genes
  gene_category <- "common DE genes"
  common_DE_genes <- as.character(annotation_1_to_2_match_with_nonDE_df[which(annotation_1_to_2_match_with_nonDE_df$category_2=="1:1_1"), "feature"])
  sub_data_df <- data_df[which(data_df[[sprintf("gene_ID_%s", gsub(" ", "", annotation_1))]] %in% common_DE_genes),]
  sub_data_df <- cbind(sub_data_df, specificity=rep(gene_category, dim(sub_data_df)[1]))
  df2write <- rbind(df2write, sub_data_df)
  exon_mergelength_df2ggplot <- rbind(exon_mergelength_df2ggplot, sub_data_df)
  ### specific DE genes
  #### exp matrix 1
  gene_category <- sprintf("%s specific DE genes", annotation_1)
  specific_DE_genes_1 <- as.character(annotation_1_to_2_match_with_nonDE_df[which(annotation_1_to_2_match_with_nonDE_df$category_2=="1:1_0_1"), "feature"])
  sub_data_df <- sub_data_df[which(sub_data_df[[sprintf("gene_ID_%s", gsub(" ", "", annotation_1))]] %in% specific_DE_genes_1),]
  sub_data_df <- cbind(sub_data_df, specificity=rep(gene_category, dim(sub_data_df)[1]))
  df2write <- rbind(df2write, sub_data_df)
  exon_mergelength_df2ggplot <- rbind(exon_mergelength_df2ggplot, sub_data_df)
  #### exp matrix 2
  gene_category <- sprintf("%s specific DE genes", annotation_2)
  specific_DE_genes_2 <- as.character(annotation_2_to_1_match_with_nonDE_df[which(annotation_2_to_1_match_with_nonDE_df$category_2=="1:1_0_1"), "feature"])
  sub_data_df <- sub_data_df[which(sub_data_df[[sprintf("gene_ID_%s", gsub(" ", "", annotation_2))]] %in% specific_DE_genes_2),]
  sub_data_df <- cbind(sub_data_df, specificity=rep(gene_category, dim(sub_data_df)[1]))
  df2write <- rbind(df2write, sub_data_df)
  exon_mergelength_df2ggplot <- rbind(exon_mergelength_df2ggplot, sub_data_df)
  
  write.csv(df2write, file=sprintf("%s/%s_exon_mergelength.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  ## only 'inconsistent_stringent' and 'inconsistent_nonstringent' categories
  inconsistent_exon_mergelength_df <- exon_mergelength_df2ggplot[which(exon_mergelength_df2ggplot[[category_colname]] %in% c("inconsistent_stringent", "inconsistent_nonstringent")),]
  inconsistent_exon_mergelength_df <- droplevels(inconsistent_exon_mergelength_df)
  
  return(inconsistent_exon_mergelength_df)
}




add_pvalues2exon_mergelength_density_plot <- function(data_df, log2ratio_colname, specificity_colname, category_colname, plot_colors, annotation_1, annotation_2) {
  
  inconsistent_exon_mergelength_plot2 <- ggplot(data_df, aes_string(x=log2ratio_colname, color=specificity_colname, linetype=category_colname)) +
    geom_line(stat="density") +
    scale_color_manual(values=c(plot_colors[1], plot_colors[2], plot_colors[3])) +
    geom_vline(xintercept=0, linetype="dashed", alpha=0.5) +
    labs(title=sprintf("Exon length ratio per DE specificity", annotation_1, annotation_2), x="exon length log2 ratio", y="Density", color="Specificity", linetype="Category") +
    plot_theme +
    theme(legend.position="bottom", legend.box="vertical")
  # guides(color=guide_legend(reverse=TRUE))
  
  
  #### add p-values on plot
  specificity_vector <- category_vector <- nb_genes_vector <- test_vector <- pval_vector <- c()
  for (one_specificity in levels(data_df[[specificity_colname]])) {
    one_specificity_data_df <- data_df[which(data_df[[specificity_colname]]==one_specificity),]
    for (one_category in levels(one_specificity_data_df[[category_colname]])) {
      one_category_data_df <- one_specificity_data_df[which(one_specificity_data_df[[category_colname]]==one_category),]
      nb_genes <- dim(one_category_data_df)[1]
      # TODO: test if sample size is greater than or equal to 3 to avoid shapiro.test() function error
      if (nb_genes >= 3) {
        differences <- one_category_data_df[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1))]] - one_category_data_df[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2))]]
        #### the Shapiro-Wilk test function is restricted to sample sizes up to 5000
        normality_test_pval <- ifelse(nb_genes > 5000,  ad.test(differences)$p.value, shapiro.test(differences)$p.value)
        test_alternative <- ifelse(mean(one_category_data_df$exon_mergelength_log2ratio) < 0, "less", "greater")
        wilcoxon_test_result <- wilcox.test(one_category_data_df[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1))]], one_category_data_df[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2))]], paired=TRUE, alternative=test_alternative)
        t_test_result <- t.test(one_category_data_df[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_1))]], one_category_data_df[[sprintf("exon_mergelength_%s", gsub(" ", "", annotation_2))]], paired=TRUE, alternative=test_alternative)
        pval_vector <- c(pval_vector, c(normality_test_pval, wilcoxon_test_result$p.value, t_test_result$p.value))
      } else {
        pval_vector <- c(pval_vector, c(NA, NA, NA))
      }
      test_vector <- c(test_vector, c("normality", "Wilcoxon", "t_test"))
      nb_genes_vector <- c(nb_genes_vector, rep(nb_genes, 3))
      category_vector <- c(category_vector, rep(one_category, 3))
    }
    specificity_vector <- c(specificity_vector, rep(one_specificity, length(levels(data_df$consistency_stringency_category))*3))
  }
  data_test_df <- data.frame(specificity=specificity_vector, category=category_vector, nb_genes=nb_genes_vector, test=test_vector, pval=pval_vector)
  ##### get plot ranges
  x_scale <- layer_scales(inconsistent_exon_mergelength_plot2)$x$range$range
  x_min <- x_scale[1]
  x_max <- x_scale[2]
  x_range <- x_max - x_min
  y_scale <- layer_scales(inconsistent_exon_mergelength_plot2)$y$range$range
  y_max <- y_scale[2]
  y_range <- y_max - y_scale[1]
  
  specificity_vector <- category_vector <- x_vector <- y_vector <- label_vector <- c()
  for (i in 1:length(levels(data_test_df$specificity))) {
    one_specificity <- levels(data_test_df$specificity)[i]
    
    one_category <- "inconsistent stringent"
    pval <- data_test_df[which(data_test_df$specificity==one_specificity & data_test_df$category==one_category & data_test_df$test=="Wilcoxon"), "pval"]
    label <- sprintf("p = %.2e", pval)
    label_vector <- c(label_vector, label)
    x_vector <- c(x_vector, x_min+0.1*x_range)
    y_vector <- c(y_vector, y_max-0.05*y_range-((i-1)*0.1*y_range))
    category_vector <- c(category_vector, one_category)
    specificity_vector <- c(specificity_vector, one_specificity)
    
    one_category <- "inconsistent non stringent"
    pval <- data_test_df[which(data_test_df$specificity==one_specificity & data_test_df$category==one_category & data_test_df$test=="Wilcoxon"), "pval"]
    label <- sprintf("p = %.2e", pval)
    label_vector <- c(label_vector, label)
    x_vector <- c(x_vector, x_max-0.1*x_range)
    y_vector <- c(y_vector, y_max-0.05*y_range-((i-1)*0.1*y_range))
    category_vector <- c(category_vector, one_category)
    specificity_vector <- c(specificity_vector, one_specificity)
  }
  annotation_df <- data.frame(specificity=specificity_vector, consistency_stringency_category=category_vector, x=x_vector, y=y_vector, label=label_vector)
  inconsistent_exon_mergelength_plot2 <- inconsistent_exon_mergelength_plot2 +
    geom_text(data=annotation_df, aes(x=x, y=y, label=label, color=specificity), size=2.5, show.legend=FALSE)
  
  return(inconsistent_exon_mergelength_plot2)
}













#####################################
# Ensembl and NCBI gene ID matching #
#####################################

remove_not_requested_getBM_ids <- function(getBM_df, gene_id_col, gene_ids, out_dir, out_basename) {
  keep_ids <- as.character(getBM_df[[gene_id_col]]) %in% gene_ids
  keep_getBM_df <- getBM_df[keep_ids,]
  rejected_getBM_df <- getBM_df[! keep_ids,]
  write.csv(rejected_getBM_df, file=sprintf("%s/%s_not_requested_ids.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
  return(keep_getBM_df)
}

Ensembl_NCBI_ID_matching <- function(match_df, ID1_col, ID2_col, genes_ID1_interest, genes_ID2_interest, out_dir, out_basename) {
  
  # vectors for data frame to return
  feature_vector <- matching_vector <- category_vector <- category_vector_2 <- c()
  
  # gene IDs with a matching gene ID
  match_df_matching <- match_df[which(match_df[[ID2_col]] != ""),]
  ## gene IDs with duplicated matching gene IDs
  dup_matching_IDs <- unique(match_df_matching[[ID2_col]][duplicated(match_df_matching[[ID2_col]])])
  match_df_dup_matching <- match_df_matching[which(match_df_matching[[ID2_col]] %in% dup_matching_IDs),]
  gene_ids_with_dup_matching <- unique(match_df_dup_matching[[ID1_col]])
  for (one_gene_id_with_dup_matching in gene_ids_with_dup_matching) {
    dup_matching_gene_ids <- unique(match_df_matching[which(match_df_matching[[ID1_col]] == one_gene_id_with_dup_matching), ID2_col])
    number_dup_matching_gene_ids <- length(dup_matching_gene_ids)
    if (number_dup_matching_gene_ids > 1) {
      # N:P
      dup_matching_gene_ids_interest <- length(dup_matching_gene_ids[dup_matching_gene_ids %in% genes_ID2_interest])
      if (dup_matching_gene_ids_interest == 0) {
        # N:P_0
        one_category <- "multiple matches shared by >1 genes: no gene of interest"
        one_category_2 <- "N:P_0"
      } else {
        if (dup_matching_gene_ids_interest == 1) {
          # N:P_1
          one_category <- "multiple matches shared by >1 genes: 1 gene of interest"
          one_category_2 <- "N:P_1"
        } else {
          if (dup_matching_gene_ids_interest > 1) {
            # N:P_Q
            one_category <- "multiple matches shared by >1 genes: >1 genes of interest"
            one_category_2 <- "N:P_Q"
          }
        }
      }
      feature_vector <- c(feature_vector, rep(one_gene_id_with_dup_matching, number_dup_matching_gene_ids))
      matching_vector <- c(matching_vector, dup_matching_gene_ids)
      category_vector <- c(category_vector, rep(one_category, number_dup_matching_gene_ids))
      category_vector_2 <- c(category_vector_2, rep(one_category_2, number_dup_matching_gene_ids))
    } else {
      if (number_dup_matching_gene_ids == 1) {
        # N:1
        one_dup_matching_id <- match_df_matching[which(match_df_matching[[ID1_col]] == one_gene_id_with_dup_matching), ID2_col]
        if (one_dup_matching_id %in% genes_ID2_interest) {
          # N:1_1
          one_category <- "unique match in genes of interest shared by >1 genes"
          one_category_2 <- "N:1_1"
        } else {
          # N:1_0
          one_category <- "unique match not in genes of interest shared by >1 genes"
          one_category_2 <- "N:1_0"
        }
      }
      feature_vector <- c(feature_vector, one_gene_id_with_dup_matching)
      matching_vector <- c(matching_vector, one_dup_matching_id)
      category_vector <- c(category_vector, one_category)
      category_vector_2 <- c(category_vector_2, one_category_2)
    }
    
  }
  write.csv(match_df_dup_matching, file=sprintf("%s/%s_N_X_matching.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)

  ### remove gene IDs with duplicated matching gene IDs
  match_df_matching_no_dup_matching_IDs <- match_df_matching[which(! match_df_matching[[ID1_col]] %in% gene_ids_with_dup_matching),]
  
  ## gene IDs with multiple matching gene IDs
  dup_gene_ids <- unique(match_df_matching_no_dup_matching_IDs[[ID1_col]][duplicated(match_df_matching_no_dup_matching_IDs[[ID1_col]])])
  match_df_matching_dup <- match_df_matching_no_dup_matching_IDs[which(match_df_matching_no_dup_matching_IDs[[ID1_col]] %in% dup_gene_ids),]
  write.csv(match_df_matching_dup, file=sprintf("%s/%s_1_N_matching.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
  ### multiple matching gene IDs in genes_ID2_interest ?
  #### no --> 1:0
  #### 1 --> 1:1
  #### >1 --> true 1:N
  #one_dup_gene_id <- dup_gene_ids[1]
  for (one_dup_gene_id in dup_gene_ids) {
    dup_matching_gene_ids <- match_df_matching_dup[which(match_df_matching_dup[[ID1_col]] == one_dup_gene_id), ID2_col]
    dup_matching_gene_ids_interest <- length(dup_matching_gene_ids[dup_matching_gene_ids %in% genes_ID2_interest])
    if (dup_matching_gene_ids_interest == 0) {
      # 1:0
      one_category <- "multiple matches: no gene of interest"
      one_category_2 <- "1:N_0"
    } else {
      if (dup_matching_gene_ids_interest == 1) {
        # 1:1
        one_category <- "multiple matches: 1 gene of interest"
        one_category_2 <- "1:N_1"
      } else {
        if (dup_matching_gene_ids_interest > 1) {
          # 1:N
          one_category <- "multiple matches: >1 genes of interest"
          one_category_2 <- "1:N_P"
        }
      }
    }
    number_dup_matching_gene_ids <- length(dup_matching_gene_ids)
    feature_vector <- c(feature_vector, rep(one_dup_gene_id, number_dup_matching_gene_ids))
    matching_vector <- c(matching_vector, dup_matching_gene_ids)
    category_vector <- c(category_vector, rep(one_category, number_dup_matching_gene_ids))
    category_vector_2 <- c(category_vector_2, rep(one_category_2, number_dup_matching_gene_ids))
  }
  ### remove gene IDs with multiple matching gene IDs
  match_df_matching_no_dup <- match_df_matching_no_dup_matching_IDs[which(! match_df_matching_no_dup_matching_IDs[[ID1_col]] %in% dup_gene_ids),]
  
  ## gene with matching ID not in genes of interest
  match_df_matching_no_dup_matching_non_interest <- match_df_matching_no_dup[which(! match_df_matching_no_dup[[ID2_col]] %in% genes_ID2_interest),]
  write.csv(match_df_matching_no_dup_matching_non_interest, file=sprintf("%s/%s_1_1_non_interest_matching.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
  features <- unique(match_df_matching_no_dup_matching_non_interest[[ID1_col]])
  matching_features <- unique(match_df_matching_no_dup_matching_non_interest[[ID2_col]])
  one_category <- "unique match not in genes of interest"
  one_category_2 <- "1:1_0"
  feature_vector <- c(feature_vector, features)
  matching_vector <- c(matching_vector, matching_features)
  category_vector <- c(category_vector, rep(one_category, length(features)))
  category_vector_2 <- c(category_vector_2, rep(one_category_2, length(features)))
  
  ## only keep gene IDs with a matching gene ID in the second expression matrix
  match_df_matching_no_dup_matching_interest <- match_df_matching_no_dup[which(match_df_matching_no_dup[[ID2_col]] %in% genes_ID2_interest),]
  write.csv(match_df_matching_no_dup_matching_interest, file=sprintf("%s/%s_1_1_interest_matching.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
  features <- unique(match_df_matching_no_dup_matching_interest[[ID1_col]])
  matching_features <- unique(match_df_matching_no_dup_matching_interest[[ID2_col]])
  one_category <- "unique match in genes of interest"
  one_category_2 <- "1:1_1"
  feature_vector <- c(feature_vector, features)
  matching_vector <- c(matching_vector, matching_features)
  category_vector <- c(category_vector, rep(one_category, length(features)))
  category_vector_2 <- c(category_vector_2, rep(one_category_2, length(features)))
  
  # gene IDs without a matching gene ID (deal with these genes at the end to prevent write.csv() function from writing 'NA' for all matching feature of df2return)
  ## gene IDs with empty matching ID in match_df data frame
  match_df_no_matching_ID <- match_df[which(match_df[[ID2_col]] == ""),]
  write.csv(match_df_no_matching_ID, file=sprintf("%s/%s_1_0_matching.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
  features <- match_df_no_matching_ID[[ID1_col]]
  ## gene IDs not in match_df data frame
  features <- c(features, genes_ID1_interest[! genes_ID1_interest %in% match_df[[ID1_col]]])
  one_category <- "no match"
  one_category_2 <- "1:0"
  feature_vector <- c(feature_vector, features)
  matching_vector <- c(matching_vector, rep(NA, length(features)))
  category_vector <- c(category_vector, rep(one_category, length(features)))
  category_vector_2 <- c(category_vector_2, rep(one_category_2, length(features)))
  
  
  df2return <- data.frame(feature=feature_vector, matching_feature=matching_vector, category=category_vector, category_2=category_vector_2)
  return(df2return)
}

gene_ID_matching_add_second_match <- function(ID_matching_categories_df, categories_of_interest, matching_genes_of_interest_1, matching_genes_of_interest_2, out_dir, out_basename) {
  
  # matching with genes of first interest and genes of second interest
  ID_matching_categories_not_in_cat_interest_df <- ID_matching_categories_df[which(ID_matching_categories_df$category_2 %in% categories_of_interest),]
  in_matching_1 <- as.character(ID_matching_categories_not_in_cat_interest_df$matching_feature) %in% matching_genes_of_interest_1
  ID_matching_categories_not_in_cat_interest_df <- cbind(ID_matching_categories_not_in_cat_interest_df, in_matching_1)
  in_matching_2 <- as.character(ID_matching_categories_not_in_cat_interest_df$matching_feature) %in% matching_genes_of_interest_2
  ID_matching_categories_not_in_cat_interest_df <- cbind(ID_matching_categories_not_in_cat_interest_df, in_matching_2)
  
  # new categories according to matching with genes of first interest and genes of second interest
  category_2_new_vector <- c()
  for (one_feature in ID_matching_categories_not_in_cat_interest_df$feature) {
    one_feature_df <- ID_matching_categories_not_in_cat_interest_df[which(ID_matching_categories_not_in_cat_interest_df$feature==one_feature),]
    one_feature_category_2 <- unique(one_feature_df$category_2)
    if (one_feature_category_2 == "1:1_0" | one_feature_category_2 == "N:1_0") {
      one_feature_category_2_new <- ifelse(one_feature_df$in_matching_2, sprintf("%s_1", one_feature_category_2), sprintf("%s_0", one_feature_category_2))
    } else {
      nb_matching_2 <- sum(one_feature_df$in_matching_2)
      if (nb_matching_2 > 1) {
        if (one_feature_category_2 == "1:N_0" | one_feature_category_2 == "1:N_1") {
          matching_letter <- "P"
        } else {
          if (one_feature_category_2 == "1:N_P" | one_feature_category_2 == "N:P_0" | one_feature_category_2 == "N:P_1") {
            matching_letter <- "Q"
          } else {
            if (one_feature_category_2 == "N:P_Q") {
              matching_letter <- "R"
            }
          }
        }
        one_feature_category_2_new <- sprintf("%s_%s", one_feature_category_2, matching_letter)
      } else {
        # nb_matching_2={0;1}
        one_feature_category_2_new <- sprintf("%s_%d", one_feature_category_2, nb_matching_2)
      }
      # one_feature_category_2_new <- ifelse(nb_matching_2>1, sprintf("%s_N", one_feature_category_2), sprintf("%s_%d", one_feature_category_2, nb_matching_2))
    }
    category_2_new_vector <- c(category_2_new_vector, one_feature_category_2_new)
  }
  ID_matching_categories_not_in_cat_interest_df2write <- cbind(ID_matching_categories_not_in_cat_interest_df, category_2_nonDE=category_2_new_vector)
  write.csv(ID_matching_categories_not_in_cat_interest_df2write, file=sprintf("%s/%s_matching_categories_not_in_DE.csv", out_dir, out_basename), quote=FALSE, row.names=FALSE)
  
  # update category_2 column
  ID_matching_categories_not_in_cat_interest_df$category_2 <- as.factor(category_2_new_vector)
  ID_matching_categories_not_in_cat_interest_df <- ID_matching_categories_not_in_cat_interest_df[, ! colnames(ID_matching_categories_not_in_cat_interest_df) %in% c("in_matching_1", "in_matching_2")]
  # update matching data frame
  ID_matching_categories_with_second_match_df <- rbind(ID_matching_categories_df[which(! ID_matching_categories_df$category_2 %in% categories_of_interest),], ID_matching_categories_not_in_cat_interest_df)
  ID_matching_categories_with_second_match_df <- droplevels(ID_matching_categories_with_second_match_df)
  return(ID_matching_categories_with_second_match_df)
}


## counts per category
gene_counts_per_ID_matching_category <- function(ID_matching_categories_df, nb_total_ids) {
  # categories_2 <- ifelse(length(ID_matching_categories_df$category_2) == 1, ID_matching_categories_df$category_2, levels(ID_matching_categories_df$category_2))
  categories_2 <- levels(ID_matching_categories_df$category_2)
  categories <- value_vector <- annotation_vector <- c()
  
  annotation <- unique(as.character(ID_matching_categories_df$annotation))
  if (length(annotation) > 1) {
    print(sprintf("Warning ! More than 1 unique category_2 for category: %s", one_category))
    annotation <- "Warning: more than 1 unique annotation"
  }
  
  for (one_category_2 in categories_2) {
    ### category
    unique_category <- unique(as.character(ID_matching_categories_df[which(ID_matching_categories_df$category_2==one_category_2),"category"]))
    if (length(unique_category) == 1) {
      categories <- c(categories, unique_category)
    } else {
      print(sprintf("Warning ! More than 1 unique category for category_2: %s", one_category_2))
    }
    ### counts
    value_vector <- c(value_vector, length(unique(ID_matching_categories_df[which(ID_matching_categories_df$category_2==one_category_2),"feature"])))
    ### annotation
    annotation_vector <- c(annotation_vector, annotation)
  }
  df2return <- data.frame(annotation=annotation_vector, category=categories, category_2=categories_2, value=value_vector, pct=value_vector/nb_total_ids*100)
  return(df2return)
}

count_pct_per_category <- function(df, category_colname, feature_colname) {
  value_vector <- c()
  nb_total_ids <- length(unique(df[[feature_colname]]))
  categories <- levels(df[[category_colname]])
  for (one_category in categories) {
    #### counts
    value_vector <- c(value_vector, length(unique(df[which(df[[category_colname]]==one_category), feature_colname])))
  }
  df2return <- data.frame(category=categories, count=value_vector, pct=value_vector/nb_total_ids*100)
}

## counts per category barplots and mean expression densities
gene_ID_matching_plots <- function(df2barplot, df2expdensityplot, out_dir, out_name) {
  # plot theme
  plot_theme <- theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8))
  
  pdf(sprintf("%s/%s.pdf", out_dir, out_name))
  # counts per gene ID match category barplots
  p <- ggplot(df2barplot, aes(x=annotation, y=value, fill=category)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%.2f%%", pct)), size=4, position=position_stack(vjust=0.5)) +
    labs(title="Gene ID counts per match category", x="Annotation", y="# IDs", fill="ID category") +
    plot_theme
  print(p)
  p <- ggplot(df2barplot, aes(x=annotation, y=value, fill=category_2)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%.2f%%", pct)), size=4, position=position_stack(vjust=0.5)) +
    labs(title="Gene ID counts per match category", x="Annotation", y="# IDs", fill="ID category") +
    plot_theme +
    theme(legend.title=element_blank())
  print(p)
  # mean expression per gene ID match category density plots
  p <- ggplot(df2expdensityplot, aes(x=mean_exp, color=category_2)) +
    geom_line(stat="density") +
    facet_wrap(~annotation, scales="free") +
    labs(title="Mean expression per gene ID match category", x="log2(Mean expression)", y="Density", color="match category") +
    plot_theme +
    theme(legend.position="bottom")
  print(p)
  
  if (length(levels(df2expdensityplot$category_2)[-which(levels(df2expdensityplot$category_2) %in% c("1:0", "1:1_0", "1:1_1"))])  > 0) {
    ## select some gene ID match categories if other categories than 1:0, 1:1_0 and 1:1_1
    df2ggplot <- df2expdensityplot[which(df2expdensityplot$category_2 %in% c("1:0", "1:1_0", "1:1_1")),]
    p <- ggplot(df2ggplot, aes(x=mean_exp, color=category_2)) +
      geom_line(stat="density") +
      facet_wrap(~annotation, scales="free") +
      labs(title="Mean expression per gene ID match category", x="log2(Mean expression)", y="Density", color="match category") +
      plot_theme +
      theme(legend.position="bottom")
    print(p)
    df2ggplot <- df2expdensityplot[which(df2expdensityplot$category_2 %in% c("1:N_0", "1:N_1", "1:N_P")),]
    p <- ggplot(df2ggplot, aes(x=mean_exp, color=category_2)) +
      geom_line(stat="density") +
      facet_wrap(~annotation, scales="free") +
      labs(title="Mean expression per gene ID match category", x="log2(Mean expression)", y="Density", color="match category") +
      plot_theme +
      theme(legend.position="bottom")
    print(p)
    df2ggplot <- df2expdensityplot[which(df2expdensityplot$category_2 %in% c("N:1_0", "N:1_1", "N:P_0", "N:P_1", "N:P_Q")),]
    p <- ggplot(df2ggplot, aes(x=mean_exp, color=category_2)) +
      geom_line(stat="density") +
      facet_wrap(~annotation, scales="free") +
      labs(title="Mean expression per gene ID match category", x="log2(Mean expression)", y="Density", color="match category") +
      plot_theme +
      theme(legend.position="bottom")
    print(p)
  }
  dev.off()
}

gene_ID_matching_barplots <- function(df2barplot, out_dir, out_name) {
  
  pdf(sprintf("%s/%s.pdf", out_dir, out_name))
  p <- ggplot(df2barplot, aes(x=annotation, y=value, fill=category)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%.2f%%", pct)), size=4, position=position_stack(vjust=0.5)) +
    labs(title="Gene ID counts per matching category", x="Annotation", y="# IDs", fill="ID category") +
    theme_bw() +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=12)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8))
  print(p)
  
  one_barplot <- ggplot(df2barplot, aes(x=annotation, y=value, fill=category_2)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%.2f%%", pct)), size=4, position=position_stack(vjust=0.5)) +
    labs(title="Gene ID counts per matching category", x="Annotation", y="# IDs", fill="ID category") +
    theme_bw() +
    theme(legend.title=element_blank()) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=12)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=10))
  print(one_barplot)
  dev.off()
}

common_features_in_matching_df <- function(match_df_1_to_2, match_df_2_to_1, name_1, name_2) {
  
  # common matching categories
  common_categories <- c("1:1_1", "1:N_1_0", "1:N_1_1", "1:N_1_P", "1:N_P_0", "1:N_P_1", "1:N_P_Q", "N:1_1", "N:P_1_0", "N:P_1_1", "N:P_1_Q", "N:P_Q_0", "N:P_Q_1", "N:P_Q_R")
  
  # add genome names in colnames
  colnames(match_df_1_to_2)[colnames(match_df_1_to_2)=="feature"] <- sprintf("gene_ID_%s", name_1)
  colnames(match_df_1_to_2)[colnames(match_df_1_to_2)=="matching_feature"] <- sprintf("gene_ID_%s", name_2)
  colnames(match_df_1_to_2)[colnames(match_df_1_to_2)=="category_2"] <- sprintf("category_2_%s_%s", name_1, name_2)
  colnames(match_df_2_to_1)[colnames(match_df_2_to_1)=="feature"] <- sprintf("gene_ID_%s", name_2)
  colnames(match_df_2_to_1)[colnames(match_df_2_to_1)=="matching_feature"] <- sprintf("gene_ID_%s", name_1)
  colnames(match_df_2_to_1)[colnames(match_df_2_to_1)=="category_2"] <- sprintf("category_2_%s_%s", name_2, name_1)
  
  # merge matching id data frames
  match_df_1_to_2_4merge <- match_df_1_to_2[which(match_df_1_to_2$category_2 %in% common_categories), c(sprintf("gene_ID_%s", name_1), sprintf("gene_ID_%s", name_2), sprintf("category_2_%s_%s", name_1, name_2))]
  match_df_2_to_1_4merge <- match_df_2_to_1[which(match_df_2_to_1$category_2 %in% common_categories), c(sprintf("gene_ID_%s", name_2), sprintf("gene_ID_%s", name_1), sprintf("category_2_%s_%s", name_2, name_1))]
  common_df <- merge(match_df_1_to_2_4merge, match_df_2_to_1_4merge, by=c(sprintf("gene_ID_%s", name_1), sprintf("gene_ID_%s", name_2)))
  common_df <- droplevels(common_df)
  
  return(common_df)
}

specific_features_in_matching_df <- function(match_df, name_1, name_2) {
  
  # specific matching categories
  specific_categories <- c("1:1_0_1", "1:N_0_1", "1:N_0_P", "N:1_0_1", "N:P_0_1", "N:P_0_Q")
  
  # add genome names in colnames
  colnames(match_df)[colnames(match_df)=="feature"] <- sprintf("gene_ID_%s", name_1)
  colnames(match_df)[colnames(match_df)=="matching_feature"] <- sprintf("gene_ID_%s", name_2)
  
  # select specific matching categories
  specific_df <- match_df[which(match_df$category_2 %in% specific_categories), c(sprintf("gene_ID_%s", name_1), sprintf("gene_ID_%s", name_2), "category_2")]
  specific_df <- droplevels(specific_df)
  
  return(specific_df)
}


##############################
# Consistency and stringency #
##############################

consistency_stringency_category <- function(consistency, stringency, threshold) {
  if (consistency >=threshold) {
    category <- ifelse(stringency >= threshold, "consistent_stringent", "consistent_nonstringent")
  } else {
    category <- ifelse(stringency >= threshold, "inconsistent_stringent", "inconsistent_nonstringent")
  }
  return(category)
}

## computations
consistency_stringency_computation <- function(exp_matrix_1, exp_matrix_1_name, exp_matrix_2, exp_matrix_2_name, match_df, n) {
  
  # columns of the 2 matrices in the same order
  exp_matrix_1_col <- colnames(exp_matrix_1)
  exp_matrix_2 <- exp_matrix_2[,exp_matrix_1_col]
  
  # consistency and stringency computations
  exp_matrix_1_id_vector <- exp_matrix_2_id_vector <- category_vector <- category_2_vector <- mean_vector <- sd_vector <- exp_matrix_1_mean_log_exp_vector <- exp_matrix_2_mean_log_exp_vector <- c()
  
  # exp_matrix_1_gene_IDs <- levels(match_df$feature)
  nb_matches <- dim(match_df)[1]
  for (i in 1:nb_matches) {
    one_exp_matrix_1_gene_ID <- as.character(match_df[i, "feature"])
    one_exp_matrix_2_matching_gene_ID <- as.character(match_df[i, "matching_feature"])
    if (one_exp_matrix_2_matching_gene_ID %in% rownames(exp_matrix_2)) {
      # count ratios
      count_ratios <- as.numeric(log((exp_matrix_1[one_exp_matrix_1_gene_ID,]+n) / (exp_matrix_2[one_exp_matrix_2_matching_gene_ID,]+n), 2))
      mean_vector <- c(mean_vector, mean(count_ratios))
      sd_vector <- c(sd_vector, sd(count_ratios))
      
      # mean exprssion
      exp_matrix_1_mean_log_exp <- apply(log(exp_matrix_1[one_exp_matrix_1_gene_ID,]+1, 2), 1, mean)
      exp_matrix_1_mean_log_exp_vector <- c(exp_matrix_1_mean_log_exp_vector, exp_matrix_1_mean_log_exp)
      exp_matrix_2_mean_log_exp <- apply(log(exp_matrix_2[one_exp_matrix_2_matching_gene_ID,]+1, 2), 1, mean)
      exp_matrix_2_mean_log_exp_vector <- c(exp_matrix_2_mean_log_exp_vector, exp_matrix_2_mean_log_exp)
      
      # IDs and categories
      exp_matrix_1_id_vector <- c(exp_matrix_1_id_vector, one_exp_matrix_1_gene_ID)
      exp_matrix_2_id_vector <- c(exp_matrix_2_id_vector, one_exp_matrix_2_matching_gene_ID)
      category_vector <- c(category_vector, as.character(match_df[i, "category"]))
      category_2_vector <- c(category_2_vector, as.character(match_df[i, "category_2"]))
    }
  }
  # mean_vector and sd_vector may contain null values
  ## add the minimum non-null absolute value to null values to avoid Inf values after log10 transformation
  mean_vector[mean_vector == 0] <- min(abs(mean_vector[mean_vector != 0]), na.rm=TRUE)
  sd_vector[sd_vector == 0] <- min(sd_vector[sd_vector != 0], na.rm=TRUE)
  consistency_vector <- -log(abs(mean_vector),10)
  stringency_vector <- -log(sd_vector,10)
  df2return <- data.frame(exp_matrix_1_gene_ID=exp_matrix_1_id_vector, exp_matrix_2_gene_ID=exp_matrix_2_id_vector, category=category_vector, category_2=category_2_vector, mean_log2ratio=mean_vector, sd_log2ratio=sd_vector, consistency=consistency_vector, stringency=stringency_vector, exp_matrix_1_mean_exp=exp_matrix_1_mean_log_exp_vector, exp_matrix_2_mean_exp=exp_matrix_2_mean_log_exp_vector)
  ## put expression matrices names in column names
  colnames(df2return)[colnames(df2return)=="exp_matrix_1_gene_ID"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
  colnames(df2return)[colnames(df2return)=="exp_matrix_2_gene_ID"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
  colnames(df2return)[colnames(df2return)=="exp_matrix_1_mean_exp"] <- sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name))
  colnames(df2return)[colnames(df2return)=="exp_matrix_2_mean_exp"] <- sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name))
  
  # add consistency/stringency categories
  consistency_stringency_category <- c()
  consistency_stringency_category_2 <- c()
  consistency_stringency_category_3 <- c()
  for (i in 1:dim(df2return)[1]) {
    one_consistency <- df2return[i, "consistency"]
    one_stringency <- df2return[i, "stringency"]
    ## threshold: 1
    threshold <- 1
    one_category <- consistency_stringency_category(one_consistency, one_stringency, threshold)
    consistency_stringency_category <- c(consistency_stringency_category, one_category)
    ## threshold: 5% difference of expression, i.e. consistency/stringency threshold: -log10(log2(1.1))=1.152493
    threshold <- -log(log(1.05, 2), 10)
    one_category <- consistency_stringency_category(one_consistency, one_stringency, threshold)
    consistency_stringency_category_2 <- c(consistency_stringency_category_2, one_category)
    ## threshold: 10% difference of expression, i.e. consistency/stringency threshold: -log10(log2(1.1))=0.8616862
    threshold <- -log(log(1.1, 2), 10)
    one_category <- consistency_stringency_category(one_consistency, one_stringency, threshold)
    consistency_stringency_category_3 <- c(consistency_stringency_category_3, one_category)
  }
  df2return <- cbind(df2return, consistency_stringency_category, consistency_stringency_category_2, consistency_stringency_category_3)
  df2return$consistency_stringency_category <- as.factor(consistency_stringency_category)
  df2return$consistency_stringency_category_2 <- as.factor(consistency_stringency_category_2)
  df2return$consistency_stringency_category_3 <- as.factor(consistency_stringency_category_3)
  return(df2return)
}

consistency_stringency_computation_2 <- function(exp_matrix_1, exp_matrix_1_name, exp_matrix_2, exp_matrix_2_name, match_df, n) {
  
  # columns of the 2 matrices in the same order
  exp_matrix_1_col <- colnames(exp_matrix_1)
  exp_matrix_2 <- exp_matrix_2[,exp_matrix_1_col]
  
  # consistency and stringency computations
  exp_matrix_1_id_vector <- exp_matrix_2_id_vector <- category_vector <- category_2_vector <- mean_vector <- sd_vector <- mixed_signs_vector <- exp_matrix_1_mean_log_exp_vector <- exp_matrix_2_mean_log_exp_vector <- c()
  
  # exp_matrix_1_gene_IDs <- levels(match_df$feature)
  nb_matches <- dim(match_df)[1]
  for (i in 1:nb_matches) {
    one_exp_matrix_1_gene_ID <- as.character(match_df[i, "feature"])
    one_exp_matrix_2_matching_gene_ID <- as.character(match_df[i, "matching_feature"])
    if (one_exp_matrix_2_matching_gene_ID %in% rownames(exp_matrix_2)) {
      # count ratios
      count_ratios <- as.numeric(log((exp_matrix_1[one_exp_matrix_1_gene_ID,]+n) / (exp_matrix_2[one_exp_matrix_2_matching_gene_ID,]+n), 2))
      mixed_signs_vector <- c(mixed_signs_vector, ifelse(sum(count_ratios > 0) != 0 & sum(count_ratios < 0) != 0, 1, 0))
      mean_vector <- c(mean_vector, mean(abs(count_ratios)))
      sd_vector <- c(sd_vector, sd(abs(count_ratios)))
      
      # mean exprssion
      exp_matrix_1_mean_log_exp <- apply(log(exp_matrix_1[one_exp_matrix_1_gene_ID,]+1, 2), 1, mean)
      exp_matrix_1_mean_log_exp_vector <- c(exp_matrix_1_mean_log_exp_vector, exp_matrix_1_mean_log_exp)
      exp_matrix_2_mean_log_exp <- apply(log(exp_matrix_2[one_exp_matrix_2_matching_gene_ID,]+1, 2), 1, mean)
      exp_matrix_2_mean_log_exp_vector <- c(exp_matrix_2_mean_log_exp_vector, exp_matrix_2_mean_log_exp)
      
      # IDs and categories
      exp_matrix_1_id_vector <- c(exp_matrix_1_id_vector, one_exp_matrix_1_gene_ID)
      exp_matrix_2_id_vector <- c(exp_matrix_2_id_vector, one_exp_matrix_2_matching_gene_ID)
      category_vector <- c(category_vector, as.character(match_df[i, "category"]))
      category_2_vector <- c(category_2_vector, as.character(match_df[i, "category_2"]))
    }
  }
  # mean_vector and sd_vector may contain null values
  ## add the minimum non-null absolute value to null values to avoid Inf values after log10 transformation
  mean_vector[mean_vector == 0] <- min(mean_vector[mean_vector != 0], na.rm=TRUE)
  sd_vector[sd_vector == 0] <- min(sd_vector[sd_vector != 0], na.rm=TRUE)
  consistency_vector <- -log(mean_vector,10)
  stringency_vector <- -log(sd_vector,10)
  df2return <- data.frame(exp_matrix_1_gene_ID=exp_matrix_1_id_vector, exp_matrix_2_gene_ID=exp_matrix_2_id_vector, category=category_vector, category_2=category_2_vector, mean_log2ratio=mean_vector, sd_log2ratio=sd_vector, mixed_signs_log2ratio=mixed_signs_vector, consistency=consistency_vector, stringency=stringency_vector, exp_matrix_1_mean_exp=exp_matrix_1_mean_log_exp_vector, exp_matrix_2_mean_exp=exp_matrix_2_mean_log_exp_vector)
  ## put expression matrices names in column names
  colnames(df2return)[colnames(df2return)=="exp_matrix_1_gene_ID"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
  colnames(df2return)[colnames(df2return)=="exp_matrix_2_gene_ID"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
  colnames(df2return)[colnames(df2return)=="exp_matrix_1_mean_exp"] <- sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name))
  colnames(df2return)[colnames(df2return)=="exp_matrix_2_mean_exp"] <- sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name))
  
  # add consistency/stringency categories
  consistency_stringency_category <- c()
  consistency_stringency_category_2 <- c()
  consistency_stringency_category_3 <- c()
  for (i in 1:dim(df2return)[1]) {
    one_consistency <- df2return[i, "consistency"]
    one_stringency <- df2return[i, "stringency"]
    ## threshold: 1
    threshold <- 1
    one_category <- consistency_stringency_category(one_consistency, one_stringency, threshold)
    consistency_stringency_category <- c(consistency_stringency_category, one_category)
    ## threshold: 5% difference of expression, i.e. consistency/stringency threshold: -log10(log2(1.1))=1.152493
    threshold <- -log(log(1.05, 2), 10)
    one_category <- consistency_stringency_category(one_consistency, one_stringency, threshold)
    consistency_stringency_category_2 <- c(consistency_stringency_category_2, one_category)
    ## threshold: 10% difference of expression, i.e. consistency/stringency threshold: -log10(log2(1.1))=0.8616862
    threshold <- -log(log(1.1, 2), 10)
    one_category <- consistency_stringency_category(one_consistency, one_stringency, threshold)
    consistency_stringency_category_3 <- c(consistency_stringency_category_3, one_category)
  }
  df2return <- cbind(df2return, consistency_stringency_category, consistency_stringency_category_2, consistency_stringency_category_3)
  df2return$consistency_stringency_category <- as.factor(consistency_stringency_category)
  df2return$consistency_stringency_category_2 <- as.factor(consistency_stringency_category_2)
  df2return$consistency_stringency_category_3 <- as.factor(consistency_stringency_category_3)
  return(df2return)
}

## plots
consistency_stringency_plot <- function(data_df, exp_matrix_1_name, exp_matrix_2_name, gene_category, nb_genes, out_dir, out_name) {
  
  # plot theme
  plot_theme <- theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8))
  
  # nb_genes <- dim(data_df)[1]
  
  # mean expression density plot
  mean_exp_df2ggplot <- data.frame(exp=c(data_df[, sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name))], data_df[, sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name))]), genome=c(rep(exp_matrix_1_name, dim(data_df)[1]), rep(exp_matrix_2_name, dim(data_df)[1])))
  density_plot_title <- sprintf("%s: %d\nMean log2 expression", gene_category, nb_genes)
  percentiles_df <- density_percentiles_plot(mean_exp_df2ggplot, "exp", "genome", density_plot_title)
  write.csv(percentiles_df, file=sprintf("%s/%s_mean_expression_percentiles.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  # consistency/stringency per gene ID match category
  p <- ggplot(data_df, aes(x=consistency, color=category_2)) +
    geom_line(stat="density") +
    labs(title=sprintf("%s: %d\nConsistency per gene ID matching category", gene_category, nb_genes), x="Consistency", y="Density", color="ID category") +
    plot_theme
  print(p)
  p <- ggplot(data_df, aes(x=stringency, color=category_2)) +
    geom_line(stat="density") +
    labs(title=sprintf("%s: %d\nStringency per gene ID matching category", gene_category, nb_genes), x="Stringency", y="Density", color="ID category") +
    plot_theme
  print(p)
  
  # consistency/stringency thresholds: 1, -log10(log2(1.1))
  threshold_comp_df <- data.frame()
  for (threshold in c(1, -log10(log2(1.05)), -log10(log2(1.1)))) {
    ## count genes according consistency/stringency categories
    if (threshold == 1) {
      category_name <- "consistency_stringency_category"
      threshold_to_print <- sprintf("%d", threshold)
      threshold_in_path <- sprintf("%d", threshold)
    } else {
      category_name <- ifelse(threshold == -log10(log2(1.05)), "consistency_stringency_category_2", "consistency_stringency_category_3")
      threshold_to_print <- sprintf("%f", threshold)
      threshold_in_path <- sub("[.]", "_", sprintf("%.2f", threshold))
    }
    nb_genes_consist_string <- unname(table(data_df[, category_name])["consistent_stringent"])
    nb_genes_consist_nonstring <- unname(table(data_df[, category_name])["consistent_nonstringent"])
    nb_genes_inconsist_string <- unname(table(data_df[, category_name])["inconsistent_stringent"])
    nb_genes_inconsist_nonstring <- unname(table(data_df[, category_name])["inconsistent_nonstringent"])
    plot_title <- sprintf("inconsistent & stringent: %.1f%%, consistent & stringent: %.1f%%\ninconsistent & non stringent: %.1f%%, consistent & non stringent: %.1f%%", nb_genes_inconsist_string/nb_genes*100, nb_genes_consist_string/nb_genes*100, nb_genes_inconsist_nonstring/nb_genes*100, nb_genes_consist_nonstring/nb_genes*100)
    
    ## consistency/stringency density plots
    if (nb_genes > 1) {
      consistency_plot_title <- sprintf("%s: %d\nConsistency and stringency - threshold: %s\n%s", gene_category, nb_genes, threshold_to_print, plot_title)
      consistency_stringency_scat_density_2d_plots(data_df, threshold, consistency_plot_title)
      
      consistency_max_threshold <- as.integer(floor(max(data_df$consistency))-2) # use of as.integer() function to avoid Error in sprintf: invalid format '%d'; use format %f, %e, %g or %a for numeric objects
      data_df2plot <- data_df[which(data_df$consistency<consistency_max_threshold),]
      nb_genes2plot <- dim(data_df2plot)[1]
      if (nb_genes2plot > 1) {
        consistency_plot_title <- sprintf("%s: %d\nConsistency and stringency (consistency < %d) - threshold: %s\n%s", gene_category, nb_genes, consistency_max_threshold, threshold_to_print, plot_title)
        consistency_stringency_scat_density_2d_plots(data_df2plot, threshold, consistency_plot_title)
      }
    }
    
    ## mean expression per consistency/stringency category
    density_plot_title <- sprintf("%s: %d\n%s mean expression per consistency/stringency category", gene_category, nb_genes, exp_matrix_1_name)
    percentiles_df_1 <- density_percentiles_plot(data_df, sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name)), category_name, density_plot_title)
    density_plot_title <- sprintf("%s: %d\n%s mean expression per consistency/stringency category", gene_category, nb_genes, exp_matrix_2_name)
    percentiles_df_2 <- density_percentiles_plot(data_df, sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name)), category_name, density_plot_title)
    percentiles_df <- merge(percentiles_df_1, percentiles_df_2, by="percentile")
    write.csv(percentiles_df, file=sprintf("%s/%s_mean_expression_categories_thres_%s.csv", out_dir, out_name, threshold_in_path), quote=FALSE, row.names=FALSE)
    
    # data frame for threshold comparison
    one_threshold_df <- data_df[, c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name)), sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name)), category_name)]
    colnames(one_threshold_df)[colnames(one_threshold_df)==category_name] <- "category"
    threshold <- rep(threshold, dim(one_threshold_df)[1])
    one_threshold_df <- cbind(one_threshold_df, threshold)
    threshold_comp_df <- rbind(threshold_comp_df, one_threshold_df)
  }
  threshold_comp_df$threshold <- as.factor(threshold_comp_df$threshold)
  # mean expression per consistency/stringency category and threshold
  for (one_category in levels(threshold_comp_df$category)) {
    plot_title <- sprintf("%s: %d\n%s mean expression per consistency/stringency threshold\nconsistency/stringency category: %s", gene_category, nb_genes, exp_matrix_1_name, gsub("_", " ", one_category))
    percentiles_df_1 <- density_percentiles_plot(threshold_comp_df[which(threshold_comp_df$category==one_category),], sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name)), "threshold", plot_title)
    plot_title <- sprintf("%s: %d\n%s mean expression per consistency/stringency threshold\nconsistency/stringency category: %s", gene_category, nb_genes, exp_matrix_2_name, gsub("_", " ", one_category))
    percentiles_df_2 <- density_percentiles_plot(threshold_comp_df[which(threshold_comp_df$category==one_category),], sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name)), "threshold", plot_title)
  }
}

consistency_stringency_scat_density_2d_plots <- function(data_df, threshold, plot_title) {
  # plot theme
  plot_theme <- theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8), strip.text=element_text(size=10))
  
  p <- ggplot(data_df, aes(x=consistency, y=stringency)) +
    geom_point(alpha=0.5) +
    geom_vline(xintercept=threshold, linetype="dashed", alpha=0.5) + 
    geom_hline(yintercept=threshold, linetype="dashed", alpha=0.5) + 
    labs(title=plot_title, x="Consistency", y="Stringency") +
    plot_theme
  print(p)
  p <- ggplot(data_df, aes(x=consistency, y=stringency)) +
    geom_density_2d_filled(bins=9) +
    scale_fill_brewer(palette="Reds") +
    geom_vline(xintercept=threshold, linetype="dashed", alpha=0.3) + 
    geom_hline(yintercept=threshold, linetype="dashed", alpha=0.3) + 
    labs(title=plot_title, x="Consistency", y="Stringency") +
    plot_theme +
    theme(legend.position="bottom")
  print(p)
}

consistency_stringency_plots_for_DE_analysis <- function(consistency_df, common_df, specific_1_df, specific_2_df, exp_matrix_1_name, exp_matrix_2_name, out_dir, out_name) {
  write.csv(consistency_df, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  exp_matrix_1_ID_colname <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
  exp_matrix_2_ID_colname <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
  
  pdf(sprintf("%s/%s.pdf", out_dir, out_name))
  # common DE genes
  gene_category <- "common DE genes"
  # common_DE_genes <- as.character(match_df_1_to_2[which(match_df_1_to_2$category_2=="1:1_1"), "feature"])
  common_DE_genes_exp_matrix_1_name <- unique(as.character(common_df[, exp_matrix_1_ID_colname]))
  common_DE_genes_exp_matrix_2_name <- unique(as.character(common_df[, exp_matrix_2_ID_colname]))
  nb_genes <- dim(common_df)[1]
  data_df <- consistency_df[which(consistency_df[[exp_matrix_1_ID_colname]] %in% common_DE_genes_exp_matrix_1_name & consistency_df[[exp_matrix_2_ID_colname]] %in% common_DE_genes_exp_matrix_2_name),]
  if (nb_genes > 0) {
    consistency_stringency_plot(data_df, exp_matrix_1_name, exp_matrix_2_name, gene_category, nb_genes, out_dir, out_name)
  }
  
  # specific DE genes
  ## exp matrix 1
  gene_category <- sprintf("%s specific DE genes", exp_matrix_1_name)
  # specific_DE_genes_1 <- as.character(match_df_1_to_2[which(match_df_1_to_2$category_2=="1:1_0_1"), "feature"])
  nb_genes <- length(unique(as.character(specific_1_df[, exp_matrix_1_ID_colname])))
  specific_DE_genes_1_exp_matrix_1_name <- unique(as.character(specific_1_df[, exp_matrix_1_ID_colname]))
  specific_DE_genes_1_exp_matrix_2_name <- unique(as.character(specific_1_df[, exp_matrix_2_ID_colname]))
  data_df <- consistency_df[which(consistency_df[[exp_matrix_1_ID_colname]] %in% specific_DE_genes_1_exp_matrix_1_name & consistency_df[[exp_matrix_2_ID_colname]] %in% specific_DE_genes_1_exp_matrix_2_name),]
  if (nb_genes > 0) {
    consistency_stringency_plot(data_df, exp_matrix_1_name, exp_matrix_2_name, gene_category, nb_genes, out_dir, out_name)
  }
  
  ## exp matrix 2
  gene_category <- sprintf("%s specific DE genes", exp_matrix_2_name)
  # specifc_DE_genes_2 <- as.character(match_df_2_to_1[which(match_df_2_to_1$category_2=="1:1_0_1"), "feature"])
  nb_genes <- length(unique(as.character(specific_2_df[, exp_matrix_1_ID_colname])))
  specific_DE_genes_1_exp_matrix_2_name <- unique(as.character(specific_2_df[, exp_matrix_1_ID_colname]))
  specific_DE_genes_2_exp_matrix_2_name <- unique(as.character(specific_2_df[, exp_matrix_2_ID_colname]))
  data_df <- consistency_df[which(consistency_df[[exp_matrix_1_ID_colname]] %in% specific_DE_genes_1_exp_matrix_2_name & consistency_df[[exp_matrix_2_ID_colname]] %in% specific_DE_genes_2_exp_matrix_2_name),]
  
  if (nb_genes > 0) {
    consistency_stringency_plot(data_df, exp_matrix_1_name, exp_matrix_2_name, gene_category, nb_genes, out_dir, out_name)
  }
  dev.off()
}

GTF_gene_stats <- function(gtf_GRanges) {
  stats2return <- list()
  # gene biotype
  gene_biotype <- unique(gtf_GRanges$gene_biotype)
  if (length(gene_biotype[! is.na(gene_biotype)]) == 1) {
    gene_biotype <- gene_biotype[! is.na(gene_biotype)] 
  } else {
    if (length(gene_biotype[! is.na(gene_biotype)]) > 1) {
      print("Warning: gene has more than one biotype")
      gene_biotype <- NA
    } else {
      if (length(gene_biotype[! is.na(gene_biotype)]) == 0) {
        print("Warning: gene has no biotype")
        gene_biotype <- NA
      }
    }
  }
  stats2return[["gene_biotype"]] <- gene_biotype
  # gene length: gtf_GRanges must only contain 1 gene
  if (length(gtf_GRanges[which(gtf_GRanges$type=="gene"),]) == 0) {
    stats2return[["gene_length"]] <- 0
  } else {
    stats2return[["gene_length"]] <- width(gtf_GRanges[which(gtf_GRanges$type=="gene"),])
  }
  # transcripts: merge annotations with overlapping ranges
  gtf_transcript_GRanges <- reduce(gtf_GRanges[which(gtf_GRanges$type=="transcript"),])
  stats2return[["transcript_nb"]] <- length(gtf_transcript_GRanges)
  stats2return[["transcript_mergelength"]] <- sum(width(gtf_transcript_GRanges))
  # exons: merge annotations with overlapping ranges
  gtf_exon_GRanges <- reduce(gtf_GRanges[which(gtf_GRanges$type=="exon"),])
  stats2return[["exon_nb"]] <- length(gtf_exon_GRanges)
  stats2return[["exon_mergelength"]] <- sum(width(gtf_exon_GRanges))
  stats2return[["exon_maxlength"]] <- max(width(gtf_exon_GRanges))
  stats2return[["exon_meanlength"]] <- mean(width(gtf_exon_GRanges))
  return(stats2return)
}

GTF_gene_stats_comparison_df <- function(GRanges_1, GRanges_2, match_df) {
  # match_df data frame must contain 2 columns
  ## the gene IDs in the first column must correspond with gene IDs found in GRanges_1
  ## the gene IDs in the second column must correspond with gene IDs found in GRanges_2
  
  # initialize vectors
  gene_biotype_1_vector <- gene_length_1_vector <- nb_transcript_1_vector <- transcript_mergelength_1_vector <- exon_nb_1_vector <- exon_mergelength_1_vector <- exon_maxlength_1_vector <- exon_meanlength_1_vector <- c()
  gene_biotype_2_vector <- gene_length_2_vector <- nb_transcript_2_vector <- transcript_mergelength_2_vector <- exon_nb_2_vector <- exon_mergelength_2_vector <- exon_maxlength_2_vector <- exon_meanlength_2_vector <- c()

  # get statistics for each gene
  nb_matches <- dim(match_df)[1]
  for (i in 1:nb_matches) {
    gene_ID <- as.character(match_df[i, 1])
    matching_gene_ID <- as.character(match_df[i, 2])
    ## gene stats in GTF file 1
    gtf_1_gene_stats <- GTF_gene_stats(GRanges_1[which(GRanges_1$gene_id==gene_ID),])
    gene_biotype_1_vector <- c(gene_biotype_1_vector, gtf_1_gene_stats[["gene_biotype"]])
    gene_length_1_vector <- c(gene_length_1_vector, gtf_1_gene_stats[["gene_length"]])
    nb_transcript_1_vector <- c(nb_transcript_1_vector, gtf_1_gene_stats[["transcript_nb"]])
    transcript_mergelength_1_vector <- c(transcript_mergelength_1_vector, gtf_1_gene_stats[["transcript_mergelength"]])
    exon_nb_1_vector <- c(exon_nb_1_vector, gtf_1_gene_stats[["exon_nb"]])
    exon_mergelength_1_vector <- c(exon_mergelength_1_vector, gtf_1_gene_stats[["exon_mergelength"]])
    exon_maxlength_1_vector <- c(exon_maxlength_1_vector, gtf_1_gene_stats[["exon_maxlength"]])
    exon_meanlength_1_vector <- c(exon_meanlength_1_vector, gtf_1_gene_stats[["exon_meanlength"]])
    ## gene stats in GTF file 2
    # print(sprintf("gene id: %s", matching_gene_ID))
    gtf_2_gene_stats <- GTF_gene_stats(GRanges_2[which(GRanges_2$gene_id==matching_gene_ID),])
    gene_biotype_2_vector <- c(gene_biotype_2_vector, gtf_2_gene_stats[["gene_biotype"]])
    # print(gtf_2_gene_stats[["gene_biotype"]])
    if (length(gtf_2_gene_stats[["gene_biotype"]]) == 0) {
      print(sprintf("gene id: %s", matching_gene_ID))
    }
    gene_length_2_vector <- c(gene_length_2_vector, gtf_2_gene_stats[["gene_length"]])
    # print(gtf_2_gene_stats[["gene_length"]])
    if (length(gtf_2_gene_stats[["gene_length"]]) == 0) {
      print(sprintf("gene id: %s", matching_gene_ID))
    }
    nb_transcript_2_vector <- c(nb_transcript_2_vector, gtf_2_gene_stats[["transcript_nb"]])
    transcript_mergelength_2_vector <- c(transcript_mergelength_2_vector, gtf_2_gene_stats[["transcript_mergelength"]])
    exon_nb_2_vector <- c(exon_nb_2_vector, gtf_2_gene_stats[["exon_nb"]])
    exon_mergelength_2_vector <- c(exon_mergelength_2_vector, gtf_2_gene_stats[["exon_mergelength"]])
    exon_maxlength_2_vector <- c(exon_maxlength_2_vector, gtf_2_gene_stats[["exon_maxlength"]])
    exon_meanlength_2_vector <- c(exon_meanlength_2_vector, gtf_2_gene_stats[["exon_meanlength"]])
    if (i %% 1000 == 0) {
      print(sprintf("%d genes read", i))
    }
  }
  
  # data frame to return
  df2return <- data.frame(feature=as.character(match_df[, 1]), matching_feature=as.character(match_df[, 2]),
                          gene_biotype_1=gene_biotype_1_vector, gene_biotype_2=gene_biotype_2_vector,
                          gene_length_1=gene_length_1_vector, gene_length_2=gene_length_2_vector,
                          nb_transcript_1=nb_transcript_1_vector, nb_transcript_2=nb_transcript_2_vector,
                          transcript_mergelength_1=transcript_mergelength_1_vector, transcript_mergelength_2=transcript_mergelength_2_vector,
                          exon_nb_1=exon_nb_1_vector, exon_nb_2=exon_nb_2_vector,
                          exon_mergelength_1=exon_mergelength_1_vector, exon_mergelength_2=exon_mergelength_2_vector,
                          exon_maxlength_1=exon_maxlength_1_vector, exon_maxlength_2=exon_maxlength_2_vector,
                          exon_meanlength_1=exon_meanlength_1_vector, exon_meanlength_2=exon_meanlength_2_vector)
  return(df2return)
}

consistency_stringency_gene_stats_for_DE_analysis <- function(consistency_df, common_df, specific_1_df, specific_2_df, gtf_1, gtf_2, exp_matrix_1_name, exp_matrix_2_name) {
  # only keep gene ID and consistency/stringency category columns
  # consistency_df <- consistency_df[, c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name)), "consistency_stringency_category")]
  df2return <- data.frame()
  specificity_vector <- c()
  # read GTF files
  print("read GTF files")
  gtf_1_GRanges <- import(gtf_1, format="gtf")
  gtf_2_GRanges <- import(gtf_2, format="gtf")
  
  # common DE genes
  print(sprintf("%s and %s common DE genes", exp_matrix_1_name, exp_matrix_2_name))
  # common_DE_genes_match_df <- match_df_1_to_2[which(match_df_1_to_2$category_2=="1:1_1"),]
  nb_genes <- dim(common_df)[1]
  if (nb_genes > 0) {
    common_DE_genes_gene_stats_df <- GTF_gene_stats_comparison_df(gtf_1_GRanges, gtf_2_GRanges, common_df)
    ## put expression matrix names in column names of gene statsistics data frame
    colnames(common_DE_genes_gene_stats_df)[colnames(common_DE_genes_gene_stats_df)=="feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
    colnames(common_DE_genes_gene_stats_df)[colnames(common_DE_genes_gene_stats_df)=="matching_feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
    for (one_stat in c("gene_length", "nb_transcript", "transcript_mergelength", "exon_nb", "exon_mergelength", "exon_maxlength", "exon_meanlength")) {
      colnames(common_DE_genes_gene_stats_df)[colnames(common_DE_genes_gene_stats_df)==sprintf("%s_1", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_1_name))
      colnames(common_DE_genes_gene_stats_df)[colnames(common_DE_genes_gene_stats_df)==sprintf("%s_2", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_2_name))
    }
    ## merge consistency/stringency and gene statsistics data frames
    consistency_genes_stats_df <- merge(consistency_df, common_DE_genes_gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
    df2return <- rbind(df2return, consistency_genes_stats_df)
    specificity_vector <- c(specificity_vector, rep("common", dim(consistency_genes_stats_df)[1]))
  }
  
  # specific DE genes
  ## exp matrix 1
  print(sprintf("%s specific DE genes", exp_matrix_1_name))
  # specific_DE_genes_1_match_df <- match_df_1_to_2[which(match_df_1_to_2$category_2=="1:1_0_1"),]
  nb_genes <- length(unique(as.character(specific_1_df[, sprintf("gene_id_%s", gsub(" ", "", exp_matrix_1_name))])))
  if (nb_genes > 0) {
    specific_DE_genes_1_gene_stats_df <- GTF_gene_stats_comparison_df(gtf_1_GRanges, gtf_2_GRanges, specific_1_df)
    ### put expression matrix names in column names of gene statsistics data frame
    colnames(specific_DE_genes_1_gene_stats_df)[colnames(specific_DE_genes_1_gene_stats_df)=="feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
    colnames(specific_DE_genes_1_gene_stats_df)[colnames(specific_DE_genes_1_gene_stats_df)=="matching_feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
    for (one_stat in c("gene_length", "nb_transcript", "transcript_mergelength", "exon_nb", "exon_mergelength", "exon_maxlength", "exon_meanlength")) {
      colnames(specific_DE_genes_1_gene_stats_df)[colnames(specific_DE_genes_1_gene_stats_df)==sprintf("%s_1", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_1_name))
      colnames(specific_DE_genes_1_gene_stats_df)[colnames(specific_DE_genes_1_gene_stats_df)==sprintf("%s_2", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_2_name))
    }
    ### merge consistency/stringency and gene statsistics data frames
    consistency_genes_stats_1_df <- merge(consistency_df, specific_DE_genes_1_gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
    df2return <- rbind(df2return, consistency_genes_stats_1_df)
    specificity_vector <- c(specificity_vector, rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_1_name)), dim(consistency_genes_stats_1_df)[1]))
  }
  
  ## exp matrix 2
  print(sprintf("%s specific DE genes", exp_matrix_2_name))
  # specific_DE_genes_2_match_df <- match_df_2_to_1[which(match_df_2_to_1$category_2=="1:1_0_1"),]
  nb_genes <- length(unique(as.character(specific_2_df[, sprintf("gene_id_%s", gsub(" ", "", exp_matrix_2_name))])))
  if (nb_genes > 0) {
    specific_DE_genes_2_gene_stats_df <- GTF_gene_stats_comparison_df(gtf_2_GRanges, gtf_1_GRanges, specific_2_df)
    ### put expression matrix names in column names of gene statsistics data frame
    colnames(specific_DE_genes_2_gene_stats_df)[colnames(specific_DE_genes_2_gene_stats_df)=="feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
    colnames(specific_DE_genes_2_gene_stats_df)[colnames(specific_DE_genes_2_gene_stats_df)=="matching_feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
    for (one_stat in c("gene_length", "nb_transcript", "transcript_mergelength", "exon_nb", "exon_mergelength", "exon_maxlength", "exon_meanlength")) {
      colnames(specific_DE_genes_2_gene_stats_df)[colnames(specific_DE_genes_2_gene_stats_df)==sprintf("%s_1", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_2_name))
      colnames(specific_DE_genes_2_gene_stats_df)[colnames(specific_DE_genes_2_gene_stats_df)==sprintf("%s_2", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_1_name))
    }
    ### merge consistency/stringency and gene statsistics data frames
    consistency_genes_stats_2_df <- merge(consistency_df, specific_DE_genes_2_gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
    df2return <- rbind(df2return, consistency_genes_stats_2_df)
    specificity_vector <- c(specificity_vector, rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_2_name)), dim(consistency_genes_stats_2_df)[1]))
  }
  
  # data frame to return
  # df2return <- rbind(consistency_genes_stats_df, consistency_genes_stats_1_df, consistency_genes_stats_2_df)
  # specificity_vector <- c(rep("common", dim(consistency_genes_stats_df)[1]), rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_1_name)), dim(consistency_genes_stats_1_df)[1]), rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_2_name)), dim(consistency_genes_stats_2_df)[1]))
  df2return <- cbind(df2return, specificity=specificity_vector)
  # colnames(df2return)[colnames(df2return)=="specificity_vector"] <- "specificity"
  return(df2return)
}

consistency_stringency_gene_stats_for_DE_analysis_2 <- function(consistency_df, gene_stats_df, common_df, specific_1_df, specific_2_df, exp_matrix_1_name, exp_matrix_2_name) {
  # only keep gene ID and consistency/stringency category columns
  # consistency_df <- consistency_df[, c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name)), "consistency_stringency_category")]
  df2return <- data.frame()
  specificity_vector <- c()
  
  exp_matrix_1_gene_ID_colname <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
  exp_matrix_2_gene_ID_colname <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
  
  # common DE genes
  print(sprintf("%s and %s common DE genes", exp_matrix_1_name, exp_matrix_2_name))
  # common_DE_genes_match_df <- match_df_1_to_2[which(match_df_1_to_2$category_2=="1:1_1"),]
  nb_genes <- dim(common_df)[1]
  if (nb_genes > 0) {
    ## merge consistency/stringency and gene statsistics data frames
    common_DE_genes_gene_stats_df <- gene_stats_df[which(gene_stats_df[[exp_matrix_1_gene_ID_colname]] %in% common_df[[exp_matrix_1_gene_ID_colname]] & gene_stats_df[[exp_matrix_2_gene_ID_colname]] %in% common_df[[exp_matrix_2_gene_ID_colname]]),]
    consistency_genes_stats_df <- merge(consistency_df, common_DE_genes_gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
    df2return <- rbind(df2return, consistency_genes_stats_df)
    specificity_vector <- c(specificity_vector, rep("common", dim(consistency_genes_stats_df)[1]))
  }
  
  # specific DE genes
  ## exp matrix 1
  print(sprintf("%s specific DE genes", exp_matrix_1_name))
  # specific_DE_genes_1_match_df <- match_df_1_to_2[which(match_df_1_to_2$category_2=="1:1_0_1"),]
  nb_genes <- length(unique(as.character(specific_1_df[, sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))])))
  if (nb_genes > 0) {
    ### merge consistency/stringency and gene statsistics data frames
    specific_DE_genes_1_gene_stats_df <- gene_stats_df[which(gene_stats_df[[exp_matrix_1_gene_ID_colname]] %in% specific_1_df[[exp_matrix_1_gene_ID_colname]] & gene_stats_df[[exp_matrix_2_gene_ID_colname]] %in% specific_1_df[[exp_matrix_2_gene_ID_colname]]),]
    consistency_genes_stats_1_df <- merge(consistency_df, specific_DE_genes_1_gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
    df2return <- rbind(df2return, consistency_genes_stats_1_df)
    specificity_vector <- c(specificity_vector, rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_1_name)), dim(consistency_genes_stats_1_df)[1]))
  }
  
  ## exp matrix 2
  print(sprintf("%s specific DE genes", exp_matrix_2_name))
  # specific_DE_genes_2_match_df <- match_df_2_to_1[which(match_df_2_to_1$category_2=="1:1_0_1"),]
  nb_genes <- length(unique(as.character(specific_2_df[, sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))])))
  if (nb_genes > 0) {
    ### merge consistency/stringency and gene statsistics data frames
    specific_DE_genes_2_gene_stats_df <- gene_stats_df[which(gene_stats_df[[exp_matrix_1_gene_ID_colname]] %in% specific_2_df[[exp_matrix_1_gene_ID_colname]] & gene_stats_df[[exp_matrix_2_gene_ID_colname]] %in% specific_2_df[[exp_matrix_2_gene_ID_colname]]),]
    consistency_genes_stats_2_df <- merge(consistency_df, specific_DE_genes_2_gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
    df2return <- rbind(df2return, consistency_genes_stats_2_df)
    specificity_vector <- c(specificity_vector, rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_2_name)), dim(consistency_genes_stats_2_df)[1]))
  }
  
  # data frame to return
  # df2return <- rbind(consistency_genes_stats_df, consistency_genes_stats_1_df, consistency_genes_stats_2_df)
  # specificity_vector <- c(rep("common", dim(consistency_genes_stats_df)[1]), rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_1_name)), dim(consistency_genes_stats_1_df)[1]), rep(sprintf("%s_specific", gsub(" ", "", exp_matrix_2_name)), dim(consistency_genes_stats_2_df)[1]))
  df2return <- cbind(df2return, specificity=specificity_vector)
  # colnames(df2return)[colnames(df2return)=="specificity_vector"] <- "specificity"
  return(df2return)
}


consistency_stringency_gene_stats_df4plots <- function(data_df, exp_matrix_1_name, exp_matrix_2_name, category_name) {
  # build data frame for ggplot
  ## select columns
  df2ggplot <- data_df[, c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("gene_length_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_length_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("exon_nb_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("exon_nb_%s", gsub(" ", "", exp_matrix_2_name)),
                           sprintf("exon_mergelength_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("exon_mergelength_%s", gsub(" ", "", exp_matrix_2_name)), 
                           category_name)]
  ## add log2 ratios
  for (stat_basename in c("gene_length", "exon_nb", "exon_mergelength")) {
    log2ratio <- log(df2ggplot[, sprintf("%s_%s", stat_basename, gsub(" ", "", exp_matrix_1_name))] / df2ggplot[, sprintf("%s_%s", stat_basename, gsub(" ", "", exp_matrix_2_name))], 2)
    df2ggplot <- cbind(df2ggplot, log2ratio)
    colnames(df2ggplot)[colnames(df2ggplot)=="log2ratio"] <- sprintf("%s_%s_%s_log2ratio", stat_basename, gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
  }
  return(df2ggplot)
}

consistency_stringency_one_gene_stat_plots <- function(data_df, stat, name_1, name_2, category_name, title) {
  # plot theme
  plot_theme <- theme_bw() +
    theme(legend.position="bottom") +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8))
  # plot titles
  ## add nb genes per consistency/stringency category
  nb_genes <- dim(data_df)[1]
  nb_genes_consist_string <- unname(table(data_df[, category_name])["consistent_stringent"])
  nb_genes_consist_nonstring <- unname(table(data_df[, category_name])["consistent_nonstringent"])
  nb_genes_inconsist_string <- unname(table(data_df[, category_name])["inconsistent_stringent"])
  nb_genes_inconsist_nonstring <- unname(table(data_df[, category_name])["inconsistent_nonstringent"])
  nb_genes_title <- sprintf("inconsistent & stringent: %d (%.1f%%), consistent & stringent: %d (%.1f%%)\ninconsistent & non stringent: %d (%.1f%%), consistent & non stringent: %d (%.1f%%)", nb_genes_inconsist_string, nb_genes_inconsist_string/nb_genes*100, nb_genes_consist_string, nb_genes_consist_string/nb_genes*100, nb_genes_inconsist_nonstring, nb_genes_inconsist_nonstring/nb_genes*100, nb_genes_consist_nonstring, nb_genes_consist_nonstring/nb_genes*100)
  plot_title <- sprintf("%s\n%s", title, nb_genes_title)
  nb_genes_title_2 <- sprintf("inconsistent & stringent: %d (%.1f%%), inconsistent & non stringent: %d (%.1f%%)", nb_genes_inconsist_string, nb_genes_inconsist_string/nb_genes*100, nb_genes_inconsist_nonstring, nb_genes_inconsist_nonstring/nb_genes*100)
  plot_title_2 <- sprintf("%s\n%s", title, nb_genes_title_2)
  
  
  # plots
  p <- ggplot(data_df, aes_string(x=sprintf("%s_%s", stat, gsub(" ", "", name_1)), y=sprintf("%s_%s", stat, gsub(" ", "", name_2)), color=sprintf("%s", category_name))) +
    geom_point() +
    stat_ellipse() +
    geom_abline(intercept=0, slope=1, linetype="dashed", alpha=0.3) +
    labs(title=sprintf("%s\n%s", plot_title, gsub("_", " ", stat)), x=name_1, y=name_2) +
    plot_theme +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
  print(p)
  p <- ggplot(data_df, aes_string(x=sprintf("%s_%s_%s_log2ratio", stat, gsub(" ", "", name_1), gsub(" ", "", name_2)), color=sprintf("%s", category_name))) +
    geom_density() +
    geom_vline(xintercept=0, linetype="dashed", alpha=0.3) +
    labs(title=sprintf("%s\n%s/%s %s log2 ratio", plot_title, name_1, name_2, gsub("_", " ", stat)), x="log2 ratio") +
    plot_theme +
    guides(color=guide_legend(nrow=2, byrow=TRUE))
  print(p)
  
  ## only 'inconsistent_stringent' and 'inconsistent_nonstringent' categories
  inconsistent_data_df <- data_df[which(data_df$consistency_stringency_category %in% c("inconsistent_stringent", "inconsistent_nonstringent")),]
  p <- ggplot(inconsistent_data_df, aes_string(x=sprintf("%s_%s_%s_log2ratio", stat, gsub(" ", "", name_1), gsub(" ", "", name_2)), color=sprintf("%s", category_name))) +
    geom_density() +
    geom_vline(xintercept=0, linetype="dashed", alpha=0.3) +
    labs(title=sprintf("%s\n%s/%s %s log2 ratio", plot_title_2, name_1, name_2, gsub("_", " ", stat)), x="log2 ratio") +
    plot_theme
  print(p)
  ## mean expression by consistency/stringency categories for 'inconsistent_stringent' and 'inconsistent_nonstringent' genes scatter plots
  p <- ggplot(inconsistent_data_df, aes_string(x=sprintf("%s_%s_%s_log2ratio", stat, gsub(" ", "", name_1), gsub(" ", "", name_2)), y=sprintf("mean_exp_%s", gsub(" ", "", name_1)), color=sprintf("%s", category_name))) +
    geom_point() +
    stat_ellipse() +
    labs(title=sprintf("%s\n%s/%s %s log2 ratio and %s mean expression", plot_title_2, name_1, name_2, gsub("_", " ", stat), name_1), x=sprintf("%s log2 ratio", gsub("_", " ", stat)), y=sprintf("%s mean log2 expression", gsub(" ", "", name_1))) +
    plot_theme
  print(p)
  p <- ggplot(inconsistent_data_df, aes_string(x=sprintf("%s_%s_%s_log2ratio", stat, gsub(" ", "", name_1), gsub(" ", "", name_2)), y=sprintf("mean_exp_%s", gsub(" ", "", name_2)), color=sprintf("%s", category_name))) +
    geom_point() +
    stat_ellipse() +
    labs(title=sprintf("%s\n%s/%s %s log2 ratio and %s mean expression", plot_title_2, name_1, name_2, gsub("_", " ", stat), name_2), x=sprintf("%s log2 ratio", gsub("_", " ", stat)), y=sprintf("%s mean log2 expression", gsub(" ", "", name_2))) +
    plot_theme
  print(p)
}

consistency_stringency_gene_stats_plots <- function(data_df, exp_matrix_1_name, exp_matrix_2_name, category_name, title) {
  for (stat_basename in c("gene_length", "exon_nb", "exon_mergelength")) {
    consistency_stringency_one_gene_stat_plots(data_df, stat_basename, exp_matrix_1_name, exp_matrix_2_name, category_name, title)
  }
}

consistency_stringency_gene_stats_plots_for_DE_analysis <- function(data_df, exp_matrix_1_name, exp_matrix_2_name, category_name, plot_title, out_dir, out_name) {
  # data frame for ggplot
  ## select columns
  df2ggplot <- data_df[, c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("mean_exp_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("gene_length_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_length_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("exon_nb_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("exon_nb_%s", gsub(" ", "", exp_matrix_2_name)), 
                           sprintf("exon_mergelength_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("exon_mergelength_%s", gsub(" ", "", exp_matrix_2_name)), 
                           category_name, "specificity"
                           )]
  ## add log2 ratios
  for (stat_basename in c("gene_length", "exon_nb", "exon_mergelength")) {
    log2ratio <- log(df2ggplot[, sprintf("%s_%s", stat_basename, gsub(" ", "", exp_matrix_1_name))] / df2ggplot[, sprintf("%s_%s", stat_basename, gsub(" ", "", exp_matrix_2_name))], 2)
    df2ggplot <- cbind(df2ggplot, log2ratio)
    colnames(df2ggplot)[colnames(df2ggplot)=="log2ratio"] <- sprintf("%s_%s_%s_log2ratio", stat_basename, gsub(" ", "", exp_matrix_1_name), gsub(" ", "", exp_matrix_2_name))
  }
  write.csv(df2ggplot, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  # plots
  pdf(sprintf("%s/%s.pdf", out_dir, out_name))
  ## common DE genes
  common_df2ggplot <- df2ggplot[which(df2ggplot$specificity=="common"),]
  plot_title_base <- sprintf("%s\n%d common DE genes", plot_title, dim(common_df2ggplot)[1])
  consistency_stringency_gene_stats_plots(common_df2ggplot, exp_matrix_1_name, exp_matrix_2_name, category_name, plot_title_base)
  ## specific DE genes
  ### exp matrix 1
  specific_df2ggplot <- df2ggplot[which(df2ggplot$specificity==sprintf("%s_specific", gsub(" ", "", exp_matrix_1_name))),]
  plot_title_base <- sprintf("%s\n%d %s DE genes", plot_title, dim(specific_df2ggplot)[1], exp_matrix_1_name)
  consistency_stringency_gene_stats_plots(specific_df2ggplot, exp_matrix_1_name, exp_matrix_2_name, category_name, plot_title_base)
  ### exp matrix 2
  specific_df2ggplot <- df2ggplot[which(df2ggplot$specificity==sprintf("%s_specific", gsub(" ", "", exp_matrix_2_name))),]
  plot_title_base <- sprintf("%s\n%d %s DE genes", plot_title, dim(specific_df2ggplot)[1], exp_matrix_2_name)
  consistency_stringency_gene_stats_plots(specific_df2ggplot, exp_matrix_1_name, exp_matrix_2_name, category_name, plot_title_base)
  dev.off()
}

# comparisons
consistency_stringency_pipeline_ensembl_ncbi_gene_ids <- function(ensembl_matrix, ncbi_matrix, ensembl_gtf, ncbi_gtf, ensembl_matrix_name, ncbi_matrix_name, ensembl_ID_col, ncbi_ID_col, consistency_n, comparison_name, dataset_name, out_dir) {
  comparison_basename <- sprintf("%s_%s", dataset_name, comparison_name)
  # Genes of interest
  ensembl_genes_of_interest <- rownames(ensembl_matrix)
  ncbi_genes_of_interest <- rownames(ncbi_matrix)
  
  # Ensembl IDs and NCBI gene names
  print("Ensembl IDs and NCBI gene names")
  gene_ID_matching_output_dir <- sprintf("%s/gene_ID_matching", out_dir)
  gene_ID_matching_output_dir <- sprintf("%s/gene_ID_matching", out_dir)
  if (! dir.exists(gene_ID_matching_output_dir)) {
    dir.create(gene_ID_matching_output_dir, recursive=TRUE, mode="0775")
  }
  # ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)
  ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", ensemblRedirect = FALSE)
  ## get matching gene IDs
  ### get Entrez accessions corresponding to Ensembl gene IDs
  gene.data.ensembl <- getBM(attributes=c(ensembl_ID_col, ncbi_ID_col), filters = ensembl_ID_col, values = ensembl_genes_of_interest, mart = ensembl)
  output_basename <- sprintf("%s_%s_useMart_getBM", dataset_name, gsub(" ", "", ensembl_matrix_name))
  #### getBM() function may include not requested IDs: remove them if any
  gene.data.ensembl.keep <- remove_not_requested_getBM_ids(gene.data.ensembl, ensembl_ID_col, ensembl_genes_of_interest, gene_ID_matching_output_dir, output_basename)
  write.csv(gene.data.ensembl.keep, file=sprintf("%s/%s.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
  ### get Ensembl gene IDs corresponding to Entrez accessions
  # gene.data.ncbi <- getBM(attributes=c('ensembl_gene_id', 'entrezgene_accession'), filters = 'entrezgene_accession', values = ncbi_genes_of_interest, mart = ensembl)
  gene.data.ncbi <- getBM(attributes=c(ncbi_ID_col, ensembl_ID_col), filters = ncbi_ID_col, values = ncbi_genes_of_interest, mart = ensembl)
  output_basename <- sprintf("%s_%s_useMart_getBM", dataset_name, gsub(" ", "", ncbi_matrix_name))
  #### getBM() function may include not requested IDs: remove them if any
  gene.data.ncbi.keep <- remove_not_requested_getBM_ids(gene.data.ncbi, ncbi_ID_col, ncbi_genes_of_interest, gene_ID_matching_output_dir, output_basename)
  write.csv(gene.data.ncbi.keep, file=sprintf("%s/%s.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
  
  ## get matching gene IDs in genes of interest
  ### ensembl to entrez gene ID matching
  print("Ensembl to Entrez gene ID matching")
  gene_IDs_col_1 <- ensembl_ID_col
  gene_IDs_col_2 <- ncbi_ID_col
  output_basename <- sprintf("%s_%s_useMart_getBM", dataset_name, gsub(" ", "", ensembl_matrix_name))
  ensembl2entrez_df <- Ensembl_NCBI_ID_matching(gene.data.ensembl.keep, gene_IDs_col_1, gene_IDs_col_2,  ensembl_genes_of_interest, ncbi_genes_of_interest, gene_ID_matching_output_dir, output_basename)
  #### add mean expression
  ensembl_matrix_mean_exp_log <- log(apply(ensembl_matrix, 1, mean) + 1, 2)
  mean_exp <- ensembl_matrix_mean_exp_log[as.character(ensembl2entrez_df$feature)]
  ensembl2entrez_df <- cbind(ensembl2entrez_df, mean_exp)
  #### add annotation column
  annotation <- rep(ensembl_matrix_name, length(mean_exp))
  ensembl2entrez_df <- cbind(ensembl2entrez_df, annotation)
  write.csv(ensembl2entrez_df, file=sprintf("%s/%s_matching_categories.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
  #### stats data frame
  nb_all_ensembl_ids <- length(ensembl_genes_of_interest)
  one_df2barplot <- gene_counts_per_ID_matching_category(ensembl2entrez_df, nb_all_ensembl_ids)
  annotation_IDs_df2barplot <- one_df2barplot
  annotation_IDs_df2expdensityplot <- ensembl2entrez_df
  
  ## entrez to ensembl gene ID matching
  print("Entrez to Ensembl gene ID matching")
  gene_IDs_col_1 <- ncbi_ID_col
  gene_IDs_col_2 <- ensembl_ID_col
  output_basename <- sprintf("%s_%s_useMart_getBM", dataset_name, gsub(" ", "", ncbi_matrix_name))
  entrez2ensembl_df <- Ensembl_NCBI_ID_matching(gene.data.ncbi.keep, gene_IDs_col_1, gene_IDs_col_2, ncbi_genes_of_interest, ensembl_genes_of_interest, gene_ID_matching_output_dir, output_basename)
  #### add mean expression
  ncbi_matrix_mean_exp_log <- log(apply(ncbi_matrix, 1, mean) + 1, 2)
  mean_exp <- ncbi_matrix_mean_exp_log[as.character(entrez2ensembl_df$feature)]
  entrez2ensembl_df <- cbind(entrez2ensembl_df, mean_exp)
  #### add annotation column
  annotation <- rep(ncbi_matrix_name, length(mean_exp))
  entrez2ensembl_df <- cbind(entrez2ensembl_df, annotation)
  write.csv(entrez2ensembl_df, file=sprintf("%s/%s_matching_categories.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
  #### stats data frame
  nb_all_entrez_accessions <- length(ncbi_genes_of_interest)
  one_df2barplot <- gene_counts_per_ID_matching_category(entrez2ensembl_df, nb_all_entrez_accessions)
  
  annotation_IDs_df2barplot <- rbind(annotation_IDs_df2barplot, one_df2barplot)
  write.csv(annotation_IDs_df2barplot, file=sprintf("%s/%s_gene_ID_matching_barplot_data.csv", out_dir, comparison_basename), quote=FALSE, row.names=FALSE)
  annotation_IDs_df2expdensityplot <- rbind(annotation_IDs_df2expdensityplot, entrez2ensembl_df)
  write.csv(annotation_IDs_df2expdensityplot, file=sprintf("%s/%s_gene_ID_matching_density_data.csv", out_dir, comparison_basename), quote=FALSE, row.names=FALSE)
  
  ## sanity check
  print("sanity check: verify that there are the same number of ID matchings for both matching request (Ensembl to Entrez and Entrez to Ensembl)")
  ensembl2entrez_matching_number <- sum(annotation_IDs_df2barplot[which(annotation_IDs_df2barplot$annotation==ensembl_matrix_name & annotation_IDs_df2barplot$category_2 %in% c("1:1_1", "1:N_1", "N:1_1", "N:P_1")), "value"]) + sum(gene.data.ensembl.keep[which(gene.data.ensembl.keep[[ensembl_ID_col]] %in% as.character(ensembl2entrez_df[which(ensembl2entrez_df$category_2=="1:N_P"),"feature"])), ncbi_ID_col] %in% ncbi_genes_of_interest) + sum(gene.data.ensembl.keep[which(gene.data.ensembl.keep[[ensembl_ID_col]] %in% as.character(ensembl2entrez_df[which(ensembl2entrez_df$category_2=="N:P_Q"),"feature"])), ncbi_ID_col] %in% ncbi_genes_of_interest)
  print(sprintf("number of Ensembl to Entrez ID matchings: %d", ensembl2entrez_matching_number))
  entrez2ensembl_matching_number <- sum(annotation_IDs_df2barplot[which(annotation_IDs_df2barplot$annotation==ncbi_matrix_name & annotation_IDs_df2barplot$category_2 %in% c("1:1_1", "1:N_1", "N:1_1", "N:P_1")), "value"]) + sum(gene.data.ncbi.keep[which(gene.data.ncbi.keep[[ncbi_ID_col]] %in% as.character(entrez2ensembl_df[which(entrez2ensembl_df$category_2=="1:N_P"),"feature"])), ensembl_ID_col] %in% ensembl_genes_of_interest) + sum(gene.data.ncbi.keep[which(gene.data.ncbi.keep[[ncbi_ID_col]] %in% as.character(entrez2ensembl_df[which(entrez2ensembl_df$category_2=="N:P_Q"),"feature"])), ensembl_ID_col] %in% ensembl_genes_of_interest)
  print(sprintf("number of Entrez to Ensembl ID matchings: %d", entrez2ensembl_matching_number))

  ## gene IDs matching categories barplots and mean expression densities
  print("Gene ID matching plots")
  output_name <- sprintf("%s_gene_ID_matching", comparison_basename)
  gene_ID_matching_plots(annotation_IDs_df2barplot, annotation_IDs_df2expdensityplot, out_dir, output_name)
  
  # consistency and stringency
  print("Consistency and stringency")
  consistency_stringency_output_basename <- sprintf("%s_consistency_stringency", comparison_basename)
  ## computations
  print("computations")
  consistency_stringency_df <- consistency_stringency_computation(ensembl_matrix, ensembl_matrix_name, ncbi_matrix, ncbi_matrix_name, ensembl2entrez_df, consistency_n)
  # consistency_stringency_df <- consistency_stringency_computation_2(ensembl_matrix, ensembl_matrix_name, ncbi_matrix, ncbi_matrix_name, ensembl2entrez_df, consistency_n)
  write.csv(consistency_stringency_df, file=sprintf("%s/%s.csv", out_dir, consistency_stringency_output_basename), quote=FALSE, row.names=FALSE)
  ## plots
  print("plots")
  pdf(sprintf("%s/%s.pdf", out_dir, consistency_stringency_output_basename))
  nb_genes <- dim(consistency_stringency_df)[1]
  consistency_stringency_plot(consistency_stringency_df, ensembl_matrix_name, ncbi_matrix_name, "expressed genes", nb_genes, out_dir, consistency_stringency_output_basename)
  dev.off()
  
  # gene statistics
  print("Gene statistics")
  ## read GTF files
  print("read GTF files")
  ensembl_gtf_GRanges <- import(ensembl_gtf, format="gtf")
  ncbi_gtf_GRanges <- import(ncbi_gtf, format="gtf")
  ## get gene stats for consistency/stringency genes
  consistency_stringency_genes <- as.character(consistency_stringency_df[[sprintf("gene_ID_%s", gsub(" ", "", ensembl_matrix_name))]])
  consistency_stringency_genes_match_df <- consistency_stringency_df[which(consistency_stringency_df[[sprintf("gene_ID_%s", gsub(" ", "", ensembl_matrix_name))]] %in% consistency_stringency_genes), c(sprintf("gene_ID_%s", gsub(" ", "", ensembl_matrix_name)), sprintf("gene_ID_%s", gsub(" ", "", ncbi_matrix_name)))]
  gene_stats_df <- GTF_gene_stats_comparison_df(ensembl_gtf_GRanges, ncbi_gtf_GRanges, consistency_stringency_genes_match_df)
  ### put expression matrices names in column names
  colnames(gene_stats_df)[colnames(gene_stats_df)=="feature"] <- sprintf("gene_ID_%s", gsub(" ", "", ensembl_matrix_name))
  colnames(gene_stats_df)[colnames(gene_stats_df)=="matching_feature"] <- sprintf("gene_ID_%s", gsub(" ", "", ncbi_matrix_name))
  for (one_stat in c("gene_biotype", "gene_length", "nb_transcript", "transcript_mergelength", "exon_nb", "exon_mergelength", "exon_maxlength", "exon_meanlength")) {
    colnames(gene_stats_df)[colnames(gene_stats_df)==sprintf("%s_1", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", ensembl_matrix_name))
    colnames(gene_stats_df)[colnames(gene_stats_df)==sprintf("%s_2", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", ncbi_matrix_name))
  }
  consistency_stringency_gene_stats_output_basename <- sprintf("%s_gene_stats", comparison_basename)
  write.csv(gene_stats_df, file=sprintf("%s/%s.csv", out_dir, consistency_stringency_gene_stats_output_basename), quote=FALSE, row.names=FALSE)
  ## gene statistics per consistency/stringency category plots
  print("gene statistics per consistency/stringency category plots")
  consistency_stringency_gene_stats_plots_output_dir <- sprintf("%s/consistency_stringency_gene_stats_plots", out_dir)
  if (! dir.exists(consistency_stringency_gene_stats_plots_output_dir)) {
    dir.create(consistency_stringency_gene_stats_plots_output_dir, recursive=TRUE, mode="0775")
  }
  ### merge consistency/stringency and gene stats data frames by gene IDs of both expression matrix
  consistency_stringency_genes_stats_df <- merge(consistency_stringency_df, gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", ensembl_matrix_name)), sprintf("gene_ID_%s", gsub(" ", "", ncbi_matrix_name))))
  for (threshold in c(1, -log10(log2(1.05)), -log10(log2(1.1)))) {
    #### consistency/stringency categories
    if (threshold == 1) {
      consistency_stringency_category_name <- "consistency_stringency_category"
      threshold_in_path <- sprintf("%d", threshold)
    } else {
      consistency_stringency_category_name <- ifelse(threshold == -log10(log2(1.05)), "consistency_stringency_category_2", "consistency_stringency_category_3")
      threshold_in_path <- sub("[.]", "_", sprintf("%.2f", threshold))
    }
    
    
    consistency_stringency_gene_stats_plots_output_basename <- sprintf("%s_gene_stats_plots_thres_%s", consistency_stringency_output_basename, threshold_in_path)
    pdf(sprintf("%s/%s.pdf", consistency_stringency_gene_stats_plots_output_dir, consistency_stringency_gene_stats_plots_output_basename))
    ### gene IDs per consistency/stringency category and gene biotype
    ensembl2entrez_count_df <- data.frame()
    entrez2ensembl_count_df <- data.frame()
    for (one_category in levels(consistency_stringency_genes_stats_df[[consistency_stringency_category_name]])) {
      consistency_stringency_genes_stats_one_category_df <- consistency_stringency_genes_stats_df[which(consistency_stringency_genes_stats_df[[consistency_stringency_category_name]] == one_category),]
      consistency_stringency_genes_stats_one_category_df <- droplevels(consistency_stringency_genes_stats_one_category_df)
      #### ensembl2entrez
      consistency_stringency_genes_stats_one_category_count_df <- count_pct_per_category(consistency_stringency_genes_stats_one_category_df, sprintf("gene_biotype_%s", gsub(" ", "", ensembl_matrix_name)), sprintf("gene_ID_%s", gsub(" ", "", ensembl_matrix_name)))
      colnames(consistency_stringency_genes_stats_one_category_count_df)[colnames(consistency_stringency_genes_stats_one_category_count_df) == "category"] <- "gene_biotype"
      consistency_stringency_category <- rep(one_category, dim(consistency_stringency_genes_stats_one_category_count_df)[1])
      consistency_stringency_genes_stats_one_category_count_df <- cbind(consistency_stringency_category, consistency_stringency_genes_stats_one_category_count_df)
      ensembl2entrez_count_df <- rbind(ensembl2entrez_count_df, consistency_stringency_genes_stats_one_category_count_df)
      #### entrez2ensembl
      consistency_stringency_genes_stats_one_category_count_df <- count_pct_per_category(consistency_stringency_genes_stats_one_category_df, sprintf("gene_biotype_%s", gsub(" ", "", ncbi_matrix_name)), sprintf("gene_ID_%s", gsub(" ", "", ncbi_matrix_name)))
      colnames(consistency_stringency_genes_stats_one_category_count_df)[colnames(consistency_stringency_genes_stats_one_category_count_df) == "category"] <- "gene_biotype"
      consistency_stringency_category <- rep(one_category, dim(consistency_stringency_genes_stats_one_category_count_df)[1])
      consistency_stringency_genes_stats_one_category_count_df <- cbind(consistency_stringency_category, consistency_stringency_genes_stats_one_category_count_df)
      entrez2ensembl_count_df <- rbind(entrez2ensembl_count_df, consistency_stringency_genes_stats_one_category_count_df)
    }
    ##### add annotation column
    annotation <- rep(ensembl_matrix_name, dim(ensembl2entrez_count_df)[1])
    ensembl2entrez_count_df <- cbind(annotation, ensembl2entrez_count_df)
    annotation <- rep(ncbi_matrix_name, dim(entrez2ensembl_count_df)[1])
    entrez2ensembl_count_df <- cbind(annotation, entrez2ensembl_count_df)
    
    ## gene IDs per consistency/stringency category and gene biotype plot
    consistency_stringency_category_biotype_df <- rbind(ensembl2entrez_count_df, entrez2ensembl_count_df)
    write.csv(consistency_stringency_category_biotype_df, file=sprintf("%s/%s_consistency_stringency_gene_biotype_counts.csv", consistency_stringency_gene_stats_plots_output_dir, consistency_stringency_gene_stats_plots_output_basename), quote=FALSE, row.names=FALSE)
    p <- ggplot(consistency_stringency_category_biotype_df, aes(x=annotation, y=count, fill=gene_biotype)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=sprintf("%.2f%%", pct)), size=2, position=position_stack(vjust=0.5)) +
      facet_wrap(~consistency_stringency_category, scales="free") +
      labs(title="Gene ID counts per consistency/stringency category and gene biotype", x="Annotation", y="# IDs", fill="Gene biotype") +
      theme_bw() +
      theme(legend.position="bottom") +
      theme(panel.border=element_rect(color="grey50")) +
      theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=5))
    print(p)
    
    #### build data frame for plots
    df2ggplot <- consistency_stringency_gene_stats_df4plots(consistency_stringency_genes_stats_df, ensembl_matrix_name, ncbi_matrix_name, consistency_stringency_category_name)
    write.csv(df2ggplot, file=sprintf("%s/%s.csv", consistency_stringency_gene_stats_plots_output_dir, consistency_stringency_gene_stats_plots_output_basename), quote=FALSE, row.names=FALSE)
    #### plots
    nb_genes <- dim(df2ggplot)[1]
    plot_title_base <- sprintf("%s and %s - %d common genes", ensembl_matrix_name, ncbi_matrix_name, nb_genes)
    consistency_stringency_gene_stats_plots(df2ggplot, ensembl_matrix_name, ncbi_matrix_name, consistency_stringency_category_name, plot_title_base)
    dev.off()
  }
}

same_gene_id_type_matching <- function(genes_ID1_interest, genes_ID2_interest, all_genes_ID2) {
  # vectors for data frame to return
  feature_vector <- matching_vector <- category_vector <- category_vector_2 <- c()
  
  # ## 1:0 --> gene IDs not in the 2nd expression matrix
  # ncbi_matrix_1_specific_features <- genes_ID1_interest[! genes_ID1_interest %in% all_genes_ID2]
  # one_category <- "specific"
  # one_category_2 <- "1:0"
  # feature_vector <- c(feature_vector, ncbi_matrix_1_specific_features)
  # matching_vector <- c(matching_vector, rep(NA, length(ncbi_matrix_1_specific_features)))
  # category_vector <- c(category_vector, rep(one_category, length(ncbi_matrix_1_specific_features)))
  # category_vector_2 <- c(category_vector_2, rep(one_category_2, length(ncbi_matrix_1_specific_features)))
  # ## 1:1_0 --> gene IDs in the 2nd expression matrix but not in genes of interest
  # common_features <- genes_ID1_interest[genes_ID1_interest %in% all_genes_ID2]
  # common_features_not_in_interest <- common_features[! common_features %in% genes_ID2_interest]
  # one_category <- "common not in genes of interest"
  # one_category_2 <- "1:1_0"
  # feature_vector <- c(feature_vector, common_features_not_in_interest)
  # matching_vector <- c(matching_vector, common_features_not_in_interest)
  # category_vector <- c(category_vector, rep(one_category, length(common_features_not_in_interest)))
  # category_vector_2 <- c(category_vector_2, rep(one_category_2, length(common_features_not_in_interest)))
  # ## 1:1_1 --> gene IDs in the 2nd expression matrix and in genes of interest
  # common_features_in_interest <- common_features[common_features %in% genes_ID2_interest]
  # one_category <- "common in genes of interest"
  # one_category_2 <- "1:1_1"
  # feature_vector <- c(feature_vector, common_features_in_interest)
  # matching_vector <- c(matching_vector, common_features_in_interest)
  # category_vector <- c(category_vector, rep(one_category, length(common_features_in_interest)))
  # category_vector_2 <- c(category_vector_2, rep(one_category_2, length(common_features_in_interest)))
  
  
  ## specific gene ID
  specific_features <- genes_ID1_interest[! genes_ID1_interest %in% genes_ID2_interest]
  nb_specific_features <- length(specific_features)
  one_category <- "unique match not in genes of interest"
  one_category_2 <- "1:1_0"
  feature_vector <- c(feature_vector, specific_features)
  matching_vector <- c(matching_vector, specific_features)
  category_vector <- c(category_vector, rep(one_category, nb_specific_features))
  category_vector_2 <- c(category_vector_2, rep(one_category_2, nb_specific_features))
  
  ## common gene id
  common_features <- genes_ID1_interest[genes_ID1_interest %in% genes_ID2_interest]
  nb_common_features <- length(common_features)
  one_category <- "unique match in genes of interest"
  one_category_2 <- "1:1_1"
  feature_vector <- c(feature_vector, common_features)
  matching_vector <- c(matching_vector, common_features)
  category_vector <- c(category_vector, rep(one_category, nb_common_features))
  category_vector_2 <- c(category_vector_2, rep(one_category_2, nb_common_features))
  
  df2return <- data.frame(feature=feature_vector, matching_feature=matching_vector, category=category_vector, category_2=category_vector_2)
  return(df2return)
}

consistency_stringency_pipeline_same_gene_id_type <- function(exp_matrix_1, exp_matrix_2, gtf_1, gtf_2, exp_matrix_1_name, exp_matrix_2_name, consistency_n, comparison_name, dataset_name, out_dir) {
  comparison_basename <- sprintf("%s_%s", dataset_name, comparison_name)
  # matching gene IDs
  print("Gene ID matching")
  print(sprintf("%s to %s gene ID matching", exp_matrix_1_name, exp_matrix_2_name))
  gene_ID_matching_output_dir <- sprintf("%s/gene_ID_matching", out_dir)
  if (! dir.exists(gene_ID_matching_output_dir)) {
    dir.create(gene_ID_matching_output_dir, recursive=TRUE, mode="0775")
  }
  genes_ID1_interest <- rownames(exp_matrix_1)
  genes_ID2_interest <- rownames(exp_matrix_2)
  output_basename <- sprintf("%s_%s", dataset_name, gsub(" ", "", exp_matrix_1_name))
  exp_matrix_1_to_2_match_df <- same_gene_id_type_matching(genes_ID1_interest, genes_ID2_interest, rownames(exp_matrix_2))
  #### add mean expression
  exp_matrix_1_mean_exp_log <- log(apply(exp_matrix_1, 1, mean) + 1, 2)
  mean_exp <- exp_matrix_1_mean_exp_log[as.character(exp_matrix_1_to_2_match_df$feature)]
  exp_matrix_1_to_2_match_df <- cbind(exp_matrix_1_to_2_match_df, mean_exp)
  #### add annotation column
  annotation <- rep(exp_matrix_1_name, length(mean_exp))
  exp_matrix_1_to_2_match_df <- cbind(exp_matrix_1_to_2_match_df, annotation)
  write.csv(exp_matrix_1_to_2_match_df, file=sprintf("%s/%s_matching_categories.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
  #### stats data frame
  nb_all_exp_matrix_1_ids <- length(genes_ID1_interest)
  one_df2barplot <- gene_counts_per_ID_matching_category(exp_matrix_1_to_2_match_df, nb_all_exp_matrix_1_ids)
  annotation_IDs_df2barplot <- one_df2barplot
  annotation_IDs_df2expdensityplot <- exp_matrix_1_to_2_match_df
  
  print(sprintf("%s to %s gene ID matching", exp_matrix_2_name, exp_matrix_1_name))
  genes_ID1_interest <- rownames(exp_matrix_2)
  genes_ID2_interest <- rownames(exp_matrix_1)
  output_basename <- sprintf("%s_%s", dataset_name, gsub(" ", "", exp_matrix_2_name))
  exp_matrix_2_to_1_match_df <- same_gene_id_type_matching(genes_ID1_interest, genes_ID2_interest, rownames(exp_matrix_1))
  #### add mean expression
  exp_matrix_2_mean_exp_log <- log(apply(exp_matrix_2, 1, mean) + 1, 2)
  mean_exp <- exp_matrix_2_mean_exp_log[as.character(exp_matrix_2_to_1_match_df$feature)]
  exp_matrix_2_to_1_match_df <- cbind(exp_matrix_2_to_1_match_df, mean_exp)
  #### add annotation column
  annotation <- rep(exp_matrix_2_name, length(mean_exp))
  exp_matrix_2_to_1_match_df <- cbind(exp_matrix_2_to_1_match_df, annotation)
  write.csv(exp_matrix_2_to_1_match_df, file=sprintf("%s/%s_matching_categories.csv", gene_ID_matching_output_dir, output_basename), quote=FALSE, row.names=FALSE)
  #### stats data frame
  nb_all_exp_matrix_2_ids <- length(genes_ID1_interest)
  one_df2barplot <- gene_counts_per_ID_matching_category(exp_matrix_2_to_1_match_df, nb_all_exp_matrix_2_ids)
  
  annotation_IDs_df2barplot <- rbind(annotation_IDs_df2barplot, one_df2barplot)
  write.csv(annotation_IDs_df2barplot, file=sprintf("%s/%s_gene_ID_matching_barplot_data.csv", out_dir, comparison_basename), quote=FALSE, row.names=FALSE)
  annotation_IDs_df2expdensityplot <- rbind(annotation_IDs_df2expdensityplot, exp_matrix_2_to_1_match_df)
  write.csv(annotation_IDs_df2expdensityplot, file=sprintf("%s/%s_gene_ID_matching_density_data.csv", out_dir, comparison_basename), quote=FALSE, row.names=FALSE)
  
  ## gene IDs matching categories barplots and mean expression densities
  print("plots")
  output_name <- sprintf("%s_gene_ID_matching", comparison_basename)
  gene_ID_matching_plots(annotation_IDs_df2barplot, annotation_IDs_df2expdensityplot, out_dir, output_name)
  
  # consistency and stringency
  print("Consistency and stringency")
  consistency_stringency_output_basename <- sprintf("%s_consistency_stringency", comparison_basename)
  ## computations
  print("computations")
  consistency_stringency_df <- consistency_stringency_computation(exp_matrix_1, exp_matrix_1_name, exp_matrix_2, exp_matrix_2_name, exp_matrix_1_to_2_match_df, consistency_n)
  # consistency_stringency_df <- consistency_stringency_computation_2(exp_matrix_1, exp_matrix_1_name, exp_matrix_2, exp_matrix_2_name, exp_matrix_1_to_2_match_df, consistency_n)
  write.csv(consistency_stringency_df, file=sprintf("%s/%s.csv", out_dir, consistency_stringency_output_basename), quote=FALSE, row.names=FALSE)
  ## plots
  print("plots")
  pdf(sprintf("%s/%s.pdf", out_dir, consistency_stringency_output_basename))
  nb_genes <- dim(consistency_stringency_df)
  consistency_stringency_plot(consistency_stringency_df, exp_matrix_1_name, exp_matrix_2_name, "expressed genes", nb_genes, out_dir, consistency_stringency_output_basename)
  dev.off()
  
  # gene statistics
  print("Gene statistics")
  ## read GTF files
  print("read GTF files")
  gtf_1_GRanges <- import(gtf_1, format="gtf")
  gtf_2_GRanges <- import(gtf_2, format="gtf")
  ## get gene stats for consistency/stringency genes
  ### get gene ID matching for consistency/stringency genes
  consistency_stringency_genes <- as.character(consistency_stringency_df[[sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))]])
  consistency_stringency_genes_match_df <- exp_matrix_1_to_2_match_df[which(exp_matrix_1_to_2_match_df$feature %in% consistency_stringency_genes),]
  gene_stats_df <- GTF_gene_stats_comparison_df(gtf_1_GRanges, gtf_2_GRanges, consistency_stringency_genes_match_df)
  ### put expression matrices names in column names
  colnames(gene_stats_df)[colnames(gene_stats_df)=="feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name))
  colnames(gene_stats_df)[colnames(gene_stats_df)=="matching_feature"] <- sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))
  for (one_stat in c("gene_biotype", "gene_length", "nb_transcript", "transcript_mergelength", "exon_nb", "exon_mergelength", "exon_maxlength", "exon_meanlength")) {
    colnames(gene_stats_df)[colnames(gene_stats_df)==sprintf("%s_1", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_1_name))
    colnames(gene_stats_df)[colnames(gene_stats_df)==sprintf("%s_2", one_stat)] <- sprintf("%s_%s", one_stat, gsub(" ", "", exp_matrix_2_name))
  }
  consistency_stringency_gene_stats_output_basename <- sprintf("%s_gene_stats", comparison_basename)
  write.csv(gene_stats_df, file=sprintf("%s/%s.csv", out_dir, consistency_stringency_gene_stats_output_basename), quote=FALSE, row.names=TRUE)
  ## gene statistics per consistency/stringency category plots
  print("gene statistics per consistency/stringency category plots")
  consistency_stringency_gene_stats_plots_output_dir <- sprintf("%s/consistency_stringency_gene_stats_plots", out_dir)
  if (! dir.exists(consistency_stringency_gene_stats_plots_output_dir)) {
    dir.create(consistency_stringency_gene_stats_plots_output_dir, recursive=TRUE, mode="0775")
  }
  ### merge consistency/stringency and gene stats data frames
  # consistency_stringency_genes_stats_df <- merge(consistency_stringency_df, gene_stats_df, by.x=sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), by.y="row.names")
  consistency_stringency_genes_stats_df <- merge(consistency_stringency_df, gene_stats_df, by=c(sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name))))
  for (threshold in c(1, -log10(log2(1.05)), -log10(log2(1.1)))) {
    #### consistency/stringency categories
    if (threshold == 1) {
      consistency_stringency_category_name <- "consistency_stringency_category"
      threshold_in_path <- sprintf("%d", threshold)
    } else {
      consistency_stringency_category_name <- ifelse(threshold == -log10(log2(1.05)), "consistency_stringency_category_2", "consistency_stringency_category_3")
      threshold_in_path <- sub("[.]", "_", sprintf("%.2f", threshold))
    }
    
    consistency_stringency_gene_stats_plots_output_basename <- sprintf("%s_gene_stats_plots_thres_%s", consistency_stringency_output_basename, threshold_in_path)
    pdf(sprintf("%s/%s.pdf", consistency_stringency_gene_stats_plots_output_dir, consistency_stringency_gene_stats_plots_output_basename))
    ### gene IDs per consistency/stringency category and gene biotype
    exp_matrix_1_to_2_count_df <- data.frame()
    exp_matrix_2_to_1_count_df <- data.frame()
    for (one_category in levels(consistency_stringency_genes_stats_df[[consistency_stringency_category_name]])) {
      consistency_stringency_genes_stats_one_category_df <- consistency_stringency_genes_stats_df[which(consistency_stringency_genes_stats_df[[consistency_stringency_category_name]] == one_category),]
      consistency_stringency_genes_stats_one_category_df <- droplevels(consistency_stringency_genes_stats_one_category_df)
      #### ensembl2entrez
      consistency_stringency_genes_stats_one_category_count_df <- count_pct_per_category(consistency_stringency_genes_stats_one_category_df, sprintf("gene_biotype_%s", gsub(" ", "", exp_matrix_1_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_1_name)))
      colnames(consistency_stringency_genes_stats_one_category_count_df)[colnames(consistency_stringency_genes_stats_one_category_count_df) == "category"] <- "gene_biotype"
      consistency_stringency_category <- rep(one_category, dim(consistency_stringency_genes_stats_one_category_count_df)[1])
      consistency_stringency_genes_stats_one_category_count_df <- cbind(consistency_stringency_category, consistency_stringency_genes_stats_one_category_count_df)
      exp_matrix_1_to_2_count_df <- rbind(exp_matrix_1_to_2_count_df, consistency_stringency_genes_stats_one_category_count_df)
      #### entrez2ensembl
      consistency_stringency_genes_stats_one_category_count_df <- count_pct_per_category(consistency_stringency_genes_stats_one_category_df, sprintf("gene_biotype_%s", gsub(" ", "", exp_matrix_2_name)), sprintf("gene_ID_%s", gsub(" ", "", exp_matrix_2_name)))
      colnames(consistency_stringency_genes_stats_one_category_count_df)[colnames(consistency_stringency_genes_stats_one_category_count_df) == "category"] <- "gene_biotype"
      consistency_stringency_category <- rep(one_category, dim(consistency_stringency_genes_stats_one_category_count_df)[1])
      consistency_stringency_genes_stats_one_category_count_df <- cbind(consistency_stringency_category, consistency_stringency_genes_stats_one_category_count_df)
      exp_matrix_2_to_1_count_df <- rbind(exp_matrix_2_to_1_count_df, consistency_stringency_genes_stats_one_category_count_df)
    }
    ##### add annotation column
    annotation <- rep(exp_matrix_1_name, dim(exp_matrix_1_to_2_count_df)[1])
    exp_matrix_1_to_2_count_df <- cbind(annotation, exp_matrix_1_to_2_count_df)
    annotation <- rep(exp_matrix_2_name, dim(exp_matrix_2_to_1_count_df)[1])
    exp_matrix_2_to_1_count_df <- cbind(annotation, exp_matrix_2_to_1_count_df)
    
    ## gene IDs per consistency/stringency category and gene biotype plot
    consistency_stringency_category_biotype_df <- rbind(exp_matrix_1_to_2_count_df, exp_matrix_2_to_1_count_df)
    write.csv(consistency_stringency_category_biotype_df, file=sprintf("%s/%s_consistency_stringency_gene_biotype_counts.csv", consistency_stringency_gene_stats_plots_output_dir, consistency_stringency_gene_stats_plots_output_basename), quote=FALSE, row.names=FALSE)
    p <- ggplot(consistency_stringency_category_biotype_df, aes(x=annotation, y=count, fill=gene_biotype)) +
      geom_bar(stat="identity") +
      geom_text(aes(label=sprintf("%.2f%%", pct)), size=2, position=position_stack(vjust=0.5)) +
      facet_wrap(~consistency_stringency_category, scales="free") +
      labs(title="Gene ID counts per consistency/stringency category and gene biotype", x="Annotation", y="# IDs", fill="Gene biotype") +
      theme_bw() +
      theme(legend.position="bottom") +
      theme(panel.border=element_rect(color="grey50")) +
      theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=5))
    print(p)
    
    #### build data frame for plots
    df2ggplot <- consistency_stringency_gene_stats_df4plots(consistency_stringency_genes_stats_df, exp_matrix_1_name, exp_matrix_2_name, consistency_stringency_category_name)
    write.csv(df2ggplot, file=sprintf("%s/%s.csv", consistency_stringency_gene_stats_plots_output_dir, consistency_stringency_gene_stats_plots_output_basename), quote=FALSE, row.names=FALSE)
    #### plots
    nb_genes <- dim(df2ggplot)[1]
    plot_title_base <- sprintf("%s and %s - %d common genes", exp_matrix_1_name, exp_matrix_2_name, nb_genes)
    consistency_stringency_gene_stats_plots(df2ggplot, exp_matrix_1_name, exp_matrix_2_name, consistency_stringency_category_name, plot_title_base)
    dev.off()
  }
}


###############
# DE analysis #
###############

fit_NB_GLM <- function(exp_mat, condition, design, out_dir, out_name) {
  
  y <- DGEList(counts=exp_mat, group=condition)
  
  # filter lowly expressed genes
  keep_genes <- filterByExpr(y)
  ## previously: keep_genes_old <- rowMeans(cpm(y)) >= 1
  y <- y[keep_genes, , keep.lib.sizes=FALSE]
  
  # normalization
  y <- calcNormFactors(y)
  
  # data exploration
  ## muldidimensional scaling plot of distances between gene expression profiles
  pdf(sprintf("%s/%s_MDS.pdf", out_dir, out_name))
  plotMDS(y, col=as.numeric(condition))
  dev.off()
  
  # write design matrix
  write.csv(design, file=sprintf("%s/%s_design.csv", out_dir, out_name), quote=FALSE)
  
  # estimate dispersion
  y <- estimateDisp(y, design)
  ## dispersion estimates on a BCV plot
  pdf(sprintf("%s/%s_BCV.pdf", out_dir, out_name))
  plotBCV(y)
  dev.off()
  
  # fit model
  ## The quasi-likelihood method is highly recommended for differential expression analyses of bulk RNA-seq data as it gives stricter error rate control by accounting for the uncertainty in dispersion estimation. (edgeR user's guide)
  fit <- glmQLFit(y, design)
  write.csv(fit$coefficients, file=sprintf("%s/%s_coef.csv", out_dir, out_name), quote=FALSE)
  write.csv(fit$samples, file=sprintf("%s/%s_samples.csv", out_dir, out_name), quote=FALSE)
  
  return(fit)
}

all_combinations_mutliple_factors_contrasts <- function(condition, intercept, append, first) {
  contrast_vector <- c()
  contrast_name_vector <- c()
  nb_levels <- length(levels(condition))
  buid_1vsall_contrasts <- ifelse(nb_levels > 2, TRUE, FALSE)
  contrast_name_append <- ifelse(is.null(append), "", sprintf("%s.", gsub(" ", "", append)))
  for(i in 1:nb_levels) {
    level_i <- levels(condition)[i]
    level_i_contrast <- ifelse(is.null(append), make.names(level_i), ifelse(first, sprintf("%s.%s", make.names(append), make.names(level_i)), sprintf("%s.%s", make.names(level_i), make.names(append))))
    if (buid_1vsall_contrasts) {
      if (i == 1) {
        # the level is not required to make contrast (reference level) if intercept is used
        one_1vsall_contrast <- ifelse(intercept, "(", sprintf("%s-(", level_i_contrast))
      } else {
        one_1vsall_contrast <- sprintf("%s-(", level_i_contrast)
      }
      one_1vsall_contrast_name <- sprintf("%svs", gsub(" ", "", level_i))
      first_level4comparison <- TRUE
    }
    for(j in 1:nb_levels) {
      if (j != i) {
        level_j <- levels(condition)[j]
        level_j_contrast <- ifelse(is.null(append), make.names(level_j), ifelse(first, sprintf("%s.%s", make.names(append), make.names(level_j)), sprintf("%s.%s", make.names(level_j), make.names(append))))
        # 1 vs 1 contrast
        if (j != 1) {
          if (i == 1) {
            # the level is not required to make contrast (reference level)
            one_contrast <- ifelse(intercept, level_j_contrast, sprintf("%s-%s", level_j_contrast, level_i_contrast))
            contrast_vector <- c(contrast_vector, one_contrast)
            contrast_name_vector <- c(contrast_name_vector, sprintf("%s%svs%s", contrast_name_append, gsub(" ", "", level_j), gsub(" ", "", level_i)))
          } else {
            if (j > i) {
              contrast_vector <- c(contrast_vector, sprintf("%s-%s", level_j_contrast, level_i_contrast))
              contrast_name_vector <- c(contrast_name_vector, sprintf("%s%svs%s", contrast_name_append, gsub(" ", "", level_j), gsub(" ", "", level_i)))
            }
          }
        }
        # 1 vs all contrast
        if (buid_1vsall_contrasts) {
          if (j == 1) {
            if (! intercept) {
              ## if j == 1 and no intercept is used, then the level is required to make contrast (reference level)
              if (first_level4comparison) {
                one_1vsall_contrast <- sprintf("%s%s", one_1vsall_contrast, level_j_contrast)
                first_level4comparison <- FALSE
              } else {
                one_1vsall_contrast <- sprintf("%s+%s", one_1vsall_contrast, level_j_contrast)
              }
            }
          } else {
            if (first_level4comparison) {
              one_1vsall_contrast <- sprintf("%s%s", one_1vsall_contrast, level_j_contrast)
              first_level4comparison <- FALSE
            } else {
              one_1vsall_contrast <- sprintf("%s+%s", one_1vsall_contrast, level_j_contrast)
            }
            
          }
          one_1vsall_contrast_name <- sprintf("%s%s", one_1vsall_contrast_name, gsub(" ", "", level_j))
        }
      }
    }
    if (buid_1vsall_contrasts) {
      one_1vsall_contrast <- sprintf("%s)/%d", one_1vsall_contrast, nb_levels-1)
      contrast_vector <- c(contrast_vector, one_1vsall_contrast)
      contrast_name_vector <- c(contrast_name_vector, sprintf("%s%s", contrast_name_append, one_1vsall_contrast_name))
    }
  }
  contrast_list <- list(contrasts=contrast_vector, names=contrast_name_vector)
  return(contrast_list)
}






build_1vs1_1vsall_contrasts <- function(design, condition) {
  contrast_vector <- c()
  contrast_name_vector <- c()
  
  
  
  
  
  
  nb_levels <- length(levels(condition))
  buid_1vsall_contrasts <- ifelse(nb_levels > 2, TRUE, FALSE)
  for(i in 1:nb_levels) {
    level_i <- levels(condition)[i]
    if (buid_1vsall_contrasts) {
      if (i == 1) {
        # the level is not required to make contrast (reference level)
        one_1vsall_contrast <- "("
      } else {
        one_1vsall_contrast <- sprintf("%s-(", make.names(level_i))
      }
      one_1vsall_contrast_name <- sprintf("%svs", gsub(" ", "", level_i))
      first_level4comparison <- TRUE
    }
    for(j in 1:nb_levels) {
      if (j != i) {
        level_j <- levels(condition)[j]
        if (j != 1) {
          # 1 vs 1 contrast
          if (i == 1) {
            # the level is not required to make contrast (reference level)
            contrast_vector <- c(contrast_vector, make.names(level_j))
            contrast_name_vector <- c(contrast_name_vector, sprintf("%svs%s", gsub(" ", "", level_j), gsub(" ", "", level_i)))
          } else {
            if (j > i) {
              contrast_vector <- c(contrast_vector, sprintf("%s-%s", make.names(level_j), make.names(level_i)))
              contrast_name_vector <- c(contrast_name_vector, sprintf("%svs%s", gsub(" ", "", level_j), gsub(" ", "", level_i)))
            }
          }
          
          # 1 vs all contrast
          if (buid_1vsall_contrasts) {
            ## if j == 1, then the level is not required to make contrast (reference level)
            if (first_level4comparison) {
              one_1vsall_contrast <- sprintf("%s%s", one_1vsall_contrast, make.names(level_j))
              first_level4comparison <- FALSE
            } else {
              one_1vsall_contrast <- sprintf("%s+%s", one_1vsall_contrast, make.names(level_j))
            }
          }
        }
        if (buid_1vsall_contrasts) {
          one_1vsall_contrast_name <- sprintf("%s%s", one_1vsall_contrast_name, gsub(" ", "", level_j))
        }
      }
    }
    if (buid_1vsall_contrasts) {
      one_1vsall_contrast <- sprintf("%s)/%d", one_1vsall_contrast, nb_levels-1)
      contrast_vector <- c(contrast_vector, one_1vsall_contrast)
      contrast_name_vector <- c(contrast_name_vector, one_1vsall_contrast_name)
    }
  }
  
  
  contrasts <- makeContrasts(contrasts=contrast_vector, levels=design)
  colnames(contrasts) <- contrast_name_vector
  return(contrasts)
}


DE_tests <- function(fit, con, nb_features, out_dir, out_name) {
  
  print ("Testing for DE genes")
  # the quasi-likelihood F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It provides more robust and reliable error rate control when the number of replicates is small. (edgeR user's guide)
  
  # fold-change: 1
  fc_threshold <- 1
  print (sprintf("FC threshold: %.2f", fc_threshold))
  fc_out_name <- sprintf("%s_FC_%d", out_name, fc_threshold)
  qlf <- glmQLFTest(fit, contrast=con)
  results <- topTags(qlf, n=nb_features)
  write.csv(results$table, file=sprintf("%s/%s_pval.csv", out_dir, fc_out_name), quote=FALSE)
  print(summary(decideTests(qlf)))
  ## plot logFC against log-counts per million
  pdf(sprintf("%s/%s_MD.pdf", out_dir, fc_out_name))
  plotMD(qlf)
  abline(h=0, col="blue")
  dev.off()
  
  fc_thresholds <- c(1.1, 1.25, 1.5, 2)
  for (one_fc_threshold in fc_thresholds) {
    print (sprintf("FC threshold: %.2f", one_fc_threshold))
    log2_fc_threshold <- log2(one_fc_threshold)
    fc_out_name <- sprintf("%s_FC_%s", out_name, sub("\\.", "_", one_fc_threshold))
    tr <- glmTreat(fit, contrast=con, lfc=log2_fc_threshold)
    results <- topTags(tr, n=nb_features)
    write.csv(results$table, file=sprintf("%s/%s_pval.csv", out_dir, fc_out_name), quote=FALSE)
    print(summary(decideTests(tr)))
    ## plot logFC against log-counts per million
    pdf(sprintf("%s/%s_MD.pdf", out_dir, fc_out_name))
    plotMD(tr)
    abline(h=c(-log2_fc_threshold, log2_fc_threshold), col="blue")
    dev.off()
  }
}

DE_test <- function(fit, con, fc, out_dir, out_name) {
  log2_fc <- log2(fc)
  if (fc == 1) {
    # fold-change: 1
    # the quasi-likelihood F-test is preferred as it reflects the uncertainty in estimating the dispersion for each gene. It provides more robust and reliable error rate control when the number of replicates is small. (edgeR user's guide)
    dge <- glmQLFTest(fit, contrast=con)
    ablines <- log2_fc
  } else {
    if (fc > 1) {
      # fold-change > 1
      dge <- glmTreat(fit, contrast=con, lfc=log2_fc)
      ablines <- c(-log2_fc, log2_fc)
    }
  }
  # plot logFC against log-counts per million
  pdf(sprintf("%s/%s_MD.pdf", out_dir, out_name))
  plotMD(dge)
  abline(h=ablines, col="blue")
  dev.off()
  
  return(dge)
}

DE_analysis_pipeline <- function(exp_mat, condition, design, contrasts, fc_thresholds, annotationdbi_id_type, biomart_id_type, out_dir, out_name) {
  fit <- fit_NB_GLM(exp_mat, condition, design, out_dir, out_name)
  # write contrasts
  write.csv(contrasts, file=sprintf("%s/%s_contrasts.csv", out_dir, out_name), quote=FALSE, row.names=TRUE)
  for (one_contrast in colnames(contrasts)) {
    print(sprintf("contrast: %s", one_contrast))
    contrast_output_dir <- sprintf("%s/%s", out_dir, gsub("\\.", "_", one_contrast))
    contrast_basename <- sprintf("%s_%s", out_name, gsub("\\.", "_", one_contrast))
    for (one_fc_threshold in fc_thresholds) {
      print(sprintf("FC threshold: %s", one_fc_threshold))
      fc_output_dir <- sprintf("%s/FC_%s", contrast_output_dir, sub("\\.", "_", one_fc_threshold))
      if (! dir.exists(fc_output_dir)) {
        dir.create(fc_output_dir, recursive=TRUE, mode="0775")
      }
      # DE tests
      print ("Testing for DE genes")
      fc_output_name <- sprintf("%s_FC_%s", contrast_basename, sub("\\.", "_", one_fc_threshold))
      dge <- DE_test(fit, contrasts[,one_contrast], one_fc_threshold, fc_output_dir, fc_output_name)
      print(summary(decideTests(dge)))
      de_results <- topTags(dge, n=dim(dge$table)[1])
      write.csv(de_results$table, file=sprintf("%s/%s_pval.csv", fc_output_dir, fc_output_name), quote=FALSE)
      de_genes <- rownames(de_results[which(de_results$table$FDR < 0.05),])
      if (length(de_genes) > 0) {
        # GO term enrichment analysis
        go_output_dir <- sprintf("%s/GO_analysis", fc_output_dir)
        go_output_name <- sprintf("%s_GO", fc_output_name)
        if (! dir.exists(go_output_dir)) {
          dir.create(go_output_dir, recursive=TRUE, mode="0775")
        }
        print ("GO term enrichment analysis")
        go <- GO_term_analysis(dge, annotationdbi_id_type, biomart_id_type, go_output_dir, go_output_name)
        write.table(go, file=sprintf("%s/%s_analysis.tsv", go_output_dir, go_output_name), sep="\t", quote=FALSE, row.names=TRUE)
        ## sort and order enriched go terms
        ont <- "BP"
        sorting <- "Up"
        go_BP_up <- topGO(go, ontology=ont, sort=sorting, number=Inf)
        write.table(go_BP_up, file=sprintf("%s/%s_analysis_%s_sort_%s.tsv", go_output_dir, go_output_name, ont, sorting), sep="\t", quote=FALSE, row.names=TRUE)
        sorting <- "Down"
        go_BP_down <- topGO(go, ontology=ont, sort=sorting, number=Inf)
        write.table(go_BP_down, file=sprintf("%s/%s_analysis_%s_sort_%s.tsv", go_output_dir, go_output_name, ont, sorting), sep="\t", quote=FALSE, row.names=TRUE)
      }
    }
  }
}

DE_pval_scatter_plot <- function(de_results_df_1, de_results_df_2, common_df, specific_1_df, specific_2_df, name_1, name_2, out_dir, out_name) {
  # plot theme
  plot_theme <- theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=10)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8), strip.text=element_text(size=10))
  
  gene_ID_colname_1 <- sprintf("gene_ID_%s", gsub(" ", "", name_1))
  gene_ID_colname_2 <- sprintf("gene_ID_%s", gsub(" ", "", name_2))
  
  # DE gene p-value scatter plot: common and specific genes
  ## common DE genes
  gene_category <- "common"
  nb_common <- dim(common_df)[1]
  exp_matrix_1_common_DE_pval <- de_results_df_1[as.character(common_df[[gene_ID_colname_1]]), "FDR"]
  exp_matrix_1_common_FC <- de_results_df_1[as.character(common_df[[gene_ID_colname_1]]), "logFC"]
  exp_matrix_2_common_DE_pval <- de_results_df_2[as.character(common_df[[gene_ID_colname_2]]), "FDR"]
  exp_matrix_2_common_FC <- de_results_df_2[as.character(common_df[[gene_ID_colname_2]]), "logFC"]
  de_pval_df2ggplot <- data.frame(gene_id_1=as.character(common_df[[gene_ID_colname_1]]), gene_id_2=as.character(common_df[[gene_ID_colname_2]]), de_pval_1=exp_matrix_1_common_DE_pval, fc_1=exp_matrix_1_common_FC, de_pval_2=exp_matrix_2_common_DE_pval, fc_2=exp_matrix_2_common_FC, category=rep(gene_category, length(exp_matrix_1_common_DE_pval)))
  ## specific DE genes
  ### genome 1
  gene_category <- sprintf("%s specific", name_1)
  nb_exp_matrix_1_specific <- length(unique(as.character(specific_1_df[, gene_ID_colname_1])))
  exp_matrix_1_specific_DE_pval_1 <- de_results_df_1[as.character(specific_1_df[[gene_ID_colname_1]]), "FDR"]
  exp_matrix_1_specific_FC_1 <- de_results_df_1[as.character(specific_1_df[[gene_ID_colname_1]]), "logFC"]
  exp_matrix_1_specific_DE_pval_2 <- de_results_df_2[as.character(specific_1_df[[gene_ID_colname_2]]), "FDR"]
  exp_matrix_1_specific_FC_2 <- de_results_df_2[as.character(specific_1_df[[gene_ID_colname_2]]), "logFC"]
  de_pval_df2ggplot <- rbind(de_pval_df2ggplot, data.frame(gene_id_1=as.character(specific_1_df[[gene_ID_colname_1]]), gene_id_2=as.character(specific_1_df[[gene_ID_colname_2]]), de_pval_1=exp_matrix_1_specific_DE_pval_1, fc_1=exp_matrix_1_specific_FC_1, de_pval_2=exp_matrix_1_specific_DE_pval_2, fc_2=exp_matrix_1_specific_FC_2, category=rep(gene_category, length(exp_matrix_1_specific_DE_pval_1))))
  ### genome 2
  gene_category <- sprintf("%s specific", name_2)
  nb_exp_matrix_2_specific <- length(unique(as.character(specific_2_df[, gene_ID_colname_2])))
  exp_matrix_2_specific_DE_pval_1 <- de_results_df_1[as.character(specific_2_df[[gene_ID_colname_1]]), "FDR"]
  exp_matrix_2_specific_FC_1 <- de_results_df_1[as.character(specific_2_df[[gene_ID_colname_1]]), "logFC"]
  exp_matrix_2_specific_DE_pval_2 <- de_results_df_2[as.character(specific_2_df[[gene_ID_colname_2]]), "FDR"]
  exp_matrix_2_specific_FC_2 <- de_results_df_2[as.character(specific_2_df[[gene_ID_colname_2]]), "logFC"]
  de_pval_df2ggplot <- rbind(de_pval_df2ggplot, data.frame(gene_id_1=as.character(specific_2_df[[gene_ID_colname_1]]), gene_id_2=as.character(specific_2_df[[gene_ID_colname_2]]), de_pval_1=exp_matrix_2_specific_DE_pval_1, fc_1=exp_matrix_2_specific_FC_1, de_pval_2=exp_matrix_2_specific_DE_pval_2, fc_2=exp_matrix_2_specific_FC_2, category=rep(gene_category, length(exp_matrix_2_specific_DE_pval_1))))
  ## sign of -log10(de_pval) depending on the sign of FC
  de_pval_df2ggplot <- cbind(de_pval_df2ggplot, signed_de_pval_1=-log10(de_pval_df2ggplot$de_pval_1)*sign(de_pval_df2ggplot$fc_1), signed_de_pval_2=-log10(de_pval_df2ggplot$de_pval_2)*sign(de_pval_df2ggplot$fc_2))
  
  # output_name <- sprintf("%s_%s_%s_%s_FC_%s_DE_pval_per_specificity", dataset, comparison, design_output_name, one_contrast_output_name, sub("\\.", "_", one_fc_threshold))
  pdf(sprintf("%s/%s.pdf", out_dir, out_name), width=9, height=9)
  p <- ggplot(de_pval_df2ggplot, aes(x=-log10(de_pval_1) , y=-log10(de_pval_2), color=category)) +
    geom_point() +
    geom_vline(xintercept=0, alpha=0.5) + 
    geom_hline(yintercept=0, alpha=0.5) +
    geom_vline(xintercept=-log10(0.05), linetype="dashed", alpha=0.5) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
    labs(title=sprintf("DE p-values per specificity categories: %s and %s\n#common DE genes: %d\n#%s specific DE genes: %d\n#%s specific DE genes: %d", name_1, name_2, nb_common, name_1, nb_exp_matrix_1_specific, name_2, nb_exp_matrix_2_specific), x=sprintf("-log10(DE p-value): %s", name_1), y=sprintf("-log10(DE p-value): %s", name_2), color="specificity") +
    plot_theme +
    theme(legend.position="bottom")
  print(p)
  
  p <- ggplot(de_pval_df2ggplot, aes(x=signed_de_pval_1 , y=signed_de_pval_2, color=category)) +
    geom_point() +
    geom_vline(xintercept=0, alpha=0.5) + 
    geom_hline(yintercept=0, alpha=0.5) +
    geom_vline(xintercept=-log10(0.05), linetype="dashed", alpha=0.5) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", alpha=0.5) +
    geom_vline(xintercept=log10(0.05), linetype="dashed", alpha=0.5) + 
    geom_hline(yintercept=log10(0.05), linetype="dashed", alpha=0.5) +
    labs(title=sprintf("DE p-values per specificity categories: %s and %s\n#common DE genes: %d\n#%s specific DE genes: %d\n#%s specific DE genes: %d", name_1, name_2, nb_common, exp_matrix_1_name, nb_exp_matrix_1_specific, name_2, nb_exp_matrix_2_specific), x=sprintf("-log10(DE p-value) * sign(FC): %s", name_1), y=sprintf("-log10(DE p-value) * sign(FC): %s", name_2), color="specificity") +
    plot_theme +
    theme(legend.position="bottom")
  print(p)
  dev.off()
  colnames(de_pval_df2ggplot)[colnames(de_pval_df2ggplot)=="gene_id_1"] <- gene_ID_colname_1
  colnames(de_pval_df2ggplot)[colnames(de_pval_df2ggplot)=="gene_id_2"] <- gene_ID_colname_2
  colnames(de_pval_df2ggplot)[colnames(de_pval_df2ggplot)=="de_pval_1"] <- sprintf("de_pval_%s", gsub(" ", "", name_1))
  colnames(de_pval_df2ggplot)[colnames(de_pval_df2ggplot)=="de_pval_2"] <- sprintf("de_pval_%s", gsub(" ", "", name_2))
  write.csv(de_pval_df2ggplot, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
}





###############
# GO analysis #
###############

ID2Entrez_conversion_retained_and_discarded_ids <- function(ids, ids_df, id_type_col, entrez_col) {
  # IDs not in ids_df
  ids_not_conversion <- ids[! ids %in% ids_df[[id_type_col]]]
  ids_not_ids_df <- data.frame(ids=ids_not_conversion, entrez=rep(NA, length(ids_not_conversion)))
  colnames(ids_not_ids_df) <- c(id_type_col, entrez_col)
  # remove IDs without Entrez ID (NA values)
  ids_no_Entrez_df <- ids_df[which(is.na(ids_df[[entrez_col]])),]
  ids_no_Entrez_df <- rbind(ids_no_Entrez_df, ids_not_ids_df)
  ids_df_no_Entrez_nb <- dim(ids_no_Entrez_df)[1]
  # IDs with Entrez ID (NA values)
  ids_Entrez_df <- ids_df[which(! is.na(ids_df[[entrez_col]])),]
  # remove duplicated IDs
  duplicated_ids_df <- ids_Entrez_df[which(ids_Entrez_df[[entrez_col]] %in% ids_Entrez_df[[entrez_col]][duplicated(ids_Entrez_df[[entrez_col]])] | ids_Entrez_df[[id_type_col]] %in% ids_Entrez_df[[id_type_col]][duplicated(ids_Entrez_df[[id_type_col]])]),]
  duplicated_Entrez_ID_nb <- length(ids_Entrez_df[which(ids_Entrez_df[[entrez_col]] %in% ids_Entrez_df[[entrez_col]][duplicated(ids_Entrez_df[[entrez_col]])]), id_type_col])
  duplicated_ID_nb <- length(unique(ids_Entrez_df[which(ids_Entrez_df[[id_type_col]] %in% ids_Entrez_df[[id_type_col]][duplicated(ids_Entrez_df[[id_type_col]])]), id_type_col]))
  discarded_ids_df <- rbind(duplicated_ids_df, ids_no_Entrez_df)
  unique_Entrez_df <- ids_Entrez_df[which(! ids_Entrez_df[[entrez_col]] %in% ids_Entrez_df[[entrez_col]][duplicated(ids_Entrez_df[[entrez_col]])] & ! ids_Entrez_df[[id_type_col]] %in% ids_Entrez_df[[id_type_col]][duplicated(ids_Entrez_df[[id_type_col]])]),]
  unique_ID_Entrez_nb <- length(ids_Entrez_df[which(! ids_Entrez_df[[entrez_col]] %in% ids_Entrez_df[[entrez_col]][duplicated(ids_Entrez_df[[entrez_col]])] & ! ids_Entrez_df[[id_type_col]] %in% ids_Entrez_df[[id_type_col]][duplicated(ids_Entrez_df[[id_type_col]])]), id_type_col])
  # counts
  counts_df <- data.frame(category=c("unique ID-Entrez match", "multiple IDs with identical Entrez ID", "multiple Entrez IDs", "no Entrez ID"), count=c(unique_ID_Entrez_nb, duplicated_Entrez_ID_nb, duplicated_ID_nb, ids_df_no_Entrez_nb))
  
  return(list(retained=unique_Entrez_df, discarded=discarded_ids_df, counts=counts_df))
}

ID2Entrez_conversion <- function(ids, method, annotationdbi_id_type, biomart_id_type, diff_method, out_dir, out_name) {
  # objective: convert Ensembl or RefSeq IDs into Entrez IDs
  # both conversions using AnnotationDdi+org.Rn.eg.db and biomaRt packages are performed
  # parameters:
  ## main_id_type: type of IDs, which enables to identify the conversion method
  ## secondary_id_type: enables to identify the conversion to use for IDs with different Entrez IDs according to both methods
  
  plot_theme <- theme_bw() +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=8)) +
    theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=7), strip.text=element_text(size=8))
  
  # AnnotationDbi and org.Rn.eg.db conversion
  library(org.Rn.eg.db)
  ids_df <- select(org.Rn.eg.db, keys=ids, keytype=annotationdbi_id_type, column="ENTREZID")
  annotationdbi_conversion <- ID2Entrez_conversion_retained_and_discarded_ids(ids, ids_df, annotationdbi_id_type, "ENTREZID")
  annotationdbi_conversion_ids_df <- annotationdbi_conversion$retained
  write.csv(annotationdbi_conversion_ids_df, file=sprintf("%s/%s_AnnotationDbi_retained.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  annotationdbi_conversion_discarded_df <- annotationdbi_conversion$discarded
  write.csv(annotationdbi_conversion_discarded_df, file=sprintf("%s/%s_AnnotationDbi_discarded.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  # biomaRt conversion
  library(biomaRt)
  # ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", host="uswest.ensembl.org", ensemblRedirect = FALSE)
  ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", ensemblRedirect = FALSE)
  ids_df <- getBM(attributes=c(biomart_id_type, "entrezgene_id"), filters=biomart_id_type, values=ids, mart=ensembl)
  biomart_conversion <- ID2Entrez_conversion_retained_and_discarded_ids(ids, ids_df, biomart_id_type, "entrezgene_id")
  biomart_conversion_ids_df <- biomart_conversion$retained
  write.csv(biomart_conversion_ids_df, file=sprintf("%s/%s_biomaRt_retained.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  biomart_conversion_discarded_df <- biomart_conversion$discarded
  write.csv(biomart_conversion_discarded_df, file=sprintf("%s/%s_biomaRt_discarded.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  # IDs with different Entrez IDs according to AnnotationDbi and biomaRt
  conversion_merge_df <- merge(annotationdbi_conversion_ids_df, biomart_conversion_ids_df, by.x=annotationdbi_id_type, by.y=biomart_id_type)
  ids_different_Entrez <- conversion_merge_df[which(conversion_merge_df$ENTREZ != as.integer(conversion_merge_df$entrezgene_id)), annotationdbi_id_type]
  if (method == "AnnotationDbi") {
    conversion_ids_df <- annotationdbi_conversion_ids_df
    id_col <- annotationdbi_id_type
    entrez_col <- "ENTREZID"
    discarded_ids_df <- annotationdbi_conversion_discarded_df
    counts_df <- annotationdbi_conversion$counts
    ### add a column for conversion method
    conversion_ids_df <- cbind(conversion_ids_df, method=rep(method, dim(conversion_ids_df)[1]))
  } else {
    if (method == "biomaRt") {
      conversion_ids_df <- biomart_conversion_ids_df
      id_col <- biomart_id_type
      entrez_col <- "entrezgene_id"
      discarded_ids_df <- biomart_conversion_discarded_df
      counts_df <- biomart_conversion$counts
      ### add a column for conversion method
      conversion_ids_df <- cbind(conversion_ids_df, method=rep(method, dim(conversion_ids_df)[1]))
    }
  }
  ## replace Entrez IDs for IDs with different Entrez IDs according to AnnotationDbi and biomaRt
  if (diff_method != method) {
    conversion_ids_df$method <- as.character(conversion_ids_df$method) # to enable modifications in this column
    if (diff_method == "AnnotationDbi") {
      for (one_id in ids_different_Entrez) {
        conversion_ids_df[which(conversion_ids_df[[id_col]]==one_id), entrez_col] <- annotationdbi_conversion_ids_df[which(annotationdbi_conversion_ids_df[[annotationdbi_id_type]]==one_id), "ENTREZID"]
        conversion_ids_df[which(conversion_ids_df[[id_col]]==one_id), "method"] <- diff_method
      }
    } else {
      if (diff_method == "biomaRt") {
        for (one_id in ids_different_Entrez) {
          conversion_ids_df[which(conversion_ids_df[[id_col]]==one_id), entrez_col] <- biomart_conversion_ids_df[which(biomart_conversion_ids_df[[biomart_id_type]]==one_id), "entrezgene_id"]
          conversion_ids_df[which(conversion_ids_df[[id_col]]==one_id), "method"] <- diff_method
        }
      }
    }
    conversion_ids_df$method <- factor(conversion_ids_df$method)
  }
  final_conversion_counts_df <- data.frame(category=names(table(conversion_ids_df$method)), count=as.numeric(unname(table(conversion_ids_df$method))))
  final_conversion_counts_df <- rbind(final_conversion_counts_df, counts_df[which(counts_df$category %in% c("multiple IDs with identical Entrez ID", "multiple Entrez IDs", "no Entrez ID")),])
  final_conversion_counts_df <- droplevels(final_conversion_counts_df)
  final_conversion_counts_df$ID <- factor(1)
  final_conversion_counts_df$category <- factor(final_conversion_counts_df$category, levels=rev(final_conversion_counts_df$category))
  
  # plots
  df2ggplot <- cbind(annotationdbi_conversion$counts, conversion=rep("AnnotationDbi", dim(annotationdbi_conversion$counts)[1]))
  df2ggplot <- rbind(df2ggplot, cbind(biomart_conversion$counts, conversion=rep("biomaRt", dim(biomart_conversion$counts)[1])))
  pdf(sprintf("%s/%s.pdf", out_dir, out_name))
  conversion <- ggplot(df2ggplot, aes(x=conversion, y=count, fill=category)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%d", count)), size=4, position=position_stack(vjust=0.5)) +
    labs(x="Conversion", y="# IDs", fill="ID matching") +
    plot_theme +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(reverse=TRUE, nrow=length(levels(df2ggplot$category)), byrow=TRUE))
  # guides(fill=guide_legend(nrow=length(levels(df2ggplot$matching_type)), byrow=TRUE))
  
  final <- ggplot(final_conversion_counts_df, aes(x=ID, y=count, fill=category)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%d", count)), size=4, position=position_stack(vjust=0.5)) +
    labs(x="Final conversion", y="# IDs", fill="ID matching") +
    plot_theme +
    theme(legend.position="bottom") +
    guides(fill=guide_legend(reverse=TRUE, nrow=length(levels(final_conversion_counts_df$category)), byrow=TRUE)) +
    theme(axis.text.x=element_blank()) +
    theme(axis.ticks.x=element_blank())
  
  multiplot <- ggdraw() +
    draw_plot(conversion, 0, 0, 0.5, 0.98) +
    draw_plot(final, 0.5, 0, 0.5, 0.98) +
    draw_label("ID conversion into Entrez IDs", x=0.02, y=0.99, hjust=0, vjust=1, size=12)
  print(multiplot)
  dev.off()
  return(list(retained=conversion_ids_df, discarded=discarded_ids_df, counts=final_conversion_counts_df))
}

GO_term_analysis <- function(dge, annotationdbi_id_type, biomart_id_type, out_dir, out_name) {
  # convert IDs into Entrez IDs
  id_output_dir <- sprintf("%s/ID_conversion", out_dir)
  id_output_name <- sprintf("%s_Entrez_ID_conversion", out_name)
  if (! dir.exists(id_output_dir)) {
    dir.create(id_output_dir, recursive=TRUE, mode="0775")
  }
  if (annotationdbi_id_type == "SYMBOL") {
    conversion_method <- "AnnotationDbi"
    id_type <- annotationdbi_id_type
    entrez_ID_name <- "ENTREZID"
  } else {
    if (annotationdbi_id_type == "ENSEMBL") {
      conversion_method <- "biomaRt"
      id_type <- biomart_id_type
      entrez_ID_name <- "entrezgene_id"
    }
  }
  id_conversion <- ID2Entrez_conversion(rownames(dge$table), conversion_method, annotationdbi_id_type, biomart_id_type, "biomaRt", id_output_dir, id_output_name)
  write.csv(id_conversion$retained, file=sprintf("%s/%s_retained.csv", id_output_dir, id_output_name), quote=FALSE, row.names=FALSE)
  write.csv(id_conversion$discarded, file=sprintf("%s/%s_discarded.csv", id_output_dir, id_output_name), quote=FALSE, row.names=FALSE)
  write.csv(id_conversion$counts, file=sprintf("%s/%s_counts.csv", id_output_dir, id_output_name), quote=FALSE, row.names=FALSE)
  
  # # get corresponding Entrez IDs
  # if (id_type == "SYMBOL") {
  #   ## NCBI gene accessions: use AnnotationDbi and org.Rn.eg.db packages to get Entrez IDs
  #   entrez_ID_name <- "ENTREZID"
  #   ids_df <- select(org.Rn.eg.db, keys=rownames(dge$table), keytype=id_type, column=entrez_ID_name)
  #   ### remove duplicated Entrez IDs (first duplicate is kept) and NA values
  #   gene_ids <- ids_df[!duplicated(ids_df[[entrez_ID_name]]) & !is.na(ids_df[[entrez_ID_name]]),]
  #   write.csv(gene_ids, file=sprintf("%s/%s_genes.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  #   discarded_genes <- ids_df[duplicated(ids_df[[entrez_ID_name]]) | is.na(ids_df[[entrez_ID_name]]),]
  #   write.csv(discarded_genes, file=sprintf("%s/%s_discarded_genes.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  # } else {
  #   if (id_type == "ensembl_gene_id") {
  #     ## Ensembl gene IDs: use biomaRt package to get Entrez IDs
  #     entrez_ID_name <- "entrezgene_id"
  #     ensembl <- useMart("ensembl", dataset="rnorvegicus_gene_ensembl", ensemblRedirect = FALSE)
  #     ids_df <- getBM(attributes=c(id_type, "entrezgene_accession", entrez_ID_name), filters=id_type, values=rownames(dge$table), mart=ensembl)
  #     ### remove duplicated Entrez IDs (first duplicate is kept) and NA values
  #     gene_ids <- ids_df[!duplicated(ids_df[[entrez_ID_name]]) & !is.na(ids_df[[entrez_ID_name]]),]
  #     discarded_genes <- ids_df[duplicated(ids_df[[entrez_ID_name]]) | is.na(ids_df[[entrez_ID_name]]),]
  #   }
  # }
  # write.csv(gene_ids, file=sprintf("%s/%s_genes.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  # write.csv(discarded_genes, file=sprintf("%s/%s_discarded_genes.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  ## only retain genes with Entrez ID
  dge <- dge[id_conversion$retained[[id_type]],]
  dge$genes <- id_conversion$retained
  write.csv(dge$genes, file=sprintf("%s/%s_gene_IDs.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  write.csv(topTags(dge, n=dim(dge$table)[1])$table, file=sprintf("%s/%s_pval_entrez_ID.csv", out_dir, out_name), quote=FALSE)
  # GO term enrichment analysis
  go <- goana(dge, geneid=entrez_ID_name, species="Rn")
  # corrections for multiple testing
  go <- cbind(go, P.Up.BH=p.adjust(go$P.Up, method="BH"), P.Up.BY=p.adjust(go$P.Up, method="BY"), P.Down.BH=p.adjust(go$P.Down, method="BH"), P.Down.BY=p.adjust(go$P.Down, method="BY"))
  return(go)
}

feature_matching <- function(features_1, features_2) {
  # vectors for data frame to return
  feature_vector <- category_vector <- c()
  ## common GO IDs
  common_features <- features_1[features_1 %in% features_2]
  one_category <- "common"
  feature_vector <- c(feature_vector, common_features)
  category_vector <- c(category_vector, rep(one_category, length(common_features)))
  ## specific GO IDs in features_1
  features_1_specific <- features_1[! features_1 %in% features_2]
  one_category <- "specific"
  feature_vector <- c(feature_vector, features_1_specific)
  category_vector <- c(category_vector, rep(one_category, length(features_1_specific)))
  # return data frame   
  df2return <- data.frame(feature=feature_vector, category=category_vector)
  return(df2return)
}

feature_matching_plot <- function(match_1_df, match_2_df, name_1, name_2, plot_title, out_dir, out_name) {
  # count GO IDs per category
  categories <- levels(match_1_df$category)
  value_vector <- annotation_vector <- c()
  for (one_category in categories) {
    value_vector <- c(value_vector, length(match_1_df[which(match_1_df$category==one_category),"feature"]))
    annotation_vector <- c(annotation_vector, name_1)
  }
  one_df2barplot <- data.frame(annotation=annotation_vector, category=categories, value=value_vector, pct=(value_vector/length(match_1_df$feature))*100)
  df2barplot <- one_df2barplot
  
  categories <- levels(match_2_df$category)
  value_vector <- annotation_vector <- c()
  for (one_category in categories) {
    value_vector <- c(value_vector, length(match_2_df[which(match_2_df$category==one_category),"feature"]))
    annotation_vector <- c(annotation_vector, name_2)
  }
  one_df2barplot <- data.frame(annotation=annotation_vector, category=categories, value=value_vector, pct=(value_vector/length(match_2_df$feature))*100)
  df2barplot <- rbind(df2barplot, one_df2barplot)
  write.csv(df2barplot, file=sprintf("%s/%s.csv", out_dir, out_name), quote=FALSE, row.names=FALSE)
  
  # plot
  pdf(sprintf("%s/%s.pdf", out_dir, out_name))
  p <- ggplot(df2barplot, aes(x=annotation, y=value, fill=category)) +
    geom_bar(stat="identity") +
    geom_text(aes(label=sprintf("%.2f%%", pct)), size=4, position=position_stack(vjust=0.5)) +
    labs(title=plot_title, x="Annotation", y="# IDs", fill="ID category") +
    theme_bw() +
    theme(axis.title.x=element_text(size=12), axis.title.y=element_text(size=12), legend.title=element_text(size=12)) +
    theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=8))
  print(p)
  dev.off()
}

get_GO_analysis_specific_results <- function(match_df, go_df, ordering_col) {
  # get GO terms with 'specific' category
  specific_terms <- as.character(match_df[which(match_df$category=="specific"), "feature"])
  go_specific_terms <- go_df[specific_terms, ]
  return(go_specific_terms[order(go_specific_terms[[ordering_col]]),])
}

GO2EntrezID <- function(GO_ID) {
  # get gene IDs associated to a GO ID
  library(org.Rn.eg.db)
  x <- org.Rn.egGO2ALLEGS
  Rkeys(x) <- GO_ID
  EG <- mappedLkeys(x)
  return(EG)
}

GO_terms_per_Entrez_IDs <- function(entrez_ids, evidence_df) {
  library(org.Rn.eg.db)
  x <- org.Rn.egGO
  # Get the entrez gene identifiers that are mapped to a GO ID
  # mapped_genes <- mappedkeys(x)
  # Convert to a list
  # xx <- as.list(x[c(annotationdbi_specific_Entrez_ID, different_matching_annotationdbi_Entrez_ID)])
  xx <- as.list(x[entrez_ids[entrez_ids %in% keys(x)]])
  GO_ID_paste_vector <- GO_ID_nb_vector <- c()
  entrez_ID_vector <- GO_ID_vector <- evidence_vector <- c()
  if(length(xx) > 0) {
    for (entrez_ID in names(xx)) {
      got <- xx[[entrez_ID]]
      if (! is.na(got)) {
        GO_IDs <- unname(unlist(lapply(got, function(x) { return(x$GOID) })))
        evidence <- unname(unlist(lapply(got, function(x) { return(x$Evidence) })))
        GO_ID_paste_vector <- c(GO_ID_paste_vector, paste(GO_IDs, collapse="/"))
        GO_ID_nb_vector <- c(GO_ID_nb_vector, length(GO_IDs))
        GO_ID_vector <- c(GO_ID_vector, GO_IDs)
        evidence_vector <- c(evidence_vector, evidence)
        entrez_ID_vector <- c(entrez_ID_vector, rep(entrez_ID, length(GO_IDs)))
      } else {
        GO_ID_paste_vector <- c(GO_ID_paste_vector, NA)
        GO_ID_nb_vector <- c(GO_ID_nb_vector, 0)
        GO_ID_vector <- c(GO_ID_vector, NA)
        evidence_vector <- c(evidence_vector, NA)
        entrez_ID_vector <- c(entrez_ID_vector, entrez_ID)
      }
    }
  }
  entrez2GO_df <- data.frame(Entrez_ID=entrez_ID_vector, GO_ID=GO_ID_vector, evidence=evidence_vector)
  ## add IDs not in keys(x)
  ids_not_in_keys <- entrez_ids[! entrez_ids %in% keys(x)]
  ids_not_in_keys_nb <- length(ids_not_in_keys)
  if (ids_not_in_keys_nb > 0) {
    entrez2GO_df_not_in_keys <- data.frame(Entrez_ID=ids_not_in_keys, GO_ID=rep(NA, ids_not_in_keys_nb), evidence=rep(NA, ids_not_in_keys_nb))
    entrez2GO_df <- rbind(entrez2GO_df, entrez2GO_df_not_in_keys)
  }
  return(entrez2GO_df)
}

