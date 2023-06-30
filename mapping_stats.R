library(ggplot2)
library(ggforce)
library(RColorBrewer)

##############
# Parameters #
##############
work_dir <- "./"
# dataset ID
# dataset <- "GSE101798"
# dataset <- "GSE136913"
# dataset <- "GSE140420"
dataset <- "GSE179101"

# plot colors
plot_colors <- c("Ensembl 104" = brewer.pal(12, "Paired")[5], 
                 "Ensembl 105" = brewer.pal(12, "Paired")[6], 
                 "RefSeq 106" = brewer.pal(12, "Paired")[1], 
                 "RefSeq 108" = brewer.pal(12, "Paired")[2])

#############
# Constants #
#############
input_dir <- sprintf("%s/20-Expression_analysis/input/%s", work_dir, dataset)
output_dir <- sprintf("%s/20-Expression_analysis/output/00-Mapping_stats/%s", work_dir, dataset)
mapping_dir <- sprintf("%s/projects/20-R_norvegicus_ref_genome_eval/10-Read_mapping/output", work_dir)

# reference genomes
genomes <- c("NCBIRefSeq106_Ensembl104GTF", "NCBIRefSeq106_NCBIRefSeq106GTF", "NCBIRefSeq108_NCBIRefSeq108GTF", "NCBIRefSeq108_Ensembl105GTF")


sample_vector <- genome_vector <- stat_vector <- value_vector <- c()
stat_names_vector <- c("total reads", "mapped reads", "total splices", "annotated splices", "multi mapped reads", "multi mapped reads too many", "unmapped reads mismatches", "unmapped reads short", "unmapped reads other", "chimeric reads", "alignment not unique unannotated reads", "ambiguous unannotated reads", "no feature unannotated reads", "too low aQual unannotated reads", "annotated reads")
for (one_ref_genome in genomes) {
  print(one_ref_genome)
  
  # load metadata table
  load(sprintf("%s/%s_%s_exp_matrix.RData", input_dir, dataset, one_ref_genome))
  samples <- rownames(metadata_df)

  for (one_sample in samples) {
    # read mapping stats file
    mapping_stats_file <- sprintf("%s/%s/%s/%s/%s_%s_Log.final.out", mapping_dir, dataset, one_sample, one_ref_genome, one_sample, one_ref_genome)
    mapping_stats_df <- read.table(mapping_stats_file, sep="\t", quote="", row.names=NULL, col.names=c("stat", "value"), fill=TRUE)
    # get statistics
    total_reads <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of input reads", mapping_stats_df$stat)), "value"]))
    mapped_reads <- as.numeric(as.character(mapping_stats_df[which(grepl("Uniquely mapped reads number", mapping_stats_df$stat)), "value"]))
    sj_total <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of splices: Total", mapping_stats_df$stat)), "value"]))
    sj_annotated <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of splices: Annotated \\(sjdb\\)", mapping_stats_df$stat)), "value"]))
    multi_mapped_reads <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of reads mapped to multiple loci", mapping_stats_df$stat)), "value"]))
    multi_mapped_reads_too_many <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of reads mapped to too many loci", mapping_stats_df$stat)), "value"]))
    unmapped_reads_mismatches <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of reads unmapped: too many mismatches", mapping_stats_df$stat)), "value"]))
    unmapped_reads_short <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of reads unmapped: too short", mapping_stats_df$stat)), "value"]))
    unmapped_reads_other <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of reads unmapped: other", mapping_stats_df$stat)), "value"]))
    chimeric_reads <- as.numeric(as.character(mapping_stats_df[which(grepl("Number of chimeric reads", mapping_stats_df$stat)), "value"]))
    ## annotation
    alignment_not_unique_unannotated_reads <- exp_matrix["__alignment_not_unique", one_sample]
    ambiguous_unannotated_reads <- exp_matrix["__ambiguous", one_sample]
    no_feature_unannotated_reads <- exp_matrix["__no_feature", one_sample]
    too_low_aQual_unannotated_reads <- exp_matrix["__too_low_aQual", one_sample]
    annotated_reads <- sum(exp_matrix[! grepl("^_", rownames(exp_matrix)), one_sample])
    
    # fill vector for all samples statistics data frame
    one_sample_values <- c(total_reads, mapped_reads, sj_total, sj_annotated, multi_mapped_reads, multi_mapped_reads_too_many, unmapped_reads_mismatches, unmapped_reads_short, unmapped_reads_other, chimeric_reads, alignment_not_unique_unannotated_reads, ambiguous_unannotated_reads, no_feature_unannotated_reads, too_low_aQual_unannotated_reads, annotated_reads)
    stat_vector <- c(stat_vector, stat_names_vector)
    value_vector <- c(value_vector, one_sample_values)
    sample_vector <- c(sample_vector, rep(one_sample, length(one_sample_values)))
    genome_vector <- c(genome_vector, rep(one_ref_genome, length(one_sample_values)))
  }
}
all_samples_mapping_stats_df <- data.frame(sample=sample_vector, genome=genome_vector, stat=stat_vector, value=value_vector)
write.csv(all_samples_mapping_stats_df, file=sprintf("%s/%s_mapping_annotation_statistics.csv", output_dir, dataset), quote=FALSE, row.names=FALSE)


# compute percentages
sample_vector <- genome_vector <- stat_vector <- pct_vector <- c() 
for (one_ref_genome in levels(all_samples_mapping_stats_df$genome)) {
  print(one_ref_genome)
  one_ref_genome_stats_df <- all_samples_mapping_stats_df[which(all_samples_mapping_stats_df$genome==one_ref_genome),]
  one_ref_genome_stats_df <- droplevels(one_ref_genome_stats_df) # if there are different number of samples per reference genome
  for (one_sample in levels(one_ref_genome_stats_df$sample)) {
    one_sample_stats_df <- one_ref_genome_stats_df[which(one_ref_genome_stats_df$sample==one_sample),]
    for (one_stat in stat_names_vector) {
      pct <- one_sample_stats_df[which(one_sample_stats_df$stat==one_stat), "value"] / one_sample_stats_df[which(one_sample_stats_df$stat=="total reads"), "value"] * 100
      pct_vector <- c(pct_vector, pct)
    }
    stat_vector <- c(stat_vector, stat_names_vector)
    sample_vector <- c(sample_vector, rep(one_sample, length(stat_names_vector)))
    genome_vector <- c(genome_vector, rep(one_ref_genome, length(stat_names_vector)))
  }
}
all_samples_mapping_pct_df <- data.frame(sample=sample_vector, genome=genome_vector, stat=stat_vector, pct=pct_vector)
write.csv(all_samples_mapping_pct_df, file=sprintf("%s/%s_mapping_annotation_statistics_pct.csv", output_dir, dataset), quote=FALSE, row.names=FALSE)


# plots
## sanity check barplots
sanity_check_barplot <- function(df, input_stats, output_stats, title_string) {
  
  for (one_ref_genome in levels(df$genome)) {
    print(one_ref_genome)
    one_ref_genome_df <- df[which(df$genome==one_ref_genome),]
    
    sample_vector <- category_vector <- stat_vector <- value_vector <- c()
    labeller_vector <- c()
    for (one_sample in levels(one_ref_genome_df$sample)) {
      one_sample_df <- one_ref_genome_df[which(one_ref_genome_df$sample==one_sample),]
      category_vector <- c(category_vector, c(rep("input categories", length(input_stats)), rep("output categories", length(output_stats))))
      value_vector <- c(value_vector, one_sample_df[which(one_sample_df$stat %in% input_stats), "value"])
      stat_vector <- c(stat_vector, input_stats)
      value_vector <- c(value_vector, one_sample_df[which(one_sample_df$stat %in% output_stats), "value"])
      stat_vector <- c(stat_vector, output_stats)
      sample_vector <- c(sample_vector, rep(one_sample, length(input_stats) + length(output_stats)))
      
      input_categories_sum <- sum(one_sample_df[which(one_sample_df$stat %in% input_stats), "value"])
      output_categories_sum <- sum(one_sample_df[which(one_sample_df$stat %in% output_stats), "value"])
      sample_label <- sprintf("%s\ninput sum: %d\noutput sum: %d", one_sample, input_categories_sum, output_categories_sum)
      labeller_vector <- c(labeller_vector, sample_label)
    }
    df2barplot <- data.frame(sample=sample_vector, category=category_vector, stat=stat_vector, value=value_vector)
    names(labeller_vector) <- levels(df$sample)
    
    for (i in 1:ceiling(length(levels(df$sample))/10)) {
      p <- ggplot(df2barplot, aes(x=category, y=value, fill=stat)) +
        facet_wrap_paginate(~sample, scales="free_y", labeller=as_labeller(labeller_vector), nrow=2, ncol=5, page=i) +
        geom_bar(stat="identity") +
        geom_text(aes(label=sprintf("%d", value)), size=2, position=position_stack(vjust=0.5)) +
        labs(title=sprintf("%s - %s", title_string, one_ref_genome), x="Read category", y="Count", fill="Read category") +
        theme_bw() +
        theme(legend.position="bottom") +
        theme(plot.title=element_text(size=12), axis.title.x=element_text(size=10), axis.title.y=element_text(size=10), legend.title=element_text(size=10)) +
        theme(axis.text.x=element_text(size=8), axis.text.y=element_text(size=8), legend.text=element_text(size=8), strip.text.x=element_text(size=8)) +
        theme(axis.text.x=element_text(angle=15, hjust=1))
      print(p)
    }
  }
}
## mapping stats boxplots
mapping_stats_plot <- function(stats_df, pct_df, stat_names) {
  
  for (one_stat in stat_names) {
    ## values
    p <- ggplot(stats_df[which(stats_df$stat==one_stat),], aes(x=genome, y=value, fill=genome)) +
      geom_boxplot() +
      labs(title=sprintf("%s", one_stat), x="Reference genome", y="Value", fill="Reference genome")+
      theme_bw() +
      theme(panel.border=element_rect(color="grey50")) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=8)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=6))
    print(p)
    
    p <- ggplot(stats_df[which(stats_df$stat==one_stat),], aes(x=sample, y=value, fill=genome)) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs(title=sprintf("%s", one_stat), x="Sample", y="Value", fill="Reference genome")+
      theme_bw() +
      theme(panel.border=element_rect(color="grey50")) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=8)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=6))
    print(p)
    
    ## percentages
    p <- ggplot(pct_df[which(pct_df$stat==one_stat),], aes(x=genome, y=pct, fill=genome)) +
      geom_boxplot() +
      labs(title=sprintf("Percentage of %s", one_stat), x="Reference genome", y="Percentage", fill="Reference genome")+
      theme_bw() +
      theme(panel.border=element_rect(color="grey50")) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=8)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=6))
    print(p)
    
    p <- ggplot(pct_df[which(pct_df$stat==one_stat),], aes(x=sample, y=pct, fill=genome)) +
      geom_bar(stat="identity", position=position_dodge()) +
      labs(title=sprintf("Percentage of %s", one_stat), x="Sample", y="Percentage", fill="Reference genome")+
      theme_bw() +
      theme(panel.border=element_rect(color="grey50")) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) +
      theme(axis.title.x=element_text(size=14), axis.title.y=element_text(size=14), legend.title=element_text(size=8)) +
      theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), legend.text=element_text(size=6))
    print(p)
  }
  
  # Change factor levels for plot legend labels and facet titles
  df2ggplot <- pct_df[which(pct_df$stat %in% c("mapped reads", "multi mapped reads", "unmapped reads short", "annotated reads")),]
  df2ggplot <- droplevels(df2ggplot)
  levels(df2ggplot$genome) <- list("Ensembl 104" = "NCBIRefSeq106_Ensembl104GTF",
                                   "RefSeq 106" = "NCBIRefSeq106_NCBIRefSeq106GTF",
                                   "RefSeq 108" = "NCBIRefSeq108_NCBIRefSeq108GTF",
                                   "Ensembl 105" = "NCBIRefSeq108_Ensembl105GTF")
  levels(df2ggplot$stat) <- list("Mapped reads" = "mapped reads",
                                 "Multi-mapped reads" = "multi mapped reads",
                                 "Unmapped reads" = "unmapped reads short",
                                 "Annotated reads" = "annotated reads")
  
  p <- ggplot(df2ggplot, aes(x=genome, y=pct, fill=genome)) +
    geom_boxplot() +
    facet_wrap(~stat, nrow=1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y="Percentage", fill="Reference genome") +
    theme_bw() +
    theme(legend.position = "top") +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # no ticks on x-axis
    theme(axis.title.y=element_text(size=14), legend.title=element_text(size=14)) +
    theme(axis.text.y=element_text(size=12), legend.text=element_text(size=12))
  print(p)
  
  # Change factor levels for plot legend labels and facet titles
  df2ggplot <- pct_df[which(pct_df$stat %in% c("mapped reads", "multi mapped reads", "annotated reads")),]
  df2ggplot <- droplevels(df2ggplot)
  levels(df2ggplot$genome) <- list("Ensembl 104" = "NCBIRefSeq106_Ensembl104GTF",
                                   "RefSeq 106" = "NCBIRefSeq106_NCBIRefSeq106GTF",
                                   "RefSeq 108" = "NCBIRefSeq108_NCBIRefSeq108GTF",
                                   "Ensembl 105" = "NCBIRefSeq108_Ensembl105GTF")
  levels(df2ggplot$stat) <- list("Mapped reads" = "mapped reads",
                                 "Multi-mapped reads" = "multi mapped reads",
                                 "Annotated reads" = "annotated reads")
  
  p <- ggplot(df2ggplot, aes(x=genome, y=pct, fill=genome)) +
    geom_boxplot() +
    facet_wrap(~stat, nrow=1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y="Percentage", fill="Reference genome") +
    theme_bw() +
    theme(legend.position = "top") +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # no ticks on x-axis
    theme(axis.title.y=element_text(size=14), legend.title=element_text(size=14)) +
    theme(axis.text.y=element_text(size=12), legend.text=element_text(size=12))
  print(p)
  
  ## mapping stats
  df2ggplot <- pct_df[which(pct_df$stat %in% c("total reads", "mapped reads", "multi mapped reads", "multi mapped reads too many", "unmapped reads short", "unmapped reads other")),]
  df2ggplot <- droplevels(df2ggplot)
  levels(df2ggplot$genome) <- list("Ensembl 104" = "NCBIRefSeq106_Ensembl104GTF",
                                   "RefSeq 106" = "NCBIRefSeq106_NCBIRefSeq106GTF",
                                   "RefSeq 108" = "NCBIRefSeq108_NCBIRefSeq108GTF",
                                   "Ensembl 105" = "NCBIRefSeq108_Ensembl105GTF")
  levels(df2ggplot$stat) <- list("Total reads" = "total reads",
                                 "Mapped reads" = "mapped reads", 
                                 "Multi-mapped \nreads" = "multi mapped reads", 
                                 "Multi-mapped \nreads too many" = "multi mapped reads too many",
                                 "Unmapped \nreads short" = "unmapped reads short", 
                                 "Unmapped \nreads other" = "unmapped reads other")
  
  p <- ggplot(df2ggplot, aes(x=genome, y=pct, fill=genome)) +
    geom_boxplot() +
    facet_wrap(~stat, nrow=1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y="Percentage", fill="Reference genome") +
    theme_bw() +
    theme(legend.position = "top") +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # no ticks on x-axis
    theme(axis.title.y=element_text(size=14), legend.title=element_text(size=14)) +
    theme(axis.text.y=element_text(size=12), legend.text=element_text(size=12))
  print(p)
  
  ## annotation statistics
  df2ggplot <- pct_df[which(pct_df$stat %in% c("mapped reads", "multi mapped reads", "alignment not unique unannotated reads", "ambiguous unannotated reads", "no feature unannotated reads", "annotated reads")),]
  df2ggplot <- droplevels(df2ggplot)
  levels(df2ggplot$genome) <- list("Ensembl 104" = "NCBIRefSeq106_Ensembl104GTF",
                                   "RefSeq 106" = "NCBIRefSeq106_NCBIRefSeq106GTF",
                                   "RefSeq 108" = "NCBIRefSeq108_NCBIRefSeq108GTF",
                                   "Ensembl 105" = "NCBIRefSeq108_Ensembl105GTF")
  levels(df2ggplot$stat) <- list("Mapped reads" = "mapped reads",
                                 "Multi-mapped \nreads" = "multi mapped reads",
                                 "Annotated reads" = "annotated reads", 
                                 "Unannotated \naln not unique" = "alignment not unique unannotated reads",
                                 "Unannotated \nambiguous" = "ambiguous unannotated reads",
                                 "Unannotated \nno feature" = "no feature unannotated reads")
  
  p <- ggplot(df2ggplot, aes(x=genome, y=pct, fill=genome)) +
    geom_boxplot() +
    facet_wrap(~stat, nrow=1) +
    coord_cartesian(ylim = c(0, 100)) +
    labs(y="Percentage", fill="Reference genome") +
    theme_bw() +
    theme(legend.position = "top") +
    theme(panel.border=element_rect(color="grey50")) +
    theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # no ticks on x-axis
    theme(axis.title.y=element_text(size=14), legend.title=element_text(size=14)) +
    theme(axis.text.y=element_text(size=12), legend.text=element_text(size=12))
  print(p)
}

## all samples
### sanity check barplots
pdf(sprintf("%s/%s_mapping_annotation_sanity_checks_all_samples.pdf", output_dir, dataset), width=9)
print("sanity check barplots")
#### read mapping
print("read mapping")
input_stat_names_vector <- c("total reads")
output_stat_names_vector <- c("mapped reads", "multi mapped reads", "multi mapped reads too many", "unmapped reads mismatches", "unmapped reads short", "unmapped reads other", "chimeric reads")
sanity_check_barplot(all_samples_mapping_stats_df, input_stat_names_vector, output_stat_names_vector, "Read mapping")
#### annotation
print("annotation")
input_stat_names_vector <- c("mapped reads", "multi mapped reads")
output_stat_names_vector <- c("alignment not unique unannotated reads", "ambiguous unannotated reads", "no feature unannotated reads", "too low aQual unannotated reads", "annotated reads")
sanity_check_barplot(all_samples_mapping_stats_df, input_stat_names_vector, output_stat_names_vector, "Annotation")
dev.off()

### mapping stats boxplots
print("mapping stats boxplots")
pdf(sprintf("%s/%s_mapping_annotation_statistics_all_samples.pdf", output_dir, dataset))
print("all samples")
mapping_stats_plot(all_samples_mapping_stats_df, all_samples_mapping_pct_df, stat_names_vector)
dev.off()
#### per region
if (dataset == "GSE101798" | dataset == "GSE140420" | dataset == "GSE179101") {
  if (dataset == "GSE101798") {
    region_col_name <- "brain_tissue"
  } else {
    if (dataset == "GSE140420") {
      region_col_name <- "Region"
    } else {
      if (dataset == "GSE179101") {
        region_col_name <- "Subfield"
      }
    }
  }
  for (one_region in levels(metadata_df[[region_col_name]])) {
    print(one_region)
    pdf(sprintf("%s/%s_mapping_annotation_statistics_%s_samples.pdf", output_dir, dataset, one_region))
    one_region_samples <- rownames(metadata_df[which(metadata_df[[region_col_name]]==one_region),])
    one_region_samples_mapping_stats_df <- all_samples_mapping_stats_df[which(all_samples_mapping_stats_df$sample %in% one_region_samples),]
    one_region_samples_mapping_pct_df <- all_samples_mapping_pct_df[which(all_samples_mapping_pct_df$sample %in% one_region_samples),]
    mapping_stats_plot(one_region_samples_mapping_stats_df, one_region_samples_mapping_pct_df, stat_names_vector)
    dev.off()
  }
}




# mapping statistics figure
df2ggplot <- all_samples_mapping_pct_df[which(all_samples_mapping_pct_df$stat %in% c("mapped reads", "multi mapped reads", "annotated reads")),]
df2ggplot <- droplevels(df2ggplot)
levels(df2ggplot$genome) <- list("Ensembl 104" = "NCBIRefSeq106_Ensembl104GTF",
                                 "Ensembl 105" = "NCBIRefSeq108_Ensembl105GTF",
                                 "RefSeq 106" = "NCBIRefSeq106_NCBIRefSeq106GTF",
                                 "RefSeq 108" = "NCBIRefSeq108_NCBIRefSeq108GTF")
levels(df2ggplot$stat) <- list("Mapped reads" = "mapped reads",
                               "Multi-mapped reads" = "multi mapped reads",
                               "Annotated reads" = "annotated reads")
df2ggplot_mean <- aggregate(df2ggplot$pct, by=list(df2ggplot$genome, df2ggplot$stat), mean)
colnames(df2ggplot_mean) <- c("genome", "stat", "mean")
write.csv(df2ggplot_mean, file=sprintf("%s/%s_mapping_annotation_statistics_all_samples_figure_mean.csv", output_dir, dataset), quote=FALSE, row.names=FALSE)
pdf(sprintf("%s/%s_mapping_annotation_statistics_all_samples_figure.pdf", output_dir, dataset), width=9)
levels(df2ggplot$genome) <- list("Ensembl 104" = "NCBIRefSeq106_Ensembl104GTF",
                                 "RefSeq 106" = "NCBIRefSeq106_NCBIRefSeq106GTF", 
                                 "RefSeq 108" = "NCBIRefSeq108_NCBIRefSeq108GTF",
                                 "Ensembl 105" = "NCBIRefSeq108_Ensembl105GTF")
p <- ggplot(df2ggplot, aes(x=genome, y=pct, fill=genome)) +
  geom_boxplot() +
  scale_fill_manual(values=plot_colors) +
  facet_wrap(~stat, nrow=1) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(y="Percentage", fill="Reference\ngenome") +
  theme_bw() +
  theme(legend.position = "top") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE)) +
  theme(panel.border=element_rect(color="grey50")) +
  theme(axis.title.x=element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank()) + # no ticks on x-axis
  theme(axis.title.y=element_text(size=20), legend.title=element_text(size=18)) +
  theme(axis.text.y=element_text(size=18), legend.text=element_text(size=16), strip.text.x=element_text(size=18))
print(p)

dev.off()
