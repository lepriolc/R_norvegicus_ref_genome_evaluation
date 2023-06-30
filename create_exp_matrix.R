##############
# Parameters #
##############
work_dir <- "./"
# dataset ID
# dataset <- "GSE101798"
# dataset <- "GSE136913"
# dataset <- "GSE140420"
dataset <- "GSE179101"

#############
# Constants #
#############
raw_data_dir <- sprintf("%s/00-Raw_data", work_dir)
input_dir <- sprintf("%s/10-Read_mapping/output/%s", work_dir, dataset)
output_dir <- sprintf("%s/20-Expression_analysis/input/%s", work_dir, dataset)
sra_run_table_file <- sprintf("%s/RNAseq/%s/SraRunTable.txt", raw_data_dir, dataset)
sra_result_file <- sprintf("%s/RNAseq/%s/sra_result.csv", raw_data_dir, dataset)

# reference genomes
ref_genomes <- c("NCBIRefSeq106_Ensembl104GTF", "NCBIRefSeq106_NCBIRefSeq106GTF", "NCBIRefSeq108_NCBIRefSeq108GTF", "NCBIRefSeq108_Ensembl105GTF")

############
# Metadata #
############
# SRA run table
sra_run_table_df <- read.csv(file=sra_run_table_file)

# SRA results
sra_result_df <- read.csv(file=sra_result_file)

# build metadata table
## GSE101798
build_metadata_df_GSE101798 <- function(sra_run_df) {
  meta_df <- sra_run_df[,c("Run", "AGE", "brain_tissue")]
  rownames(meta_df) <- as.character(meta_df$Run)
  return(meta_df)
}
## GSE136913
build_metadata_df_GSE136913 <- function(sra_run_df, sra_res_df) {
  sra_run_sub_df <- sra_run_df[,c("Run", "Sample.Name", "Experiment", "Instrument")]
  experiment_title <- sra_res_df$Experiment.Title
  condition <- ifelse(grepl("kainate", sapply(strsplit(as.character(experiment_title), ":"), "[", 3)), "kainate", "control")
  litter <- sapply(strsplit(sapply(strsplit(as.character(experiment_title), ":"), "[", 3), " ", 2), "[", 3)
  sra_res_sub_df <- data.frame(Experiment=sra_res_df$Experiment.Accession, Condition=condition, Litter=litter)
  meta_df <- merge(sra_run_sub_df, sra_res_sub_df, by="Experiment")
  rownames(meta_df) <- as.character(meta_df$Run)
  return(meta_df)
}
## GSE140420
build_metadata_df_GSE140420 <- function(sra_run_df, sra_res_df) {
  region <- sapply(strsplit(as.character(sra_run_df$tissue), " "), "[", 1)
  sra_run_sub_df <- data.frame(Run=sra_run_df$Run, Age=sra_run_df$AGE, Region=region, Experiment=sra_run_df$Experiment, Instrument=sra_run_df$Instrument)
  # get what seems to be individual information (needs confirmation)
  individual <- substr(sapply(strsplit(sapply(strsplit(as.character(sra_res_df$Experiment.Title), ";"), "[", 1), ": "), "[", 2), 1, 4)
  sra_res_sub_df <- data.frame(Experiment=sra_res_df$Experiment.Accession, Individual=individual)
  meta_df <- merge(sra_run_sub_df, sra_res_sub_df, by="Experiment")
  rownames(meta_df) <- as.character(sra_run_sub_df$Run)
  return(meta_df)
}

build_metadata_df_GSE179101 <- function(sra_run_df, sra_res_df) {
  sra_run_sub_df <- sra_run_df[,c("Run", "Experiment")]
  subfield <- sapply(strsplit(sapply(strsplit(as.character(sra_res_df$Experiment.Title), ": "), "[", 2), " subfield"), "[", 1)
  sra_res_sub_df <- data.frame(Experiment=sra_res_df$Experiment.Accession, Subfield=subfield)
  meta_df <- merge(sra_run_sub_df, sra_res_sub_df, by.x="Experiment")
  rownames(meta_df) <- as.character(meta_df$Run)
  return(meta_df)
}

if (dataset == "GSE101798") {
  dataset_metadata_df <- build_metadata_df_GSE101798(sra_run_table_df)
} else {
  if (dataset == "GSE136913") {
    dataset_metadata_df <- build_metadata_df_GSE136913(sra_run_table_df, sra_result_df)
  } else {
    if (dataset == "GSE140420") {
      dataset_metadata_df <- build_metadata_df_GSE140420(sra_run_table_df, sra_result_df)
    } else {
      if (dataset == "GSE179101") {
        dataset_metadata_df <- build_metadata_df_GSE179101(sra_run_table_df, sra_result_df)
      }
    }
  }
}


samples <- as.character(sra_run_table_df$Run)
for (one_ref_genome in ref_genomes) {
  print(one_ref_genome)
  # expression matrix
  ## initialize data frame
  one_sample <- samples[1]
  ## count filenames
  counts_file <- sprintf("%s/%s/%s/%s_%s_Aligned.sortedByCoord.out.counts.txt", input_dir, one_sample, one_ref_genome, one_sample, one_ref_genome)
  exp_matrix <- read.table(counts_file, header=FALSE, sep="\t", quote="", row.names=1)
  colnames(exp_matrix) <- one_sample
  
  ## add counts from the other samples
  for (one_sample in samples[-1]) {
    ### count filenames
    counts_file <- sprintf("%s/%s/%s/%s_%s_Aligned.sortedByCoord.out.counts.txt", input_dir, one_sample, one_ref_genome, one_sample, one_ref_genome)
    if (file.exists(counts_file)) {
      print(one_sample)
      one_sample_counts_df <- read.table(counts_file, header=FALSE, sep="\t", quote="", row.names=1)
      colnames(one_sample_counts_df) <- one_sample
      
      exp_matrix <- merge(exp_matrix, one_sample_counts_df, by="row.names", all.x=TRUE)
      rownames(exp_matrix) <- exp_matrix$Row.names
      exp_matrix <- subset(exp_matrix, select=-Row.names)  
    }
  }
  
  # metadata table
  metadata_df <- dataset_metadata_df[colnames(exp_matrix),]
  
  # save expression matrix and metadata tablein a RData file
  save(exp_matrix, metadata_df, file=sprintf("%s/%s_%s_exp_matrix.RData", output_dir, dataset, one_ref_genome))
}

