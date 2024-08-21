#!/usr/bin/env Rscript

# Suppress TF messages
Sys.setenv(TF_CPP_MIN_LOG_LEVEL = 3)

# Load libraries
suppressWarnings(suppressPackageStartupMessages({
  library(deepG)
  library(magrittr)
  library(microseq)
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(zoo)
  library(keras)
  library(reticulate)
  library(futile.logger)
}))


# Get the base directory of the Conda environment
conda_prefix <- Sys.getenv("CONDA_PREFIX")

#invisible(source("utils.r"))
#invisible(source("setup_logger.r"))

# Source the R scripts
invisible(source(file.path(conda_prefix, "bin", "utils.r")))
invisible(source(file.path(conda_prefix, "bin", "setup_logger.r")))

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
command_args <- list()
for (i in seq(1, length(args), 2)) {
  command_args[[gsub("--", "", args[i])]] <- args[i + 1]
}

# Check if the output directory exists, if not, create it
output_directory <- command_args$output
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

custom_log("INFO", "Checking input file", "1/8")
fasta_data <- readFasta(command_args$input)
num_fasta_entries <- nrow(fasta_data)

custom_log("INFO", paste("Number of FASTA entries in the file:", num_fasta_entries), "1/8")

if (num_fasta_entries == 0) {
  custom_log("ERROR", "Input FASTA file is empty. Please provide a non-empty FASTA file.", "1/8")
  stop("Empty input file.")
}

custom_log("INFO", "Loading genus model", "2/8")

# Suppress the h5py warning message
suppressWarnings({
  model_genus <- keras::load_model_hdf5(command_args$model_genus, compile = FALSE)
})

custom_log("INFO", "Performing genus predictions", "3/8")
temp_file_genus <- tempfile(fileext = ".h5")

pred <- predict_model(
  output_format = "by_entry_one_file",
  model = model_genus,
  layer_name = "dense_3",
  path_input = command_args$input,
  round_digits = 3,
  step = as.numeric(command_args$step_size),
  batch_size = as.numeric(command_args$batch_size),
  verbose = FALSE,
  return_states = FALSE,
  padding = "standard",
  mode = "label",
  format = "fasta",
  filename = temp_file_genus,
  return_int = length(model_genus$input_shape) == 2
)

genus_labels <- readRDS(command_args$labels)
prediction_df <- load_prediction(h5_path = temp_file_genus, get_sample_position = TRUE, verbose = FALSE)

number_of_processed_contigs <- length(prediction_df)
contig_names <- fasta_data$Header

all_contig_states <- lapply(1:number_of_processed_contigs, function(i) {
  contig_states <- as.data.frame(prediction_df[[i]]$states)
  contig_states$sample_end_position <- prediction_df[[i]]$sample_end_position
  contig_states$contig <- rep(i, nrow(contig_states))
  contig_states
})

# Combine all data frames into one
final_df <- do.call(rbind, all_contig_states)
colnames(final_df) <- c(genus_labels, "position", "contig_ID")
final_df$contig_name <- contig_names[final_df$contig]

# Summarize data for each contig_name
contig_summary <- final_df %>%
  group_by(contig_name) %>%
  summarise(across(1:(which(colnames(.) == "position") - 1), mean, .names = "mean_{.col}"), .groups = "drop")

# Calculate the virus with the highest mean probability for each contig
contig_summary <- contig_summary %>%
  rowwise() %>%
  mutate(max_index = which.max(c_across(starts_with("mean_"))),
         max_virus = names(.)[max_index + 1],  # +1 to adjust for the contig_name column
         max_probability = max(c_across(starts_with("mean_")))) %>%
  ungroup()
contig_summary <- as.data.frame(contig_summary)

df <- data.frame(contig_index = seq_along(contig_summary$contig_name),
                 max_virus = substring(contig_summary$max_virus, 6),
                 prop = round(contig_summary$max_probability, digits = 2) * 100)

cat("Contig Summary:\n")
for (i in 1:nrow(df)) {
  cat(sprintf("Contig %d: Virus - %s, Probability - %d%%\n", df$contig_index[i], df$max_virus[i], df$prop[i]))
}

write.csv(final_df, paste0(output_directory, "/genus_output.csv"), row.names = FALSE, quote = FALSE)
write.csv(contig_summary, paste0(output_directory, "/genus_summary_output.csv"), row.names = FALSE, quote = FALSE)
saveRDS(final_df, paste0(output_directory, "/genus_output.rds"))
saveRDS(contig_summary, paste0(output_directory, "/genus_summary_output.rds"))