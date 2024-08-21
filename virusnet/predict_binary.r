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
# Source the R scripts
invisible(source(file.path(conda_prefix, "bin", "utils.r")))
invisible(source(file.path(conda_prefix, "bin", "setup_logger.r")))

# Function to safely write FASTA if not empty
safe_write_fasta <- function(fasta_subset, file_path) {
  if (!is.null(fasta_subset) && nrow(fasta_subset) > 0) {
    writeFasta(fasta_subset, file_path)
    custom_log("INFO", paste("FASTA data written to:", file_path), "3/8")
  } else {
    custom_log("WARN", "No contigs found", "3/8")
  }
}

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

custom_log("INFO", "Loading binary model", "2/8")

suppressWarnings({
  model_binary <- keras::load_model_hdf5(command_args$model_binary, custom_objects = custom_objects)
})

custom_log("INFO", "Performing predictions", "3/8")
temp_file_binary <- tempfile(fileext = ".h5")

custom_log("INFO", "Using non-metagenomic mode", "3/8")
pred <- predict_model(
  output_format = "by_entry_one_file",
  model = model_binary,
  layer_name = "dense_26",
  path_input = command_args$input,
  round_digits = 3,
  step = as.numeric(command_args$step_size),
  batch_size = as.numeric(command_args$batch_size),
  verbose = FALSE,
  return_states = FALSE,
  padding = "standard",
  mode = "label",
  format = "fasta",
  filename = temp_file_binary,
  return_int = length(model_binary$input_shape) == 2
)

# Interpret the output file
prediction_df <- load_prediction(h5_path = temp_file_binary, get_sample_position = TRUE, verbose = FALSE)
number_of_processed_contigs <- length(prediction_df)
contig_names <- fasta_data$Header

all_contig_states <- lapply(1:number_of_processed_contigs, function(i) {
  contig_states <- as.data.frame(prediction_df[[i]]$states)
  contig_states$sample_end_position <- prediction_df[[i]]$sample_end_position
  colnames(contig_states) <- c("is_non_virus", "is_virus", "position")
  contig_states$contig <- rep(i, nrow(contig_states))
  contig_states
})

# Combine all data frames into one
final_df <- do.call(rbind, all_contig_states)
final_df$contig_name <- contig_names[final_df$contig]

# Summarize data for each contig_name
contig_summary <- final_df %>%
  group_by(contig_name) %>%
  summarise(
    mean_is_virus = mean(is_virus),
    median_is_virus = median(is_virus),
    iqr_is_virus = IQR(is_virus),
    sd_is_virus = sd(is_virus),
    number_of_entries = n(),
    is_virus_binary = mean(is_virus) >= 0.5
  ) %>%
  ungroup()

contig_summary_df <- as.data.frame(contig_summary)

# Define file paths for the summarized and non-summarized data
summarized_output_path <- file.path(output_directory, "binary_results_summarized.csv")
non_summarized_output_path <- file.path(output_directory, "binary_results.csv")

# Write the summarized and non-summarized data frames to CSV
write.csv(contig_summary_df, summarized_output_path, row.names = FALSE, quote = FALSE)
write.csv(final_df, non_summarized_output_path, row.names = FALSE, quote = FALSE)

custom_log("INFO", paste("Summarized results written to:", summarized_output_path), "3/8")
custom_log("INFO", paste("Non-summarized results written to:", non_summarized_output_path), "3/8")

# Open a PDF device to save the plots
pdf(paste0(output_directory, "/binary_results.pdf"), width = 20, height = 10)

# Generate a plot for each unique contig
unique_contigs <- unique(final_df$contig_name)
for (contig_name in unique_contigs) {
  contig_data <- final_df[final_df$contig_name == contig_name, ]
  p <- ggplot(contig_data, aes(x = position, y = is_virus)) +
    geom_line() +
    geom_point(size = 1, color = "green") +
    ggtitle(paste("Contig:", contig_name)) +
    xlab("Position") +
    ylab("Is Non-Virus Probability") +
    theme_classic() +
    ylim(0, 1) +
    geom_hline(yintercept = 0.5, linetype = 2)
  print(p)
}

invisible(dev.off())

# Filter contig IDs based on is_virus_binary
viral_contigs <- contig_summary_df %>%
  filter(is_virus_binary == TRUE) %>%
  pull(contig_name)

non_viral_contigs <- contig_summary_df %>%
  filter(is_virus_binary == FALSE) %>%
  pull(contig_name)

# Output summary using the logging system
custom_log("INFO", paste("Number of contigs classified as viral:", length(viral_contigs)), "3/8")
custom_log("INFO", paste("Number of contigs classified as non-viral:", length(non_viral_contigs)), "3/8")

# Subset and write FASTA data for viral and non-viral contigs
viral_fasta <- if (length(viral_contigs) > 0) fasta_data[fasta_data$Header %in% viral_contigs, ] else NULL
safe_write_fasta(viral_fasta, file.path(output_directory, "viral_contigs.fasta"))

non_viral_fasta <- if (length(non_viral_contigs) > 0) fasta_data[fasta_data$Header %in% non_viral_contigs, ] else NULL
safe_write_fasta(non_viral_fasta, file.path(output_directory, "non_viral_contigs.fasta"))

on.exit({
  unlink(temp_file_binary)
}, add = TRUE)