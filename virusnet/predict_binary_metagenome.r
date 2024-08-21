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

custom_log("INFO", "Using metagenomic mode", "3/8")
result <- tryCatch({
  sink("/dev/null")
  pred <- predict_model(
    output_format = "one_pred_per_entry",
    model = model_binary,
    layer_name = "dense_26",
    path_input = command_args$input,
    round_digits = 3,
    step = as.numeric(command_args$step_size),
    batch_size = as.numeric(command_args$batch_size),
    verbose = TRUE,
    return_states = FALSE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = temp_file_binary,
    return_int = length(model_binary$input_shape) == 2
  )
  sink()
  success <<- TRUE
  custom_log("INFO", "Successful prediction", "3/8")
}, error = function(e) {
  custom_log("WARN", paste("Error during prediction: ", e$message), "3/8")
})

if (!success) {
  custom_log("ERROR", "Failed to predict. Please check your input data and model.", "3/8")
  stop("Prediction failed.")
}

# Interpret the output file
prediction_df <- load_prediction(h5_path = temp_file_binary, get_sample_position = FALSE, verbose = FALSE)
prediction_df <- as.data.frame(prediction_df$states)

colnames(prediction_df) <- c("non_virus", "virus")
prediction_df$contig_name <- fasta_data$Header

# Output summary using the logging system
custom_log("INFO", paste("Number of contigs classified as viral:", length(which(prediction_df$virus >= 0.5))), "3/8")
custom_log("INFO", paste("Number of contigs classified as non-viral:", length(which(!prediction_df$virus >= 0.5))), "3/8")
non_summarized_output_path <- file.path(output_directory, "binary_results.csv")
write.csv(prediction_df, non_summarized_output_path, row.names = FALSE, quote = FALSE)

# Subset and write FASTA data for viral and non-viral contigs
viral_contigs <- prediction_df$contig_name[prediction_df$virus >= 0.5]
non_viral_contigs <- prediction_df$contig_name[prediction_df$virus < 0.5]

viral_fasta <- if (length(viral_contigs) > 0) fasta_data[fasta_data$Header %in% viral_contigs, ] else NULL
safe_write_fasta(viral_fasta, file.path(output_directory, "viral_contigs.fasta"))

non_viral_fasta <- if (length(non_viral_contigs) > 0) fasta_data[fasta_data$Header %in% non_viral_contigs, ] else NULL
safe_write_fasta(non_viral_fasta, file.path(output_directory, "non_viral_contigs.fasta"))

on.exit({
  unlink(temp_file_binary)
}, add = TRUE)