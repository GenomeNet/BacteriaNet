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


# generator args
maxlen_fragment <- 10000
maxlen <- 2000000
total_seq_len <- maxlen
step <- floor(maxlen/4)
samples_per_target <- maxlen/maxlen_fragment
sample_by_file_size <- FALSE
random_sampling <- TRUE
batch_size <- 10
concat_seq <- ""
vocabulary <- c("A", "C", "G", "T")
proportion_per_seq <- NULL
seed <- 33

pred_list <- vector('list', 1)


to_time_dist <- function(x, samples_per_target) {
  x_dim <- dim(x)
  x_dim_td <- c(x_dim[1], samples_per_target, x_dim[2]/samples_per_target, x_dim[3])
  x_td <- keras::k_reshape(x, shape = x_dim_td)
  keras::k_eval(x_td)
}

custom_log("INFO", "Checking input file", "1/8")

fasta_df <- readFasta(command_args$input)
nt_seq <- paste(fasta_df$Sequence, collapse = "")

 unpadded_seq_len <- nchar(nt_seq)
  
  if (unpadded_seq_len < total_seq_len) {
    nt_seq <- strrep(nt_seq, ceiling(total_seq_len / unpadded_seq_len))
    nt_seq <- substr(nt_seq, 1, total_seq_len)
  }



num_fasta_entries <- nrow(fasta_df)

custom_log("INFO", paste("Number of FASTA entries in the file:", num_fasta_entries), "1/8")

if (num_fasta_entries == 0) {
  custom_log("ERROR", "Input FASTA file is empty. Please provide a non-empty FASTA file.", "1/8")
  stop("Empty input file.")
}

custom_log("INFO", "Loading binary model", "2/8")

suppressWarnings({
  model_phenotypes <- keras::load_model_hdf5(command_args$model_phenotypes, custom_objects = custom_objects)
})

print(model_phenotypes)

custom_log("INFO", "Performing predictions", "3/8")
temp_file_binary <- tempfile(fileext = ".h5")



start_ind <- seq(1, nchar(nt_seq) - total_seq_len + 1, by = step)

x <- seq_encoding_label(maxlen = total_seq_len,
                        vocabulary = vocabulary,
                        start_ind = start_ind,
                        char_sequence = nt_seq)

x <- to_time_dist(x, samples_per_target = samples_per_target)
y_pred <- model_phenotypes$predict(x)

# Process predictions for a single file
target_split <- list(
  cellsize = c("rgStart_len", "rgEnd_len", "rgStart_wid", "rgEnd_wid"),
  cellshape = c("is_cell_shape_rod.shaped", "is_cell_shape_other", "is_cell_shape_coccus.shaped", "is_cell_shape_vibrio.shaped", "is_cell_shape_filament.shaped", "is_cell_shape_sphere.shaped", "is_cell_shape_ovoid.shaped", "is_cell_shape_pleomorphic.shaped", "is_cell_shape_spiral.shaped", "is_cell_shape_curved.shaped", "is_cell_shape_oval.shaped"),
  flagellum = c("is_flagellum_arrangement_monotrichous", "is_flagellum_arrangement_monotrichous_polar", "is_flagellum_arrangement_polar", "is_flagellum_arrangement_peritrichous", "is_flagellum_arrangement_lophotrichous", "is_flagellum_arrangement_gliding"),
  gram = c("is_gram_stain_positive", "is_gram_stain_variable", "is_gram_stain_negative"),
  motility = c("is_motile"),
  biosafety = c("biosafety_level"),
  pathogenicity_human = c("pathogenicity_human"),
  pathogenicity_animal = c("pathogenicity_animal"),
  pathogenicity_plant = c("pathogenicity_plant"),
  oxygen_growth = c("aerobe", "anaerobe"),
  oxygen_facultative = c("facultative.aerobe", "facultative.anaerobe"),
  oxygen_obligate = c("obligate.aerobe", "obligate.anaerobe"),
  oxygen_microaerophile = c("microaerophile"),
  spore = c("ability_spore")
)
target_names <- names(target_split)

column_means <- function(mat) {
  apply(mat, 2, mean)
}


get_pred <- function(l, target_split) {
  res_list <- list()
  
  label <- 'cellsize'
  res_list[[label]] <- l[[label]]
  
  label <- 'cellshape'
  m <- l[[label]]
  pred_index <- which.max(m)
  pred <- target_split[[label]][pred_index] %>% stringr::str_remove('is_cell_shape_')
  res_list[[label]] <- pred
  
  label <- 'flagellum'
  m <- l[[label]]
  pred_index <- which.max(m)
  pred <- target_split[[label]][pred_index] %>% stringr::str_remove('is_flagellum_arrangement_')
  res_list[[label]] <- pred
  
  label <- 'gram'
  m <- l[[label]]
  pred_index <- which.max(m)
  pred <- target_split[[label]][pred_index] %>% stringr::str_remove('is_gram_stain_')
  res_list[[label]] <- pred
  
  label <- 'motility'
  m <- l[[label]]
  pred <- ifelse(m > 0.5, TRUE, FALSE)
  res_list[['is_motile']] <- pred
  
  label <- 'biosafety'
  m <- l[[label]]
  if (m < 1.15) {
    pred <- 1
  } else if (m > 2.1) {
    pred <- 3
  } else {
    pred <- 2
  }
  res_list[[label]] <- pred
  
  label <- 'pathogenicity_human'
  m <- l[[label]]
  res_list[[label]] <- ifelse(m > 1.4, TRUE, FALSE)
  
  label <- 'pathogenicity_animal'
  m <- l[[label]]
  res_list[[label]] <- ifelse(m > 1.4, TRUE, FALSE)
  
  label <- 'pathogenicity_plant'
  m <- l[[label]]
  res_list[[label]] <- ifelse(m > 1, TRUE, FALSE)
  
  label <- 'oxygen_growth'
  m <- l[[label]]
  pred_index <- which.max(m)
  pred <- target_split[[label]][pred_index] 
  res_list[[label]] <- pred
  
  label <- 'oxygen_facultative'
  m <- l[[label]]
  pred_index <- which.max(m)
  pred <- target_split[[label]][pred_index] 
  res_list[[label]] <- pred
  
  label <- 'oxygen_obligate'
  m <- l[[label]]
  pred_index <- which.max(m)
  pred <- target_split[[label]][pred_index] 
  res_list[[label]] <- pred
  
  label <- 'oxygen_microaerophile'
  m <- l[[label]]
  pred <- ifelse(m > 0.5, TRUE, FALSE)
  res_list[['microaerophile']] <- pred
  
  label <- 'spore'
  m <- l[[label]]
  pred <- ifelse(m > 0.5, TRUE, FALSE)
  res_list[['ability_spore']] <- pred
  
  res <- do.call(c, res_list)
  names(res)[1:4] <- target_split[['cellsize']]
  return(res)
}



# Aggregate predictions
l_sub <- lapply(y_pred, column_means)
names(l_sub) <- target_names
predictions <- get_pred(l = l_sub, target_split)

# Convert predictions to data frame and save
df <- as.data.frame(t(predictions))

print(df)