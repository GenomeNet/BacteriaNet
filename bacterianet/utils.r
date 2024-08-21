# Functions

predict_and_write <- function(input_path, step, batch_size, output_path, model, genus_labels, tmp_file) {
  pred <- predict_model(
    output_format = "one_seq",
    model = model,
    layer_name = "dense_3",
    path_input = input_path,
    round_digits = 4,
    step = step,
    batch_size = batch_size,
    verbose = FALSE,
    return_states = TRUE,
    padding = "standard",
    mode = "label",
    format = "fasta",
    filename = tmp_file
  )

  df <- data.frame(pred$states)
  names(df) <- genus_labels
  write.table(df, file = output_path, sep = "\t",
  row.names = FALSE, quote = FALSE)

  # Validate the prediction DataFrame
  if (is.null(df) || nrow(df) == 0) {
    stop("Prediction failed. The resulting data frame is empty.")
  }
  return(df)
}

print_top_predictions <- function(df) {
  # aggregate
  agg <- colMeans(df)
  agg_o <- agg[order(agg, decreasing = TRUE)]
  message("Top 5 predictions of the sample:")
  for (i in 1:5) {
    message(paste0("Predicted FASTA sample as ",
    names(agg_o[i]), " (" , round(agg_o[i] * 100, digits = 1), "%)"))
  }
}

layer_pos_embedding <- keras::new_layer_class(
  "layer_pos_embedding",
  
  initialize = function(maxlen=100, vocabulary_size=4, embed_dim=64, ...) {
    super$initialize(...)
    if (embed_dim != 0) {
      self$token_emb <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(vocabulary_size),
                                                              output_dim = as.integer(embed_dim))
      self$position_embeddings <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(maxlen),
                                                                        output_dim = as.integer(embed_dim))
    } else {
      self$position_embeddings <- tensorflow::tf$keras$layers$Embedding(input_dim = as.integer(maxlen),
                                                                        output_dim = as.integer(vocabulary_size))
    }
    self$embed_dim <- as.integer(embed_dim)
    self$maxlen <- as.integer(maxlen)
    self$vocabulary_size <- as.integer(vocabulary_size)
  },
  
  call = function(inputs) {
    positions <- tensorflow::tf$range(self$maxlen, dtype = "int32") 
    embedded_positions <- self$position_embeddings(positions)
    if (self$embed_dim != 0) inputs <- self$token_emb(inputs)
    inputs + embedded_positions
  },
  
  get_config = function() {
    config <- super$get_config()
    config$maxlen <- self$maxlen
    config$vocabulary_size <- self$vocabulary_size
    config$embed_dim <- self$embed_dim
    config
  }
)

layer_pos_sinusoid <- keras::new_layer_class(
  "layer_pos_sinusoid",
  initialize = function(maxlen, vocabulary_size, n, embed_dim, ...) {
    super$initialize(...)
    self$maxlen <- as.integer(maxlen)
    self$vocabulary_size <- vocabulary_size
    self$n <- as.integer(n)
    self$pe_matrix <- positional_encoding(seq_len = maxlen,
                                          d_model = ifelse(embed_dim == 0,
                                                           as.integer(vocabulary_size),
                                                           as.integer(embed_dim)),  
                                          n = n)
    
    if (embed_dim != 0) {
      self$token_emb <- tensorflow::tf$keras$layers$Embedding(input_dim = vocabulary_size, output_dim = as.integer(embed_dim))
    }
    self$embed_dim <- as.integer(embed_dim)
    
    # self$position_embedding <- tensorflow::tf$keras$layers$Embedding(
    #   input_dim = as.integer(maxlen),
    #   output_dim = ifelse(is.null(embed_dim),
    #                       as.integer(vocabulary_size),
    #                       as.integer(embed_dim)), 
    #   weights = list(pe_matrix),
    #   trainable = FALSE,
    #   name="position_embedding"
    # )(tensorflow::tf$range(start=0, limit=as.integer(maxlen), delta=1))
    
  },
  
  call = function(inputs) {
    if (self$embed_dim != 0) {
      inputs <- self$token_emb(inputs)
    } 
    inputs + self$pe_matrix
  },
  
  get_config = function() {
    config <- super$get_config()
    config$maxlen <- self$maxlen
    config$vocabulary_size <- self$vocabulary_size
    config$n <- self$n
    config$embed_dim <- self$embed_dim
    config$pe_matrix <- self$pe_matrix
    config
  }
)

layer_transformer_block <- keras::new_layer_class(
  "layer_transformer_block",
  initialize = function(num_heads=2, head_size=4, dropout_rate=0, ff_dim=64L, vocabulary_size=4, embed_dim=64, ...) {
    super$initialize(...)
    self$num_heads <- num_heads
    self$head_size <- head_size
    self$dropout_rate <- dropout_rate
    self$ff_dim <- ff_dim
    self$embed_dim <- as.integer(embed_dim)
    self$vocabulary_size <- vocabulary_size
    self$att <- tensorflow::tf$keras$layers$MultiHeadAttention(num_heads=as.integer(num_heads),
                                                               key_dim=as.integer(head_size))
    
    self$ffn <- keras::keras_model_sequential() %>% keras::layer_dense(units=as.integer(ff_dim), activation="relu") %>%
      keras::layer_dense(units=ifelse(embed_dim == 0, as.integer(vocabulary_size), as.integer(embed_dim)))
    
    self$layernorm1 <- keras::layer_layer_normalization(epsilon=1e-6)
    self$layernorm2 <- keras::layer_layer_normalization(epsilon=1e-6)
    self$dropout1 <- keras::layer_dropout(rate=dropout_rate)
    self$dropout2 <- keras::layer_dropout(rate=dropout_rate)
  },
  
  call = function(inputs) {
    attn_output <- self$att(inputs, inputs, inputs)
    attn_output <- self$dropout1(attn_output)
    out1 <- self$layernorm1(inputs + attn_output)
    ffn_output <- self$ffn(out1)
    ffn_output <- self$dropout2(ffn_output)
    seq_output <- self$layernorm2(out1 + ffn_output)
    return(seq_output)
  },
  
  get_config = function() {
    config <- super$get_config()
    config$num_heads <- self$num_heads
    config$head_size <- self$head_size
    config$dropout_rate <- self$dropout_rate
    config$ff_dim <- self$ff_dim
    config$vocabulary_size <- self$vocabulary_size
    config$embed_dim <- self$embed_dim
    config
  }
)

layer_aggregate_time_dist <- keras::new_layer_class(
  "layer_aggregate_time_dist",
  
  initialize = function(method, ...) {
    super$initialize(...)
    self$method <- method
  },
  
  call = function(inputs, mask = NULL) {
    out <- list()
    if ("sum" %in% self$method) {
      out <- c(out, tensorflow::tf$math$reduce_sum(inputs, axis = 1L))
    }
    if ("mean" %in% self$method) {
      out <- c(out, tensorflow::tf$math$reduce_mean(inputs, axis = 1L))
    }
    if ("max" %in% self$method) {
      out <- c(out, tensorflow::tf$math$reduce_max(inputs, axis = 1L))
    }
    
    if (length(out) > 1) {
      out <- tensorflow::tf$concat(out, axis = -1L)
    } else {
      out <- out[[1]]
    }
    
    out
  },
  
  get_config = function() {
    config <- super$get_config()
    config$method <- self$method
    config
  }
)

custom_objects <- list(
  "layer_pos_embedding" = layer_pos_embedding,
  "layer_pos_sinusoid" = layer_pos_sinusoid,
  "layer_transformer_block" = layer_transformer_block,
  "layer_aggregate_time_dist" = layer_aggregate_time_dist
)