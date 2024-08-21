
# Configure the logger
flog.threshold(INFO)

# Define colors and formatting
reset <- "\033[0m"
bold <- "\033[1m"
red <- "\033[31m"
green <- "\033[32m"
yellow <- "\033[33m"
blue <- "\033[34m"


# Custom logging function
custom_log <- function(level, msg, step, info = "") {
  # Trim leading and trailing whitespaces from the log message
  msg <- trimws(msg)
  info <- trimws(info)
  
  # Format the log message
  formatted_msg <- paste0(
    "[", bold, level, reset, "] ",
    bold, "Step ", step, reset, " - ",
    bold, msg, reset,
    " ", info, "\n"
  )
  
  # Print the formatted log message
  switch(level,
    INFO = cat(green, formatted_msg, reset),
    WARN = cat(yellow, formatted_msg, reset),
    ERROR = cat(red, formatted_msg, reset),
    cat(formatted_msg)
  )
}
