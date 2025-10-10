# Test Runner for FCG Analysis Scripts
#
# This script runs analysis scripts and checks for expected outputs

library(tidyverse)

# Set testing mode
Sys.setenv(TESTING_MODE = "TRUE")

# Function to run a single script and log results
test_single_script <- function(script_path,
                                expected_output = NULL,
                                timeout_min = 10,
                                description = "") {

  script_name <- basename(script_path)
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  log_file <- paste0("tests/logs/", gsub("\\.R$", "", script_name), "_", timestamp, ".log")

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("Testing:", script_name, "\n")
  if (description != "") cat("Description:", description, "\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Check if script exists
  if (!file.exists(script_path)) {
    cat("ERROR: Script not found:", script_path, "\n")
    write(paste("ERROR: Script not found:", script_path), log_file)
    return(list(status = "FAIL", reason = "Script not found"))
  }

  # Create log file
  if (!dir.exists("tests/logs")) {
    dir.create("tests/logs", recursive = TRUE)
  }

  # Capture start time
  start_time <- Sys.time()

  # Run the script and capture output
  result <- tryCatch({
    # Redirect output to log file
    sink(log_file, split = TRUE)
    cat("Script:", script_path, "\n")
    cat("Start time:", as.character(start_time), "\n")
    cat(rep("-", 80), "\n\n", sep = "")

    # Source the script
    source(script_path, local = new.env())

    cat("\n", rep("-", 80), "\n", sep = "")
    cat("Completed successfully\n")
    sink()

    list(status = "PASS", reason = "Script executed without errors")

  }, error = function(e) {
    cat("\n", rep("-", 80), "\n", sep = "")
    cat("ERROR:", e$message, "\n")
    sink()

    list(status = "FAIL", reason = paste("Error:", e$message))
  }, warning = function(w) {
    cat("\nWARNING:", w$message, "\n")
    list(status = "WARN", reason = paste("Warning:", w$message))
  })

  # Calculate elapsed time
  end_time <- Sys.time()
  elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

  # Check for expected outputs
  output_check <- "N/A"
  if (!is.null(expected_output)) {
    if (grepl("\\*", expected_output)) {
      # Pattern matching for outputs
      pattern <- gsub("\\*", ".*", expected_output)
      matching_files <- list.files(dirname(expected_output),
                                    pattern = basename(pattern),
                                    full.names = TRUE)
      if (length(matching_files) > 0) {
        output_check <- "FOUND"
        cat("Expected outputs found:", length(matching_files), "files\n")
      } else {
        output_check <- "MISSING"
        cat("WARNING: Expected outputs not found\n")
      }
    } else {
      # Exact file check
      if (file.exists(expected_output)) {
        output_check <- "FOUND"
        cat("Expected output found:", expected_output, "\n")
      } else {
        output_check <- "MISSING"
        cat("WARNING: Expected output not found:", expected_output, "\n")
      }
    }
  }

  # Print summary
  cat("\n")
  cat("Status:", result$status, "\n")
  cat("Elapsed time:", round(elapsed, 2), "minutes\n")
  cat("Output check:", output_check, "\n")
  cat("Log file:", log_file, "\n")

  if (result$status == "FAIL") {
    cat("Reason:", result$reason, "\n")
  }

  return(list(
    script = script_name,
    status = result$status,
    elapsed_min = elapsed,
    output_check = output_check,
    log_file = log_file,
    reason = result$reason
  ))
}

# Function to run all tests from config
run_all_tests <- function(config_file = "tests/test_config.csv") {

  if (!file.exists(config_file)) {
    cat("No test configuration found at:", config_file, "\n")
    cat("Creating example configuration...\n")
    create_example_config()
    return(invisible(NULL))
  }

  # Read config
  config <- read_csv(config_file, show_col_types = FALSE)

  # Run each test
  results <- list()
  for (i in 1:nrow(config)) {
    row <- config[i, ]
    result <- test_single_script(
      script_path = row$script_path,
      expected_output = row$expected_output,
      timeout_min = row$timeout_min,
      description = ifelse(is.na(row$description), "", row$description)
    )
    results[[i]] <- result
  }

  # Create summary
  results_df <- bind_rows(results)

  cat("\n", rep("=", 80), "\n", sep = "")
  cat("TEST SUMMARY\n")
  cat(rep("=", 80), "\n\n", sep = "")

  print(results_df %>% select(script, status, elapsed_min, output_check))

  cat("\nTotal tests:", nrow(results_df), "\n")
  cat("Passed:", sum(results_df$status == "PASS"), "\n")
  cat("Failed:", sum(results_df$status == "FAIL"), "\n")
  cat("Warnings:", sum(results_df$status == "WARN"), "\n")

  # Save results
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_file <- paste0("tests/logs/test_results_", timestamp, ".csv")
  write_csv(results_df, results_file)
  cat("\nResults saved to:", results_file, "\n")

  return(results_df)
}

# Function to create example test configuration
create_example_config <- function() {
  config <- tibble(
    script_path = c(
      "scripts/lipid_structure_visualisation.R",
      "scripts/lipidomic_bubbles.R"
    ),
    test_type = c("visualization", "visualization"),
    expected_output = c(
      "tests/outputs/lipid_structure_vis_*.csv",
      "tests/outputs/lipidomic_bubbles_*.csv"
    ),
    timeout_min = c(5, 5),
    description = c(
      "Lipid structural analysis visualization",
      "Lipidomic bubble plots"
    )
  )

  write_csv(config, "tests/test_config.csv")
  cat("Example configuration created at: tests/test_config.csv\n")
  cat("Edit this file to add more tests.\n")

  return(config)
}

# Main execution
if (!interactive()) {
  cat("FCG Analysis Test Runner\n")
  cat("========================\n\n")

  # Check if config exists, if not create example
  if (!file.exists("tests/test_config.csv")) {
    create_example_config()
  } else {
    # Run all tests
    run_all_tests()
  }
}
