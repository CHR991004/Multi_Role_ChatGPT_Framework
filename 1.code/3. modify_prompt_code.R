# Required libraries
library(rentrez)
library(XML)
library(GEOquery)
library(httr)
library(jsonlite)
library(stringi)
library(dplyr)
library(purrr)

# -------------------------
# 1. Functions Definition
# -------------------------

# 1.1. Prompt Generation for Classification Logic
generate_classification_logic_prompt <- function(ID, gse_data, relevant_columns, original_prompt) {
  prompt <- paste0(
    "I am conducting a bioinformatics analysis and need assistance with automated sample grouping for differential expression analysis. Here are the details:\n",
    "gse_id: ", ID,
    "\nabstract: ", gse_data@experimentData@abstract,
    "\nExperimental design: ", gse_data@experimentData@other[["overall_design"]]
  )
  
  samples_info <- apply(relevant_columns, 1, function(x) paste(x, collapse = "\t"))
  prompt <- paste(prompt, paste("\nSample info: ", samples_info, collapse = "\n"), "\n")
  prompt <- paste(prompt, "All requirements you need to follow:Generate a tab-separated table based on the provided sample information. The table should strictly contain the following columns only/: gse_id, cell_type, ctrl_ids, pert_ids, type, and pert_name. Each row must represent a unique perturbation. Use '\t' for column separation and '|' to separate multiple GSM ids within the same cell. Please ensure the response contains only the table, with no introductory, explanatory, or additional text.",
                  original_prompt)
  
  return(prompt)
}

# 1.2. ChatGPT API Interaction
get_chatGPT_answer <- function(content, api_key, max_retries = 100) {
  attempt <- 1
  
  repeat {
    tryCatch({
      response <- httr::POST(
        url = "https://api.openai.com/v1/chat/completions",
        httr::add_headers(Authorization = paste("Bearer", api_key)),
        httr::content_type_json(),
        encode = "json",
        body = list(
          model = "gpt-4-1106-preview",
          messages = list(list(role = "user", content = content))
        )
      )
      
      response_text <- rawToChar(response$content)
      response_json <- jsonlite::fromJSON(response_text)
      chatGPT_answer <- response_json[["choices"]][["message"]][["content"]]
      return(chatGPT_answer)
    }, error = function(e) {
      if (grepl("schannel: failed to receive handshake, SSL/TLS connection failed", e$message)) {
        if (attempt <= max_retries) {
          Sys.sleep(5)  # Wait for 5 seconds before retrying
          attempt <- attempt + 1
        } else {
          stop("Reached maximum number of retries. Last error: ", e$message)
        }
      } else {
        stop("An error occurred: ", e$message)
      }
    })
  }
}

# 1.3. Quality Assurance Prompt Generation
generate_quality_assurance_prompt <- function(original_prompt, original_output, target_output) {
  qc_prompt <- paste0(readLines("qc_prompt.txt"), collapse = "\t")
  prompt <- paste0(
    "(Required: Directly provide me with the prompt content without any description or introduction)As a quality assurance engineer, I am reviewing the output of a ChatGPT-based classification system for bioinformatics data. Below are the details:\n\n",
    "Fixed Original Prompt:\n", "All requirements you need to follow:Generate a tab-separated table based on the provided sample information. The table should strictly contain the following columns only/: gse_id, cell_type, ctrl_ids, pert_ids, type, and pert_name. Each row must represent a unique perturbation. Use '\t' for column separation and '|' to separate multiple GSM ids within the same cell. Please ensure the response contains only the table, with no introductory, explanatory, or additional text.",
    "Modifiable Original Prompt:\n", original_prompt, "\n\n",
    "gse_id: ", ID,
    "\nabstract: ", gse_data@experimentData@abstract,
    "\nExperimental design: ", gse_data@experimentData@other[["overall_design"]],
    "Original Output (before any prompt modification):\n", original_output, "\n\n",
    "Target Output:\n", target_output, "\n\n",
    qc_prompt
  )
  return(prompt)
}

# 1.4. Process Classification Results into a Dataframe
process_classification_results <- function(group_results, expected_column_count) {
  df <- try(read.table(text = group_results, header = FALSE, sep = "\t", stringsAsFactors = FALSE, fill = FALSE, quote = ""), silent = TRUE)
  if (!inherits(df, "try-error") && ncol(df) == expected_column_count) {
    return(df)
  } else {
    cat("Conversion to data frame failed.\n")
    return(NULL)
  }
}

# -------------------------
# 2. Workflow
# -------------------------

# 2.1 Data Preparation

# Data Loading (uncomment if needed)
# load("gse_list.Rdata")
# load("train_list.Rdata")

max_attempts <- 3
expected_column_count <- 6  # Define expected column count based on your data structure
api_key <- "sk-yj0SfCYBRHGETByP9e1a3b136bA9*********************"
train_data <- openxlsx::read.xlsx("GEO_train.xlsx")

# Calculate proportions of each sample type
prop <- table(train_data$type) / nrow(train_data)
# Calculate sample size per type
sample_size <- round(prop * 40)
# Set random seed for reproducibility
set.seed(123)
# Stratified sampling using map function
selected_sample <- train_data %>%
  split(.$type) %>%
  map2(sample_size, ~slice_sample(.x, n = .y)) %>%
  bind_rows() %>%
  sample_n(nrow(.))  
train_data <- train_data[train_data$geo_id %in% selected_sample$geo_id,]

final_results <- list()
optimized_prompts <- list()

# Initialize dataframe for prompt outputs
prompts_outputs_df <- data.frame(ID = character(), original_prompt = character(), qa_prompt = character(), qa_output = character(), stringsAsFactors = FALSE)

# Load original_prompt from file
original_prompt <- paste0(readLines("2. supple_data/requirement.txt"), collapse = "\t")
out_of_length <- c()

# Processing loop for each unique GEO ID
for (i in 1:length(unique(train_data$geo_id))) {
  print(i)
  ID <- unique(train_data$geo_id)[i]
  gse_data <- gse_list[[ID]]
  if (is.null(gse_data)) {
    print("don't have series data", gse_data)
    next
  }
  pheno_data <- pData(phenoData(gse_data))
  relevant_columns <- pheno_data[, grepl("title|geo_accession|characteristics_", names(pheno_data))]
  
  # Subset training data
  subset_data <- train_data[train_data$geo_id == ID, ]
  
  # Building target output string
  output_str <- paste(colnames(subset_data), collapse = "\t")
  data_rows <- apply(subset_data, 1, function(row) paste(row, collapse = "\t"))
  target_output <- paste(output_str, paste(data_rows, collapse = "\n"), sep = "\n")
  
  satisfied <- FALSE
  attempt <- 1
  
  while (!satisfied && attempt <= max_attempts) {
    prompt <- generate_classification_logic_prompt(ID, gse_data, relevant_columns, original_prompt)
    if (stri_count(prompt, regex = "\\S+") > 2000) {
      out_of_length <- c(out_of_length, ID)
      break
    }
    original_output <- get_chatGPT_answer(prompt, api_key)
    
    while (length(original_output) == 0) {
      Sys.sleep(5)  # Wait before retrying
      original_output <- get_chatGPT_answer(prompt, api_key)
    }
    qa_prompt <- generate_quality_assurance_prompt(original_prompt, original_output, target_output)
    qa_output <- get_chatGPT_answer(qa_prompt, api_key)
    while (length(qa_output) == 0) {
      Sys.sleep(5)  # Wait before retrying
      qa_output <- get_chatGPT_answer(qa_prompt, api_key)
    }
    if (!grepl("NO_CHANGE_REQUIRED", qa_output)) {
      prompts_outputs_df <- rbind(prompts_outputs_df, data.frame(ID, original_prompt, qa_prompt, qa_output))
      original_prompt <- qa_output # Update original prompt based on feedback
    } else {
      df <- process_classification_results(original_output, expected_column_count)
      if (!is.null(df)) {
        final_results[[ID]] <- df
        optimized_prompts[[ID]] <- original_prompt # Save optimized prompt
        satisfied <- TRUE
      } else {
        cat("Attempt", attempt, "to process GSE", ID, "failed. Retrying...\n")
      }
    }
    
    attempt <- attempt + 1
  }
  
  # Handling unsatisfied conditions after maximum attempts
  if (!satisfied) {
    final_results[[ID]] <- df # df might be NULL if all attempts were unsuccessful
    optimized_prompts[[ID]] <- original_prompt
    openxlsx::write.xlsx(prompts_outputs_df, file ="final_prompts_outputs.xlsx", overwrite = TRUE)
    cat("Max attempts reached for GSE", ID, ". Current result saved.\n")
  }
  
  # Save data backup
  save_file_name <- paste0("GEO_data_backup_", format(Sys.time(), "%H_%M_%S"), ".RData")
  save(final_results, optimized_prompts, file = save_file_name)
  cat("Data saved in", save_file_name, "\n")
}

# Combine final results into a dataframe
final_results_df <- do.call(rbind, final_results)
colnames(final_results_df) <- c('gse_id', 'cell_type', 'ctrl_ids', 'pert_ids', 'type', 'pert_name')

# Write final results to Excel file
openxlsx::write.xlsx(final_results_df, file = "final_results_df.xlsx")
