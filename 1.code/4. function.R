library(rentrez)
library(XML)
library(GEOquery)
library(httr)
library(jsonlite)
library(stringi)
library(dplyr)
library(purrr)
library(stringr)

# 1. Functions Definition
# 1.1. Prompt Generation for Classification Logic
generate_classification_logic_prompt <- function(ID, gse_data, relevant_columns, original_prompt) {
  prompt <- paste0(
    "I am conducting a bioinformatics analysis and need assistance with automated sample grouping for differential expression analysis. Here are the details:\n",
    "gse_id: ",ID,
    "\nabstract: ",gse_data@experimentData@abstract,
    "\nExperimental design: ",gse_data@experimentData@other[["overall_design"]]
  )
  
  samples_info <- apply(relevant_columns, 1, function(x) paste(x, collapse = "\t"))
  prompt <- paste(prompt, paste("\nSample info: ",samples_info, collapse = "\n"), "\n")
  prompt <- paste(prompt, "All requirements you need to follow:Generate a tab-separated table based on the provided sample information. The table must strictly contain the following columns only/: gse_id, cell_type, ctrl_ids, pert_ids, type(According to the description, categorize into gene/disease/drug/other. If it's definitely not possible to categorize, only then classify as 'other'.), and pert_name. Each row must represent a unique perturbation. Use '\t' for column separation and '|' to separate multiple GSM ids within the same cell. Please ensure the response contains only the table, with no introductory, explanatory, or additional text.",
                  original_prompt)
  
  return(prompt)
}
# 1.2. ChatGPT API Interaction
get_chatGPT_answer <- function(content, api_key, max_retries = 100) {
  attempt <- 1
  errors <- list()
  
  repeat {
    tryCatch({
      response <- httr::POST(
        url = "https://newapi.ikungpt.com/v1/chat/completions",
        httr::add_headers(Authorization = paste("Bearer", api_key)),
        httr::content_type_json(),
        encode = "json",
        body = list(
          model = "gpt-4",
          messages = list(list(role = "user", content = content))
        ),
        timeout(60) # 设置超时时间为60秒
      )
      
      response_text <- rawToChar(response$content)
      response_json <- jsonlite::fromJSON(response_text)
      chatGPT_answer <- response_json[["choices"]][["message"]][["content"]]
      return(chatGPT_answer)
    }, error = function(e) {
      errors[[length(errors) + 1]] <- e$message
      if (attempt <= max_retries) {
        Sys.sleep(5) # 等待5秒
        attempt <- attempt + 1
      } else {
        return(list("error" = "Reached maximum number of retries", "errors" = errors))
      }
    })
  }
}


# 1.3. Quality Assurance Prompt Generation
generate_quality_assurance_prompt <- function(original_prompt, original_output, target_output) {
  qc_prompt <- paste0(readLines("qc_prompt.txt"),collapse = "\t")
  prompt <- paste0(
    "(Required: Directly provide me with the prompt content without any description or introduction)As a quality assurance engineer, I am reviewing the output of a ChatGPT-based classification system for bioinformatics data. Below are the details:\n\n",
    "Fixed Original Prompt:\n","All requirements you need to follow:Generate a tab-separated table based on the provided sample information. The table should strictly contain the following columns only/: gse_id, cell_type, ctrl_ids, pert_ids, type, and pert_name. Each row must represent a unique perturbation. Use '\t' for column separation and '|' to separate multiple GSM ids within the same cell. Please ensure the response contains only the table, with no introductory, explanatory, or additional text.",
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
process_classification_results <- function(group_results) {
  group_results <- gsub("```", "", group_results)
  df <- try(read.table(text = group_results, header = F, sep = "\t", stringsAsFactors = FALSE, fill = FALSE, quote = ""), silent = TRUE)
  if (!inherits(df, "try-error") && ncol(df) == 6) {
    return(df)
  } else {
    cat("Conversion to data frame failed.\n")
    return(NULL)
  }
}

# 1.5 all_run
# Function to process GSE data
process_gse_data <- function(gse_list, train_data, original_prompt) {
  for (i in seq_along(unique(train_data$geo_id))) {
    cat("Processing GSE ID: ", i, "\n")
    ID <- unique(train_data$geo_id)[i]
    gse_data <- gse_list[[ID]]
    if (is.null(gse_data)) {
      cat("No series data for GSE ID: ", ID, "\n")
      next
    }
    
    pheno_data <- pData(phenoData(gse_data))
    relevant_columns <- pheno_data[, grepl("description|title|geo_accession|characteristics_", names(pheno_data))]
    # subset_data <- train_data[train_data$geo_id == ID, ]
    output_str <- paste(colnames(relevant_columns), collapse = "\t")
    data_rows <- apply(relevant_columns, 1, function(row) paste(row, collapse = "\t"))
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
      
      while (is.null(original_output) || nchar(original_output) == 0) {
        Sys.sleep(5)
        original_output <- get_chatGPT_answer(prompt, api_key)
      }
      
      df <- process_classification_results(original_output)
      if (!is.null(df)) {
        final_results[[ID]] <- df
        optimized_prompts[[ID]] <- original_prompt
        satisfied <- TRUE
      } else {
        cat("Attempt ", attempt, " for GSE ", ID, " failed. Retrying...\n")
      }
      attempt <- attempt + 1
    }
    
    if (!satisfied) {
      final_results[[ID]] <- df
      optimized_prompts[[ID]] <- original_prompt
      cat("Max attempts reached for GSE ", ID, ". Result saved.\n")
    }
    
    if (i %% 5 == 0) {
      save_file_name <- paste0("GEO_data_backup_", i, ".RData")
      save(final_results, optimized_prompts,out_of_length, file = save_file_name)
      cat("Data saved in ", save_file_name, "\n")
    }
  }
  
  return(list(final_results = final_results, optimized_prompts = optimized_prompts, out_of_length = out_of_length))
}

is_true <- function(expr) {
  return(!is.null(expr) && expr)
}

perform_quality_checks <- function(final_results_df, gse_list, api_key, save_path="./") {
  check_results <- vector("list", length(unique(final_results_df$gse_id)))
  
  for (row in seq_len(length(unique(final_results_df$gse_id)))) {
    print(row)
    gse_id <- final_results_df$gse_id[row]
    group_info <- final_results_df[final_results_df$gse_id%in%gse_id,]
    
    final_string <- group_info %>%
      mutate(across(everything(), as.character)) %>%  # 确保所有列都是字符型
      {
        col_names <- paste(names(.), collapse = "\t")
        data_rows <- apply(., 1, paste, collapse = "\t")
        paste(col_names, paste(data_rows, collapse = "\n"), sep = "\n")
      }
    
    gse_data <- gse_list[[gse_id]]
    pheno_data <- pData(phenoData(gse_data))
    GSM_id <- unique(unlist(str_match_all(as.character(group_info),pattern = "GSM[0-9]+")))
    pheno_data <- pheno_data[pheno_data$geo_accession%in%GSM_id,]
    
    
    relevant_columns <- pheno_data[, grepl("description|title|geo_accession|characteristics_", names(pheno_data))]
    prompt <- build_quality_check_prompt(gse_id, final_string,  gse_data, relevant_columns)
    
    if(str_length(prompt)>15000|stri_count(prompt, regex = "\\S+") > 2000){
      relevant_columns <- pheno_data[, grepl("title|geo_accession|characteristics_", names(pheno_data))]
      prompt <- build_quality_check_prompt(gse_id, final_string,  gse_data, relevant_columns)
    }
    # 持续尝试获取非空结果
    response <- NULL
    while (is_true(is.null(response) || length(response) == 0) || 
           is_true(is.data.frame(df) && nrow(df) != nrow(group_info)) || 
           is_true(is.data.frame(df) && ncol(df) != 2)){
      df <- data.frame()
      -
      response <- get_chatGPT_answer(content = prompt, api_key = api_key)
      df <- try(read.table(text = response, header = FALSE,quote = "'", sep = "\t", stringsAsFactors = FALSE, fill = FALSE), silent = TRUE)
      
      if (!is.null(response) && length(response) > 0 && is.data.frame(df)) {
        check_results[[row]] <- df
      } else {
        print("Retrying to get response...")
      }
    }
    
    # 每5行保存一次
    if (row %% 1 == 0 || row == length(unique(final_results_df$gse_id))) {
      save_file_path <- file.path(save_path, paste0("check_results_", row, ".RData"))
      save(check_results, file = save_file_path)
      print(paste("Data saved to", save_file_path))
    }
  }
  
  return(check_results)
}


build_quality_check_prompt <- function(gse_id, final_string, gse_data, relevant_columns) {

  # Extract relevant information for the prompt
  abstract <- gse_data@experimentData@abstract
  # design <- gse_data@experimentData@other[["overall_design"]]
  sample_characteristics <- paste(apply(relevant_columns, 1, paste, collapse="\t"), collapse="\n")%>%
    str_replace_all(.,pattern = "X.",replacement = "")%>%
  str_replace_all(.,pattern = "\\\"",replacement = "")

  # Check if abstract is empty and substitute if necessary
  if (is.null(gse_data@experimentData@abstract) || gse_data@experimentData@abstract == "") {
    abstract <- "No description"
  } else {
    abstract <- gse_data@experimentData@abstract
  }
  
  # Check if design is empty and substitute if necessary
  # if (is.null(gse_data@experimentData@other[["overall_design"]]) || gse_data@experimentData@other[["overall_design"]] == "") {
  #   design <- "No description"
  # } else {
  #   design <- gse_data@experimentData@other[["overall_design"]]
  # }
  
  # Construct the prompt
  prompt <- sprintf(
    "I need assistance to check if the grouping is correct.\nHere are my grouping details(Note:You should use the form of a table to understand the information I have given you,'\t' is a tab and '\n' is a line break,The entries in the Control IDs (ctrl_ids) and Treatment IDs (pert_ids) are separated by vertical bars (|).Examine each row independently, without being influenced by multiple rows.): (%s)\n---\n",
    final_string
  )
  verification_conditions <- readLines("2. supple_data/verify_prompt.txt")
  prompt <- paste0(prompt,sprintf("The information you can refer to is as following:\n---\nExperimental Abstract: %s\n---\nCharacteristics of each GSM (One GSM per row):\n%s\n---\n",
                                                           abstract,# design,
                                                           sample_characteristics), verification_conditions)
  return(prompt)
}

identify_reanalysis_candidates <- function(final_results_df) {
  # Identify rows where quality check failed (contains "×")
  reanalysis_GSE <- final_results_df$gse_id[grepl("×", final_results_df$qc_results)]
  
  return(reanalysis_GSE)
}

create_results_dataframe <- function(results) {
  # Convert the list of results into a dataframe
  results_df <- do.call(rbind, results)
  
  # Ensure column names are set correctly
  colnames(results_df) <- c('gse_id', 'cell_type', 'ctrl_ids', 'pert_ids', 'type', 'pert_name')
  
  # Example: Filter out rows that don't match a specific pattern in 'gse_id'
  results_df <- results_df[grep(pattern = "GSE[0-9]+", x = results_df$gse_id), ]
  
  return(results_df)
}

merge_results_fun <- function(original_results_df, reanalysis_results_df) {
  # Identify the GSE IDs in the reanalysis results
  reanalyzed_ids <- unique(reanalysis_results_df$gse_id)
  
  # Remove rows from the original results that have been reanalyzed
  original_results_df <- original_results_df[!original_results_df$gse_id %in% reanalyzed_ids,]
  
  # Merging the original and reanalysis results
  merged_results_df <- rbind(original_results_df, reanalysis_results_df)
  
  return(merged_results_df)
}

reanalyze_data <- function(reanalysis_GSE, gse_list, train_data, original_prompt, api_key, max_reanalysis_attempts = 5) {
  final_results_reanalysis <- list()
  
  for (ID in reanalysis_GSE) {
    attempt_count <- 0
    satisfied <- FALSE
    
    while (!satisfied && attempt_count < max_reanalysis_attempts) {
      # Reprocess the GSE data
      prompt <- generate_classification_logic_prompt(ID, gse_list[[ID]], relevant_columns, original_prompt)
      original_output <- get_chatGPT_answer(prompt, api_key)
      
      # Check for successful processing and handle accordingly
      df <- process_classification_results(original_output)
      if (!is.null(df)) {
        final_results_reanalysis[[ID]] <- df
        satisfied <- TRUE
      } else {
        attempt_count <- attempt_count + 1
      }
    }
    
    # Handle cases where max attempts were reached but no satisfactory result was obtained
    if (!satisfied) {
      df <- matrix(ncol = 6,dimnames = list(NULL,c('gse_id', 'cell_type', 'ctrl_ids', 'pert_ids', 'type', 'pert_name')))
      final_results_reanalysis[[ID]] <- df
    }
  }
  
  return(final_results_reanalysis)
}

final_quality_checks <- function(merged_results_df, gse_list, api_key) {
  final_check_results <- vector("list", nrow(merged_results_df))
  
  for (row in seq_len(nrow(merged_results_df))) {
    gse_id <- merged_results_df$gse_id[row]
    gse_data <- gse_list[[gse_id]]
    pheno_data <- pData(phenoData(gse_data))
    relevant_columns <- pheno_data[, grepl("description|title|geo_accession|characteristics_", names(pheno_data))]
    
    # 这里构建的提示可能与初次检查不同，可能需要反映重新分析后的情况
    prompt <- build_quality_check_prompt(gse_id, merged_results_df, row, gse_data, relevant_columns)
    
    final_check_results[[row]] <- get_chatGPT_answer(content = prompt, api_key = api_key)
  }
  
  return(final_check_results)
}

get_gse_info <- function(gse_id, gse_list) {
  phenoData <- gse_list[[gse_id]]@phenoData@data
  experimentData <- gse_list[[gse_id]]@experimentData@other
  
  organism <- unique(phenoData[["organism_ch1"]])
  GPL <- unique(phenoData[["platform_id"]])
  Exp_type <- unique(experimentData[["type"]])
  
  # 如果存在多个值，使用逗号和空格连接它们
  info <- c(
    organism = paste(organism, collapse = ", "),
    GPL = paste(GPL, collapse = ", "),
    Exp_type = paste(Exp_type, collapse = ", ")
  )
  
  return(info)
}
