library(rentrez)
library(XML)
library(GEOquery)
library(httr)
library(jsonlite)
library(stringi)
library(dplyr)
library(purrr)

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
        url = "https://api.ikungpt.com/v1/chat/completions",
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
          Sys.sleep(5) # 等待5秒
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
  df <- try(read.table(text = group_results, header = F, sep = "\t", stringsAsFactors = FALSE, fill = FALSE, quote = ""), silent = TRUE)
  if (!inherits(df, "try-error") && ncol(df) == expected_column_count) {
    return(df)
  } else {
    cat("Conversion to data frame failed.\n")
    return(NULL)
  }
}
# 2. workflow
## 2.1 data prepare

# load("gse_list.Rdata")
# load("train_list.Rdata")
max_attempts <- 3
expected_column_count <- 6 # 根据您的数据结构定义期望的列数
api_key <- "sk-yj0SfCYBRHGETByP9e1a3b136bA94684AdF24f7cBe519578"
train_data <- openxlsx::read.xlsx("GEO_train.xlsx")

# 计算每个类型的样本比例
prop <- table(train_data$type) / nrow(train_data)
# 按比例计算每个类型应抽取的数量
sample_size <- round(prop * 40)
# 设置随机种子
set.seed(123)
# 使用map函数进行分层抽样
selected_sample <- train_data %>%
  split(.$type) %>%
  map2(sample_size, ~slice_sample(.x, n = .y)) %>%
  bind_rows() %>%
  sample_n(nrow(.))  
train_data <- train_data[train_data$geo_id%in%selected_sample$geo_id,]

final_results <- list()
optimized_prompts <- list()

# 初始化变量
prompts_outputs_df <- data.frame(ID = character(), original_prompt = character(), qa_prompt = character(), qa_output = character(), stringsAsFactors = FALSE)

# 设置初始的 original_prompt
# original_prompt <- "Return the sample grouping results in a tab-separated table format.\tGSM samples need to be divided into two groups: control and test, with each row representing one type of perturbation.\tIf all samples do not meet the criteria for comparison or no biological duplications, only GSE_ID will be returned\tRecord the corresponding tissue, cell line, and time for each perturbation in the rows. Use 'NA' to display if there is no information.\tApart from the table, no additional text should be returned.\tDo not use '|' as a delimiter and do not use Markdown syntax. Use '\t' as the separator.\tThe output must be a table ('\t' separate) with columns: geo_id(GSE_id), cell_type, ctrl_id(control GSM id), pert_id(perturbation GSM id), type(only can in disease/drug/gene), pert_name(Detailed perturbation name).\tIf the GSM with the same processing conditions is put into the same cells, separated by '|'\tExtremely important: No comparison is required between the control and control groups.\ttest_name only can be one perturbation factor"
# original_prompt <- "\n\n1. While generating the table, if the treatment involves a combination of factors (like LPS and Simvastatin), do not group them under a single row with combined 'type' and 'pert_name'. Instead, separate them into individual rows for each condition (treatment or control).\n\n2. Clarify that the 'cell_type' should include more specific terms where possible, such as 'Lung (from C57BL/6J)' instead of merely 'mouse lung' or 'Lung Tissue'. \n\n3. The 'type' for LPS should be 'disease' and not 'other'. \n\n4. For 'pert_name', use a more specific term like 'Acute Lung Injury' where applicable instead of simply 'LPS'. \n\nAfter these modifications, the revised prompt should be as follows:\n\nReturn the sample grouping results in a tab-separated table format. GSM samples should be divided into groups representing different perturbations. These groups should include 'control' and 'test'. If the samples do not meet the criteria for comparison or no biological duplicates exist, return only the GSE_ID. Record the specific tissue or cell line from the original sample, along with the time for each perturbation. If the cell line is from a tissue source, include this (for example, 'Lung (from C57BL/6J)' not 'mouse lung'). If no information is available, use 'NA'. Exclude any additional text. Use '|' to separate multiple GSM ids in the same cell, and use '\t' as the other separator. Extremly important: The output should be a table ('\t' separated) with columns: GSE_id, cell_type, ctrl_ids, pert_ids, type (only 'disease', 'drug', 'gene', or 'other'), and a concise, individual perturbation name (pert_name). Group samples with similar processing conditions in the same cells. Exclude comparisons between control groups. The 'pert_name' should only include a single perturbation factor and should not be a combined perturbation. You don't need any extraneous statements but a table at the beginning of your answer to me."
original_prompt <- paste0(readLines("requirement_revision.txt"),collapse = "\t")
out_of_length <- c()
for (i in 1:length(unique(train_data$geo_id))) {
  print(i)
  ID <-unique(train_data$geo_id)[i]
  gse_data <- gse_list[[ID]]
  if (is.null(gse_data)) {
    print("don't have series data",gse_data)
    next
  }
  pheno_data <- pData(phenoData(gse_data))
  relevant_columns <- pheno_data[, grepl("title|geo_accession|characteristics_", names(pheno_data))]
  
  # 获取训练数据子集
  subset_data <- train_data[train_data$geo_id == ID, ]
  
  # 目标输出的字符串构建
  output_str <- paste(colnames(subset_data), collapse = "\t")
  data_rows <- apply(subset_data, 1, function(row) paste(row, collapse = "\t"))
  target_output <- paste(output_str, paste(data_rows, collapse = "\n"), sep = "\n")
  
  satisfied <- FALSE
  attempt <- 1
  
  while (!satisfied && attempt <= max_attempts) {
    prompt <- generate_classification_logic_prompt(ID, gse_data, relevant_columns, original_prompt)
    if (stri_count(prompt,regex = "\\S+")>2000) {
      out_of_length <- c(out_of_length,ID)
      break
    }
    original_output <- get_chatGPT_answer(prompt, api_key)
    # original_output <- substr(original_output, regexpr(":", qa_output) + 1, nchar(qa_output))
    
    while (length(original_output) == 0) {
      Sys.sleep(5)  # 等待一段时间后重试
      original_output <- get_chatGPT_answer(prompt, api_key)
    }
    qa_prompt <- generate_quality_assurance_prompt(original_prompt, original_output, target_output)
    qa_output <- get_chatGPT_answer(qa_prompt, api_key)
    while (length(qa_output) == 0) {
      Sys.sleep(5)  # 等待一段时间后重试
      qa_output <- get_chatGPT_answer(qa_prompt, api_key)
    }
    if (!grepl("NO_CHANGE_REQUIRED", qa_output)) {
      # qa_output <- substr(qa_output, regexpr(":", qa_output) + 1, nchar(qa_output))
      prompts_outputs_df <- rbind(prompts_outputs_df, data.frame(ID, original_prompt, qa_prompt, qa_output))
      original_prompt <- qa_output # 根据反馈修改原始prompt
    } else {
      df <- process_classification_results(original_output)
      if (!is.null(df)) {
        final_results[[ID]] <- df
        optimized_prompts[[ID]] <- original_prompt # 保存优化后的prompt
        satisfied <- TRUE
      } else {
        cat("Attempt", attempt, "to process GSE", ID, "failed. Retrying...\n")
      }
    }
    
    attempt <- attempt + 1
  }
  # 如果在最大尝试次数后仍未满足条件，则保存当前结果
  if (!satisfied) {
    final_results[[ID]] <- df # df可能是NULL，如果在所有尝试中都未成功
    optimized_prompts[[ID]] <- original_prompt
    openxlsx::write.xlsx(prompts_outputs_df, file ="final_prompts_outputs.xlsx",overwrite = )
    cat("Max attempts reached for GSE", ID, ". Current result saved.\n")
  }

  save_file_name <- paste0("GEO_data_backup_", format(Sys.time(), "%H_%M_%S"), ".RData")
  save(final_results, optimized_prompts, file = save_file_name)
  cat("Data saved in", save_file_name, "\n")
}
final_results_df <- do.call(rbind,final_results)
colnames(final_results_df)<-c('gse_id',
                              'cell_type',
                              'ctrl_ids',
                              'pert_ids',
                              'type',
                              'pert_name')
openxlsx::write.xlsx(final_results_df,file = "finalresults_df.xlsx")
