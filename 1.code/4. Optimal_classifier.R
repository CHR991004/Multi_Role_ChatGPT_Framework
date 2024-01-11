# Load required libraries and source files
source("1. code/4. function.R")
#-------
original_prompt <- "Additional Requirements:1.If the patient is on medication, the type is the drug
2.single ctrl_ids and pert_ids must be IDs starting with GSM
3.If the sample description is vague, you can refer to abstract
If the type is a drug, the pert_name should specify the drug name, dosage, and duration of action (if applicable,Format is drug(dose, time))."

load("2. supple_data/gse_list.Rdata")
load("2. supple_data/1128-1367.Rdata")

colnames(results_df_ingselist)[1] <- "geo_id"

# 计算每个 list 元素的列数
ncols <- lapply(gse_list, ncol)

# 找出列数小于或等于 50 的元素的位置
indices <- which(unlist(ncols) <= 50)

# 选择对应的 names
selected_names <- names(gse_list)[indices]

results_df_ingselist <- results_df_ingselist[results_df_ingselist$geo_id%in%selected_names,]

# Randomly select 5 for subsequent grouping.
grouped_df <- results_df_ingselist %>% group_by(drug)
sampled_df <- grouped_df %>% 
  sample_n(size = 5, replace = T) %>% 
  ungroup()%>% distinct()

gse_list <- gse_list[sampled_df$geo_id]
gc()
train_data <- sampled_df[1211:1841,]


# Initialize variables
final_results <- list()
optimized_prompts <- list()
out_of_length <- c()
max_attempts <- 6
api_key <- "sk-KkzHmlYm97TjgLXw63DdDcEc01C14562A383370eF07397Ce" # Replace with actual API key
# Main execution
result_data <- process_gse_data(gse_list, train_data, original_prompt)
final_results_df <- create_results_dataframe(result_data$final_results)
openxlsx::write.xlsx(x = final_results_df,file = "last_achieve.xlsx")


load("2. supple_data/paper_use_list.Rdata")
drug_results <- openxlsx::read.xlsx(xlsxFile = "3. Intermediate files/qc_failed_drug_2.xlsx",sheet = 1)
drug_results$ctrl_count <- sapply(strsplit(drug_results$ctrl_ids, "\\|"), length)
drug_results$pert_count <- sapply(strsplit(drug_results$pert_ids, "\\|"), length)
drug_results <- drug_results[drug_results$ctrl_count >= 3 & drug_results$pert_count >= 3, 1:6]
gse_list <- gse_list[unique(drug_results$gse_id)]


# save(gse_list,drug_results,file = "paper_use_list.Rdata")

#"GSE229146",
# Quality Checks
drug_results2 <- drug_results[drug_results$gse_id%in%unique(drug_results$gse_id)[1:length(unique(drug_results$gse_id))],]
check_results <- perform_quality_checks(final_results_df = drug_results2, gse_list = gse_list, api_key = api_key,save_path = "H:/chat_project/drugGPT/4. achieve/qc2/")

qc_results <- do.call(what = rbind,args = check_results)
# 使用stringr提取√或×或?字符

qc_results$V1 <- unlist(formatted_qc)
openxlsx::write.xlsx(x = final_results_df,file = "final_results_qced.xlsx")

final_results_df$organism <- NA
final_results_df$GPL <- NA
final_results_df$Exp_type <- NA
# 遍历 drug_results 中的每个 GSE ID
for(i in 1:nrow(final_results_df)) {
  gse_id <- final_results_df$gse_id[i]
  
  # 获取每个 GSE ID 的信息
  info <- get_gse_info(gse_id, gse_list)
  
  # 将信息添加到 final_results
  final_results_df$organism[i] <- info["organism"]
  final_results_df$GPL[i] <- info["GPL"]
  final_results_df$Exp_type[i] <- info["Exp_type"]
}

openxlsx::write.xlsx(x = final_results_df,file = "final_results_qced_ann.xlsx")

# Reanalysis
reanalysis_GSE <- identify_reanalysis_candidates(final_results_df)
final_results_reanalysis <- reanalyze_data(reanalysis_GSE, gse_list, train_data, original_prompt, api_key, max_reanalysis_attempts = 5)
final_results_df_reanalysis <- create_results_dataframe(final_results_reanalysis)

# Merging the results
merged_results <- merge_results_fun(final_results_df, final_results_df_reanalysis)

# Final Quality Checks
final_checked_results <- final_quality_checks(merged_results)

# Exporting Data
export_file_path <- "final_results.csv"
export_data(final_checked_results, export_file_path)

# Generating Report
report <- generate_report(final_checked_results)

