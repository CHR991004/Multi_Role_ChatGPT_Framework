rm(list = ls())
root_dir <- "E:/0.HR/1_tasks/past/0.auto_array"#需修改
setwd(root_dir) 
source(paste(getwd(),"/code/function.R",sep = ""),encoding = "UTF-8")
Sys.setlocale(category = "LC_ALL",locale = "English_United States.1252")
options( 'download.file.method.GEOquery' = 'libcurl' )
readr::local_edition(1)

require_CRAN("httr")
require_CRAN("jsonlite")
require_CRAN("stringr")
require_CRAN("uuid")
require_CRAN("request")
require_CRAN("cli")
require_bioc("GEOquery")
require_CRAN("openxlsx")


file_dir<-paste(root_dir,"/file",sep = "")
input_file_name <- "D2D_group_check.xlsx"
setwd(file_dir)
input<-openxlsx::read.xlsx(input_file_name,sheet = 1)
ser_dir <- paste(getwd(),"/series", sep = "")
GPL_dir <- paste(getwd(),"/GPL", sep = "")
setwd(ser_dir)
error_list<-c()

for (i in 1:nrow(input)) {
  GSE<-input[i,1]
  control_GSM<-input[i,"control_GSM"]
  control_GSM <- strsplit(control_GSM,split = ";")[[1]]
  test_GSM<-input[i,"test_GSM"]
  test_GSM <- strsplit(test_GSM,split = ";")[[1]]
  filename<-grep(paste0(GSE,"_"),x = list.files(path = ser_dir),value = T)
  destfile = file.path(ser_dir, filename)
  GSE_data <- GEOquery:::parseGSEMatrix(destfile, destdir = ser_dir, 
                                        AnnotGPL = F, getGPL = F)$eset
  sam_data_matrix <- GSE_data@phenoData@data # all sample data
  
  organism_name <- grep(x = names(sam_data_matrix),pattern = "organism",fixed = T,value = T)
  EX_type_name <- GSE_data@experimentData@other[["type"]]
  source_name <- grep(x = names(sam_data_matrix),pattern = "source_name",fixed = T,value = T)
  title_name <- grep(names(sam_data_matrix),pattern = "title",value = T,fixed = T)
  sam_char <- grep(names(GSE_data@phenoData@data),pattern = ":",value = T,fixed = T)
  abstract_char <- GSE_data@experimentData@abstract
  overall_design_char <- GSE_data@experimentData@other[["overall_design"]]
  all_choice <- c(source_name,title_name,sam_char)
  col_data <- sam_data_matrix[,all_choice]
  control_GSM_data <- col_data[control_GSM,]
  test_GSM_data <- col_data[test_GSM,]
  View(control_GSM_data)
  control_check <- readline("Control group ok? 1 for ok and 0 for reanalysis: ")
  if (control_check == 0) {
    error_list<-i
  }
  View(test_GSM_data)
  test_check <- readline("Test group ok? 1 for ok and 0 for reanalysis: ")
  if (test_check == 0) {
    error_list<-i
  }
  setwd(file_dir)
  if (length(error_list>0)) {
    write.table(error_list,"error_GSE_check_result.txt",
                sep = "\t",col.names = F,row.names = F,quote = F,append = T)
  }
  
}