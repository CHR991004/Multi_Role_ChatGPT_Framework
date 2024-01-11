setwd("H:/chat_project/drugGPT/6. auto_array/")
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

root_dir <- "E:/0.HR/1_tasks/past/0.auto_array"
input_file_name <- "D2D.xlsx"
ncol <- 1
start_number <- 1
end_number <- 10
output_file_name <- "D2D_group.xlsx"



results<-GEO_auto(root_dir = root_dir,
                  input_file_name = input_file_name,
                  ncol = ncol,
                  start_number = start_number,
                  end_number = end_number,
                  output_file_name = output_file_name,
                  translate = T)
