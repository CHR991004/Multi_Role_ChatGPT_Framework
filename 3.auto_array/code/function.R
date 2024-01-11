
# 
require_CRAN <- function(package_name){
  if(require(package_name,character.only = T)){
    print(paste(package_name," is loaded correctly",sep = ""))
  }else {
    print(paste("trying to install ",package_name,sep = ""))
    install.packages(package_name)
    if(require(package_name,character.only = T)){
      print(paste(package_name, " installed and loaded",sep = ""))
    } else {
      stop(paste("could not install ",package_name,sep = ""))
    }
  }
}

require_bioc <- function(package_name){
  if(require(package_name,character.only = T)){
    print(paste(package_name," is loaded correctly",sep = ""))
  }else {
    print(paste("trying to install ",package_name,sep = ""))
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    
    BiocManager::install(package_name)
    if(require(package_name,character.only = T)){
      print(paste(package_name, " installed and loaded",sep = ""))
    } else {
      stop(paste("could not install ",package_name,sep = ""))
    }
  }
}
youdao_translate<-function(abstract_char){
  abstract_char <- enc2utf8(abstract_char)
  Encoding(abstract_char)<-"UTF-8"
  if (nchar(abstract_char) > 20 ) {
    input <- paste(str_sub(abstract_char,1,10),
                   nchar(abstract_char),
                   str_sub(abstract_char,(str_length(abstract_char)-9),
                           str_length(abstract_char)),
                   sep = "")
  }else{
    input <- abstract_char
  }
  url <- "https://openapi.youdao.com/api"
  appKey <- "1651b49eba901xxx" # 
  secret_key <- "CugoLiczjbkY7r4L6lzu8htr5xxxxxx" # 
  salt <- UUIDgenerate(use.time = F)
  curtime <- as.character(as.integer(Sys.time()))
  sign_char <- paste(appKey,input,salt,curtime,secret_key,sep = "") 
  Encoding(sign_char)<-"UTF-8"
  sign <- hash_sha256(sign_char)
  res = POST(url = url,
             body = list(q = abstract_char,
                         from = "en",
                         to ="zh-CHS",
                         appKey = appKey,
                         salt = salt,
                         sign = sign,
                         signType = "v3",
                         curtime = curtime,
                         strict = "true"
             ),
             encode = "form"
  )
  
  
  r2 <- rawToChar(res$content) 
  r2 <- str_conv(r2,"UTF-8")
  r3 <- fromJSON(r2)
  output <- r3[["translation"]]
  return(output)
}
GEO_auto <- function(root_dir, input_file_name, ncol,
                     start_number, end_number,
                     output_file_name,translate){

  file_dir<-paste(root_dir,"/file",sep = "")
  setwd(file_dir)
  environment_check <- try(setwd(paste(getwd(),"/series", sep = "")),silent = T)
  if ("try-error" %in% class(environment_check)) {
    dir.create(paste(getwd(),"/series", sep = ""))
    ser_dir <- paste(getwd(),"/series", sep = "")
    dir.create(paste(getwd(),"/GPL", sep = ""))
    GPL_dir <- paste(getwd(),"/GPL", sep = "")
    setwd(file_dir)
    print("Environment fix successfully!")
  }else{
    setwd(file_dir)
    ser_dir <- paste(getwd(),"/series", sep = "")
    GPL_dir <- paste(getwd(),"/GPL", sep = "")
    print("Environment check successfully!")
  }
  setwd(file_dir)
  xlsx_file <- openxlsx::read.xlsx(input_file_name,sheet = 1)
  for (n in xlsx_file[,ncol][start_number:end_number]) { #delete
    repeat{
      group_definition_matrix<-matrix(ncol = 2)
      GSE_list <- getGEO(GEO = n,destdir = ser_dir,getGPL = F)
      sam_data_matrix <- GSE_list[[1]]@phenoData@data # all sample data
      
      # names 
      organism_name <- grep(x = names(sam_data_matrix),pattern = "organism",fixed = T,value = T)
      EX_type_name <- GSE_list[[1]]@experimentData@other[["type"]]
      source_name <- grep(x = names(sam_data_matrix),pattern = "source_name",fixed = T,value = T)
      title_name <- grep(names(sam_data_matrix),pattern = "title",value = T,fixed = T)
      sam_char <- grep(names(GSE_list[[1]]@phenoData@data),pattern = ":",value = T,fixed = T)
      abstract_char <- GSE_list[[1]]@experimentData@abstract
      overall_design_char <- GSE_list[[1]]@experimentData@other[["overall_design"]]
      
      # read translation
      print("Abstract:")
      print(abstract_char)
      if (translate == T) {
        print(youdao_translate(abstract_char))
      }
      readline("Press ENTER to continue!")
      print("Overall design:")
      print(overall_design_char)
      if (translate == T) {
        print(youdao_translate(overall_design_char))
      }
      readline("Press ENTER to continue!")
      
      # choose characteristic
      all_choice <- c(source_name,title_name,sam_char)
      for (i in all_choice) {
        col_data <- sam_data_matrix[,i]
        choice_value <- unique(col_data)
        cat(paste("    The ", i ," levels is:\n", 
                  paste(choice_value,collapse = "\n"),sep = ""))
        choice_1 <- readline(
          "press 'the key words' to choose it, and press 0 to uncheck.")
        if (choice_1 !=0) {
          repeat_choice<- readline(
            "Please enter discarded fields(If no, enter 0):")
          if (repeat_choice != 0) {
            repeat_GSM_1 <- sam_data_matrix[grep(pattern = repeat_choice,
                                                x = col_data,fixed = T),
                                           "geo_accession"]
            group_GSM_1 <- sam_data_matrix[grep(pattern = choice_1,
                                                x = col_data,fixed = T),
                                           "geo_accession"]
            group_GSM_1 <- setdiff(group_GSM_1,repeat_GSM_1)
          }else{
            group_GSM_1 <- sam_data_matrix[grep(pattern = choice_1,
                                                x = col_data,fixed = T),
                                           "geo_accession"]
          }
          
          group_type_1 <- readline("Control: 0, test:1, Then press enter:")
          choice_2 <- readline("Please input the other key word:")
          
          repeat_choice_2<- readline(
            "Please enter discarded fields(If no, enter 0):")
          if (repeat_choice_2 != 0) {
            repeat_GSM_2 <- sam_data_matrix[grep(pattern = repeat_choice_2,
                                                 x = col_data,fixed = T),
                                           "geo_accession"]
            group_GSM_2 <- sam_data_matrix[grep(pattern = choice_2,
                                                x = col_data,fixed = T),
                                           "geo_accession"]
            group_GSM_2 <- setdiff(group_GSM_2,repeat_GSM_2)
          }else{
            group_GSM_2 <- sam_data_matrix[grep(pattern = choice_2,
                                                x = col_data,fixed = T),
                                           "geo_accession"]
          }
          
          if (group_type_1 == 0) {
            group_definition_matrix <- rbind(group_definition_matrix,
                                             cbind(group_GSM_1,
                                                   rep("control",
                                                       length(group_GSM_1))))
            group_definition_matrix <- rbind(group_definition_matrix,
                                             cbind(group_GSM_2,
                                                   rep("test",
                                                       length(group_GSM_2))))
          }else if(group_type_1 == 1){
            group_definition_matrix <- rbind(group_definition_matrix,
                                             cbind(group_GSM_1,
                                                   rep("test",
                                                       length(group_GSM_1))))
            group_definition_matrix <- rbind(group_definition_matrix,
                                             cbind(group_GSM_2,
                                                   rep("control",
                                                       length(group_GSM_2))))
          }
        }else{
          next
        }
      }
      
      group_definition_matrix<-na.omit(group_definition_matrix)
      control_GSM <- group_definition_matrix[group_definition_matrix[,2]=="control",1]
      control_list <- paste(control_GSM, collapse = ";")
      test_GSM <-group_definition_matrix[group_definition_matrix[,2]=="test",1]
      test_list <- paste(test_GSM,collapse = ";")
      control_GPL <- sam_data_matrix %>%
        .[.[,"geo_accession"]%in%control_GSM,"platform_id"]
      test_GPL <- sam_data_matrix %>%
        .[.[,"geo_accession"]%in%test_GSM,"platform_id"]
      GPL_list <- paste(unique(c(control_GPL,test_GPL)),collapse = ";")
      
      control_organism <- sam_data_matrix %>%
        .[.[,"geo_accession"]%in%control_GSM,organism_name]
      test_organism <- sam_data_matrix %>%
        .[.[,"geo_accession"]%in%test_GSM,organism_name]
      organism_list <- paste(unique(c(control_organism,test_organism)),collapse = ";")

      
      disturbance <- readline("Input the Omics disturbance:")
      time_col <- readline("Input the time:")
      perturbation_type_col <- readline("Disease or drug,0 for disease,1 for drug:")
      if (perturbation_type_col == 0) {
        perturbation_type <- "disease"
        drug_dose_col <- "NA"
      }else{
        perturbation_type <- "drug"
        drug_dose_col <- readline("Input the drug dose:")
      }
      
      tissue_col <- readline("Input the tissue:")
      cell_line_col <- readline("Input the cell line:")
      one_row <- c(n,organism_list,EX_type_name,GPL_list,disturbance,
                   perturbation_type,drug_dose_col,time_col,
                   tissue_col,cell_line_col,control_list,test_list)

      output_matrix <- matrix(data = one_row, ncol = 12, 
                              dimnames = list(NULL,c("GSE","organism",
                                                     "ex_type","GPL",
                                                     "perturbation","perturbation_type",
                                                     "drug_dose","time","tissue",
                                                     "cell_line","control_GSM",
                                                     "test_GSM")))
      output_matrix <- as.data.frame(output_matrix)
      output_check <- try(read.xlsx(xlsxFile = output_file_name,sheet = 1,
                                    check.names = F),silent = T)
      if ("try-error" %in% class(output_check)) {
        openxlsx::write.xlsx(output_matrix,file = output_file_name,
                             asTable = F,colnames = T)
      }else{
        output_matrix_final <- rbind(output_check,output_matrix)
        openxlsx::write.xlsx(output_matrix_final, file = output_file_name,
                             asTable = F)
      }
      repeat_logi<- readline(
        "Do you want to continue analyzing this GSE? 1 to continue and 0 to next.")
      if (repeat_logi==1) {
        print("Analysis once again!")
      }else{
        break
      }
    }
  }
}

getDirListing <- function(url) {
  # Takes a URL and returns a character vector of filenames
  a <- xml2::read_html(url)
  fnames = grep('^G',xml_text(xml_find_all(a,'//a/@href')),value=TRUE)
  return(fnames)
}
