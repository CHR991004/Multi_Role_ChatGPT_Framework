library(limma)
library(impute)
library(stringr)
library(biomaRt)
library(openxlsx)
library(tidyr)
rm(list = ls())

logFoldChange=0
adjustP=0.05

infor_matrix <- read.xlsx("./file/final_1215_array.xlsx")
matrix_ID <- infor_matrix$drug_id
output_name <- infor_matrix$drug_id

for (i in 1:length(matrix_ID)) {
  HBEXP_ID <- matrix_ID[i]
  infor_location <- i
  Herb_name <- infor_matrix$pert_name[infor_location]
  con_group <- str_split(string = infor_matrix$ctrl_ids[infor_location] ,pattern = "[|]",simplify = T)[1,]
  con_group <-  con_group[con_group!=""]
  drug_group <- str_split(string = infor_matrix$pert_ids[infor_location] ,pattern = "[|]",simplify = T)[1,] 
  drug_group <-  drug_group[drug_group!=""]
  if ((length(con_group)==1)|(length(drug_group)==1)) {
    next
  }
  setwd("H:\\chat_project\\drugGPT\\6. auto_array\\file\\genematrix/")
  rt=read.table(paste0(HBEXP_ID,"_gene_matrix.txt"),sep="\t",header=T,check.names = F,fill = T)
  rt <- rt[!is.na(rt$genes),]
  # rt <- na.omit(rt)
  rt<-rt[,colnames(rt)!="UniGene symbol"]
  con_group <- con_group[con_group%in%colnames(rt)]
  drug_group <- drug_group[drug_group%in%colnames(rt)]
  rt = as.matrix(rt)
  rownames(rt) = rt[,1]
  exp =rt[,2:ncol(rt)]
  dimnames = list(rownames(exp),colnames(exp))
  exp = matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  exp <- exp[,c(drug_group,con_group)]
  exp <- na.omit(exp)
  exp <- subset(exp, rowMeans(exp) > 0)
  exp <- exp[apply(exp, 1, var) > 0, ]
  #impute missing expression data
  mat=impute.knn(exp)
  rt=mat$data
  rt=avereps(rt)
  #normalize
  setwd("H:\\chat_project\\drugGPT\\6. auto_array\\file\\diffsig/")
  rt=normalizeBetweenArrays(as.matrix(rt))
  
  qx <- as.numeric(quantile(rt, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0)
  if (LogC) { rt[which(rt <= 0)] <- NaN
  rt <- log2(rt) }
  #differential
  class <- c(rep("ITP",length(drug_group)),rep("nom",length(con_group)))    #??Ҫ?޸?
  design <- model.matrix(~0+factor(class))
  colnames(design) <- c("ITP","nom")
  fit <- lmFit(rt,design)
  cont.matrix<-makeContrasts(ITP-nom,levels=design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  allDiff=topTable(fit2,adjust='fdr',number=200000)
  allDiff <- cbind(genes=rownames(allDiff),allDiff)
  allDiff <- na.omit(allDiff)
  #write table
  logFoldChange = mean(abs(allDiff$logFC))+sd(abs(allDiff$logFC))
  diffSig <- allDiff[with(allDiff, (abs(logFC)>logFoldChange & P.Value < adjustP )), ]
  diffSig<-separate_rows(data = diffSig,genes,sep = "///")
  write.table(diffSig,file=paste(output_name[i],".xls",sep = ""),sep="\t",quote=F,row.names = F)
  setwd("H:\\chat_project\\drugGPT\\6. auto_array\\file\\diffall")
  write.table(allDiff,file=paste(output_name[i],"_alldiff.xls",sep = ""),sep="\t",quote=F,row.names = F)
  
}
