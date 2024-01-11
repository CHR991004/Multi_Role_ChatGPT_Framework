# BiocManager::install("DESeq2")
rm(list = ls())

library(DESeq2)
library(limma)
library(edgeR)
library(ggplot2)
library(dplyr)
library(tinyarray)

logFC_cutoff = 0.5 # Filter the minimal logFC
P_adj = F 
adj_Pvalue = 0.05 # if P_adj = F, the adj_Pvalue just means P_value

infor_matrix <- read.table("E:\\0.HR\\1_tasks\\2021\\herb\\RNA-seq\\RNA-seq.txt",sep = "\t",
																											header = T,check.names = F)
matrix_now <- list.files("E:\\0.HR\\1_tasks\\2021\\herb\\RNA-seq\\gene_matrix")
matrix_ID <- str_split(matrix_now,".txt",simplify = T)[,1]
setwd("E:\\0.HR\\1_tasks\\2021\\herb\\RNA-seq\\gene_matrix")

for (i in 1:length(matrix_ID)) {
	HBEXP_ID <- matrix_ID[i]
	infor_location <- as.numeric(grep(pattern = HBEXP_ID,x = infor_matrix$GSE_id))
	Herb_name <- infor_matrix$`Herb/ingredient_name`[infor_location]
	con_group <- str_split(string = infor_matrix$Control_samples[infor_location] ,pattern = "; ",simplify = T)
	con_group <-  con_group[con_group!=""]
	drug_group <- str_split(string = infor_matrix$Treatment_samples[infor_location] ,pattern = "; ",simplify = T)[1,] 
	drug_group <-  drug_group[drug_group!=""]
	setwd("E:\\0.HR\\1_tasks\\2021\\herb\\RNA-seq\\gene_matrix")
	rt <- read.table(matrix_now[i],sep = "\t",header = T,check.names = F)
	colnames(rt)[1]<-"genes"
	rt <- aggregate(.~genes,mean,data=rt)
	rownames(rt) <- rt[,1]
	rt <- rt[,-1]
	exprSet <- round(rt)
	
	group_list <- factor(c(rep("N",length(con_group)), rep("T",length(drug_group))))
	colData <- data.frame(row.names=colnames(exprSet), group_list=group_list)
	dds <- DESeqDataSetFromMatrix(countData = exprSet,
																															colData = colData,
																															design = ~ group_list)
	
	dds <- DESeq(dds)
	res <- results(dds, contrast=c("group_list","T","N"))
	# DEG2=as.data.frame(res2)
	# DEG2 = na.omit(DEG2)
	DEG=as.data.frame(res)
	DEG = na.omit(DEG)
	
	# logFC_cutoff = 1
	# P_adj = F
	# adj_Pvalue =0.05
	
	if (P_adj == T) {
		down = DEG[(DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff),]
		up = DEG[(DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff),]
	}else{
		down = DEG[(DEG$pvalue < 0.05)&(DEG$log2FoldChange < -logFC_cutoff),]
		up = DEG[(DEG$pvalue < 0.05)&(DEG$log2FoldChange > logFC_cutoff),]
	}
	sigdiff <- rbind(up,down)
	sigdiff <- cbind(rownames(sigdiff),sigdiff)
	colnames(sigdiff)[c(1,3)]<-c("genes","logFC")
	DEG2 <- cbind(rownames(DEG),DEG)
	colnames(DEG2)[1] <- "genes"
	sigdiff <- sigdiff[,c(1,3,2,4,5,6,7)]
	DEG2 <- DEG2[,c(1,3,2,4,5,6,7)]
	setwd("E:\\0.HR\\1_tasks\\2021\\herb\\RNA-seq\\diff")
	write.table(sigdiff,"diffSig.xls",sep = "\t",row.names = F,col.names = T)
}
