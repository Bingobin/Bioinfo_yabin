##Part Library: Loading R package
library(ggplot2)
library(cowplot) 
#library(customLayout)
#library(dplyr)
library(tidyr)
library(forcats)
library(UpSetR)
library(RColorBrewer)
library(DESeq2)
library(survival)
library(survminer)
library(survivalROC)
library(glmnet)
library(MAGeCKFlute)
#library(org.Hs.eg.db)
library(clusterProfiler)
library(WGCNA)

blank <- theme(panel.border = element_blank(), 
               panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.background = element_blank()
               )
#+ theme(panel.border = element_rect(fill = NA, linetype = 1, size = 1), axis.line = element_blank())


##Part Subplement: Build Function

###01.Function of DESeq2 DEG analysis
DESeq2_DEG_analysis <-function(count.matrix,design.df){
  #SampleID Group
  dds <- DESeqDataSetFromMatrix(countData = count.matrix, colData = design.df, design =~ Group)
  dds <- DESeq(dds)
  res <- results(dds, alpha = 0.05)
  inter<-resultsNames(dds)
  resOrdered = res[order(res$padj,na.last = TRUE),]
  #  plotMA(resOrdered, alpha = 0.05)
  deseq2_result <- data.frame(resOrdered)
  deseq2_result <- na.omit(deseq2_result)
  gene_id_anno <- read.table("/Volumes/TCGA-2/BEAT_AML_RNAseq_count/03.merge_and_combat/gene_id_trans.txt", header = TRUE, sep = "\t",stringsAsFactors=FALSE,row.names = 1)
  deseq2_result <- merge(gene_id_anno, deseq2_result, by.x=0, by.y=0, all = F)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  result <- list()
  result[["inter"]] <- inter
  #  result[["mrna"]] <- deseq2_result[deseq2_result$Type=="protein_coding",]
  result[["result"]] <- deseq2_result
  result[["vsd"]] <- vsd
  return(result)
}


###02.Function of extracting up and donw genes form DESeq2 DEG result

DESeq2_DEG_extract <- function(result.df,aml_design,plot){
  deseq2_result <- result.df$result
  result.df[["up"]] <- deseq2_result[deseq2_result$log2FoldChange > log2(1.5) & deseq2_result$padj < 0.05,]
  result.df[["down"]] <- deseq2_result[deseq2_result$log2FoldChange < -log2(1.5) & deseq2_result$padj < 0.05,]
  if(plot){
    #col_anno=data.frame(row.names = rownames(aml_design),Group=aml_design$Group,WT1=aml_design$WT1_tpm)
    col_anno <- data.frame(row.names = aml_design$SampleID,Group=aml_design$Group)
    col_anno_color <- list(Group = c(CD34="#F15A24",APL="#0071BC"))
    EX_data <- assay(result.df$vsd[c(result.df$up$Row.names,result.df$down$Row.names)])
    EX_data <- log2(EX_data+1)
    EX_data.mean <- matrix(rep(apply(EX_data,1,mean),ncol(EX_data)),ncol=ncol(EX_data))
    EX_data.sd <- matrix(rep(apply(EX_data,1,sd),ncol(EX_data)),ncol=ncol(EX_data))
    EX_data <- (EX_data - EX_data.mean) / EX_data.sd
#    EX_data[EX_data > 4] =  2
#    EX_data[EX_data < -4] =  -2
    EX_data[EX_data > 2] =  2
    EX_data[EX_data < -2] =  -2
    result.df[["plot"]] <- pheatmap::pheatmap(EX_data,scale="none",
                                              show_colnames = T,
                                              show_rownames = F,
                                              annotation_col=col_anno,
                                              annotation_colors = col_anno_color,
                                              cluster_cols  = T,
                                              cluster_rows = T,
                                              clustering_method = "complete", 
                                              color = colorRampPalette(rev(brewer.pal(n = 11, name ="PRGn")))(100)
    )
  }
  return(result.df)
}


###03.Function of risk groping based on multivariable cox model

RiskScore_multivar_cox <- function(sample.list, gene.list, feature.list,tpm.matrix, clinic.df, cox.df = NULL) {
  #  sample.list <- aml_all_wt1.clinc.rm[aml_all_wt1.clinc.rm$Source.x == "BEAT",]$SampleID
  #  gene.list <- gene_sig.list
  #  tpm.matrix <- aml_all_tpm.matrix
  #  clinic.df <- aml_all_wt1.clinc
  #  feature.list <- c("WT1","TET2","IDH1","IDH2","ageAtDiagnosis")
  
  clinic.df$ageAtDiagnosis <- as.numeric(clinic.df$ageAtDiagnosis)
  sample.list <- gsub("-","\\.", sample.list)  
  col_index <- match(sample.list,colnames(tpm.matrix))
  row_index <- match(gene.list, tpm.matrix$Symbol)
  risk.tpm <- tpm.matrix[row_index,col_index,]
  rownames(risk.tpm) <- tpm.matrix$Symbol[row_index]
  risk.tpm <- t(risk.tpm)
  rownames(risk.tpm) <- gsub("\\.","-",rownames(risk.tpm))
  risk.cox <- merge(risk.tpm, clinic.df, by.x=0,by.y =2, all=FALSE)
  risk.cox[risk.cox$vitalStatus  == "Death",]$vitalStatus <- 1
  risk.cox[risk.cox$vitalStatus  == "Alive",]$vitalStatus <- 0
  risk.cox$vitalStatus <- as.numeric(risk.cox$vitalStatus)
  risk.cox$overallSurvival <- as.numeric(risk.cox$overallSurvival)
  risk.cox <- risk.cox %>% filter(overallSurvival > 0)
  #  survpre_mut <- data.table::as.data.table(risk.cox[,c("WT1","TET2","IDH1","IDH2","ageAtDiagnosis")])
  #  colnames(survpre_mut) <- c("WT1_mut","TET2_mut", "IDH1_mut", "IDH2_mut","Age")
  #  survpre_mut$Age <- as.numeric(survpre_mut$Age)
  survpre_mut <- data.table::as.data.table(risk.cox[,feature.list,drop = FALSE])
  survpre_sig <- risk.cox[,gene.list,drop = FALSE]
  survpre_os <-  risk.cox[,c("overallSurvival","vitalStatus")]
  colnames(survpre_os) <- c("Time","Status")
  surv_model.df <- cbind(survpre_sig,survpre_mut)
  surv_model.df <- as.data.frame(row.names = risk.cox$Row.names, cbind(surv_model.df,survpre_os))
  surv_model.df <- na.omit(surv_model.df)
  survobj <- Surv(surv_model.df$Time,surv_model.df$Status)
  survpre <- as.matrix(surv_model.df[,c(gene.list,feature.list),drop = FALSE])
  model <- coxph( survobj ~ survpre)
  #  if(is.null(cox.df)){
  cox.df.sum <-summary(model)
  cox.df <- as.data.frame(cbind(cox.df.sum$coefficients, cox.df.sum$conf.int))
  cox.df$Feature <- gsub("survpre","",rownames(cox.df))
  cox.df <- data.table::as.data.table(cox.df[,c(10,1,5,6,8,9)])
  colnames(cox.df) <- c("Feature","Coefficient","Pvalue","Hazard.ratio","HR.Lower.95","HR.Upper.95")
  cox.df <- na.omit(cox.df)
  #  }
  ls_df <- lapply(1:nrow(surv_model.df), function(x){
    sample <- rownames(surv_model.df)[x]
    score <-sum(sapply(1:nrow(cox.df),function(i){
      coef = cox.df$Coefficient[i]
      tpm = surv_model.df[x,colnames(surv_model.df) == cox.df$Feature[i]]
      coef * tpm
    }))
    data.frame(Sample = sample, Risk_score=score)
  })
  rs.df <- do.call(rbind, ls_df)
  rs.df$RS_norm<-sapply(1:nrow(rs.df),function(i){
    (rs.df$Risk_score[i] - mean(rs.df$Risk_score))/sd(rs.df$Risk_score)
  })
  rs.df$Risk_Group <- "High_risk"
  rs.df[rs.df$RS_norm < median(rs.df$RS_norm),]$Risk_Group <- "Low_risk"
  rs.df$Sample <- gsub("\\.","-",rs.df$Sample)
  rownames(surv_model.df) <- gsub("\\.","-",rownames(surv_model.df))
  rs.df.clinc <- merge(rs.df,surv_model.df, by.x = 1, by.y =0, all = FALSE)
  rs.df.clinc <- merge(rs.df.clinc, aml_all_wt1.clinc, by.x = 1, by.y =2, all = FALSE)
  result.ls <- list()
  result.ls[["HR"]] <- cox.df
  result.ls[["Clinc"]] <- rs.df.clinc
  return(result.ls)
}


###03.Function of survial of risk groping based on univ cox and TPM

Risk_grouping_cox<-function(gene.list,design.df){
  #gene.list<-risk_gene.list
  #design.df<-Risk_model.beat.df
  
  ls_df <- lapply(match(gene.list,colnames(design.df)), function(i){
    pretty_aml <- data.frame(time   = design.df$Time, 
                             status = design.df$Status, 
                             TPM    = design.df[,i])
    Gene_label <- colnames(design.df)[i]
    model <- coxph( Surv(time, status) ~ TPM, data = pretty_aml)
    cox_result <- summary(model)
    data.frame(row.names = Gene_label, cox_result$conf.int,PV=cox_result$coefficients[5])
  })
  hr.df <- do.call(rbind, ls_df)
  ls_df <- lapply(1:nrow(design.df), function(x){
    sample <- design.df$SampleID[x]
    score <-sum(sapply(1:length(gene.list),function(i){
      coef = log(hr.df[rownames(hr.df) == gene.list[i],]$exp.coef.) 
      tpm = design.df[x,colnames(design.df)==gene.list[i]]
      coef * tpm
    }))
    data.frame(Sample = sample, Risk_score=score)
  })
  rs.df <- do.call(rbind, ls_df)
  rs.df$RS_norm<-sapply(1:nrow(rs.df),function(i){
    (rs.df$Risk_score[i] - mean(rs.df$Risk_score))/sd(rs.df$Risk_score)
  })
  rs.df$Risk_Group <- "High_risk"
  rs.df[rs.df$RS_norm < median(rs.df$RS_norm),]$Risk_Group <- "Low_risk"
  rs.df$Sample <- gsub("\\.","-",rs.df$Sample)
  design.df$SampleID <- gsub("\\.","-",design.df$SampleID)
  rs.df.clinc <- merge(rs.df,design.df, by.x = 1, by.y =1, all = FALSE)
  return(rs.df.clinc)
}



###04.Function of getting risk model dataframe
Risk_model_df<-function(gene.list,design.df,tpm.df){
  #gene.list <- risk_gene.list
  #design.df <- aml_design_TCGA
  #tpm.df <- aml_combind_tpm.df
  index_col<-match(design.df$SampleID,colnames(tpm.df))
  index_row <- match(gene.list,tpm.df$Symbol)
  Risk.df <- tpm.df[index_row,index_col]
  rownames(Risk.df) <- gene.list
  Risk.df <- t(Risk.df)
  Risk.df <- merge(Risk.df,design.df,by.x=0,by.y=1)
  Risk.df<-Risk.df %>% filter(!(is.na(Status) | Status == "Unknown" | Time==0))
  colnames(Risk.df)[1] <- "SampleID"
  Risk.df$Time <- as.numeric(as.character(Risk.df$Time))
  Risk.df$Status <- as.character(Risk.df$Status)
  Risk.df[Risk.df$Status == "Death",]$Status <- 1
  Risk.df[Risk.df$Status == "Alive",]$Status <- 0
  Risk.df$Status <- as.numeric(Risk.df$Status)
  return(Risk.df)
}



###05.Function of survival-related plot basd on risk group

plot_survial_risk<-function(clinic_rs.df){
  #clinic_rs.df<-rs.df.clinc
  blank <- theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank())
  clinic_rs.df$rank <- rank(clinic_rs.df$Risk_score)
  x1<-ggplot(clinic_rs.df,aes(x=rank,y=Risk_score)) + 
    geom_point(aes(color=Risk_Group)) + 
    blank + 
    theme(legend.justification = c("right", "top"), legend.position = c(1,1),legend.title = element_blank()) + 
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(), axis.text.y = element_text(angle=90,hjust=0.5))  +
    #    scale_x_continuous(breaks=c(dim(clinic_rs.df)[1]*0.25, dim(clinic_rs.df)[1]*0.75),labels=c("Low_risk","High_risk")) +
    geom_vline(xintercept = dim(clinic_rs.df)[1]*0.5,linetype=5) +
    scale_color_manual(values=c(c("#FF7F00","#1F78B4"))) +
    ylab("Risk score") + 
    theme(panel.border = element_rect(linetype = 1, size = 1, fill = NA), axis.line = element_blank(),axis.text.x = element_blank()) + 
    theme(legend.background = element_rect(fill=NA), legend.key = element_rect(fill=NA))
  
  x2<-ggplot(clinic_rs.df,aes(x=rank,y=Time)) + 
    geom_point(aes(color=vitalStatus,shape=vitalStatus)) + 
    blank + 
    theme(legend.justification = c("right", "top"), legend.position = c(1,1),legend.title = element_blank()) + 
    theme(axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.text.y = element_text(angle=90,hjust=0.5)) +
    scale_x_continuous(breaks=c(dim(clinic_rs.df)[1]*0.25, dim(clinic_rs.df)[1]*0.75),labels=c("Low_risk","High_risk")) +
    geom_vline(xintercept = dim(clinic_rs.df)[1]*0.5,linetype=5) +
    scale_color_manual(values=c("grey40","black")) +
    scale_shape_manual(values=c(1,2)) +
    ylab("Survial days") + 
    theme(panel.border = element_rect(linetype = 1, size = 1, fill = NA), axis.line = element_blank())+ 
    theme(legend.background = element_rect(fill=NA),legend.key = element_rect(fill=NA))
  
  return(plot_grid(x1,x2,ncol=1))
}


###06.Function of Riks model ROC plot 

plot_roc_curve<-function(risk.df){
  #  risk.df <- rs.df.clinc
  #  risk.df$lp <- a
  tmp.df.1 <- survivalROC(Stime = risk.df$Time,
                          status = risk.df$Status,
                          marker = risk.df$Risk_score,
                          predict.time = 365*1,
                          span = 0.1*nrow(risk.df)^(-0.30)
                          #                 method="KM"
  )
  tmp.df.3 <- survivalROC(Stime = risk.df$Time,
                          status = risk.df$Status,
                          marker = risk.df$Risk_score,
                          predict.time = 365*3,
                          span = 0.1*nrow(risk.df)^(-0.30)
                          #                 method="KM"
  )
  tmp.df.5 <- survivalROC(Stime = risk.df$Time,
                          status = risk.df$Status,
                          marker = risk.df$Risk_score,
                          predict.time = 365*5,
                          span = 0.1*nrow(risk.df)^(-0.30)
                          #                 method="KM"
  )
  #  str(tmp.df)
  auc1 <- paste0("1 years", "(AUC=", round(tmp.df.1$AUC,4), ")")
  auc3 <- paste0("3 years", "(AUC=", round(tmp.df.3$AUC,4), ")")
  auc5 <- paste0("5 years", "(AUC=", round(tmp.df.5$AUC,4), ")")
  tmp.df.plot.1 <- data.frame(FP=sort(tmp.df.1$FP),TP=sort(tmp.df.1$TP), PT=rep(auc1,length(tmp.df.1$TP)))
  tmp.df.plot.3 <- data.frame(FP=sort(tmp.df.3$FP),TP=sort(tmp.df.3$TP), PT=rep(auc3,length(tmp.df.3$TP)))
  tmp.df.plot.5 <- data.frame(FP=sort(tmp.df.5$FP),TP=sort(tmp.df.5$TP), PT=rep(auc5,length(tmp.df.5$TP)))
  #  tmp.df.plot.1 <- data.frame(FP=tmp.df.1$FP,TP=tmp.df.1$TP, PT=rep(auc1,length(tmp.df.1$TP)))
  #  tmp.df.plot.3 <- data.frame(FP=tmp.df.3$FP,TP=tmp.df.3$TP, PT=rep(auc3,length(tmp.df.3$TP)))
  #  tmp.df.plot.5 <- data.frame(FP=tmp.df.5$FP,TP=tmp.df.5$TP, PT=rep(auc5,length(tmp.df.5$TP)))
  tmp.df.plot <- rbind(tmp.df.plot.1, tmp.df.plot.3)
  tmp.df.plot <- rbind(tmp.df.plot, tmp.df.plot.5)
  
  
  ggplot(tmp.df.plot,aes(x=FP,y=TP)) + 
    geom_line(aes(color=PT),size = 0.5) + blank + 
    geom_abline(slope = 1,linetype=2, size = 0.5) +
    xlab("1-Specificity") + ylab("Sensitivity") + 
    scale_color_manual(values = brewer.pal(9,"Set1")[c(2,4,5)]) + 
    theme(legend.key = element_blank(),legend.title = element_blank(), 
          legend.position = c(.95, .05), legend.justification = c("right", "bottom")) + 
    theme(panel.border = element_rect(linetype = 1, size = 0.8, fill = NA), 
          axis.line = element_blank())
}


###07.Build function of Volcano Plot

volcano_plot_Deseq2<- function(deseq2_result.df,target_gene.list,gene.list, pv = 0.05, fc = 1.5){
  #  deseq2_result.df <- WT1_MUT.deg$result
  #  target_gene.list <- WT1_status.gene$Symbol
  #  gene.list <- risk_gene.list
  #  pv <- 0.05
  #  fc <- 1.5
  gg <- deseq2_result.df
  gg <- gg[match(target_gene.list,gg$Symbol),]
  gg <-na.omit(gg)
  index <- match(gene.list,gg$Symbol)
  index <- na.omit(index)
  gg$group = "no"
  gg[gg$log2FoldChange > log2(fc) & gg$padj < pv,]$group <- "up"
  gg[gg$log2FoldChange < -log2(fc) & gg$padj < pv,]$group <- "down"
  gg$color <- gg$group
  gg$color[index] <- "black"
  #mycolour = c("grey", "#B30000", "#08519C", "black")
  mycolour = c("grey", "#810F7C", "#006D2C", "black")
  names(mycolour) = c("no", "up", "down", "black")
  gg$label <- ""
  gg$label[index] <- gg$Symbol[index]
  gg[gg$group == "no",]$label <- ""
  p <- ggplot(gg, aes(x = log2FoldChange, y = -log10(padj)))
  p <- p + geom_point(aes(fill=group), shape = 21, alpha = 0.6, show.legend = FALSE)
  p <- p + geom_point(aes(colour = color), shape = 21, alpha = 0.6, show.legend = FALSE)
  p = p + scale_color_manual(values = mycolour)
  p = p + scale_fill_manual(values = mycolour)
  p = p + geom_hline(yintercept = -log10(pv), linetype = "dotted")
  p = p + geom_vline(xintercept = c(-log2(fc), log2(fc)), linetype = "dotted")
  p = p + ggrepel::geom_text_repel(aes(label=label,color=group),show.legend = FALSE,fontface = "bold", size = 2.5,box.padding = unit(0.8, "lines"),point.padding = unit(0.3, "lines"), segment.size = 0.3)
  p = p + blank
  p = p + theme(panel.border = element_rect(fill = NA, linetype = 1,size = 1 ), axis.line = element_blank())
  #  p= p + xlim(-4,4) + ylim(0,10)
  return(p)
}


###08.Build function of enrich result combind  
enrich_combind <- function(gene,pvc,qvc){
  wp2gene <- read.gmt.wp("~/bin/gmt/wikipathways-20190510-gmt-Homo_sapiens.gmt")
#  wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
  wpid2gene <- wp2gene %>% dplyr::select(wpid, gene)
  wpid2name <- wp2gene %>% dplyr::select(wpid, name)
  ego_bp <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = pvc, qvalueCutoff  = qvc, readable = TRUE)
  ego_cc <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "CC", pAdjustMethod = "BH", pvalueCutoff  = pvc, qvalueCutoff  = qvc, readable = TRUE)
  ego_mf <- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID', ont = "MF", pAdjustMethod = "BH", pvalueCutoff  = pvc, qvalueCutoff  = qvc, readable = TRUE)
  ewp <- enricher(gene, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 0.5,qvalueCutoff  = qvc)
  ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  ekg<- enrichKEGG(gene, organism = "hsa",pvalueCutoff = 0.5,qvalueCutoff  = qvc)
  ekg <- setReadable(ekg, org.Hs.eg.db, keyType = "ENTREZID")
  tmp <- list()
  tmp[["ego_bp"]] <- ego_bp@result
  tmp[["ego_cc"]] <- ego_cc@result
  tmp[["ego_mf"]] <- ego_mf@result
  tmp[["ewp"]] <- ewp@result
  tmp[["ekg"]] <- ekg@result
  return(tmp)
}

###09. Build function of enricher_plot
enricher_plot <- function(enricher,bp,cc,mf,wp,kg,value){
  df <- data.frame(Description=enricher$ego_bp[bp,2],Pvalue=enricher$ego_bp[bp,value],Type=rep("GO_BP",length(bp)))
  df.tmp <- data.frame(Description=enricher$ego_cc[cc,2],Pvalue=enricher$ego_cc[cc,value],Type=rep("GO_CC",length(cc)))
  df <- rbind(df,df.tmp)
  df.tmp <- data.frame(Description=enricher$ego_mf[mf,2],Pvalue=enricher$ego_mf[mf,value],Type=rep("GO_MF",length(mf)))
  df <- rbind(df,df.tmp)
  df.tmp <- data.frame(Description=enricher$ewp[wp,2],Pvalue=enricher$ewp[wp,value],Type=rep("WikiPath",length(wp)))
  df <- rbind(df,df.tmp)
  df.tmp <- data.frame(Description=enricher$ekg[kg,2],Pvalue=enricher$ekg[kg,value],Type=rep("KEGG",length(kg)))
  df <- rbind(df,df.tmp)
  #ggplot(df,aes(x=Description,y=-log(Pvalue),fill=Type)) + geom_bar(stat = "identity",show.legend = TRUE) + coord_flip()
  p<-ggplot(df,aes(-log(Pvalue),fct_reorder(Description, -log(Pvalue)))) +geom_segment(aes(xend=0, yend = Description,color=Type),linetype = 2,show.legend = FALSE) + geom_point(aes(color=Type),size=5,show.legend = FALSE) + scale_color_manual(values = brewer.pal(5,"Set2")) + facet_grid(Type~.,scales = 'free',space = 'free_y', switch = "x") + blank + ylab("")
  return(p)
}

###10. Build function of paired gene expression correlation
ccor_paired_gene <- function(EXdata,Gene1,Gene2){
  #  EXdata<-aml_all_tpm.matrix
  #  Gene1 <- "TET2"
  #  Gene2 <- "TP53"
  ##########################
  #           Symbol Length           Type    BA2000    BA2003      BA2004
  #TSPAN6     TSPAN6   4535 protein_coding   0.00000  1.016602   0.3353875
  #TNMD         TNMD   1610 protein_coding   0.00000  0.000000   0.0000000
  ########################
  corr.df <- data.frame(Gene1=log2(as.numeric(EXdata[Gene1,2:ncol(EXdata)])+1),Gene2=log2(as.numeric(EXdata[Gene2,2:ncol(EXdata)])+1))
  corr.result <- cor.test(~Gene1 + Gene2,corr.df,method="pearson")
  text<-paste("k=",round(lm( Gene1~Gene2, corr.df)$coefficients[2],4),",","r=",round(corr.result$estimate,4), sep = "")
  pv <- corr.result$p.value
  if (pv < 0.05 & pv > 0.01){sig="*"}else if(pv < 0.01 & pv > 0.001){sig="**"}else if(pv < 0.001 & pv > 0.0001){sig="***"}else if(pv < 0.0001){sig="****"}else{sig="na"}
  blank <- theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), panel.background = element_blank())
  ggplot(corr.df,aes(y=Gene2,x=Gene1)) + 
    geom_point(size = 1, alpha = 0.5, color = "grey") + 
    geom_smooth(method = "lm", se = TRUE,formula = y ~ x) + ggtitle(text)  +  
    annotate("text",label=sig, 
             y=max(corr.df$Gene2,na.rm = TRUE),
             x=median(corr.df[corr.df$Gene1 != 0,]$Gene1,na.rm = TRUE), size=7)  + 
    xlab(Gene1) + ylab(Gene2) + blank
}

###11. Build function of H3K27ac enhancer rank plot
enhancer_rank_plot <- function(label, sample, top, rank.df){
###############################################  
  #  Sample    Signal Type Rank CLOSEST_GENE
  #  CHH 555644.95   SE    1         ETV6
  #  CHH 429203.43   SE    2       FNDC3B
  #  CHH 235318.70   SE    3        RREB1
  #  CHH 221833.55   SE    4      COL23A1
  #  CHH 198674.49   SE    5   BZRAP1-AS1
  #  CHH 185526.61   SE    6        CELF2
##############################################
  #  sig_gene <- "WT1"
  #  sample <- "CHH"
  #  top <- 5
  enhancer_signal <- rank.df %>% filter(Sample == sample)
  p<-ggplot(enhancer_signal, aes(x = Rank/1000, y= Signal)) + 
    geom_point(aes(color=Type), alpha=1,size=0.5)  + 
    xlab("Rank/1000") + ylab("Enhancer signal") + ggtitle(paste0(sample," H3K27ac")) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    blank +
    scale_color_manual(values=c(brewer.pal(n = 9, name ="Set1")[1], brewer.pal(n = 9, name ="Set1")[9])) +
    geom_vline(xintercept=nrow(enhancer_signal[enhancer_signal$Type == "SE", ])/1000,lty=4,col="black",lwd=0.8)
  
  
  p <- p + ggrepel::geom_label_repel(data = enhancer_signal[enhancer_signal$CLOSEST_GENE %in% label | enhancer_signal$Rank %in% 1:top,], 
                                     aes(x = Rank/1000, y = Signal, label = CLOSEST_GENE, fill = Type), 
                                     fontface = "bold", size = 2.5, 
                                     box.padding = unit(1, "lines"), 
                                     segment.color = brewer.pal(n = 9, name ="Set1")[2], 
                                     point.padding = unit(0.3, "lines"), 
                                     segment.size = 0.3, show.legend = FALSE, 
                                     label.r = 0.5, nudge_x = 0, nudge_y = 0)  + 
    scale_fill_manual(values=c(brewer.pal(n = 9, name ="Set3")[4], brewer.pal(n = 9, name ="Set3")[9]))
  p <- p + theme(legend.title = element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"))
  p <- p + theme(panel.border = element_rect(fill=NA,linetype = 1,size = 1), axis.line = element_blank())
  return(p)
}

###12. Build function of sig gene expression by mut stat(WT1 TET2 IDH1 IDH2)
expression_boxplot_byMut <- function(gene, risk.df, nm.df) {
  #  gene <- "CHRNE"
  #  risk.df <- rs.df.clinc.beat
  #  nm.df <- aml_all_tpm.matrix.nm
  index1 <- match(gene,nm.df$Symbol)
  index2 <- match(gene,colnames(risk.df))
  index3 <- match(c("WT1.x","TET2.x", "IDH1.x", "IDH2.x"),colnames(risk.df))
  group<-rowSums(risk.df[,index3])
  group[group > 0] <-  "MUT"
  group[group == 0] <-  "WT"
  gg1 <- data.frame(Group = rep("NM",length(index1)), TPM=as.numeric(nm.df[index1,4:ncol(nm.df)][1,]), Risk = rep("NM",length(index1)))
  gg2 <- data.frame(Group = group, TPM=as.numeric(risk.df[,index2]), Risk = risk.df[,4])
  
  gg <- rbind(gg1,gg2)
  gg$Group <- factor(gg$Group, levels = c("MUT","WT","NM"))
  gg$Risk <- factor(gg$Risk, levels = c("High_risk","Low_risk","NM"))
  ggplot(gg,aes(x=Risk,y=log2(TPM),color=Risk,fill=Risk)) + 
    geom_boxplot(outlier.shape = 21,show.legend = FALSE,alpha = 0.65) + 
    #  geom_jitter(size=0.5,show.legend = FALSE) + 
    theme(panel.border = element_rect(fill=NA,linetype = 1),
          panel.background = element_blank()) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5), plot.title= element_text(hjust = 0.5)) + 
    scale_color_manual(values = c("#FF7F00","#1F78B4","grey50")) + 
    scale_fill_manual(values = c("#FF7F00","#1F78B4","grey50")) + 
#    stat_compare_means(show.legend = FALSE, comparisons = list(c("MUT","WT"),c("WT","NM"),c("MUT","NM"))) +
    stat_compare_means(show.legend = FALSE, comparisons = list(c("High_risk","Low_risk"),c("Low_risk","NM"),c("High_risk","NM"))) +
    ggtitle(paste0(gene, " Expression")) + xlab("")
}
