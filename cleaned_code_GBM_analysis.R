

#load all the libraries
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readxl)
library(WGCNA)
library(tximport)
library(AnnotationDbi)
library(GO.db)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v75)
library(ggfortify)
library(factoextra)
library(cluster)
library(ggpubr)
library(DESeq2)
library(sva)
library(limma)
library(plotrix)
library(GenomicFeatures)
library(DRIMSeq)
library(DEXSeq)
library(stageR)
library(fgsea)
library(msigdbr)
library(dendextend)
library(qusage)
library(GSVA)
library(maftools)
library(preprocessCore)
library(EBSeq)
library(ggsci)
library(GSEABase)
library(GSVAdata)
library(Biobase)
library(genefilter)
library(RColorBrewer)
library(RankProd)
library(ppcor)
library(edgeR)
library(survival)
library(survminer)
library(biomaRt)
# library(rats)
library(factoextra)
library(NbClust)
library(umap)
library(pals)
library(nph)
library(ggbiplot)

library(reshape2)

clinical_data = read.table('clin_data.txt',sep='\t')
gene_abundances = read.table('output_tpm_all_samples.txt',sep='\t')

colnames(gene_abundances) <- gsub('X','',colnames(gene_abundances))
#these are the batches
names1 <- as.character(melt(outer(c('1','2','3','4','5','7','8','9','10'),c('A','B'),FUN='paste0'))[,3])
names2 <- as.character(melt(outer(c('31','32','33','34','35','36','37','38','39','40','41','42','43','44','45','46','47','48','49','50'),c('A','B'),FUN='paste0'))[,3])
names3 <- c(as.character(melt(outer(c('11','12','13','14','17','19','24','25','26','27','27','30'),c('A','B'),FUN='paste0'))[,3]), c('20A','21A','22A','23A'))
batches_vec <- c(rep(1,times=length(names1)),rep(2,times=length(names2)),rep(3,times=length(names3)))
names(batches_vec) <- c(names1,names2,names3)


#perform adjustment of the data 
gene_abundances_filtered <- gene_abundances
gene_abundances_filtered <- gene_abundances_filtered[which(rowSums(gene_abundances_filtered > 0)  > 1),] #at least 2 samples with at leat 1 gene nonzero
gene_abundances_filtered_adjusted <- removeBatchEffect(gene_abundances_filtered, batch=batches_vec[colnames(gene_abundances_filtered)])

#molecular subtype


library('biomaRt')
mart <- useEnsembl('ensembl', mirror = 'useast')
mart <- useDataset("hsapiens_gene_ensembl", mart=mart)

verhaak_sigs <- list()
for(fName in list.files('verhaak_sigs/')){
	if(fName!='neural.txt'){
		verhaak_sigs[[substr(fName,start=1,stop=nchar(fName)-4)]] <- read.table(file=paste0('verhaak_sigs/',fName),sep='\n',stringsAsFactors=F,header=T,quote='')
		verhaak_sigs[[substr(fName,start=1,stop=nchar(fName)-4)]]  <- getBM(filters= "hgnc_symbol", attributes= c("ensembl_gene_id","hgnc_symbol"),values=verhaak_sigs[[substr(fName,start=1,stop=nchar(fName)-4)]] ,mart= mart)
	}
}

library(GSVA)
output_ssgsea <- gsva(gene_abundances_filtered_adjusted,list(verhaak_sigs[[1]][,1],verhaak_sigs[[2]][,1],verhaak_sigs[[3]][,1] ),method='ssgsea') #,verhaak_sigs[[4]][,1]
heatmap(output_ssgsea)
molec_subtype <- names(verhaak_sigs)[apply(output_ssgsea,2,function(z){which.max(z)})]
names(molec_subtype) <- colnames(output_ssgsea)

#-------survival curves----------
#need to make overall survival curve
library(survival)
library(survminer)
library(ggplot2)
library(survival)
library(survminer)
library(ggplot2)

egfr_vec <- c(0,1,-1)
names(egfr_vec) <- c('Not Amplified','Amplified','Uninformative')

data_to_plot <- (clinical_data[prim_tums,])
data_to_plot$OS <- as.numeric(as.character(data_to_plot$OS))
data_to_plot$Expired <- as.numeric(as.character(data_to_plot$Expired))
data_to_plot$EGFR <- as.character(clinical_data[prim_tums,'EGFR'])
data_to_plot$EGFR[which(data_to_plot$EGFR == 'Uninformative')] = NA
data_to_plot$EGFR[which(data_to_plot$EGFR == 'Amplified')] = 'Amp'
data_to_plot$EGFR[which(data_to_plot$EGFR == 'Not Amplified')] = 'Non Amp'

fit <- survfit(Surv(OS,Expired) ~ EGFR, data = data_to_plot)

dev.new()
ggsurvplot(fit, conf.int = F,  pval=T  , risk.table = TRUE,xlim = c(0, 2500),  break.time.by = 500,#legend.labs=c("Sex=1", "Sex=2"),
          palette='npj',data=data_to_plot,) + labs(title = "Overall Survival, Primary EGFR Status")# ggtheme = theme_minimal()
# dev.copy(pdf,'overall_survival_curves_EGFR_primary.pdf',height=5,width=6)
# dev.off()
#--------------------------------

# -------check correlations between primary and recurrent samples
all_cors_gene_expr_primary_rec <- c()
all_cors_gene_expr_primary_rec_p <- c()
all_gene_expr <- c()
for(k in prim_tums){
	tum_name = gsub('A','',k)
	tmp_val = cor.test(as.numeric(gene_abundances_filtered_adjusted[,paste0(tum_name,'A')]),as.numeric(gene_abundances_filtered_adjusted[,paste0(tum_name,'B')]),method='spearman')
	all_cors_gene_expr_primary_rec <- c(all_cors_gene_expr_primary_rec,tmp_val$estimate)
	all_cors_gene_expr_primary_rec_p <- c(all_cors_gene_expr_primary_rec_p,tmp_val$p.value)
	all_gene_expr <- rbind(all_gene_expr,cbind(prim=gene_abundances_filtered_adjusted[,paste0(tum_name,'A')],rec=gene_abundances_filtered_adjusted[,paste0(tum_name,'B')]))
}
print(min(all_cors_gene_expr_primary_rec))
print(max(all_cors_gene_expr_primary_rec))
print(min(all_cors_gene_expr_primary_rec_p))
print(max(all_cors_gene_expr_primary_rec_p))

dev.new()
plot(density(all_cors_gene_expr_primary_rec),xlab='',main='',ylab='',frame.plot=F)
abline(v=0.7,lty=2)
dev.copy(pdf,file='cor_coeffs_samples.pdf',width=4,height=3.5)
dev.off()

names(all_cors_gene_expr_primary_rec) <- prim_tums

low_cor = prim_tums[which(all_cors_gene_expr_primary_rec <= 0.7627424)]
high_cor = prim_tums[which(all_cors_gene_expr_primary_rec >= 0.8561980)]

#test age
wilcox.test(as.numeric(as.character(clinical_data[low_cor,'Age'])), as.numeric(as.character(clinical_data[high_cor,'Age'])))

#test OS
wilcox.test(as.numeric(as.character(clinical_data[low_cor,'OS'])), as.numeric(as.character(clinical_data[high_cor,'OS'])))
cor.test(all_cors_gene_expr_primary_rec,as.numeric(as.character(clinical_data[prim_tums,'OS'])),method='spearman')

#test Surgery1ToSurgery2
wilcox.test(as.numeric(as.character(clinical_data[low_cor,'Surgery1ToSurgery2'])), as.numeric(as.character(clinical_data[high_cor,'Surgery1ToSurgery2'])))

#test EGFR
#fisher's test
egfr_amp = prim_tums[which(clinical_data[prim_tums,'EGFR'] == 'Amplified')]
wilcox.test(all_cors_gene_expr_primary_rec[egfr_amp], all_cors_gene_expr_primary_rec[setdiff(prim_tums,egfr_amp)])

tab_11 = length(intersect(low_cor,egfr_amp))
tab_12 = length(low_cor) - tab_11
tab_21 = length(egfr_amp) - tab_11
tab_22 = length(prim_tums) - length(low_cor) - length(egfr_amp) + tab_11

fisher.test(cbind(c(tab_11,tab_12),c(tab_21,tab_22)))

#test sex
male = prim_tums[which(clinical_data[prim_tums,'Sex'] == 'M')]

tab_11 = length(intersect(low_cor,male))
tab_12 = length(low_cor) - tab_11
tab_21 = length(male) - tab_11
tab_22 = length(prim_tums) - length(low_cor) - length(male) + tab_11

fisher.test(cbind(c(tab_11,tab_12),c(tab_21,tab_22)))

#test Max Diam
wilcox.test(as.numeric(as.character(clinical_data[low_cor,'Max Diam'])), as.numeric(as.character(clinical_data[high_cor,'Max Diam'])))
cor.test(all_cors_gene_expr_primary_rec,as.numeric(as.character(clinical_data[prim_tums,'Max Diam'])),method='spearman')

tmp_val = cor.test(all_gene_expr[,'prim'],all_gene_expr[,'rec'],method='spearman')

for(k in colnames(clinical_data)){
	if(sum(is.na(as.numeric(clinical_data[prim_tums,k]))) == 0){


		tmp = cor.test(all_cors_gene_expr_primary_rec,as.numeric(clinical_data[prim_tums,k]),method='spearman')
		if(!is.na(tmp$p.value)){
			# if(tmp$p.value < 0.1){
				print(k)
			print(tmp$p.value)

			# }
		}
	}

}

wilcox.test(all_cors_gene_expr_primary_rec[which(clinical_data[prim_tums,'EGFR'] =='Amplified')],all_cors_gene_expr_primary_rec[which(clinical_data[prim_tums,'EGFR'] !='Amplified')])


#---------------differential expression analysis

test_differential_expression <- function(groupA,groupB,paired,EGFR=NULL){

	names_to_take <- paste0(tolower(c(groupA,groupB)),'_all') 

	txi.subset <- #note this needs to be from raw quants files from the analysis of raw data tximport(paste0('~/quants/',names_to_take,'/quant.sf'),type='salmon',tx2gene=tx2gene)

	if(paired){
		sampleTable <- data.frame(patient=gsub('A|B','',c(groupA,groupB)),condition=c(rep('A',times=length(groupA)),rep('B',times=length(groupB))))
	}else{
		sampleTable <- data.frame(condition=c(rep('A',times=length(groupA)),rep('B',times=length(groupB))))
	}
	row.names(sampleTable) <- c(groupA,groupB)

	if(!is.null(EGFR)){
		sampleTable$EGFR = as.factor(rep(c('Not Amplified','Amplified')[1*(EGFR=='Amplified') + 1],times=2))
		dds <- DESeqDataSetFromTximport(txi.subset, sampleTable, ~patient + condition + EGFR )

	}else{
		if(paired){
			dds <- DESeqDataSetFromTximport(txi.subset, sampleTable, ~patient + condition )
		
		}else{
			dds <- DESeqDataSetFromTximport(txi.subset, sampleTable, ~ condition )	
		}
	}
	dds <- dds[ rowSums(counts(dds)) > 1, ]
	dds <- DESeq(dds, betaPrior=FALSE)
	res <- DESeq2::results(object = dds,contrast = c('condition','A','B'))
	res_processed <- res[which(!is.na(res$padj)),]
	res_processed
	# names_de_genes <- rownames(res_processed)[which(res_processed$padj < 0.01)]
	# names_de_genes <- names_de_genes[order(res_processed[names_de_genes,'padj'])]
	# getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=names_de_genes ,mart= mart)
}


plot_diff_exp_genes <- function(res_processed,groupA,groupB,groupA_label, groupB_label,fName){
	names_de_genes <- rownames(res_processed)[which(res_processed$padj < 0.01)]
	names_de_genes <- names_de_genes[order(res_processed[names_de_genes,'padj'])]
	names_genes_hgnc <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=names_de_genes ,mart= mart)
	names_genes_hgnc = names_genes_hgnc[which(nchar(names_genes_hgnc[,2])> 1),]
	row.names(names_genes_hgnc) <- names_genes_hgnc[,1]
	names_genes_hgnc <- names_genes_hgnc[intersect(rownames(gene_abundances_filtered_adjusted),names_genes_hgnc[,1]),]
	melted_df = cbind(gene_abundances_filtered_adjusted[names_genes_hgnc[,1],groupA],gene_abundances_filtered_adjusted[names_genes_hgnc[,1],groupB])
	melted_df = melt(melted_df)
	melted_df = cbind(melted_df, group=c(rep(groupA_label,times=(length(groupA) * length(names_genes_hgnc[,1]) )),rep(groupB_label,times=(length(groupB) * length(names_genes_hgnc[,1]) ))))
	melted_df[,'Var1'] = names_genes_hgnc[melted_df[,'Var1'],2]
	dev.new()
	p = ggviolin(melted_df, x='group',y='value',xlab='',ylab='Gene expression (log-scale)',add = "boxplot",color='group' ,ggtheme=theme_pubr()) +  stat_compare_means(comparisons = list(c(1,2)), label = 'p.signif') + 
  		facet_grid(~Var1)+ rremove('legend')
	print(p)
			width_val = 2 * length(names_genes_hgnc[,1])

	dev.copy(pdf,file=paste0(fName,'.pdf'),height=5,width=width_val)
	dev.off()
	graphics.off()
}


groupA <- prim_tums
groupB <- rec_tums
paired=T
diff_exp_prim_vs_rec_all <- test_differential_expression(groupA,groupB,paired)


groupA <- (prim_tums[which(clinical_data[prim_tums,'EGFR']=='Amplified')])
groupB <- gsub('A','B',prim_tums[which(clinical_data[prim_tums,'EGFR']=='Amplified')])
paired=T
diff_exp_prim_vs_rec_egfr_amp <- test_differential_expression(groupA,groupB,paired)
diff_pathways_prim_vs_rec_egfr_amp <- test_differential_pathways(groupA,groupB,paired)
# res_processed = diff_exp_prim_vs_rec_egfr_amp
# names_de_genes <- rownames(res_processed)[which(res_processed$padj < 0.01)]
# names_de_genes <- names_de_genes[order(res_processed[names_de_genes,'padj'])]
# t <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=names_de_genes ,mart= mart)
# t$p_val <- res_processed[t[,1],'padj']
# t$fc <- res_processed[t[,1],'log2FoldChange']


groupA <- (prim_tums[which(clinical_data[prim_tums,'EGFR']=='Not Amplified')])
groupB <- gsub('A','B',prim_tums[which(clinical_data[prim_tums,'EGFR']=='Not Amplified')])
paired=T
diff_exp_prim_vs_rec_egfr_non_amp <- test_differential_expression(groupA,groupB,paired)
diff_pathways_prim_vs_rec_egfr_non_amp <- test_differential_pathways(groupA,groupB,paired)


groupA <- (prim_tums[which(clinical_data[prim_tums,'EGFR']=='Amplified')])
groupB <- (prim_tums[which(clinical_data[prim_tums,'EGFR']=='Not Amplified')])
paired=F
diff_exp_prim_egfr_amp_vs_non_amp <- test_differential_expression(groupA,groupB,paired)
diff_pathways_prim_egfr_amp_vs_non_amp <- test_differential_pathways(groupA,groupB,paired)


groupA <- (rec_tums[which(clinical_data[rec_tums,'EGFR']=='Amplified')])
groupB <- (rec_tums[which(clinical_data[rec_tums,'EGFR']=='Not Amplified')])
paired=F
diff_exp_rec_egfr_amp_vs_non_amp <- test_differential_expression(groupA,groupB,paired)
diff_pathways_rec_egfr_amp_vs_non_amp <- test_differential_pathways(groupA,groupB,paired)

###=====making the plot of EGFR amplification vs expression status in primary and recurrent tumours ====

prim_amplified <- prim_tums[which(clinical_data[prim_tums,'EGFR']=='Amplified')]
prim_not_amp <- prim_tums[which(clinical_data[prim_tums,'EGFR']!='Amplified')]
rec_amplified <- rec_tums[which(clinical_data[rec_tums,'EGFR']=='Amplified')]
rec_not_amp <- rec_tums[which(clinical_data[rec_tums,'EGFR']!='Amplified')]
egfr_ensg = 'ENSG00000146648'
v1 <- gene_abundances_filtered_adjusted[egfr_ensg,prim_amplified]
v2 <- gene_abundances_filtered_adjusted[egfr_ensg,prim_not_amp]
v3 <- gene_abundances_filtered_adjusted[egfr_ensg,rec_amplified]
v4 <- gene_abundances_filtered_adjusted[egfr_ensg,rec_not_amp]

df <- cbind.data.frame(expr = as.numeric(c(v1,v2,v3,v4)),group = as.factor(c(rep('Prim_Amp',times=length(v1)), rep('Prim_Non_Amp',times=length(v2)), rep('Rec_Amp',times=length(v3)), rep('Rec_Non_Amp',times=length(v4)))))
row.names(df) <- c(names(v1),names(v2),names(v3),names(v4))
df$amp <- clinical_data[rownames(df),'EGFR']
df <- df[-which(df$amp == 'Uninformative'),]
dev.new()
ggboxplot(df,x='group',y='expr',color='amp',palette='npj',add='jitter',
	add.params=list(alpha=1),xlab='',ylab='Expression',
	alpha=0.2,title='EGFR amplification vs. expression') + 
	rotate_x_text(angle=30)+   stat_compare_means(comparisons=list(c(1,2),c(3,4))) 
dev.copy(pdf,file='egfr_expr_vs_amp_status.pdf',width=4,height=5)
dev.off()
###==========================================================================================
pathway <- 'REACTOME_SIGNALING_BY_EGFR_IN_CANCER'

v1 <- output_ssgsea[pathway,prim_amplified]
v2 <- output_ssgsea[pathway,prim_not_amp]
v3 <- output_ssgsea[pathway,rec_amplified]
v4 <- output_ssgsea[pathway,rec_not_amp]

df <- cbind.data.frame(expr = as.numeric(c(v1,v2,v3,v4)),group = as.factor(c(rep('Prim_Amp',times=length(v1)), rep('Prim_Non_Amp',times=length(v2)), rep('Rec_Amp',times=length(v3)), rep('Rec_Non_Amp',times=length(v4)))))
row.names(df) <- c(names(v1),names(v2),names(v3),names(v4))
df$amp <- clinical_data[rownames(df),'EGFR']
df <- df[-which(df$amp == 'Uninformative'),]
dev.new()
ggboxplot(df,x='group',y='expr',color='amp',palette='npj',add='jitter',
	add.params=list(alpha=1),xlab='',ylab='ssGSEA score',
	alpha=0.2,title=paste0('EGFR amplification vs. pathway score')) + 
	rotate_x_text(angle=30)+   stat_compare_means(comparisons=list(c(1,2),c(3,4))) 
dev.copy(pdf,file='egfr_expr_vs_pathway_status.pdf',width=4,height=5)
dev.off()


#-------make the volcano plots----

mat = diff_exp_prim_vs_rec_egfr_amp
fName = paste0('diff_exp_prim_vs_rec_egfr_amp','_volcano.pdf')

g_list_all <- getBM(filters= "ensembl_gene_id", 
		attributes= c("ensembl_gene_id","hgnc_symbol"),
		values=rownames(mat),mart= mart)
	row.names(g_list_all) = g_list_all[,1]
	mat$HGNC = g_list_all[rownames(mat),2]
	mat$HGNC[which(is.na(mat$HGNC))] <- ''


x_max = ceiling(max(abs(mat$log2FoldChange))/5) * 5
y_max = ceiling(max(abs(-log10(mat$pvalue)))/5) * 5

bad_rows = which(mat$padj >= 0.05)
mat$HGNC[bad_rows] = ''
dev.new()
	EnhancedVolcano(mat,
                lab = mat$HGNC,
                 x = 'log2FoldChange',
                 y = 'pvalue',
                 pCutoff=0.05,    
                 pCutoffCol='padj',
                 xlim=c(-x_max,x_max),ylim=c(0,y_max),drawConnectors = T,gridlines.major=F,gridlines.minor=F,
                 labSize = 5,max.overlaps =Inf,selectLab=c('BRD2','ZIC2','TFF1','PRR3','LST1'),boxedLabels=T)
	dev.copy(pdf,file=fName,width=6.5,height=6)
	dev.off()
#--------

