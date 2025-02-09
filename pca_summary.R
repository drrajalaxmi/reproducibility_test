
##______________________________
## Rajalaxmi 
# 240825 
# mm_data_240814
##-------------------

#library
library(readxl)
library(ggfortify) 
library(dplyr)
library(tidyr)
library(openxlsx)
library(DESeq2)
library(ggplot2)
library(ggfortify)
library(apeglm)
library(biomaRt)
library(clusterProfiler)
library(stringr)
library(tibble)
library(fgsea)
#BiocManager::install("biomaRt")




BE02_case_list <-read_excel("data/BE02 case list.xlsx", sheet = "FINAL")

column_name <- stringr::str_replace(colnames(BE02_case_list), " ", "_")
column_name <- stringr::str_replace_all(column_name, "\\(|\\)", "")
colnames(BE02_case_list) <- column_name
BE02_case_list$Sample_ID <- paste0("BE",BE02_case_list$Sample_ID)
str(BE02_case_list)

 #df_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
df_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
colnames(df_count) <- stringr::str_replace(colnames(df_count), "X", "BE")
condition <- BE02_case_list %>% filter(!Sample_ID %in% c("BE65" , "BE77", "BE78",  "BE80",  "BE81",  "BE120")) %>%
                          #  filter(Grade !="normal") %>%
                       filter(Sample_ID %in% colnames(df_count))

length(condition$Sample_ID)
counts <- df_count  %>% select(condition$Sample_ID)
# colnames(counts) <- paste0(condition$Path_ID,"(",condition$Sample_ID,")",condition$Disease
#                            )
length(colnames(counts))
# counts <-read.csv("data/counts.csv", header = TRUE,stringsAsFactors=FALSE)
# rownames(counts) <- counts[[1]]
# counts <- counts[-1]
# head(counts[,1:5])



cnt_pca <- prcomp(t(counts))
#summary(cnt_pca)

plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (all Grades)") 
pca.plot <- autoplot(cnt_pca,
                     data = t(counts)) 
plotly::ggplotly(autoplot(cnt_pca, data = condition,color = "Path_ID"))
         #label.size = 3)

count_mtx <- as.matrix(counts)

Ca_cor <-round(cor(count_mtx), 2)
# heatmaply::heatmaply((Ca_cor))
# df_val <- df_mtx[!rowSums(df_mtx)==0, ]
# dim(df_val)

heatmap(Ca_cor, margins = c(5, 5))

# colnames(condition)
#autoplot(kmeans(t(counts), 5), data = t(counts))
# Inspired by https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/01c_RNAseq_count_distribution.html
mean_counts <- apply(counts, 1, mean)
variance_counts <- apply(counts, 1, var)
df <- data.frame(mean_counts, variance_counts)

ggplot(df) +
  geom_point(aes(x=mean_counts, y=variance_counts)) + 
  geom_line(aes(x=mean_counts, y=mean_counts, color="red")) +
  scale_y_log10() +
  scale_x_log10()

# write.csv(condition,"data/sample_mtx.csv")
unique(condition$Disease)
# plot_name <- "isolated HGD vs progressor HGD w PDAC"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated HGD", "Progressor HGD w PDAC")) %>%
#   mutate( Grade = factor(Grade, levels= c("Isolated HGD", "Progressor HGD w PDAC")))
str(sample_mtx)
#isolated LGD vs. progressor LGD all (combined HGD or PDAC)
# plot_name <- "isolated LGD vs. progressor LGD all (combined HGD or PDAC)"
# sample_mtx <- condition %>% filter(Disease %in% c("Isolated LGD", "Progressor LGD w HGD", "Progressor LGD w PDAC")) %>%
#   mutate( Disease = ifelse(Disease =="Isolated LGD", "Isolated LGD","Progressor LGD" ) ,
#     Disease = factor(Disease, levels= c("Isolated LGD", "Progressor LGD")))

# normal-LGD vs. (normal-HGD and normal-PDAC)

plot_name <- "normal-LGD vs. (normal-HGD and normal-PDAC)"
sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
  mutate( Disease = ifelse(Disease =="Normal-LGD", "Normal-LGD","Normal-Progressor" ) ,
          Disease = factor(Disease, levels= c("Normal-LGD","Normal-Progressor")))

plot_name <- "(normal-LGD and normal-HGD) vs. normal-PDAC"
sample_mtx <- condition %>% filter(Disease %in% c("Normal-LGD", "Normal-HGD", "Normal-PDAC")) %>%
  mutate( Disease = ifelse(Disease =="Normal-PDAC", "Normal-PDAC","Normal-Non-Progressor" ) ,
          Disease = factor(Disease, levels= c("Normal-Non-Progressor","Normal-PDAC")))
sample_mtx$Disease
# 
# plot_name <- "Normal Vs Low grade"
# sample_mtx <- condition %>% filter(Grade %in% c("normal", "low grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("normal", "low grade")))

# plot_name <- "Normal Vs Invasive"
# sample_mtx <- condition %>% filter(Grade %in% c("normal", "invasive")) %>%
#   mutate( Grade = factor(Grade, levels= c("normal", "invasive")))
#
# plot_name <- "Normal Vs High grade"
# sample_mtx <- condition %>% filter(Grade %in% c("normal", "high grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("normal", "high grade")))
# # 
# plot_name <- "Low grade Vs Invasive"
# sample_mtx <- condition %>% filter(Grade %in% c("low grade", "invasive")) %>%
#   mutate( Grade = factor(Grade, levels= c("low grade", "invasive")))
# 
# plot_name <- "Low grade Vs High grade"
# sample_mtx <- condition %>% filter(Grade %in% c("low grade", "high grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("low grade", "high grade")))
# 
# plot_name <- "Invasive Vs High grade"
# sample_mtx <- condition %>% filter(Grade %in% c("invasive", "high grade")) %>%
#   mutate( Grade = factor(Grade, levels= c("invasive", "high grade")))

# c("normal", "invasive", "low grade", "high grade")





size <- sample_mtx %>% mutate(Diseases = as.character(Disease)) %>% 
                      group_by(Disease) %>% summarise(Diseases = unique(Diseases) ,
                                                    n = length(Sample_ID),)
counts_mtx <- counts %>% dplyr::select(which(colnames(counts) %in% sample_mtx$Sample_ID)) %>% as.matrix()



#Import to DEseq2
counts.DEseq = DESeqDataSetFromMatrix(countData=round(counts_mtx), colData = sample_mtx, design = ~Disease)

dds <- DESeq(counts.DEseq)
comparison <- resultsNames(dds)[length(resultsNames(dds))] #lists the coefficients


# cnt_pca <- prcomp(t(counts_mtx))
# summary(cnt_pca)
# 
# plot.iris.pca <- plot(cnt_pca, type="l", main ="PCA (Normal Vs Low grade)") 
# pca.plot <- autoplot(cnt_pca,
#                      data = t(counts_mtx)) 
# autoplot(cnt_pca, data = sample_mtx,color = "Grade", shape = FALSE, label.size = 3)


# par(mar=c(4,4,1,1))
# plotDispEsts(dds)
# 
# rld <- rlog(dds, blind = TRUE)
# 
# #plotPCA from DEseq2 plots uses the top 500 genes:
# data = plotPCA(rld, intgroup = c("Grade"), returnData = TRUE)
# p <- ggplot(data, aes(x = PC1, y = PC2, color =  Grade))
# p <- p + geom_text(aes(label = name)) + 
#   labs(title = paste0("PCA ", plot_name))+ theme ()
# print(p)




res <- results(dds)



#Substitute the '????' with a comparison, selected from the resultsNames(dds) shown above
LFC <- lfcShrink(dds, coef = comparison, type = "apeglm")
# #Add a title to reflect your comparison 
# plotMA(LFC, main = plot_name, cex = 0.5, ylim=c(-2,2))
# #The distribution of p-values
# hist(LFC$pvalue, breaks = 50, col = 'grey', main = plot_name, xlab = 'p-value')
# 
# #The false-discovery rate distribution
# hist(LFC$padj, breaks = 50, col = 'grey', main = plot_name, xlab = 'Adjusted p-value')
# 
# #Allow for more space around the borders of the plot
# par(mar = c(5, 4, 4, 4))

#Set your log-fold-change and p-value thresholds
lfc = 2
pval = 0.05

tab_res = data.frame(logFC = res$log2FoldChange, negLogPval = -log10(res$padj))#make a data frame with the log2 fold-changes and adjusted p-values


 
tab = data.frame(logFC = LFC$log2FoldChange, negLogPval = -log10(LFC$padj))#make a data frame with the log2 fold-changes and adjusted p-values



#Genes with a fold-change greater than 2 and p-value<0.05:
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))
signGenes_res = (abs(tab_res$logFC) > lfc & tab_res$negLogPval > -log10(pval))

gene_res <- rownames(res[which(signGenes_res),])
gene <- rownames(LFC[which(signGenes),])


# gene_list <- unique(c(gene,  gene_list))



# logFC <- LFC %>% as.data.frame() %>% select(log2FoldChange) %>% rownames_to_column("enseble_ID") %>% 
#           rename(nVsh_logFC = log2FoldChange) %>% #head()
#           left_join(logFC)
# 
# 
# 
# 
# logFC_signif <- logFC %>% filter(enseble_ID %in% gene_list)
# head(logFC_signif)
# 
# write.csv(logFC_signif,  "data/signif_gene_log2FC.csv")
# write.csv(data.frame("geneID" = c(gene_list[,2])),  "data/signif_genelist_lcf.csv")


# 
# length(gene_list)
# length(gene)
# sum(!gene %in% gene_res)

# plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nSignificant genes from lcfshrink")) #replace main = with your title

# gene_ensemble <- rownames(LFC)[which(LFC$log2FoldChange < -5)]
# gene_ensemble = unlist(lapply(strsplit(gene_ensemble, '\\.'), '[[',1))
# 
# converted <- getBM(attributes=c('entrezgene_id','hgnc_symbol'), filters = 'hgnc_symbol',
#                    values = gene_ensemble, mart = ensembl)
#Colour these red
points(tab[signGenes, ], pch = 16, cex = 0.5, col = "red")

plot(tab_res, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nSignificant genes from DEseq")) #replace main = with your title
points(tab_res[signGenes_res, ], pch = 16, cex = 0.5, col = "red")



#Show the cut-off lines
abline(h = -log10(pval), col = "green3", lty = 2)
abline(v = c(-lfc, lfc), col = "blue", lty = 2)

mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.6, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.6, line = 0.5)

# plot(tab_res, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nColored for significant gene from lcfshrink")) #replace main = with your title
# points(tab_res[signGenes, ], pch = 16, cex = 0.5, col = "red")


# plot(tab, pch = 16, cex = 0.4, xlab = expression(log[2]~fold~change),
#      ylab = expression(-log[10]~pvalue), main = paste0(plot_name, "\nColored for significant gene from DEseq")) #replace main = with your title
# points(tab[signGenes_res, ], pch = 16, cex = 0.5, col = "red")




# #increased expression
# attach(as.data.frame(res))
# 
# #The total number of DEGs with an adjusted p-value<0.05
# summary(res, alpha=0.05)
# 
# 
# #The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
# a <-sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
# 
# #Decreased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange < 0) #any fold-change
# b <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange > (2)) #fold-change greater than 2
# 
# 
# #Increased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
# c <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange < -2) #fold-change greater than 2
# 
# cat("Sample Size: \n",paste0(size[1,]), "\n",paste0(size[2,]), "\n \nTotal DE genes: ", a, "\nUp regulated: ", b, "\nDown regulated: ", c)

# #increased expression
# attach(as.data.frame(LFC))
# 
# #The total number of DEGs with an adjusted p-value<0.05
# summary(LFC, alpha=0.05)
# 
# 
# #The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
# a <-sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
# 
# #Decreased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange < 0) #any fold-change
# b <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange > (2)) #fold-change greater than 2
# 
# 
# #Increased expression:
# sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
# c <-sum(!is.na(padj) & padj < 0.05 & log2FoldChange < -2) #fold-change greater than 2
# 
# cat("Sample Size: \n",paste0(size[1,]), "\n",paste0(size[2,]), "\n \nTotal DE genes: ", a, "\nUp regulated: ", b, "\nDown regulated: ", c)

attach(as.data.frame(LFC))
#At this stage it may be useful to create a copy of the results with the gene version removed from the gene name, to make it easier for you to search for the gene name etc. 
#The rownames currently appear as 'ENSG00000175197.12, ENSG00000128272.15' etc.
#To change them to 'ENSG00000175197, ENSG00000128272'
LFC.gene = as.data.frame(LFC)

#Some gene names are repeated if they are in the PAR region of the Y chromosome. Since dataframes cannot have duplicate row names, we will leave these gene names as they are and rename the rest.
whichgenes = which(!grepl('PAR', rownames(LFC.gene)))
rownames(LFC.gene)[whichgenes] = unlist(lapply(strsplit(rownames(LFC.gene)[whichgenes], '\\.'), '[[',1))

#subset the significant genes
LFC.sig = LFC.gene[padj < 0.05 & !is.na(padj) & abs(log2FoldChange) >2,]#subset the significant genes
# LFC.sig %>% arrange(log2FoldChange) %>% t
# rownames(LFC.sig)
#We can add a column with the HGNC gene names
# mart = useMart('ensembl')
# listDatasets(mart)
# listEnsembl()
# ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
# mmusculus_gene_ensembl  - mouse
# hsapiens_gene_ensembl - human
# 
# ensembl <- useEnsembl(biomart = 'genes',
#            dataset = 'hsapiens_gene_ensembl', host ="useast.ensembl.org",
#            version = 113)



# (which(str_detect(listAttributes(ensembl)$name, "symbol")))
# listAttributes(ensembl)$name[98]

converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = rownames(LFC.sig), mart = ensembl)

# converted <- converted %>% filter(ensembl_gene_id %in% gene_corr_symbol$ensembl_gene_id)

#Add gene names to the LFC.sig data-frame

LFC.sig_gene <- LFC.sig %>% rownames_to_column(var ="ensembl_gene_id") %>% left_join(converted) %>%
           filter(!is.na(entrezgene_id))
#LFC.sig$hgnc = converted[converted[,1] %in% rownames(LFC.sig),2]

#View the top 10 genes with the most significant (adjusted) p-values
# head(LFC.sig, n = 10)

#The largest fold-changes with a significant p-value
# LFC.sig[order(abs(LFC.sig$log2FoldChange), decreasing = TRUE),][1:10,] #add the [1:10,] to see the top 10 rows



#org.Mm.eg.db - mouse
#org.Hs.eg.db - Human
mf_go <- enrichGO(LFC.sig_gene$entrezgene_id, 'org.Hs.eg.db', ont= "MF", pvalueCutoff=0.01) # %>%
  # as.data.frame() 
dotplot(mf_go, title = paste0(plot_name, ": GO-MF"))


bp_go <- enrichGO(LFC.sig_gene$entrezgene_id, 'org.Hs.eg.db', ont= "BP", pvalueCutoff=0.01)  #%>%
  # as.data.frame() 

dotplot(bp_go, title = paste0(plot_name, ": GO-BP"))

cc_go <- enrichGO(LFC.sig_gene$entrezgene_id, 'org.Hs.eg.db', ont= "CC", pvalueCutoff=0.01) # %>%
  # as.data.frame() 
dotplot(cc_go, title = paste0(plot_name, ": GO-CC"))
# emapplot(enrichplot::pairwise_termsim(bp_go))
kegg <- enrichKEGG(LFC.sig_gene$entrezgene_id, 'hsa', pvalueCutoff=0.01)   # for human
# kegg <- enrichKEGG(LFC.sig_gene$entrezgene_id, 'mmu', pvalueCutoff=0.01)   # for mouse
dotplot(kegg,  title = paste0(plot_name, ": KEGG"))


#######

res$row <- rownames(res)
dim(res)
# res.sig.gene <- res[which(abs(res$log2FoldChange) >2 & res$padj >0.05),]
# dim(res.sig.gene)
# ens2symbol <- AnnotationDbi::select(org.Mm.eg.db,
#                                     key=res$row, 
#                                     columns="SYMBOL",
#                                     keytype="ENSEMBL")
# names(ens2symbol)[1] <- "row"

# dim(LFC.sig_gene)


rownames(res) = unlist(lapply(strsplit(rownames(res), '\\.'), '[[',1))

converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = rownames(res), mart = ensembl)

res2 <- res %>% as.data.frame() %>% 
  rownames_to_column(var ="ensembl_gene_id") %>% left_join(converted) %>%
  filter(!is.na(entrezgene_id), !is.na(stat)) %>%
  group_by(entrezgene_id) %>% 
  summarize(stat=mean(stat)) 

nrow(res2)
ranks <- res2$stat
names(ranks) <- res2$entrezgene_id



### downloaded from https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
# pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/mm_data_240814/mm_data_240814/data/mh.all.v2024.1.Mm.entrez.gmt")   # for mouse
# head(pathways.hallmark)
pathways.hallmark <- gmtPathways("~/Documents/Rutgers_CW/De_lab/download_resources/h.all.v2024.1.Hs.entrez.gmt") # for human

#Running fgsea algorithm:
fgseaRes <- fgseaMultilevel(pathways=pathways.hallmark, stats=ranks)

# Tidy the results:
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES)) # order by normalized enrichment score (NES)

# To see what genes are in each of these pathways:
gene.in.pathway <- pathways.hallmark %>% 
  enframe("pathway", "entrezgene_id") %>% 
  unnest(cols = c(entrezgene_id)) %>% 
  mutate(entrezgene_id = as.numeric(entrezgene_id)) %>%
  inner_join( res2, by="entrezgene_id")

fgseaResTidy$adjPvalue <- ifelse(fgseaResTidy$padj <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy %>% mutate(pathway = str_remove(pathway, "HALLMARK_")), aes(reorder(pathway, NES), NES, fill = adjPvalue)) +
  geom_col(show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Hallmark Pathways", y="Normalized Enrichment Score",
       title=paste0(plot_name) )
# ggsave(paste0("docs/",str_remove(comparison, "Grade_"),"_hallmark.png"), width = 4.5, height = 7)
# dir()
# 
# plotGseaTable(pathways.hallmark[fgseaRes$pathway[fgseaRes$padj < 0.05]], ranks, fgseaRes, 
#               gseaParam=0.5)
# 
# # "HALLMARK_WNT_BETA_CATENIN_SIGNALING"
# 
# 
# plotEnrichment(pathways.hallmark[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]],
#                ranks) + labs(title="wnt signaling")
# 
# 
# 
# ### for Dr Shen
# rownames(counts) = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))
# 
# converted <- getBM(attributes=c('ensembl_gene_id','mgi_symbol' ), filters = 'ensembl_gene_id',
#                    values = rownames(counts), mart = ensembl)
# 
# converted <- converted %>% rename( "gene_symbol" = "mgi_symbol")
# 
# counts_name <- converted %>% left_join( counts %>% rownames_to_column("ensembl_gene_id") ) 
# write.xlsx(counts_name, "data/gene_counts.xlsx")                          
# 
# 
