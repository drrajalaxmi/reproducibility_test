## ssgsea

# BiocManager::install("GSVA")
library(GSVA)
library(dplyr)

ensembl <- useEnsembl(biomart = 'genes',
           dataset = 'hsapiens_gene_ensembl', host ="useast.ensembl.org",
           version = 113)
rownames(counts) = unlist(lapply(strsplit(rownames(counts), '\\.'), '[[',1))




converted <- getBM(attributes=c('entrezgene_id','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = rownames(counts), mart = ensembl)

df_count <-read.csv("data/tpm_untrimmed.csv", header=TRUE, row.names = 1)
#df_count <-read.csv("data/counts_salmon_bam_untrimmed.csv", header=TRUE, row.names = 1)
colnames(df_count) <- stringr::str_replace(colnames(df_count), "X", "BE")
condition <- BE02_case_list %>% filter(!Sample_ID %in% c("BE65" , "BE77", "BE78",  "BE80",  "BE81",  "BE120")) %>%
  #  filter(Grade !="normal") %>%
  filter(Sample_ID %in% colnames(df_count)) %>% arrange( Disease, Category)

length(condition$Sample_ID)
counts <- df_count  %>% select(condition$Sample_ID)
 colnames(counts) <- paste0(condition$Path_ID,"(",condition$Sample_ID,")",condition$Disease
                            )


counts_entreID <- counts %>% 
  rownames_to_column(var ="ensembl_gene_id") %>% 
  left_join(converted) %>%
  filter(!is.na(entrezgene_id)) 
counts_entreID <- counts_entreID[!(duplicated(counts_entreID$entrezgene_id)),] 
rownames(counts_entreID) <- counts_entreID$entrezgene_id
dim(counts_entreID)

counts_entreID <- counts_entreID %>% select(-c("ensembl_gene_id", "entrezgene_id"))
str(pathways.hallmark)
gsvapar <- ssgseaParam(as.matrix(counts_entreID), pathways.hallmark)

ssgsea_kegg <- gsva(gsvapar )

heatmap((ssgsea_kegg))
        
heatmaply::heatmaply(t(ssgsea_kNULLheatmaply::heatmaply(t(ssgsea_kegg))
