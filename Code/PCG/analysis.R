rm(list=ls())
gc()

load("~/Desktop/PCG_SV/ACS_csDM_TOAST_fdr.Rdata")
res_acs_sig <- load("/Users/bingbing/Desktop/PCG_SV/ACS_csDM_TOAST_fdr.Rdata")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # EPIC array

annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# combine and annotate significanr CpG sites 
significant_cpgs <- res_acs_sig$CpGs
annotated_cpgs <- annotation[rownames(annotation) %in% significant_cpgs, ]

# absorb gene info
genes <- annotated_cpgs$UCSC_RefGene_Name
genes <- unlist(strsplit(genes, ";"))  # some CpG many relate to more than 1 gene
genes <- genes[genes != ""]
genes <- table(genes)  # count every gene related CpG sites
top_genes <- sort(genes, decreasing = TRUE)[1:100]  # select top 100 significant genes by most significant CpGs


### Plot - Top 100 Genes Enriched by Significant CpGs
library(ggplot2)

top_genes_df <- as.data.frame(top_genes)
colnames(top_genes_df) <- c("Gene", "Count")

ggplot(top_genes_df, aes(x = Count, y = reorder(Gene, -Count))) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(title = "Top 100 Genes Enriched by Significant CpGs",
       x = "Number of Associated CpGs",
       y = "Gene") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))



### GO result
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(org.Hs.eg.db)

# transfer Gene name into ENTREZ ID（路径分析需要）
gene_list <- names(top_genes)
head(top_genes_df)
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(gene_entrez)

# GO analysis
go_enrichment <- enrichGO(gene_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05)
View(go_enrichment)

# Save GO.csv - all GSEA infomation
write.csv(go_enrichment, file = "~/Desktop/PCG_SV/GO.csv", row.names = FALSE)
# Save top 100 gene name and gene number count
write.csv(top_genes_df, file = "~/Desktop/PCG_SV/top_100_genes.csv", row.names = FALSE)

### Plot - GO Enrichment Analysis for Top 100 Genes
dotplot(go_enrichment, showCategory = 20) +
  labs(title = "GO Enrichment Analysis for Top 100 Genes")



library(dplyr)
library(tidyr)
library(org.Hs.eg.db)

go_genes <- read.csv("~/Desktop/PCG_SV/GO.csv")
head(go_genes)

# View gene_entrez data
gene_list <- names(top_genes)
gene_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(gene_entrez)


# Match geneID in GO.csv with ENTREZID in gene_entrez
go_genes_matched <- go_genes %>%
  separate_rows(geneID, sep = "/") %>%  # 分割含有多个基因ID的行
  left_join(gene_entrez, by = c("geneID" = "ENTREZID")) %>%  # 匹配geneID和ENTREZID
  group_by(Description) %>%  # 使用Description分组
  summarize(
    geneID = paste(unique(geneID), collapse = ", "),  # 合并geneID
    SYMBOL = paste(unique(SYMBOL), collapse = ", ")  # 合并SYMBOL
  )

head(go_genes_matched)


# Merge the data, match up the Variable Description of each of the two datasets with GO_Term
merged_data1 <- merge(go_genes, go_genes_matched, by.x = "Description", by.y = "Description", all.x = TRUE)
head(merged_data1)

grouped_gene <- merged_data1 %>%
  group_by(Description) %>%
  arrange(desc(Count))  # rank by gene number

write.csv(grouped_gene,file="~/Desktop/PCG_SV/GO_with_genes.csv")

# rename gene.ID.x to geneID
colnames(grouped_gene)[colnames(grouped_gene) == "geneID.x"] <- "geneID"

# delet geneID.y
grouped_gene <- grouped_gene %>%
  select(-geneID.y)
head(grouped_gene)

write.csv(grouped_gene,file="~/Desktop/PCG_SV/GO_with_genes.csv")


