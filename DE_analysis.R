# Differential expression analysis with DESeq2

#Load libraries
library(DESeq2)

setwd("/home/d/chick/gallus_de/tm/")

directory <- "/home/d/chick/gallus_de/tm/"

# Generate sample level metadata (i.e SampleTable)

#samplesFiles is a variable which points to your htseq-count output files

sampleFiles <- c("KN_tm_1.sam.sorted.bam_gene.count.deseqfile",
                 "KN_tm_2.sam.sorted.bam_gene.count.deseqfile",
                 "KN_tm_3.sam.sorted.bam_gene.count.deseqfile",
                 "Bro_tm_1.sam.sorted.bam_gene.count.deseqfile",
                 "Bro_tm_2.sam.sorted.bam_gene.count.deseqfile",
                 "Bro_tm_3.sam.sorted.bam_gene.count.deseqfile"
)

sampleNames <- c("KN_tm_1", "KN_tm_2", "KN_tm_3", "Bro_tm_1", "Bro_tm_2", "Bro_tm_3")

sampleCondition <- c("KN_tm", "KN_tm", "KN_tm", "Bro_tm", "Bro_tm", "Bro_tm")

sampleTable <- data.frame(sampleName = sampleNames, fileName = sampleFiles, condition = sampleCondition)

View(sampleTable)

# Create DESeq2 object

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable, directory = directory, design = ~ condition)

# Run DESeq2

dds <- DESeq(ddsHTSeq)

# Generate results object

res <- results(dds, contrast=c("condition","KN_tm","Bro_tm"))
head(res)

#copy results in file

write.table(res, file = "../DEseq2/results/res.csv ",
            col.name = TRUE, 
            sep="\t", 
            row.names = TRUE, 
            quote = FALSE)


# number of genes with FDR < 0.01
sum(res_tm$padj < 0.01, na.rm=TRUE)

#total number of genes with pvalues
sum(!is.na(res_tm$padj))

#Generate heatmap
# use the log transform on the data set for heat map

rld <- rlog(dds, blind=F)
sigGenes <- subset(res,padj <0.05 & !is.na(res$padj) & abs(log2FoldChange) > 2)
sigGenes1 <- head(order(sigGenes$padj),decreasing=TRUE,50)
mat <- assay(rld)[ sigGenes1, ]
mat <- mat - rowMeans(mat)
dim(mat)

library('ComplexHeatmap')
fontsize <- 0.6
Heatmap(mat, cluster_rows = TRUE, cluster_columns = T,
        row_names_side = 'left',
        row_dend_side = 'left',
        row_names_gp = gpar(cex=fontsize),
        row_dend_width = unit(3,"cm"),
        column_names_rot = 360,
        name = "Color key")

# Generate volcano plot

library('EnhancedVolcano')
EnhancedVolcano(res_tm,lab = rownames(res_tm), x = 'log2FoldChange', y = 'pvalue', title = 'KN_tm.vs.Bro_tm',
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 2.0)
