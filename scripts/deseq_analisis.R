# Differential expression of E. faecium expressed mRNA in BH and Serum 
# Illumina paired-end reads were utilized.
# The count has originated on the annotation of our assembly with Prokka

# Install required packages
install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")

# Load necessary libraries
library("DESeq2")
library("ggplot2")
library("apeglm")
library("pheatmap")
library("dplyr")

#################### Preparing the data ######################

#Save the names of the count files in a vector
count_files <- c("trim_paired_ERR1797969_serum", 
                 "trim_paired_ERR1797970_serum",
                 "trim_paired_ERR1797971_serum",
                 "trim_paired_ERR1797972_bh",
                 "trim_paired_ERR1797973_bh",
                 "trim_paired_ERR1797974_bh")


#Function that reads and saves a count file 
read.file <- function(file) {
  file_name <- paste(file, "_htseq_counts.txt", sep="")
  result <- read.delim(file_name, col.names=c("gene", "count"), sep="\t", 
                       colClasses=c("character", "numeric"))
}


# Read and combine all data in one dataframe
all_data <- read.file(count_files[1])
for (c in 2:length(count_files)) {
  temp.data <- read.file(count_files[c])
  all_data <- cbind(all_data,temp.data$count)
}

#Rename the column based on the file they were obtained
colnames(all_data)[2:7] <-count_files
colnames(all_data)


#Create metadata
metadata <- data.frame(
  row.names=count_files,
  condition=c(rep("Serum", 3), rep("BH", 3)),
  libType=rep("paired-end", 6))

#It is important the column names of the all_data is the same as the rows in metadata
dsq_data<-all_data[2:7]
rownames(dsq_data)<-all_data$gene
deseq_data <- DESeqDataSetFromMatrix(countData = dsq_data,
                                     colData = metadata,
                                     design = ~ condition)
deseq_data

#Add the gene names in DESeq2 dataframe
featureData <- data.frame(gene=all_data[1])
mcols(deseq_data) <- DataFrame(mcols(deseq_data), featureData)
mcols(deseq_data)


#Normalize data
deseq_data <- estimateSizeFactors(deseq_data)
deseq_data <- estimateDispersions(deseq_data)

#Disperion plot
plotDispEsts(deseq_data)



##################### Differential expression #######################

# Specify reference level BH to be compared with the Serum
deseq_data$condition <- relevel(dds$condition, ref = "BH") 

#Perform DE
dds <- DESeq(deseq_data)
head(results(dds, tidy=TRUE))
res <- results(dds) #ordered by p-adj
summary(res)

#Scatter plot of the relationship between log fold change and 
#average expression of genes (mean normalized counts).
plotMA(res) #plotMA(dds) same as plotMA(res)

# Shrink long-fold change - Recommended for better accuracy
resultsNames(dds) #differences in plot is in different orders
resLFC <- lfcShrink(dds, coef="condition_Serum_vs_BH", type="apeglm", lfcThreshold=1)
summary(resLFC)
plotMA(resLFC, main = "Log2Fold change - Average Expression Plot of Shrinked LogFold Change")
abline(h=c(-1,1), col="dodgerblue", lwd=2)

# Extract significant genes based on adjusted p-value and fold change
sign_genes <- subset(res, res$padj < 0.001 & abs(log2FoldChange)>1 )
summary(sign_genes)


#Count data transformation
#Normal transformation: # gives log2(n + 1)
#This step is important to ensure propositional comparison of the count reads
ntd <- normTransform(dds)

# Extract top 20 over and under expressed genes
top_upregulated <- rownames(sign_genes[order(-sign_genes$log2FoldChange),])[1:20]
top_downregulated <- rownames(sign_genes[order(sign_genes$log2FoldChange),])[1:20]

# Combine the selected genes
select_top20 <- c(top_upregulated, top_downregulated)

# Plot heatmap
pheatmap(assay(ntd)[select_top20,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=metadata, main="Top 20 differentially expressed genes")



##################### Prokka Annotation #######################

# Prokka file dir
prokka_file <- "./E_faecium.tsv"

# Import prokka annotation file
prokka_data <- read.delim(prokka_file, header = TRUE, sep = "\t")

# Filter top upregulated
filtered_upregulated <- prokka_data %>%
  filter(locus_tag %in% top_upregulated)

# Filter top downregulated
filtered_downregulated <- prokka_data %>%
  filter(locus_tag %in% top_downregulated)

# Set output paths
outfile_upreg<- "./filtered_upregulated_annotations.txt"
outfile_downreg <- "./filtered_downregulated_annotations.txt"

# Extract significant genes in .tsv 
write.table(filtered_upregulated, file = outfile_upreg, sep = "\t", row.names = FALSE, quote = FALSE)
write.table(filtered_downregulated, file = outfile_downreg, sep = "\t", row.names = FALSE, quote = FALSE)



                    
"""
#This step was performed to check that all files have the same 
#number of genes with the same name in order (per row)

#Read all files
for (fl in count_files) {
  # Construct variable name
  var_name <- paste( fl, sep="")
  # Read file and assign to variable
  assign(var_name, read.file(fl))
}

#Check that the gene names are all the same in the same order

#Function to check if two variables have the same gene names
check_gene_names <- function(var1, var2) {
  all(var1$gene == var2$gene)
}

#Empty matrix to store the comparison results
comparison_matrix <- matrix(NA, nrow = length(count_files), ncol = length(count_files))

#Iterate over pairs of variables and compare their gene names
for (i in 1:length(count_files)) {
  for (j in 1:length(count_files)) {
    var1 <- get(count_files[i])
    var2 <- get(count_files[j])
    comparison_matrix[i, j] <- check_gene_names(var1, var2)
  }
}

# Show the comparison matrix
comparison_matrix
""" 
  


