# CRISPR Gene Effect Correlation Analysis ----
# Correlates DNA damage repair genes with all other genes in CRISPR dataset

# Load required libraries
library(data.table)
library(tidyverse)
library(gplots)
library(readxl)
library(stringr)

# Read the CRISPR gene effect matrix
crispr_data <- fread('CRISPRGeneEffect_25Q3.csv', data.table = FALSE)

# Remove columns that are entirely NA
crispr_data <- crispr_data[, colSums(is.na(crispr_data)) == 0]

# Clean column names: extract text before the first space
colnames(crispr_data) <- str_match(string = colnames(crispr_data), pattern = '(.*) ')[, 2]

# Remove the first column
crispr_data <- crispr_data[, -1]

# Read the curated list of DNA damage repair genes
dnarepair <- fread('dna_damage_repair.txt')

# Split dataset to speed up correlation analysis ----

# Only DNA damage/repair genes
crispr_dna_repair <- crispr_data[, colnames(crispr_data) %in% dnarepair$Gene]

# All other genes 
crispr_other <- crispr_data[, !colnames(crispr_data) %in% dnarepair$Gene]

# Compute Pearson correlation between DNA repair genes and all other genes
cor_dna_row <- cor(crispr_dna_repair, 
                          crispr_other, 
                          use = 'pairwise', 
                          method = 'pearson')

# Calculate stats and ranks on correlation matrix ----

# Parameters 
TOP_N        <- 15     # number of top correlated genes to consider per row
COR_THRESHOLD <- 0.2   # minimum correlation to count a DNA association

# Depmap co-correlation frequency ----
# For each DDR gene, find the top N most correlated genes then count how often each gene appears across all rows
top_n_genes <- t(apply(cor_dna_row, 1, function(x) {
  colnames(cor_dna_row)[order(x, decreasing = TRUE)[1:TOP_N]]
}))
freq_table <- table(top_n_genes)   # frequency of each gene appearing in top N
# DNA correlation summary statistics
mean_cor_dna <- apply(cor_dna_row, 2, mean)
max_cor_dna    <- apply(cor_dna_row, 2, max)
n_dna_above    <- apply(cor_dna_row, 2, function(x) sum(x > COR_THRESHOLD))
# Rank each gene within each DNA feature
cor_dna_rank <- t(apply(cor_dna_row, 1, rank))
n_genes      <- ncol(cor_dna_row)   # total number of genes
n_genes
# Genes that achieve the top rank (rank == n_genes) in at least one DNA feature
top_ranked_genes <- colnames(cor_dna_row)[
  apply(cor_dna_rank, 2, function(x) any(x == n_genes))
]
# Mean and best rank across all DNA features
mean_rank_dna <- apply(cor_dna_rank, 2, mean)
top_rank_dna  <- apply(cor_dna_rank, 2, max)

# Assemble summary data frame ----
# All vectors are aligned to gene symbols via match() for safety
df <- data.frame(Symbol = colnames(cor_dna_row), stringsAsFactors = FALSE)

safe_match <- function(symbols, vec) vec[match(symbols, names(vec))]

# DNA-based metrics
df$dna_mean_cor   <- safe_match(df$Symbol, mean_cor_dna)
df$dna_max_cor    <- safe_match(df$Symbol, max_cor_dna)
df$dna_top_rank   <- (n_genes+1)-safe_match(df$Symbol, top_rank_dna)
df$dna_n_top15    <- safe_match(df$Symbol, freq_table)

writexl::write_xlsx(df,path = 'depmap_ddr_gene_stats.xlsx')

