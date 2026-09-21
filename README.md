# CRISPR DNA Damage Repair Co-Correlation Analysis

This R script performs a co-dependency analysis using CRISPR gene effect data to examine relationships between DNA Damage Repair (DDR) genes and all other genes across cell lines.

## Requirements
### Input Files Required in Working Directory
1. `CRISPRGeneEffect_25Q3.csv` — Genome-wide CRISPR knockout screen effect matrix.

Download CRISPR DepMap data (https://depmap.org/portal/data_page/?tab=allData) CRISPRGeneEffect.csv (v 25Q3). 
This contains gene effect estimates for all models, integrated using Chronos. Copy number corrected, scaled, and screen quality corrected.

2. `dna_damage_repair.txt` — Plain text file listing curated DDR genes.

---

### Subset Partitioning

* Imports DDR List: Reads a curated list of DNA damage and repair genes from `dna_damage_repair.txt`.
* Splits Dataset: Separates the cleaned matrix into two distinct subsets to streamline correlation calculations:
* `crispr_dna_repair`: Features corresponding exclusively to the curated DDR genes.
* `crispr_other`: All remaining genes in the dataset.

### Correlation Analysis

* Computes Pearson Correlation: Calculates pairwise Pearson correlation coefficients between every DDR gene and every other gene in the dataset.

### Statistical Metrics and Ranking

* *Top-N Co-Occurrences (`TOP_N = X`):* Identifies the top X most strongly correlated genes for each DDR gene and tallies their overall frequency across the dataset.
* *Summary Statistics:* Computes the mean and maximum correlation values for each gene relative to the DNA repair feature set, alongside counting how many DNA associations exceed a set threshold (COR_THRESHOLD).
* *Ranking Metrics:* Ranks genes across features, computing each gene's mean rank and best rank across all DNA repair features.

### Output Summary

* `dna_mean_cor`: Mean correlation with DNA repair features.
* `dna_max_cor`: Maximum correlation with any DNA repair feature.
* `dna_top_rank`: Best position-adjusted rank across DNA features.
* `dna_n_top15`: Frequency of appearing within the top 15 correlated partners for DDR genes.



