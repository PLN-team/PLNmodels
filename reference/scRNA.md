# Single cell RNA-seq data

A dataset containing the counts of the 500 most varying transcripts in
the mixtures of 5 cell lines in human liver (obtained with standard 10x
scRNAseq Chromium protocol).

## Usage

``` r
scRNA
```

## Format

A data frame named 'scRNA' with 3918 rows (the cells) and 3 variables:

- counts:

  a 500 trancript by 3918 count matrix

- cell_line:

  factor, the cell line of the current row (among 5)

- total_counts:

  Total number of reads for that cell

## Source

<https://github.com/LuyiTian/sc_mixology/>
