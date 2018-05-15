# Mini RNA-Seq Workshop

Show two exmples on how to analyze RNA-Seq data, from raw reads to RNA-Seq differetial expression anlaysis.

## Project 1

This project has the simplest setup an RNA-Seq project have. It contains 2 conditions, 1 control and three replicates for each.

Working directory tree:

```
-── Giardia
    ├── Rcodes
    │   ├── Rcodes.Rproj
    │   ├── report.Rmd
    │   └── report.html
    ├── data
    └── scripts
        ├── count_table.sh
        ├── htseq_count.sh
        └── star.sh
```

- `scripts` contains the bash scripts for preparing the raw counts for RNA-Seq analysis in R.
	- `star.sh` maps RNA-Seq fastq files against the reference genome using `STAR`.
	- `htseq_count.sh` uses `htseq-count` to count the raw read counts for each genes from the produced BAM files
	- `count_table.sh` produces the `metadata.csv`
- `data` contains the files which are necessary for running this, not shared in public domain.
- `Rcodes` contains the R codes for RNA-Seq analysis using `edgeR`



## Project 2

This project is a more complicated setup. It's a time-serie RNA-Seq experiment targeting host-parasite interaction. More details see the report. 

Work directory tree:

-── Eimeria
    ├── Rcodes
    │   ├── edgeR.eimeria.Rmd
    │   └── edgeR.gallus.Rmd
    ├── data
    │   ├── counts
    │   ├── eimeria_annotation.tab
    │   ├── gallus_annotation.tab
    │   └── metadata.csv
    └── scripts
        ├── count_table.sh
        ├── htseq_count.sh
        └── star_merged.sh

- `scripts` contains the bash scripts for preparing the raw counts for RNA-Seq analysis in R.
	- `star_merged.sh` maps RNA-Seq fastq files against the combined reference genomes of the host and the parasite using `STAR`.
	- `htseq_count.sh` uses `htseq-count` to count the raw read counts for each genes from the produced BAM files
	- `count_table.sh` produces the `metadata.csv`
- `data` contains the files which are necessary for running this, not shared in public domain.
- `Rcodes` contains the R codes for RNA-Seq analysis using `edgeR`. Host and parasite genes are analyzed separately.
