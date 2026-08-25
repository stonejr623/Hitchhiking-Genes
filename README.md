# Hitchhiking-Genes
Package for the identification of meaningful expression of noncoding transcripts in RNA-seq

A large issue when analyzing data for relevant noncoding transcripts such as transposable elements, miRNA, and lncRNA in RNA-seq data is increased expression due to proximity to significant coding genes. These confounding results, or 'hitchhiking' genes, can lead to an overestimation of the involvement of these transcripts in disease or other biological processes. This package creates a quick and simple way to cross-reference your DE gene list to these transcripts.

Usage: 

   ` python find_hitchhikers.py --deseq results.csv --genome GRCh38 [optional arguments] `

Inputs:

    --deseq: CSV file result from DESeq2. Must contain a column with gene identifiers (default "gene_name") matching 
    the gtf file and columns "log2FoldChange" and "pvalue"
    --genome: Name of reference genome (ex: GRCh37, GRCh38, T2T)

Outputs:

    hitchhiker_report.txt: Reports noncoding genes with a coding one in proximity at several thresholds
