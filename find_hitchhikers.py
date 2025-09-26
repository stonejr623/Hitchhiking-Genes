#!/usr/bin/python

""" 
Main script for finding hitchhiking genes in DE gene results
Usage: 
    python find_hitchhikers.py --deseq results.csv [optional arguments]

Inputs:
    --deseq: CSV file result from DESeq2. Must contain a column with gene identifiers (default "gene_name") matching 
    the gtf file and columns "log2FoldChange" and "pvalue"
    --genome: Name of reference genome (ex: GRCh37, GRCh38, T2T)

Outputs:
    hitchhiker_report.txt: Reports noncoding genes with a coding one in proximity at several thresholds

""" 

import argparse
import os
import subprocess

#Parse arguments
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--deseq', required=True, help='DESeq2 CSV output (contains gene_id/log2FoldChange/padj)')
    p.add_argument('--genome', default='GRCh38', help='Genome build (GRCh37, GRCh38)')
    p.add_argument('--out', default='hitchhiker_report.txt', help='Output report filename')
    p.add_argument('--id-col', default='gene_name', help='column name in DESeq CSV with gene IDs')
    p.add_argument('--split', default=False, help='Write split csv files of coding and noncoding sequences (boolean)')
    return p.parse_args()


def main():
    args = parse_args()

    # Identify paths to reference GTFs
    build = args.genome
    gene_gtf = f"./references/{build}_genes.gtf"
    te_gtf = f"./references/{build}_TEs.gtf"
    print(f"Using reference {build}.")

    # Run R script to format DESeq2 outputs
    subprocess.run([
        "Rscript", "format_inputs.R",
        "--deseq", args.deseq,
        "--gene_gtf", gene_gtf,
        "--te_gtf", te_gtf,
        "--id_col", args.id_col
    ], check=True)

    # Run Python script to find hitchhikers
    subprocess.run([
        "python", "multiple_thresholds.py",
        "--out", args.out
    ], check=True)

    # Clean up temporary files if --split not requested
    if not args.split:
        for f in ["noncoding_genes.csv", "coding_genes.csv"]:
            if os.path.exists(f):
                os.remove(f)
        print("Temporary split files removed.")
    else:
        print("Annotated csv files split into coding and noncoding.")

if __name__ == "__main__":
    main()


