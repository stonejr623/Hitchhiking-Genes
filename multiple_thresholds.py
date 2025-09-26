#!/usr/bin/python

import pandas as pd
import argparse

# parse arguments
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--out', required=True)
    return p.parse_args()

# thresholds (bp) to test
thresholds = [1000, 5000, 10000]

# import HML-2 loci and DE genes
noncoding = pd.read_csv("noncoding_genes.csv")
coding = pd.read_csv("coding_genes.csv")

# make dict from DE genes (grouped by chromosome)
de_dict = {}
for _, row in coding.iterrows():
    chrom = row['chr']
    de_dict.setdefault(chrom, []).append(row)

args = parse_args()
filename = args.out

# make file for results
with open(filename, "w") as out:
    for distance in thresholds:
        matches = {}

        # check each noncoding locus
        for _, seq in noncoding.iterrows():
            chrom = seq['chr']
            start = seq['start']
            end = seq['end']
            name = seq['gene/TE']
            strand = seq['strand']
            x = 'up' if seq['log2FoldChange'] > 0 else 'down'

            nearby_genes = {}
            for gene in de_dict.get(chrom, []):
                if gene['strand'] != strand:
                    continue  # skip wrong strand
                y = 'up' if gene['log2FoldChange'] > 0 else 'down'
                if x != y:
                    continue  # skip opposite regulation

                gene_start = gene['start']
                gene_end = gene['end']

                # compute distances
                dist1 = start - gene_end   # distance from locus start to gene end
                dist2 = gene_start - end   # distance from gene start to locus end

                if dist1 < 0 and dist2 < 0:
                    # overlap case
                    nearby_genes[gene['gene/TE']] = 0
                else:
                    # upstream/downstream case
                    if 0 <= dist1 <= distance:
                        nearby_genes[gene['gene/TE']] = dist1
                    if 0 <= dist2 <= distance:
                        nearby_genes[gene['gene/TE']] = dist2

            if nearby_genes:
                matches[name] = nearby_genes

        # calculate percentage
        num_with_match = len(matches)
        total_loci = len(noncoding)
        percent = num_with_match / total_loci * 100

        # write section for this threshold
        out.write(f"=== Results for threshold {distance} bp ===\n")
        out.write("Matched transcript and nearby DE genes:\n")
        for seq, genes in matches.items():
            pairs = [f"{g} (dist={d})" for g, d in genes.items()]
            out.write(f"{seq}: {'; '.join(pairs)}\n")
        out.write(f"\nPercentage of transcripts with a DE gene within {distance} bp: {percent:.2f}%\n\n")

