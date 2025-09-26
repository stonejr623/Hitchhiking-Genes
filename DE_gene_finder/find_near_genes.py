#!/usr/bin/python 

import pandas as pd

#set distance to look at
thresholds = [0, 1000, 5000, 10000, 50000]

#import HML-2 loci and DE genes
noncoding = pd.read_csv("noncoding_genes.csv")
coding = pd.read_csv("coding_genes.csv")

#make dict from DE genes
de_dict = {}
for _, row in coding.iterrows():
    chrom = row['chr']
    de_dict.setdefault(chrom, []).append(row)

def calculate_distance(a_start, a_end, b_start, b_end):
    if a_end < b_start:
        return b_start - a_end
    if b_end < a_start:
        return a_start - b_end
    return 0

with open("matches_report.txt", "w") as out:
    for distance in thresholds:
        matches = {}
        for _, seq in noncoding.iterrows():
            chrom = seq['chr']
            start = seq['start']
            end = seq['end']
            name = seq['gene/TE']
            strand = seq['strand']
            if seq['log2FoldChange'] > 0:
                x = 'up'
            else:
                x = 'down'

            nearby_genes = {}
            for _, gene in pd.DataFrame(de_dict.get(chrom, [])).iterrows():
                if gene['strand'] != strand:
                    continue  # skip genes on wrong strand
                if (gene['log2FoldChange'] > 0) != (seq['log2FoldChange'] > 0):
                    continue  # skip opposite regulation

                d = interval_distance(start, end, gene['start'], gene['end'])
                if d <= distance:
                    nearby_genes[gene['gene/TE']] = d

            if nearby_genes:
                matches[name] = nearby_genes

        # calculate percentage
        num_with_match = len(matches)
        total_loci = len(noncoding)
        percent = num_with_match / total_loci * 100

        # write results
        out.write(f"=== Results for threshold {distance} bp ===\n")
        for seq, genes in matches.items():
            formatted = "; ".join([f"{g} ({d} bp)" for g, d in genes.items()])
            out.write(f"{seq}: {formatted}\n")
        out.write(f"\nPercentage of loci with a DE gene within {distance}: {percent:.2f}%\n\n")

#find genes within distance of loci
matches = {}
for _, seq in noncoding.iterrows():
    chrom = seq['chr']
    start = seq['start']
    end = seq['end']
    name = seq['gene/TE']
    strand = seq['strand']
    if seq['log2FoldChange'] >0:
        x = 'up'
    else:
        x = 'down'

    nearby_genes = {}
    for gene in de_dict.get(chrom, []):
        if gene['strand'] != strand:
            continue #skip genes on wrong strand
        if gene['log2FoldChange'] > 0:
            y = 'up'
        else:
            y = 'down'
        if x != y:
            continue #skip genes opposite regulated
        gene_start = gene['start']
        gene_end = gene['end']
        
        if not (end < gene_start or start > gene_end):
            distance_val = 0
            nearby_genes[gene['gene/TE']] = distance_val
        else:
            a = start - gene_end
            b = gene_start - end
            if 0 <= a <= distance:
                nearby_genes[gene['gene/TE']] = a
            if 0 <= b <= distance:
                nearby_genes[gene['gene/TE']] = b

    if nearby_genes:
        matches[name] = nearby_genes

#calculate percentage
num_with_match = len(matches)
total_loci = len(noncoding)
percent = num_with_match / total_loci * 100

#output results
print("Matched loci and nearby DE genes:")
for seq, genes in matches.items():
    details = [f"{g} ({d} bp)" for g, d in genes.items()]
    print(f"{seq}: {'; '.join(details)}")
print(f"\nPercentage of loci with a DE gene within {threshold}: {percent:.2f}%")


