#!/usr/bin/python

with open("GRCh38_Ensembl_rmsk_TE.gtf") as infile, open("reformatted.gtf", "w") as outfile:
    for line in infile:
        outfile.write(line.replace("class_id", "gene_biotype"))

