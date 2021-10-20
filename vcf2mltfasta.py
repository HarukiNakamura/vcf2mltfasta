#!/usr/bin/env python3
#! coding: utf-8

from src.my_utils import fasta_seq, parse_vcf, read_file, vcf2fasta, format_fasta

fasta: str = "testData/input.fasta"
gff: str = "testData/input.gff3"
vcf: str = "testData/input.vcf"

outfasta: str = "testData/out.fasta"

####################################
# 遺伝子情報をgffから取得するコード #
####################################
chr: str = "Chr01"
start: int = 1
end: int = 70

seq: str = fasta_seq(fasta, chr, start, end)
#print(seq)

vcfs = []
for vcf_line in read_file(vcf):
    if vcf_line.startswith("##"):
        pass
    elif vcf_line.startswith("#CHROM"):
        vcf_line.rstrip("\n|\r|\r\n")
        sample_names: list = vcf_line.split("\t")[9:]
    if vcf_line.startswith(chr):
        if start <= int(vcf_line.split("\t")[1]) <= end:
            vcf_line = vcf_line.rstrip("\n|\r|\r\n")
            vcfs.append(parse_vcf(vcf_line))
#print(sample_names)
#print(vcfs)

with open(outfasta, mode="a") as out:
    for i in range(len(sample_names)):
        new_seq: str = seq
        for info in vcfs[::-1]:
            num_alt: int = info["geno"][i]
            if num_alt:
                new_seq: str = vcf2fasta(new_seq, (info["pos"] - start), info["ref"], info["alt"][num_alt-1])
        #print(new_seq)
        out.write(format_fasta(sample_names[i], new_seq))
