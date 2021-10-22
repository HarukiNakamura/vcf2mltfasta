#!/usr/bin/env python3
#! coding: utf-8
"""
This script 
"""

import argparse
from src.my_utils import fasta_seq, parse_vcf, read_file, vcf2fasta, format_fasta, read_gff3


def main():
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "-r", "--rereference-fasta-path", type=str, action="store",
        dest="fastaFilePath", required=True, help="Path to input fasta file.")

    parser.add_argument(
        "-v", "--input-vcf-path", type=str, action="store",
        dest="vcfFilePath", required=True, help="Path to input vcf file.")

    parser.add_argument(
        "-g", "--input-gff3-path", type=str, action="store",
        dest="gff3FilePath", required=True, help="Path to input gff3 file.")
    
    parser.add_argument(
        "-t", "--target-gene-name", type=str, action="store",
        dest="geneName", required=True, help="Target gene name.")

    parser.add_argument(
        "-o", "--output-fasta-prefix", type=str, action="store",
        dest="outputPrefix", required=True, help="Prefix of output fasta file.")
    
    args = parser.parse_args()
    fasta_path: str = args.fastaFilePath
    vcf_path: str = args.vcfFilePath
    gff3_path: str = args.gff3FilePath
    gene_name: str = args.geneName
    output_prefix: str = args.outputPrefix
    ################ End of setting command line arguments ################
    outfasta: str = output_prefix + ".fasta"

    chr, start, end = read_gff3(gff3_path=gff3_path, gene=gene_name)

    seq: str = fasta_seq(fasta_path=fasta_path, chr=chr, start=start, end=end)
    #print(seq)

    vcfs = []
    for vcf_line in read_file(filepath=vcf_path):
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


if __name__=="__main__":
    main()