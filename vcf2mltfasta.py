#!/usr/bin/env python3
#! coding: utf-8
"""
This script 
"""

import argparse
import sys

from src.my_utils import fasta_seq, parse_vcf, read_file, vcf2fasta, format_fasta, read_gff3


def main():
    ################ Setting command line arguments ################
    parser=argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument(
        "-r", "--rereference-fasta-path", type=str, action="store",
        dest="fastaFilePath", required=False, help="Path to input fasta file.")

    parser.add_argument(
        "-v", "--input-vcf-path", type=str, action="store",
        dest="vcfFilePath", required=False, help="Path to input vcf file.")

    parser.add_argument(
        "-i", "--index-path", type=str, action="store",
        dest="indexPath", required=False, help="Path to FASTA and VCF index.")

    parser.add_argument(
        "-g", "--input-gff3-path", type=str, action="store",
        dest="gff3FilePath", required=True, help="Path to input gff3 file.")
    
    parser.add_argument(
        "-t", "--target-gene-name", type=str, action="store",
        dest="geneName", required=True, help="Target gene name.")

    parser.add_argument(
        "-p", "--output-fasta-prefix", type=str, action="store", default="out",
        dest="outputPrefix", required=False, help="Prefix of output fasta file. (default: out)")
    
    args = parser.parse_args()
    fasta_path: str = args.fastaFilePath
    vcf_path: str = args.vcfFilePath
    index_path: str = args.indexPath
    gff3_path: str = args.gff3FilePath
    gene_name: str = args.geneName
    output_prefix: str = args.outputPrefix

    if (fasta_path is None) and (vcf_path is None) and (index_path is None):
        print("You need to input index file or FASTA and VCF.")
        sys.exit()

    ############################## Main ##############################
    outfasta: str = output_prefix + ".fasta"

    chr, start, end = read_gff3(gff3_path=gff3_path, gene=gene_name)

    if index_path is not None:
        with open(index_path, mode="r") as index:
            fasta_path: str = next(index).rstrip("\n|\r|\r\n")
            vcf_path: str = next(index).rstrip("\n|\r|\r\n")
            for index_line in index:
                if index_line.startswith(chr):
                    index_line = index_line.rstrip("\n|\r|\r\n")
                    _, fa_bit, vcf_bit = index_line.split("\t")
                    fa_bit = int(fa_bit) + 1
                    vcf_bit = int(vcf_bit) + 1
                break
            else:
                print(f"No index information for {chr}.")
                sys.exit()

    seq: str = fasta_seq(fasta_path=fasta_path, chr=chr, start=start, end=end, fa_bit=fa_bit)
    #print(seq)

    vcfs = []
    for vcf_line in read_file(filepath=vcf_path):
        if vcf_line.startswith("#CHROM"):
            vcf_line = vcf_line.rstrip("\n|\r|\r\n")
            sample_names: list = vcf_line.split("\t")[9:]
        if vcf_line.startswith(chr):
            if start <= int(vcf_line.split("\t")[1]) <= end:
                vcf_line = vcf_line.rstrip("\n|\r|\r\n")
                vcfs.append(parse_vcf(vcf_line))
    #print(sample_names)
    #print(vcfs)

    with open(outfasta, mode="w") as out:
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