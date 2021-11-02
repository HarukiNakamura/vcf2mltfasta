#!/usr/bin/env python3
#! coding: utf-8
"""
This script 
"""

import argparse
from logging import getLogger, StreamHandler, INFO, Formatter
import sys

from src.my_utils import fasta_seq, parse_vcf, read_file, vcf2fasta, format_fasta, read_gff3


def main():
    ################ Setting of logger ################
    logger = getLogger(__name__)
    logger.setLevel(INFO)
    sh = StreamHandler()
    sh.setLevel(INFO)
    sh.setFormatter(Formatter("[%(asctime)s] %(levelname)s: %(message)s"))
    logger.addHandler(sh)
    ################ End of setting of logger ################


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
        "-k", "--keep-reference-seq", action="store_true",
        dest="keepRef", required=False, help="Output reference sequence.")

    parser.add_argument(
        "-p", "--output-fasta-prefix", type=str, action="store", default="out",
        dest="outputPrefix", required=False, help="Prefix of output fasta file. (default: out)")
    
    args = parser.parse_args()
    fasta_path: str = args.fastaFilePath
    vcf_path: str = args.vcfFilePath
    index_path: str = args.indexPath
    gff3_path: str = args.gff3FilePath
    gene_name: str = args.geneName
    keep_ref: bool = args.keepRef
    output_prefix: str = args.outputPrefix
    if (fasta_path is None) and (vcf_path is None) and (index_path is None):
        logger.error("You need to input index file or FASTA and VCF.")
        sys.exit()

    ############################## Main ##############################
    outfasta: str = output_prefix + ".fasta"

    chr, start, end = read_gff3(gff3_path=gff3_path, gene=gene_name)

    if index_path is not None:
        with open(index_path, mode="r") as index:
            fasta_path: str = next(index).rstrip("\n|\r|\r\n")
            vcf_info: list = next(index).rstrip("\n|\r|\r\n").split("\t")
            vcf_path: str = vcf_info[0]
            vcf_header_byte: int = int(vcf_info[1].split("=")[1])
            for index_line in index:
                if index_line.startswith(chr):
                    index_line = index_line.rstrip("\n|\r|\r\n")
                    _, fa_byte, vcf_chrom_byte = index_line.split("\t")
                    fa_byte = int(fa_byte)
                    vcf_chrom_byte = int(vcf_chrom_byte)
                break
            else:
                logger.error(f"No index information for {chr}.")
                sys.exit()
    else:
        fa_byte = None
        vcf_chrom_byte = None
    
    seq: str = fasta_seq(fasta_path=fasta_path, chr=chr, start=start, end=end, fa_byte=fa_byte)
    #print(seq)

    vcfs = []
    if vcf_chrom_byte is not None:
        sample_names: list = next(read_file(filepath=vcf_path, byte=vcf_header_byte)).rstrip("\n|\r|\r\n").split("\t")[9:]
    for vcf_line in read_file(filepath=vcf_path, byte=vcf_chrom_byte):
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
        if keep_ref:
            ref_name = fasta_path + " " + chr + ":" + str(start) + "-" + str(end)
            out.write(format_fasta(ref_name, seq))
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