#!/usr/bin/env python3
#! coding: utf-8
"""
This script create index for FASTA and VCF.
The output index is in the following format.
The first and second line contains FASTA and VCF file name (Full path).
From the third line,
The first column contains chromosome (contig) name.
The second and third columns contain the byte index of FASTA and VCF.
Each column is tab-delimited.
"""

import argparse
import pathlib

from src.my_utils import read_file


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
        "-p", "--output-index-prefix", type=str, action="store", default="out",
        dest="outputPrefix", required=False, help="Prefix of output index file. (default: out)")

    args = parser.parse_args()
    fasta_path: str = args.fastaFilePath
    vcf_path: str = args.vcfFilePath
    output_prefix: str = args.outputPrefix

    out_index: str =  output_prefix + ".index.txt"
    ################ End of setting command line arguments ################


    ################################ Main ################################
    fa_byte: int = 0
    fa_byte_dic: dict = {}
    print("Reading input FASTA...")
    for line in read_file(fasta_path):
        fa_byte += len(line)
        if line.startswith(">"):
            chr = line.lstrip(">").rstrip("\n|\r|\r\n")
            fa_byte_dic[chr] = str(fa_byte)
    print(f"Found {len(fa_byte_dic)} chromosomes (contigs) in FASTA.")


    vcf_byte: int = 0
    vcf_byte_dic: dict = {}
    tmp_chr: str = ""
    print("Reading input VCF...")
    for line in read_file(vcf_path):
        if line.startswith("#"):
            vcf_byte += len(line)
        else:
            chr: str = line.split("\t")[0]
            if tmp_chr != chr:
                vcf_byte_dic[chr] = str(vcf_byte)
                tmp_chr = chr
            vcf_byte += len(line)
    print(f"Found {len(vcf_byte_dic)} chromosomes (contigs) in VCF.")


    # FASTAが全ゲノムで、VCFが一部の染色体だけ、ということもあると思うあるので
    # FASTAとVCFで共通して情報のある染色体のみ、そのバイトインデックスを出力する。
    # set型を用いた方が速いが、順番が保持されないのでNG。
    with open(out_index, mode="w") as out:
        out.write(str(pathlib.Path(fasta_path).resolve()) + "\n")
        out.write(str(pathlib.Path(vcf_path).resolve()) + "\n")
        for chr in fa_byte_dic.keys():
            if chr in vcf_byte_dic.keys():
                out.write("\t".join([chr, fa_byte_dic[chr], vcf_byte_dic[chr]])+"\n")
    ############################# End of Main #############################


if __name__=="__main__":
    main()