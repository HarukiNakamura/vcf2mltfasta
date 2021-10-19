#!/usr/bin/env python3
#! coding: utf-8


fasta: str = "input.fasta"
gff: str = "input.gff3"
vcf: str = "input.vcf"


####################################
# 遺伝子情報をgffから取得するコード #
####################################
chr: str = "Chr01"
start: int = 1
end: int = 70
