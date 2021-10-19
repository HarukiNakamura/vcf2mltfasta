#!/usr/bin/env python3
#! coding: utf-8
'''
関数集
'''

import re

from typing import Generator, List



def read_file(filepath: str) -> Generator[str, None, None]:
    """
    Return a text content one line at a time.

    Parameters
    ----------
    filepath : str
        File path you want to read.

    Yields
    -------
    Generator[str, None, None]
        Text content.
    """
    with open(filepath, encoding='utf-8', mode="r") as f:
        for line in f:
            yield line.rstrip('\n|\r|\r\n')


def fasta_seq(fasta_path: str, chr: str, start: int, end: int, fasta_index_path: str=None) -> str:
    """
    Read fasta file and return target seqence.
    This function similar to "samtools faidx <fasta_path> <chr>-<start>:<end>"

    Parameters
    ----------
    fasta_path : str
        Fasta file path.
    chr : str
        Target chromosome name.
    start : int
        Start position of the target sequence.
    end : int
        End position of the target sequence
    fasta_index_path : str, optional
        Fasta index file path, by default None.

    Returns
    -------
    seq: str
        Target sequence.
    """
    c: int = 0
    seq: str = ""
    # fastaのインデックスがない場合、とりあえずターゲットの染色体まで読みすすめる。
    with open(fasta_path, mode="r") as fasta:
    #インデックスファイルがない場合該当する染色体番号に至るまで読み飛ばす。
        if fasta_index_path is None:
            for fa_line in fasta:
                if fa_line.startswith(">" + chr):
                    break
        #インデックスファイルがある場合、seekで読み飛ばす。
        else:
            for fai_line in read_file(fasta_index_path):
                if fai_line.startswith(chr):
                    bit = fai_line.split("\t")[2]
                    fasta.seek(int(bit))
        #読み飛ばしたところから始める。
        for fa_line in fasta:
            fa_line = fa_line.rstrip("\n|\r|\r\n")
            #一つのfa_lineはc+1 ~ c+len(fa_line)文字目の配列を持っていることを念頭に入れて条件分岐
            if start > c + len(fa_line):
                pass
            elif c > end:
                break
            else:
                seq += fa_line[max(start - 1 - c, 0):min(end -c , len(fa_line))]
            c += len(fa_line)
    return seq


def geno2numeric(geno: str) -> int:
    """
    Change genotype field of input VCF to numeric data.
    For example, 0/0 -> 0 and 1/1 -> 1.

    Parameters
    ----------
    geno : str
        Genotype filed of input VCF.

    Returns
    -------
    num: int
        Genotype number.
    """
    # VCFの0/0:3,0:3:9:0,9,103という表記が["0", "0"]に変換する。
    gt_list: list[str] = re.split("/|\|", geno.split(":")[0])

    # 数字に変換
    # ホモの場合
    if gt_list[0] == gt_list[1]:
        num: str = gt_list[0]
        if num == ".":
            num = "0"
        num = int(num)
    # ヘテロの場合、数字として小さい方を採用
    else:
        num: int = min(int(gt_list[0]), int(gt_list[1]))
    return num


def parse_vcf(vcf_line: str) -> dict:
    """
    Parse VCF body line.
    One data line will be conberted to like following format.
    {"pos":120, "ref":"A", "alt":["ATT", "ATTTTTT"], "geno":[0,0,1,2,0]}
    [pos, ref, list[alt], list[geno_num]]

    Parameters
    ----------
    vcf_line : str
        One data line of input VCF.

    Returns
    -------
    dict_line: dict
        {"pos":POS, "ref":REF, "alt":[ALT], "geno":[GENO]]}
    """
    dict_line = {}
    splited_line: list = vcf_line.split("\t")
    dict_line["pos"] = int(splited_line[1])
    dict_line["ref"] = splited_line[3]
    dict_line["alt"] = splited_line[4].split(",")
    dict_line["geno"] = list(map(geno2numeric, splited_line[9:]))
    return dict_line


def vcf2fasta(seq: str, index: int, ref: str, alt: str) -> str:
    """
    Change fasta based on variant information.

    Parameters
    ----------
    seq : str
        Fasta sequence.
    index : int
        Relative position of the variant.
    ref : str
        Reference sequence.
    alt : str
        Alternative sequence.

    Returns
    -------
    new_seq: str
        Modified sequence.
    """
    # 一旦リストに変換する。
    # 配列は常にリストで扱ったほうが速い？
    new_seq: list[str] = list(str)
    # SNP or insertion
    if len(ref) <= len(alt):
        new_seq[index] = alt
    # deletion
    else:
        seq[index:index+len(ref)] = alt
    new_seq: str = "".join(seq)
    return new_seq
