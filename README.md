# vcf2mltfasta

## Create Multi-FASTA from VCF

Create Multi-FASTA from Reference FASTA, VCF, and gff3.

___

### Requirements

* Python3 (3.7 or higher)

### Usage

```text
usage: vcf2mltfasta.py -r FASTAFILEPATH -v VCFFILEPATH -g GFF3FILEPATH 
                       -t GENENAME [-k] [-p OUTPUTPREFIX]


       vcf2mltfasta.py -i INDEXPATH -g GFF3FILEPATH -t GENENAME
                       [-k] [-p OUTPUTPREFIX]



optional arguments:
  -h, --help            show this help message and exit
  -r FASTAFILEPATH, --rereference-fasta-path FASTAFILEPATH
                        Path to input fasta file.
  -v VCFFILEPATH, --input-vcf-path VCFFILEPATH
                        Path to input vcf file.
  -i INDEXPATH, --index-path INDEXPATH
                        Path to FASTA and VCF index.
  -g GFF3FILEPATH, --input-gff3-path GFF3FILEPATH
                        Path to input gff3 file.
  -t GENENAME, --target-gene-name GENENAME
                        Target gene name.
  -k, --keep-reference-seq
                        Output reference sequence.
  -p OUTPUTPREFIX, --output-fasta-prefix OUTPUTPREFIX
                        Prefix of output fasta file. (default: out)
```

*Create index for FASTA and VCF if you want to speed up.*

```text
usage: create_index.py [-h] -r FASTAFILEPATH -v VCFFILEPATH [-p OUTPUTPREFIX]

optional arguments:
  -h, --help            show this help message and exit
  -r FASTAFILEPATH, --rereference-fasta-path FASTAFILEPATH
                        Path to input fasta file.
  -v VCFFILEPATH, --input-vcf-path VCFFILEPATH
                        Path to input vcf file.
  -p OUTPUTPREFIX, --output-index-prefix OUTPUTPREFIX
                        Prefix of output index file. (default: out)
```

The output index is in the following format.  
The first line contains FASTA file name (Full path).
The second line contains VCF file name (Full path) and headers byte index.
From the third line,  
The first column contains chromosome (contig) name.  The second and third columns contain the byte indices of FASTA and VCF.  Each column is tab-delimited.

```text
/path/to/your/input/FASTA
/path/to/your/input/VCF header=10
Chr01   10  120
Chr02   100 11000
    .
    .
    .
```
