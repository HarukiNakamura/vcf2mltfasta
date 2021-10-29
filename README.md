# vcf2mltfasta

## Create Multi-FASTA from VCF

Reference genome 

___

### Requirements

* Python3 (3.7 or higher)

### Usage

```text
python3 create_index.py -r <FASTA> -v <VCF>
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
The first and second line contains FASTA and VCF file name (Full path).  
From the third line,  
The first column contains chromosome (contig) name.  The second and third columns contain the byte indices of FASTA and VCF.  Each column is tab-delimited.

```text
/path/to/your/input/FASTA
/path/to/your/input/VCF
Chr01   10  120
Chr02   100 11000
    .
    .
    .
```