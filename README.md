# bedtools
This repository made for python cosplay of the bedtools
## Functionaloty
The are sort, merge, getfasta, intersect -w available. See https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html for more information
## Example of usage
python3 bedtools.py --sort -fb lamina.bed
python3 bedtools.py --merge -fb lamina.bed 
python3 bedtools.py --getfasta -fb sample.bed --fasta_file sequence.fasta
piython3 bedtools.py --intersect -fb bed_1.bed -sb bed_2.bed

