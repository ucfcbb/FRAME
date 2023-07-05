## FRAME ##

FRAME is a tool used to estimate global ancestry composition. It takes as an input the following parameters:

`--i` : input directory that contains matches between query and reference panel for each chromosome. 
One match file per chromosome.

`--r` : reference population file consisting of two columns i.e. reference samples and population label representing the population the sample belongs to. 
The two columns are tab separated.

`--q` : query sample file. 
Contains query sample IDs, one sample ID per line. 

`--n` : num sites per chrom file.
Consists of two columns that is tab separated.
First column is the chromosome number and the second column is the number of variant sites in the VCF file pertaining to the the chromosome.

`--o` : Output directory to store intermediate weights per chromosome 

## Details ##

- The IBD match file should report matches in 4 columns in the following order:
QueryHaplotype, ReferenceHaplotype, StartSite, EndSite(exclusive)

## Sample Usage ##
