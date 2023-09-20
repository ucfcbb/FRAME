## FRAME ##

FRAME is a tool used to estimate global ancestry composition. It takes as an input the following parameters:

## Installation ##

Following are the prerequisites to install FRAME:
`gcc`

First, clone the git repository in a preferred directory on your local computer:
`$git clone https://github.com/jikhashkya/FRAME.git'
`$cd FRAME`

It's best advised to run FRAME in a virtual environment using `conda`. To do so:
`conda create -n FRAME`
`conda activate FRAME`
`conda install --yes --file requirements.txt`
`bash build.sh`

Add the path to the build directory to the PYTHONPATH like so:
`export PYTHONPATH="${PYTHONPATH}":"/path/to/FRAME/build/"`
This ensures that FRAME can be run from any directory.


To test if the installation worked properly, run the following command:
``
The expected output is: 
``

## Sample Usage ##
To run FRAME, you must provide the following parameters:

`--i` : input directory that contains matches between query and reference panel for each chromosome. One match file per chromosome.

`--r` : reference population file consisting of two columns i.e. reference samples and population label representing the population the sample belongs to. 
The two columns are tab separated.

`--q` : query sample file. 
Contains query sample IDs, one sample ID per line. 

`--n` : num sites per chrom file.
Consists of two columns that is tab separated.
First column is the chromosome number and the second column is the number of variant sites in the VCF file pertaining to the the chromosome.

`--o` : Output directory to store intermediate weights per chromosome 

Sample run:
```
python run.py \
    --i \
    --i \
    --i \
    --i \
```

## Details ##

- The IBD match file should report matches in 4 columns in the following order:
QueryHaplotype, ReferenceHaplotype, StartSite, EndSite(exclusive)

