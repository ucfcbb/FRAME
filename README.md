## FRAME ##

FRAME is a tool used to estimate global ancestry composition. It takes as an input the following parameters:

## Installation ##

The tool was tested on a Linux system with following configuration:
`gcc 9.4.0`
`ubuntu 20.04.1`

First, clone the git repository in a preferred directory on your local computer:

```
$git clone https://github.com/jikhashkya/FRAME.git'
$cd FRAME
```

It's recommended to run FRAME in a virtual environment using `conda`. Create the virtual environment
like shown below:

```
conda create -n FRAME
conda activate FRAME
conda install --yes --file requirements.txt
bash build.sh
```

Add the path to the `FRAME/build` directory to the `PYTHONPATH`:

`export PYTHONPATH="${PYTHONPATH}":"/path/to/FRAME/build/"`

This ensures that FRAME can be run from any directory.
To test if the installation worked properly, run the following command:

```
python run.py --h
```
A list of options required should be shown in the terminal.

## Sample Usage ##
To run FRAME, you *MUST* provide the following parameters:

|Parameters | Description |
|---:| ---:|
|`--i` | input directory that contains matches between query and reference panel for each chromosome. One match file per chromosome|
|`--r` | reference population file consisting of two columns i.e. reference samples and population label representing the population the sample belongs to. 
The two columns are tab separated|
|`--q` | query sample file. This file contains query sample IDs, one sample ID per line|
|`--n` | Num sites per chrom file. Consists of two columns that is tab separated. First column is the chromosome number and the second column is the number of variant sites in the VCF file pertaining to the the chromosome|
|`--o` | Output directory to store intermediate results and the estimated ancestry proportions|


Sample run:
```
python path/to/FRAME/run.py \
    --i \
    --i \
    --i \
    --i \
```

## Details ##

- The IBD match file should report matches in 4 columns in the following order:
QueryHaplotype, ReferenceHaplotype, StartSite, EndSite(exclusive)

