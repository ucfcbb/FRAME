# FRAME #

## Introduction
FRAME is a tool used to estimate global ancestry composition using IBDs (Identity-By-Descent) segments. 

## Dependencies
- Linux (tested on Ubuntu 20.04.1)
- GCC (at least 9.4.0)
- Python (3.10.9) 
- Pybind11 

## Installation ##

To install the tool, first, clone the git repository in preferred directory on your local computer using the terminal:

```
git clone https://github.com/ucfcbb/FRAME.git 
cd FRAME
```

It's recommended to run FRAME in a virtual environment using `conda`. Create the virtual environment
like shown below:

```
conda env create -n FRAME --file packages.yaml
conda activate FRAME
bash build.sh
```

Add the path to the `FRAME/build` directory to the `PYTHONPATH`:

```
export PYTHONPATH="${PYTHONPATH}":"/path/to/FRAME/build/"
```

This ensures that FRAME can be run from any directory.
To test if the installation worked properly, run the following command:

```
python run.py --h
```
A list of parameters required should be shown in the terminal.

## Usage ##
To run FRAME, you **MUST** provide the following parameters:

|Parameters | Description | Requirements |
|:---:|:---:| :---:|
|`--i` | Input directory that contains the VCF files for both query and reference panel (default). When using precomputed matches(`--precomputed` flag should be used in this case), this is the directory of matches for each chromosome. | The VCF and precomputed match filename MUST start with the chromosome number. E.g. 1.query.vcf, 1.ref.vcf, 1.matches|
|`--l` | Length of the matches (in sites)| N/A |
|`--r` | Reference population file consisting of two columns| First column is the reference sample name and second column is the population label the sample belongs to. The two columns are tab separated|
|`--q` | Query sample file. | This file contains query sample IDs. One sample ID per line |
|`--n` | Number of sites per chromosome file. This file consists of two columns. | First column is the chromosome number and the second column is the number of variant sites in the VCF file pertaining to this chromosome. The columns are tab separated. |
|`--o` | Output directory to store intermediate results and the estimated ancestry proportions| N/A |

Sample run:
```
python path/to/FRAME/run.py \
    --i ../data/vcf/ \
    --l 20 \
    --r ../data/ref.pop \
    --q ../data/query.sample \
    --n ../data/sites.per.chrom \
    --o ../data/results \
```

The estimated ancestry proportions are output to the file `Frame.prop` in the specified output folder. The output folder will also contain other intermediate results if in case further analyses is needed.

## Additional Options ##
There are some optional arguments that can be used with FRAME.

|Parameters | Description | Requirements |
|:---:|:---:| :---:|
|`--B` | Bit size for Syllable-PBWT | only two values can be used: 64 or 128 (default) |
|`--precomputed` | Use this flag to use precomputed matches (possibly from tools other than Syllable-PBWT) | The match file should include 4 tab-separated columns in the following order: query haplotype id, reference haplotype id, start site, end site (exclusive). E.g.: A_0    R_1 0   13 | 
|`--c` | Chromsome number. When specified, only this chromsome file is evaluated | Only autosomes, i.e. [1,22] |


## Contact ##
- Author :Pramesh Shakya (pramesh.shakya@ucf.edu)

