# CSE185_project

This is the final project of Group 15 in CSE185 in spring 2023. We build a variant calling tool that is similar to [VarScan](https://varscan.sourceforge.net/)

## Installation instructions

Installation requires `pandas` to be installed. You can install it using pip: `pip install pandas`

Installation requires `numpy` to be installed. You can install it using pip: `pip install numpy`

## Basic Usage

usage: `Vcall.py [-h] [-mvf MINVARFREQ] MPILEUPFILE`

To run mypileup on a small test example (using files in this repo):

`python Vcall.py ./test/trio.mpileup -minvarfreq 0.2`

note: file directly save to same level as Vcall.py

## Options

positional arguments:

  `MPILEUPFILE`:           path to the input mpileup file

optional arguments:

  `-h`, `--help`:            show this help message and exit
  
  `-mvf MINVARFREQ`, `--minvarfreq MINVARFREQ`: 
                        minimum frequency threshold for variants (default: 0.2)
                        
## File format

The input file format is `mpileup`, generated by [samtools](http://www.htslib.org/) using `bam` file.

The output file will be in [`VCF`](https://samtools.github.io/hts-specs/VCFv4.2.pdf).

## Contributers:

Kaifu Yang and Jiayu Chen worked on this tool, with inspiration from [VarScan](https://varscan.sourceforge.net/).

Email for questions regarding the tool:

Kaifu Yang: kay002@ucsd.edu

Jiayu Chen: jic034@ucsd.edu
