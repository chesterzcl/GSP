# Genetic Signature EXPlorer (GSEXP)

## Description
We provide a command-line tool used for fast and efficient whole-genome scale genetic signature discovery. In brief, this tool is a population-frequency based scanner that can prioritize all-category genetic signatures exclusively presented in given population(s). This tool was implemented in a multi-threaded mode to deal with large population-scale whole genome sequencing(WGS) data. This tool runs on two required inputs including:
1. A standard Variant Calling Format(VCF) file containing genomic information of the entire sample collection.
2. A sample group label file indicating the membership of each individual sample.

Besides, three main customizable modules allow this tool to fulfill a variety of reserach needs. They are:
1. Loci annotation flags (position range, variant type, variant size, .etc.)
2. Variant quality metrics (minimum tolerable sequencing depth)
3. Population frequency constraints (frequency ranges in target and reference sample groups)

One successful application of this tool was described in :

Li, Z., Wang, Z., Chen, Z. et al. Systematically identifying genetic signatures including novel SNP-clusters, nonsense variants, frame-shift INDELs, and long STR expansions that potentially link to unknown phenotypes existing in dog breeds. BMC Genomics 24, 302 (2023)


## Key Features
1. Basic variant filtering and extraction based on annotation flags.
2. Can calcualte and display variant frequency by user-defined sample groups.
3. Profiling exclusively owned genetic signatures of a single sample group or shared genetic signatures among multiple sample groups.
4. Identification of bi-allelic genetic signatures statisfying user-defined population frequency constraints.
5. Identification of Short Tandem Repeats(STR) from multi-allelic loci. Can profiling sample groups for STR signatures with extreme expansion or contraction status.
6. Independent discovery and validation of genetic sigantures in different variant datasets. 

## Installation
We recommand using GCCv11.2 or greater for compilation of this tool.

After downloading the source code, move to the working directory via:
```
cd GSEXP
```
Then type the following to compile:
```
g++ -std=c++17 -pthread -o GSEXP utility.h ann_data.h input_param.h var_list.h pop_data.h main_analysis_module.h thread_analysis_module.h
```
Then the user can run the tool by directly invoking:
```
./GSEXP
```

After invoking the tool, sucessfully interpreted arguments will be printed out before analysis.

For a more detailed guide to use the tool, please move to the wiki section of this repo.

## Quick start

For signature discovery of single population with fixed-frequency based algorithm:
```
GSEXP \
unique \
-i Input.vcf \
-o Output.txt \
-p PopulationLabel.txt \
--popnum 1 \
--p1l 0.9 --p1u 1.0 --p2l 0.0 --p2u 0.1
```

For signature discovery of single population with machine-learning based adaptive algorithm:
```
GSEXP \
mldisc \
-i Input.vcf \
-o Output.txt \
-p PopulationLabel.txt \
--sample-frac 0.8 \
--popnum 1 \
--p1l 0.9 --p1u 1.0 --p2l 0.0 --p2u 0.1
```






