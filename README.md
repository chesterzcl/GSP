# Genetic Signature EXPlorer (GSEXP)

## Description
We provide a command-line tool used for fast and efficient whole-genome scale genetic signature discovery.

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

cd GSEXP

Type the following to compile:

'#0969DA' g++ -std=c++17 -pthread -o GSEXP utility.h ann_data.h input_param.h var_list.h pop_data.h main_analysis_module.h thread_analysis_module.h

