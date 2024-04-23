# Genetic Signature EXPlorer (GSEXP)

## Description
We provide a command-line tool used for fast and efficient whole-genome scale genetic signature discovery. In brief, this tool is a population-frequency based scanner that can prioritize all-category genetic signatures exclusively presented in given population(s). This tool was implemented in a multi-threaded mode to deal with large population-scale whole genome sequencing(WGS) data. This tool runs on two required inputs including:
1. A standard Variant Calling Format(VCF) file containing genomic information of the entire sample collection.
2. A sample group label file indicating the membership of each individual sample.

Besides, three main customizable modules allow this tool to fulfill a variety of reserach needs. They are:
1. Loci annotation flags (position range, variant type, variant size, .etc.).
2. Variant quality metrics (minimum tolerable sequencing depth).
3. Signature exclusivity level (targeted frequency ranges in target and reference sample groups).

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
A cmakelist file was configured to build relevant files.

To compile the binaries (the signature profiler and one testing program), first move to the root directory where the tool locates and simply run:
```
#Configure build files
cmake .

#Build binaries
make
```
GCC v8.0/Clang 5 or greater is required for the compilation of this tool.

To test the binaries, execute the tool on a toy dataset using a pre-configured script:
```
##Run the pre-configured bash script
bash ../tool_dir/test/Toy_driver_script.sh
```
In the output, the adaptive discovery module and fixed-frequency based module should identify 2 and 3 signatures from the toy dataset, respectively.

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








