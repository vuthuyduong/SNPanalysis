# SNPanalysis
The SNP pipeline is to generate a vcf file based on the raw data (fastq files) of the given samples. The script scripts/SNPpipeline.py requires an input as a json file (examples/SNP_data_B8441.json for example) containing all information about the genomes and the tools used for the analysis. Longshot (Edge et al. 2019)  is used for long sequencing reads (nanopore for example), and gatk (Van der Auwera & O'Connor 2020) is used for short sequencing reads (Illumina for example).

## Dependencies (see examples/SNP_data_B8441.json for example):

- [longshot](https://www.nature.com/articles/s41467-019-12493-y), used for long reads
- [gatk](https://gatk.broadinstitute.org/hc/en-us), used short reads
- [sratoolkit](https://github.com/ncbi/sra-tools), optional for downloading sra
- [picard](https://broadinstitute.github.io/picard/), optional for marking duplicates for short reads

## How to create a vcf file

scripts/SNPpipeline.py -i examples/SNP_data_B8441.json -o B8441_vcf -prefix allsnps

## References
Duong Vu (2021). https://github.com/vuthuyduong/SNPanalysis
