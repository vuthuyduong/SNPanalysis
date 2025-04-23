[![DOI](https://zenodo.org/badge/374599160.svg)](https://zenodo.org/badge/latestdoi/374599160)

# SNPanalysis
The SNP pipeline is to generate a vcf file based on the raw data (fastq files) of the given samples. The script scripts/SNPpipeline.py requires an input as a json file (examples/SNP_data_B8441.json for example) containing all information about the genomes and the tools used for the analysis. Longshot (Edge et al. 2019)  is used for long sequencing reads (nanopore for example), and gatk (Van der Auwera & O'Connor 2020) is used for short sequencing reads (Illumina for example).

## Dependencies (see examples/SNP_data_B8441.json for example):

- [longshot](https://www.nature.com/articles/s41467-019-12493-y), used for long reads
- [minimap2](https://academic.oup.com/bioinformatics/article/34/18/3094/4994778), used to map long reads to the reference genome
- [gatk](https://gatk.broadinstitute.org/hc/en-us), used short reads
- [bwa](https://bio-bwa.sourceforge.net/), used to map short reads to the reference genome 
- [sratoolkit](https://github.com/ncbi/sra-tools), optional for downloading sra
- [picard](https://broadinstitute.github.io/picard/), optional for marking duplicates for short reads
- [vcftools](https://github.com/vcftools/vcftools), for merging vcf files
- [sratoolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit), optional, for downloading fastq files if SRA accession numbers are given

## How to create a vcf file

scripts/SNPpipeline.py -i examples/SNP_data_B8441.json -o B8441_vcf -prefix allsnps

## References
Duong Vu (2023). Optimized SNP variant calling pipeline for fungal genetics population analysis. https://github.com/vuthuyduong/SNPanalysis. [DOI: 10.5281/zenodo.8046747](https://doi.org/10.5281/zenodo.8046747)
