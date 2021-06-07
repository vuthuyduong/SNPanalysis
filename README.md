# SNPanalysis
The SNP pipeline is to generate a vcf file based on the raw data (fastq files) of the given samples. The script scripts/SNPpipeline.py requires an input as a json file (examples/SNP_data_B8441.json for example) containing all information about the genomes and the tools used for the analysis. Longshot (Edge et al. 2019)  is used for long sequencing reads (nanopore for example), and gatk (Van der Auwera & O'Connor 2020) is used for short sequencing reads (Illumina for example).

## References
Edge, P., Bansal, V. Longshot enables accurate variant calling in diploid genomes from single-molecule long read sequencing. Nat Commun 10, 4660 (2019).

Van der Auwera GA & O'Connor BD. (2020). Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). O'Reilly Media.
