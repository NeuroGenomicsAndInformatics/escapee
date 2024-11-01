We are performing segregation analyses addressing three potential modes of inheritance: dominant, de novo, and recessive to identify potential protective variants within a family with autosomal dominant form of Alzheimer's disease. Variants are intersected and removed if found in other participants of the family, matching the mode of inheritance (dominant, de novo and recessive). Annotation is then performed using SnpEff, CADD and REVEL scores, and gnomAD and output as a readable text file for further variant filtration.

Code is organized into 4 bash files that lay out the analyses as run for our analysis. They are not strictly runnable scripts, but would be able to be run with minimal changes. These files include a data preparation file for preparing sample gvcfs for analysis and segregation and annotation bash files for each mode of inheritance.

This analysis utilized the following software:
* bcftools v1.16
* gatk v4.5
* SnpEff v5.1d
* SnpSift v5.1d

This analysis utilized the following reference files:
* GRCh38 reference fasta
* SnpEff database GRCh38.105
* dbNSFP v4.1a
* dbSNP b153 GRCh38
* CADD release v1.6
* REVEL database
* gnomAD v3.1.2
