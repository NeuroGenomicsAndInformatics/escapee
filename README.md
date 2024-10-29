We are performing segregation analyses addressing three potential modes of inheritance: dominant, de novo, and recessive to identify potential protective variants within a family with autosomal dominant form of Alzheimer's disease. Variants are intersected and removed if found in other participants of the family, matching the mode of inheritance (dominant, de novo and recessive). Annotation is then performed using SnpEff, CADD and REVEL scores, and gnomAD and output as a readable text file for further variant filtration.

Code is organized into 4 bash files that lay out the analyses. These files include a data preparation file for preparing sample gvcfs for analysis and segregation and annotation files for each mode of inheritance.
