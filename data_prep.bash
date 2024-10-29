# This code prepares data for intersections
# Genotype gvcfs to remove <NON_REF> values
for VCF in $(find $PWD -name "*.vcf.gz"); do 
	gatk GenotypeGVCFs \
	-R Homo_sapiens_assembly38.fasta \
	-V $VCF \
	-O ${VCF%.vcf.gz}_gted.vcf.gz
done
	
# Split multiallelics
for VCF in $(find $PWD -name "*_gted.vcf.gz"); do 
	bcftools norm \
	-m -any \
	-Oz \
	-o ${VCF%.*.*}.splitMA.vcf.gz \
	$VCF

	tabix ${VCF%.*.*}.splitMA.vcf.gz
done