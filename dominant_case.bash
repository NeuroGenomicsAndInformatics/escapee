# For the dominant case
# File list must have escapee as first file
VCFS=()
for VCF in $(cat sample_vcf_list_dom.txt); do
	VCFS+="$VCF "
done

# Intersect and remove variants found in the other samples
bcftools isec \
	-w1 \
	-c some \
	-C \
	${VCFS[@]} \
	| bgzip -c > escapee-subtracted.dom.vcf.gz
		
# Removing homozygous variants because they're not possible with this expected rare variant
bcftools view \
	-i 'GT[*]="RA"' \
	escapee-subtracted.dom.vcf.gz \
	| bgzip -c > escapee.dom.het.vcf.gz
	
# switch chr form to 1..22 X Y for pipeline
bcftools annotate \
	--rename-chrs chr_change.txt \
	escapee.dom.het.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.vcf.gz

# SnpEff annotation for predicted impact
snpeff \
	-v GRCh38.105 \
	escapee.dom.het.chrs.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.vcf.gz

# dbNSFP annotation
snpsift dbnsfp \
	-db dbNSFP4.1a.txt.gz \
	-v escapee.dom.het.chrs.snpeff.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.vcf.gz

# dbSNP annotation
snpsift annotate \
	-id dbSNP_b153_GRCh38.gz \
	-v escapee.dom.het.chrs.snpeff.dbnsfp.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.vcf.gz

# CADD score
echo '##INFO=<ID=CADD_PHRED,Number=A,Type=Float,Description="CADD PHRED score">' > CADD_hdr.txt
bcftools annotate \
	-a CADD/release_v1.6/whole_genome_SNVs.tsv.gz \
	-c CHROM,POS,REF,ALT,-,CADD_PHRED \
	-h CADD_hdr.txt \
	escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.vcf.gz

# REVEL score
echo '##INFO=<ID=REVEL_SCORE,Number=A,Type=Float,Description="REVEL score">' > revel_hdr.txt
bcftools annotate \
	-a revel_with_transcript_ids_fixed.tsv.gz \
	-c CHROM,POS,REF,ALT,REVEL_SCORE \
	-h revel_hdr.txt \
	escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.vcf.gz

# gnomAD annotation
# Need to change chromosome format back to chr1,chr2,chrX
bcftools annotate \
	--rename-chrs chr_changeback.txt \
	escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.chrs.vcf.gz

for CHR in chr{1..22} chrX chrY; do
	bcftools view \
	-r $CHR \
	escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.chrs.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.${CHR}.vcf.gz
done

parallel -j 24 \
	"snpsift annotate \
	-id \
	-info AF,AC,AN,popmax,AF_popmax,AC_popmax,AN_popmax,AF_nfe,AC_nfe,AN_nfe,AF_afr,AC_afr,AN_afr,AF_amr,AC_amr,AN_amr,AF_sas,AC_sas,AN_sas,cadd_phred \
	gnomad.genomes.v3.1.2.sites.{}.vcf.bgz \
	escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.{}.vcf.gz \
	| bgzip -c > escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.gnomad.{}.vcf.gz" \
	::: chr{1..22} chrX chrY

# Merge chr vcfs back together
VCFS=()
for VCF in $(find $PWD -name "*dom*gnomad.*.vcf.gz"); do
	VCFS+="-I $VCF "
done

gatk MergeVcfs \
	${VCFS[@]} \
	-O escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.gnomad.ALLchr.vcf.gz

# make a readable text file
zcat escapee.dom.het.chrs.snpeff.dbnsfp.dbsnp.CADD.REVEL.gnomad.ALLchr.vcf.gz \
	| snpEff_v5.1d/scripts/vcfEffOnePerLine.pl \
	| snpsift extractFields \
	-e "." - CHROM POS ID REF ALT "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "CADD_PHRED" "AF_nfe" "REVEL_SCORE" > escapee.dom.FIELDS.txt