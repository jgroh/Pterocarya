#!/bin/bash

REF='/Users/Jeff/workspace/Pterocarya/denovo_assembly/scaffolded/P.stenoptera_UCD1470/hap2/P.stenoptera_UCD1470_hap2.fasta'
BAMDIR="bams_to_P.ste_hap2"
VCFDIR="vcfs_to_P.ste_hap2"

for N in J.reg P.ste_hap1 P.ste_hap2 P.mac_hap1 P.mac_hap2 P.fra_hap1 P.fra_hap2 Cpal_2PA.1 Cpal_2PA.2 Cpal_2PG; do
		minimap2 -a $REF $N.fasta | samtools sort -O bam - > $BAMDIR/$N.bam &&
		samtools index $BAMDIR/$N.bam &&
		bcftools mpileup -Ou -f $REF -r Chr11:4358494-4371311 $BAMDIR/$N.bam | bcftools call -m --ploidy 1 -Oz -o $VCFDIR/$N.vcf.gz &&
		bcftools index $VCFDIR/$N.vcf.gz
done 

bcftools merge -m snps $(find $VCFDIR -name "*.vcf.gz") > toPsteHap2.vcf

vcftools --vcf toPsteHap2.vcf --mac 2 --chr Chr11 --from-bp 4358494 --to-bp 4371311 --recode --stdout > toPsteHap2_fltd.vcf && rm toPsteHap2.vcf

vcftools --vcf toPsteHap2_fltd.vcf --plink --out toPsteHap2
