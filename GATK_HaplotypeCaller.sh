#! /bin/bash

###Req
##Indexed reference/resource files
##Same contigs for reference/resource files
##Software : BWA-Mem, samtools, picard, GATK

#input fastq files
sample="home/alex/Desktop/GATKex1/Genomuri/na12878_r1.fq home/alex/Desktop/GATKex1/Genomuri/na12878_r2.fq"
ref="home/alex/Desktop/GATKex1/Genomuri/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
hapmap="home/alex/Desktop/GATKex1/Genomuri/hapmap.fixed.vcf"
omni="home/alex/Desktop/GATKex1/Genomuri/omni.fixed.vcf"
snp="home/alex/Desktop/GATKex1/Genomuri/snp.fixed.vcf"
gatk="java -Xmx8g -jar home/alex/gatk/gatk-package-4.6.1.0-local.jar"
picard="java -jar /home/alex/picard/picard/build/libs/picard.jar"

#Create directory



#Align

date
bwa mem -t 8 $ref $sample > aligned_reads.sam

#Convert, sort and index

date
samtools view -Sb aligned_reads.sam > aligned_reads.bam
samtools sort -o sorted.bam aligned_reads.bam
samtools index sorted.bam

#Remove Duplicates

date
picard AddOrReplaceReadGroups \
	-I sorted.bam \
	-O sortedRG.bam \
	-LB lib1 \
	-PL ILLUMINA \
	-PU unit1 \
	-SM sample1 \
	-ID 1

date
picard MarkDuplicates \
	-I sortedRG.bam \
	-O noduplicates.bam \
	-M metrics.txt \
	--REMOVE_DUPLICATES true

#Haplotype Caller

date
gatk CreateSequenceDictionary -R $ref

samtools index noduplicates.bam

gatk HaplotypeCaller \
	-R $ref \
	-I noduplicates.bam \
	-O variants.g.vcf.gz \
	-ERC GVCF

#Genotype variants

date
gatk GenotypeGVCFs \
	-R $ref \
	-V variants.g.vcf.gz \
	-O genotyped_variants.vcf.gz

#Variant recalibrator

date
gatk VariantRecalibrator \
	-R $ref \
	-V genotyped_variants.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=prior,prior=15.0 $hapmap \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 $snp \
	-O recalibrate_output.recal \	
	-an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	--tranches-file output.tranches \
	--rscript-file output.plots.R

#Filter variants

date
gatk ApplyVQSR \
	-R $ref \
	-V genotyped_variants.vcf.gz \
	-O filtered_variants.vcf.gz \
	--truth-sensitivity-filter-level 99.0 \
	--recal-file recalibrate_output.recal \
	--tranches-file outputVQSR.tranches 

