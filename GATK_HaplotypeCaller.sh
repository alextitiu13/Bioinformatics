#! /bin/bash

###Req
##Indexed reference/resource files
##Same contigs for reference/resource files
##Software : BWA-Mem, samtools, picard, GATK

#Confirmation prompt
echo "Starting GATK HaplotypeCaller pipeline"
echo "==============================================================="
read -p "This will run a complete variant calling pipeline, be sure you changed working directory, named samples inside the script, changed paths to reference/resource files and GATK/Picard software. Continue ? (y/n): " confirm
if [[ "$confirm" == [Yy] ]]; then
        echo "Starting GATK Analysis..."
else
        echo "Exiting..."
        exit 1
fi

#input fa/fq files - Use sample files from working directory

sample="na12878_r1.fq na12878_r2.fq"

#Reference genome - INDEXED

ref="/home/alex/Desktop/REF_RES/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

#Path to resource files - CORECTED CONTIGS / INDEXED

hapmap="/home/alex/Desktop/REF_RES/hapmap.fixed.vcf"
omni="/home/alex/Desktop/REF_RES/omni.fixed.vcf"
snp="/home/alex/Desktop/REF_RES/snp.fixed.vcf"
dbsnp="/home/alex/Desktop/REF_RES/dbsnp_138.fixed.vcf"

#Path to software

gatk="java -Xmx8g -jar /home/alex/gatk/gatk-package-4.6.1.0-local.jar"
picard="java -jar /home/alex/picard/picard/build/libs/picard.jar"

#Set working directory

WORK_DIR="/home/alex/Desktop/Test_Script_Final"
cd "$WORK_DIR"
echo "Changed working directory to $WORK_DIR"

#Align

echo "Starting alignment $(date)"
bwa mem -t 8 $ref $sample > aligned_reads.sam

#Convert, sort and index

echo "Starting converting, sorting, indexing files $(date)"
samtools view -Sb aligned_reads.sam > aligned_reads.bam
samtools sort -o sorted.bam aligned_reads.bam
samtools index sorted.bam

#Remove Duplicates

echo "Adding ReadGroups $(date)"
$picard AddOrReplaceReadGroups \
	-I sorted.bam \
	-O sortedRG.bam \
	-LB lib1 \
	-PL ILLUMINA \
	-PU unit1 \
	-SM sample1 \
	-ID 1

echo "Starting picard MarkDuplicates $(date)"
$picard MarkDuplicates \
	-I sortedRG.bam \
	-O noduplicates.bam \
	-M metrics.txt \
	--REMOVE_DUPLICATES true

#Haplotype Caller

echo "Creating .dict file fore reference, indexing noduplicates file and starting HaplotypeCaller $(date)"
$gatk CreateSequenceDictionary -R $ref

samtools index noduplicates.bam

$gatk HaplotypeCaller \
	-R $ref \
	-I noduplicates.bam \
	-O variants.g.vcf.gz \
	-ERC GVCF

#Genotype variants

echo "Genotyping variants $(date)"
$gatk GenotypeGVCFs \
	-R $ref \
	-V variants.g.vcf.gz \
	-O genotyped_variants.vcf.gz

#Variant recalibrator

echo "Starting Variant Recalibrator $(date)"
$gatk VariantRecalibrator -R $ref -V genotyped_variants.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 $snp \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
	-O recalibrate_output.recal \
	-an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	--tranches-file output.tranches \
	--rscript-file output.plots.R

#Filter variants

echo "Aplying VQSR filters $(date)"
$gatk ApplyVQSR \
	-R $ref \
	-V genotyped_variants.vcf.gz \
	-O filtered_variants.vcf.gz \
	--truth-sensitivity-filter-level 99.0 \
	--recal-file recalibrate_output.recal \
	--tranches-file output.tranches

##END
echo "Pipeline completed $(date), check files for validation"
