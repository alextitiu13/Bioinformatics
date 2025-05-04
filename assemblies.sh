#!/bin/bash

###Homework 1 - Clinical Genomics - Guided and De-novo assemblies

###Req
##Bacterial reference genome
##Bacterial reads
##Used with example - Neisseria Meningitidis

#Confirmation prompt
echo "Starting genome assembly"
echo "============================"
read -p "This will run two genome assemblies for chosen bacterial genomes. Continue ? (y/n): " confirm
if [[ "$confirm" == [Yy] ]]; then
	echo "Starting genome assembly.."
else
	echo "Exiting..."
	exit 1
fi

##PATHS - always change this ones to the ones working for you PC
#Use relevant SRR ID/reference genome

reference="/home/alex/Desktop/sursa_date/GCF_022869645.1_ASM2286964v1_genomic.fna"
reads1="/home/alex/Desktop/sursa_date/SRR31969561_1.fastq"
reads2="/home/alex/Desktop/sursa_date/SRR31969561_2.fastq"
SPECIES="Neiserria_Meningitidis"
#Check reference/reads
###Download necessary reference genome and reads!!!

echo "Creating directories"

WORK_DIR="/TitiuAlexandru_Tema1"
mkdir -p ~/Desktop/$WORK_DIR
cd "$WORK_DIR"
echo "Changed working directory to $WORK_DIR"

##Create subdirectiories
mkdir -p ~/Desktop/$WORK_DIR/Asamblare_ghidata/
mkdir -p ~/Desktop/$WORK_DIR/Asamblare_de_novo/

echo "Subfolders have been created"

#FASTQ copy from source to working directory (paired reads)
cp -v $reads1 ~/Desktop/$WORK_DIR/Asamblare_ghidata/
cp -v $reads2 ~/Desktop/$WORK_DIR/Asamblare_ghidata/
cp -v $reads1 ~/Desktop/$WORK_DIR/Asamblare_de_novo/
cp -v $reads2 ~/Desktop/$WORK_DIR/Asamblare_de_novo/
echo "Reads have been copied"
cp -v $reference ~/Desktop/$WORK_DIR/Asamblare_ghidata/
echo "Reference genome has been copied"

### Guided assembly
echo "Guided assembly pipeline starting now: $(date)"
cd ~/Desktop/$WORK_DIR/Asamblare_ghidata/
echo "$(pwd)"

echo "Indexing reference"

bwa index "$reference"

echo "Indexing successful $(date)"

echo "Alignment starting..."

bwa mem -t 12 $reference $reads1 $reads2 > aligned_reads.sam

echo "Alignment done $(date)"
echo "Starting samtools and bcftools"

samtools view -Sb aligned_reads.sam > aligned_reads.bam
samtools sort aligned_reads.bam > sorted_reads.bam
samtools index sorted_reads.bam
bcftools mpileup --threads 10 -f "$reference" sorted_reads.bam > mpileup.bcf
bcftools call -mv -Oz -o variants.vcf.gz mpileup.bcf
# -m = multialelic caller
# -v = variants only
# -Oz = output in compressed VCF
echo "Files dones $(date)"
bcftools index variants.vcf.gz
echo "Generating consensus"

bcftools consensus -f "$reference" variants.vcf.gz > consensus_$SPECIES.fasta
echo "Consensus generated $(date)"
echo "Generating quality control on sorted reads"
echo "---------------------------------------------------------------------"

samtools flagstat sorted_reads.bam

echo "---------------------------------------------------------------------"
echo "Guided assembly done $(date)"

###De novo assembly
echo "Starting de novo assembly pipeline $(date)"

cd ~/Desktop/$WORK_DIR/Asamblare_de_novo
pwd
echo "---------------------------------------------------"
echo "k-mer size = 55"
spades -1 $reads1 -2 $reads2 -k 55 -o spades_output

echo "Spades assembly done $(date)"
echo "Copying useful files"
cp spades_output/scaffolds.fasta ./

echo "Cleaning subfolders $(date)"

cd ~/Desktop/$WORK_DIR/Asamblare_de_novo
rm -r spades_output
rm SRR31969561_1.fastq
rm SRR31969561_2.fastq

cd ~/Desktop/$WORK_DIR/Asamblare_ghidata
rm aligned_reads.*
rm GCF_022869645.1_ASM2286964v1_genomic.fna
rm mpileup.bcf
rm SRR31969561_1.fastq
rm SRR31969561_2.fastq
rm variants.vcf.gz.csi

cd ~/Desktop/$WORK_DIR/Asamblare_de_novo
quast.py -o quast_report scaffolds.fasta
cd ~/Desktop/$WORK_DIR/Asamblare_ghidata
samtools flagstat sorted_reads.bam > qc.txt
cd ~/Desktop/$WORK_DIR
echo "Assemblies done successfuly $(date)"
###END
