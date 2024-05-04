#!/bin/bash -l

############# Setting up Uppmax #############

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 06:00:00
#SBATCH -J gen_assembly
#SBATCH --mail-type=ALL
#SBATCH --mail-user aesimaiaki@gmail.com
#SBATCH --output=%x.%j.out

######### 2.Genome assembly ############                                                                                                                                                     
#The present script is hardcoded for the Illumina 

#Load modules 
module load Flye/2.9.1
module load spades/3.15.5
module load quast/5.0.2
module load bwa/0.7.17
module load bwa-mem2
module load samtools/1.19
module load Pilon/1.24
module load augustus/3.5.0-20231223-33fc04d
module load BUSCO/5.5.0


#Set variables
data_path="/home/agsi4287/Genome_Analysis/raw_data/genomics_data"
results_path="/home/agsi4287/Genome_Analysis/results"
gen_data_path="../generated_data"
pacbio_path="$results_path/pac_assembly"
spades_path="$results_path/spades_assembly"


#2.1 Genome assembly with PacBio reads with Flye
#Note: To set the genome size value, the NCBI the reference
#genome of Enterococcus faecium was used.

flye --pacbio-raw $data_path/PacBio/*  --out-dir ${results_path}/pac_assembly --genome-size 2.9m --threads 4                                                                                

#Improving PacBio assembly with Pilon
#Pilon uses Illumina short reads to correct base errors, 
#small indels and missassemblies.

#Preparing short reads for Pilon
#Mapping Illumina paired-end reads to PacBIo
bwa index $pacbio_path/assembly.fasta
bwa mem -t 8 $pacbio_path/assembly.fasta $data_path/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz $data_path/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz > $pacbio_path/assembly_pac.sam

#Convert SAM to BAM file, sort and index BAM file
samtools view -bo $results_path/pac_assembly/assembly_pac.bam $results_path/pac_assembly/assembly_pac.sam
samtools sort $pacbio_path/assembly_pac.sam -o $pacbio_path/assembly_pac_sorted.bam
samtools index $pacbio_path/assembly_pac_sorted.bam

#Improve PacBIo assembly
java -jar $PILON_HOME/pilon.jar  --genome  $pacbio_path/assembly.fasta --frags $pacbio_path/assembly_pac_sorted.bam --output pac_pilon_assembly --outdir $pacbio_path
                                                                                                                                                                                       
                                                              
#2.2 Genome assembly with PacBio & Illuma reads with Spades                                                                                                                                  
#Spades assembles short Illumina reads and then uses the PacBio reads to scaffold and conjuct reads together                                                                                 

#Create a folder to save  all data there                                                                                                                                                     
mkdir ../generated_data

#Create a file with all PacBio data in one as spades is not friendly with seperate files
cat $data_path/PacBio/* > ../generated_data/all_pacbio.fasta

spades.py -1 $data_path/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz -2 $data_path/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz --pacbio ../generated_data/all_pacbio.fasta -o ${results_path}/spades_assembly


#2.4 Download the reference Genome and Annotation from NCBI
#https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009734005.1/

#Reference genome
wget -P $gen_data_path/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/734/005/GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.fna.gz
#Annotation
wget -P $gen_data_path/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/734/005/GCF_009734005.1_ASM973400v2/GCF_009734005.1_ASM973400v2_genomic.gff.gz


#2.5 Evaluation of Genome Assemblies

#Obtaining statistic metrics with Quast
#PacBIo Pilon Assembly
quast.py $pacbio_path/pac_pilon_assembly.fasta -r $gen_data_path/GCF_009734005.1_ASM973400v2_genomic.fna -g $gen_data_path/GCF_009734005.1_ASM973400v2_genomic.gff.gz -o $pacbio_path/pac_pilon_quast_eval
#Spades assembly
quast.py $spades_path/scaffolds.fasta -r $gen_data_path/GCF_009734005.1_ASM973400v2_genomic.fna -g $gen_data_path/GCF_009734005.1_ASM973400v2_genomic.gff.gz -o $spades_path/spades_quast_eval

#Checking if the assemblies contain well-known reserved genes with BUSCO
source $AUGUSTUS_CONFIG_COPY                                                                                                                                                                
busco -i $spades_path/scaffolds.fasta -m genome -l $BUSCO_LINEAGE_SETS/bacteria_odb10 -o ../results/spades_assembly/spades_busco_eval -f                                      
busco -i $pacbio_path/pac_pilon_assembly.fasta -m genome -l $BUSCO_LINEAGE_SETS/bacteria_odb10 -o ../results/pac_assembly/pac_pilon_busco_eval -f 

#Relative path is adviced as when full path used BUSCO creates locally a folder with all the given path folder names
#Busco hierarchical database from which i can select a specific lineage and it contains all the conserved genes of 
#the lineage and checks if those are found in my assembly. Basically, it assesses how not false the assembly is based 
#on the most conserved genes.


#Delete BUSCO folders
rm -r augustus_config
rm -r busco_downloads
