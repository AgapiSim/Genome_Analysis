#!/bin/bash -l

############# Setting up Uppmax ############# 

#SBATCH -A uppmax2024-2-7
#SBATCH -M snowy
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH -J g_assembly_eval_aes
#SBATCH --mail-type=ALL
#SBATCH --mail-user aesimaiaki@gmail.com
#SBATCH --output=%x.%j.out

# Load modules - Each script calls the necessary modules
module load bioinfo-tools

#Set variables & directories
data_path="/home/agsi4287/Genome_Analysis/raw_data/genomics_data"
results_path="/home/agsi4287/Genome_Analysis/results"




########## 1.Quality Control & Trimming / Preprocessing ############

bash ./qc_preprocessing.sh
   
#1.1 Quality Control of Illumina reads
#1.2 Trimming of reads
#1.3 Quality Control on the trimmed reads



 
######### 2.Genome assembly ############

bash ./genome_assembly.sh

#2.1 Genome assembly with PacBio reads with Flye
#2.2 Genome assembly with PacBio & Illuma reads with Spades
#Pilon improves PacBio assembly utilizing the Illumina short reads
#Spades assembles short Illumina reads and then uses the PacBio reads to scaffold and conjuct reads together
#2.3 Download the reference genome and annotation of Enterococcus faecium
#2.4 Genome Assembly Evaluation




############################ 3. Annotation  ##################################                                                                               
#Annotation of the assembly with Prokka

bash ./annotation.sh




########################### 4. Syneteny Comparison / DNA mapping to refence ########################

#The assembly is aligned to the reference genome to check the similarities of the two sequencies
#and assess the quality of our assembly. Two approaches were tested 1. Blastn & Artemis visualization
# and 2. Mummer. Mummer was used as the main approach for synteny comparison.



#1. BLASTn approach and visualization with Artemis/act (manual step locally)

#DNA alignement of the assembly to the reference genome
# Decompress the gzipped database file                                                                                                                                                     
#gzip -d ./../generated_data/GCF_009734005.1_ASM973400v2_genomic.fna.gz

# Use makeblastdb to create the BLAST database                                                                                                                                               
#makeblastdb -in ./../generated_data/GCF_009734005.1_ASM973400v2_genomic.fna -dbtype nucl -out ./../generated_data/GCF_009734005.1_ASM973400v2_genomic.db

# Run blastn using the created database                                                                                                                                                      
#blastn -query $results_path/spades_assembly/scaffolds.fasta -db ./../generated_data/GCF_009734005.1_ASM973400v2_genomic.db -out $results_path/spades_assembly/E_faecium_spades.blastn

#rm ./../generated_data/GCF_009734005.1_ASM973400v2_genomic.fna  


#2. Mummer approach

bash ./mummer_align.sh


########################### 5. RNA mapping & Count reads ##########################################

#In this part the RNA Illumina reads obtained in BH and Serum are mapped against the assembly
#and the reads per region are count

bash ./rna_map_count.sh
