#!/bin/bash -l                                                                                                                                                                               

############# Setting up Uppmax #############                                                                                                                                               
#SBATCH -A uppmax2024-2-7                                                                                                                                                                    
#SBATCH -M snowy                                                                                                                                                                             
#SBATCH -p core                                                                                                                                                                       
#SBATCH -n 4
#SBATCH -t 07:00:00         
#SBATCH -J count_reads
#SBATCH --mail-type=ALL                                                                                                                                                                      
#SBATCH --mail-user aesimaiaki@gmail.com                                                                                                                                                     
#SBATCH --output=%x.%j.out                                                                                                                                                                   

#Load modules                                                                                                                                                                              
module load bwa/0.7.17
module load bwa-mem2
module load samtools/1.19
module load htseq/2.0.2

#Set variables                                                                                                                                                                              
rna_path="/home/agsi4287/Genome_Analysis/raw_data/transcriptomics_data"
results_path="/home/agsi4287/Genome_Analysis/results"
gen_data_path="../generated_data"
pacbio_path="$results_path/pac_assembly"
spades_path="$results_path/spades_assembly"


#Create a new folder in results
mkdir $spades_path/rna_mapping

#Create a txt file with the names of the files that will help to name the new files
#Here we are interested in the paired-end reads
ls $rna_path/RNA-Seq_BH |  grep  -E  "^trim_paired" |sed 's/_pass_*.*//' | uniq > $gen_data_path/bh_file_names
ls $rna_path/RNA-Seq_Serum | grep  -E "^trim_paired" |sed 's/_pass_*.*//' | uniq > $gen_data_path/serum_file_names


#RNA mapping

#Prepare reference
bwa index $spades_path/scaffolds.fasta


#BH samples
for file in $(cat $gen_data_path/bh_file_names); do
    # Alignment to your assembly
    bwa mem $spades_path/scaffolds.fasta \
            $rna_path/RNA-Seq_BH/${file}_pass_1.fastq.gz \
            $rna_path/RNA-Seq_BH/${file}_pass_2.fastq.gz | \
    samtools view -Sb - > $spades_path/rna_mapping/${file}_bh.bam

    # Sorting and indexing
    samtools sort $spades_path/rna_mapping/${file}_bh.bam \
                 -o $spades_path/rna_mapping/${file}_bh_sorted.bam
    samtools index $spades_path/rna_mapping/${file}_bh_sorted.bam

    # Delete the unsorted BAM file                                                                                                                    
    rm $spades_path/rna_mapping/${file}_bh.bam
done

#Serum
for file in $(cat $gen_data_path/serum_file_names); do
    # Alignment
    bwa mem $spades_path/scaffolds.fasta \
            $rna_path/RNA-Seq_Serum/${file}_pass_1.fastq.gz \
            $rna_path/RNA-Seq_Serum/${file}_pass_2.fastq.gz | \
    samtools view -Sb - > $spades_path/rna_mapping/${file}_serum.bam

    # Sorting and indexing
    samtools sort $spades_path/rna_mapping/${file}_serum.bam \
                 -o $spades_path/rna_mapping/${file}_serum_sorted.bam
    samtools index $spades_path/rna_mapping/${file}_serum_sorted.bam

    # Delete the unsorted BAM file
    rm $spades_path/rna_mapping/${file}_serum.bam
done


#Count reads
#HTSeq requires a .gff file that doesn't contain sequences. This file was created manually after annotation
#as prokka incorparates the sequences in the produced .gff file.


mapping_path='/home/agsi4287/Genome_Analysis/results/spades_assembly/rna_mapping'
assembly_annot_path='/home/agsi4287/Genome_Analysis/results/spades_assembly/prokka_annot'

mkdir $spades_path/count_reads
counts_path=$spades_path/count_reads

for sample in $(ls $mapping_path| grep -E ".bam$"| sed "s/_sorted.bam//")
do 

    htseq-count --format=bam --stranded=no --type=CDS --order=pos --idattr=ID $mapping_path/${sample}_sorted.bam $assembly_annot_path/E_faecium_noseq.gff > $counts_path/${sample}_htseq_counts.txt

done



#Calculate the statistics of each .bam file
#Output on terminal 

for f in $(ls $mapping_path| grep -E ".bam$")
do

    samtools flagstat $f

done
