################### 3. Annotation #######################

#3.Annotation

#Set variables                                                                                                                                                                              
data_path="/home/agsi4287/Genome_Analysis/raw_data/genomics_data"
results_path="/home/agsi4287/Genome_Analysis/results"
gen_data_path="../generated_data"
pacbio_path="$results_path/pac_assembly"
spades_path="$results_path/spades_assembly"


#Load the module locally
module load prokka/1.45-5b58020
module load bwa/0.7.17
module load samtools/1.19
module load blast/2.15.0+


#PacBio Pilon Assembly (additional step for comparison)
#Spades assembly was selected for downstream analysis

#prokka \
 #       --outdir ../results/pac_assembly/prokka_annot --force \
  #      --prefix E_faecium --addgenes  \
   #     --genus Enterococcus --species faecium --strain E745 \
    #    --kingdom Bacteria  \
     #   $pacbio_path/assembly.fasta


#Spades assembly
prokka \
        --outdir $spades_path/prokka_annot --force \
        --prefix E_faecium --addgenes  \
        --genus Enterococcus --species faecium --strain E745 \
        --kingdom Bacteria  \
        $spades_path/scaffolds.fasta


#Remove sequence from .gff file 
sed '/^##FASTA/Q' ../results/pac_assembly/prokka_annot/E_faecium.gff > ../results/pac_assembly/prokka_annot/E_faecium _noseq.gff

