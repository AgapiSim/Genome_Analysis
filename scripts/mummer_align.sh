#Set variables                                                                                                                                                                                
data_path="/home/agsi4287/Genome_Analysis/raw_data/genomics_data"
results_path="/home/agsi4287/Genome_Analysis/results"
gen_data_path="../generated_data"
pacbio_path="$results_path/pac_assembly"
spades_path="$results_path/spades_assembly"


#Run allignment
nucmer -p $spades_path/mummer_scaf_align $gen_data_path/GCF_009734005.1_ASM973400v2_genomic.fna $spades_path/scaffolds.fasta

#Find 1-1 matches
delta-filter -1 $spades_path/mummer__scaf_align.delta

#Plot and extract graph
mummerplot --png $spades_path/mummer_scaf_align.delta -p $spades_path/mummer_spades_scaf


#Additional step for comparing the difference between the contigs.fasta and assembly.fasta files that spades produces
#Minor difference were seen in favor of scaffolds.fasta hence it was selected as the assembly file for downstream analysis
 
#Contigs Spades
#Run alignment
#nucmer -p $spades_path/mummer_contigs_align $gen_data_path/GCF_009734005.1_ASM973400v2_genomic.fna $spades_path/contigs.fasta

#Find 1-1 matches                                                                                                                                                                            
#delta-filter -1 $spades_path/mummer_contigs_align.delta

#Plot                                                                                                                                                                                        
#mummerplot --png $spades_path/mummer_contigs_align.delta -p $spades_path/mummer_spades_contigs
