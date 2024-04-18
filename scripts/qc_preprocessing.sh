########## 1. Quality Control & Trimming / Preprocessing ############                                                                                                                        #The present script is hard coded for the specific Illumina files provided to us for the project 

#Load necessary modules
module load FastQC/0.11.9
module load trimmomatic/0.39

#Define variables
data_path="/home/agsi4287/Genome_Analysis/raw_data/genomics_data"
results_path="/home/agsi4287/Genome_Analysis/results"


#1.1 Quality Control with FastQC                                                                                                                                                                
cd $data_path/Illumina

for i in $(ls)
do
        echo $i
        fastqc -o $results_path $i                                                                                                                                                           
done

cd ../../../results


#1.2 Reads trimming with trimmomatic                
                                                                                                                              
adapters_path="/sw/bioinfo/trimmomatic/0.39/snowy/adapters"

java -jar /sw/bioinfo/trimmomatic/0.39/snowy/trimmomatic-0.39.jar PE -threads 2 ${data_path}/Illumina/E745-1.L500_SZAXPI015146-56_1_clean.fq.gz ${data_path}/Illumina/E745-1.L500_SZAXPI015146-56_2_clean.fq.gz E745-1.L500_1_paired.fq.gz E745-1.L500_1_unpaired.fq.gz E745-1.L500_2_paired.fq.gz E745-1.L500_2_unpaired.fq.gz ILLUMINACLIP:$adapters_path/TruSeq3-PE.fa:2:30:10:2:TrueLEADING:20 TRAILING:20                                                                                                                                                                        

#Notes: ILLUMINACLIP the default setting were selected                                                                                                                                      
#LEADING & TRAILING define the min quality allows for the astart and end respectivelly bases.                                       
#The adapters are sequences used by the Illumina file                                                                                                                                        
#adapters file can be found in the directory whwre trimmomatic is installed                                                                                                                  
#module help trimmomatic                                                                                                                                                                     
#echo $TRIMMOMATIC_ROOT                                                                                                                                                                      
#ls /sw/bioinfo/trimmomatic/0.39/snowy/adapters/    


#1.3 Quality Control of trimmed reads with FastQC 

fastqc E745-1.L500_1_paired.fq.gz
fastqc E745-1.L500_2_paired.fq.gz

cd ../scripts
