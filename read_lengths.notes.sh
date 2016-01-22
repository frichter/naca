
# FastQC command
module load fastqc
fastqc -f bam -o fastqc YOUR_BAMFILE

#Example : 
fastqc -f bam -o fastqc \
/sc/orga/projects/chdiTrios/NYGC_25WGS_Trio_30x_2/1-01094/1-01094-02.final.bam
