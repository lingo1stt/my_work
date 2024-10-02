############ alignment demo ######################\
### enter the directory having reference genome
cd ~/bwa-0.7.15

### using BWA-MEM to conduct alignment 
bwa mem -t 6 -R '@RG\tID:L_294\tSM:L_294\tLB:294lib\tPL:ILLUMINA' bwa_rice_index ~/RILs/01.RawData/L_294/L_294_DKDN220005958-1A_HJTWTDSX3_L2_1.fq.gz ~/RILs/01.RawData/L_294/L_294_DKDN220005958-1A_HJTWTDSX3_L2_2.fq.gz -o L294_output.sam

### using samtools to transform output sam file into bam file
cd ~/samtools-1.9
samtools view -hb ~/bwa-0.7.15/L294_output.sam -o L294.bam 
samtools sort L294.bam -o L294_sorted.bam
samtools index L294_sorted.bam

### using gatk to validate the avialibility of bam file
cd ~/gatk-4.3.0.0
source /sharedata/MISC/setJava8.bsh
java -jar gatk-package-4.3.0.0-local.jar ValidateSamFile -I ~/samtools-1.9/L294_sorted.bam -M SUMMARY

