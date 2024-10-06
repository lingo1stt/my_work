###### variant calling #################

cd ~/gatk-4.3.0.0
source /sharedata/MISC/setJava8.bsh

### remove duplicated reads
java -jar gatk-package-4.3.0.0-local.jar MarkDuplicatesSpark -I ~/samtools-1.9/L294_sorted.bam -O L294_marked.bam -M L294_duplicate.matrix.txt --tmp-dir ~/tmp

### use GNU parallel to conduct haplotypecaller with multiple thread
if ! command -v parallel &> /dev/null
then
	echo "GNU Parallel is installing..."
	sudo apt-get install -y parallel
fi

### do haplotype calling using HaplotypeCaller (Sample L_294 for example)
java -jar gatk-package-4.3.0.0-local.jar HaplotypeCaller -R ~/OsNB1.0.fa -I L294_marked.bam -O L294.gvcf -ERC GVCF

### do haplotype calling using HaplotypeCaller with multiple threads
mkdir "gvcf"
find . -type f -name "*.bam" > bam_name.txt
cat bam_name.txt | parallel -j 8 "java -jar gatk-package-4.3.0.0-local.jar HaplotypeCaller -R ~/OsNB1.0.fa -I {} -O ~/gatk-4.3.0.0/gvcf/{.}.gvcf -ERC GVCF"


### Combine all gvcf file of samples together
#list all the file in the gvcf_name.txt and combine them into one sample.gvcf
find ~/gatk-4.0.0.0/gvcf -type f -name "*.gvcf" > gvcf_name.txt
java -jar gatk-package-4.3.0.0-local.jar CombineGVCFs \
       	-R ~/OsNB1.0.fa \
	-V ~/gatk-4.0.0.0/gvcf_name.txt \
	-O sample.gvcf
### genotype Gvcf
java -jar gatk-package-4.3.0.0-local.jar GenotypeGVCFs -R ~/OsNB1.0.fa -V sample.gvcf -O sample_genotype.gvcf

### choose only SNP varaints
java -jar gatk-package-4.3.0.0-local.jar SelectVariants -R ~/OsNB1.0.fa -V sample_genotype.gvcf -select-type SNP -O sample_snp.vcf

### filter out low quality SNP using GATK hard-filtering criteria
java -jar gatk-package-4.3.0.0-local.jar VariantFiltration -R ~/OsNB1.0.fa -V sample_snp.vcf -filter "QD < 2.0" --filter-name "QD2"     -filter "QUAL < 30.0" --filter-name "QUAL30"     -filter "SOR > 3.0" --filter-name "SOR3"     -filter "FS > 50.0" --filter-name "FS60"     -filter "MQ < 40.0" --filter-name "MQ40"     -filter "MQRankSum < -12.5" --filter-name "MQRankSum-4"     -filter "ReadPosRankSum < -4.0" --filter-name "ReadPosRankSum-4" -O sample_filtered.vcf

