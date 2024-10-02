###### variant calling #################

cd ~/gatk-4.3.0.0
source /sharedata/MISC/setJava8.bsh

### remove duplicated reads
java -jar gatk-package-4.3.0.0-local.jar MarkDuplicatesSpark -I ~/samtools-1.9/L294_sorted.bam -O L294_marked.bam -M L294_duplicate.matrix.txt --tmp-dir ~/tmp

### do haplotype calling using HaplotypeCaller (Sample L_294 for example)
java -jar gatk-package-4.3.0.0-local.jar HaplotypeCaller -R ~/OsNB1.0.fa -I L294_marked.bam -O L294.gvcf -ERC GVCF

### Combine all gvcf file of samples together
java -jar gatk-package-4.3.0.0-local.jar CombineGVCFs -R ~/OsNB1.0.fa -V L190.gvcf -V L210.gvcf -V L214.gvcf -V L238.gvcf -V L244.gvcf -V L245.gvcf -V L265.gvcf -V L282.gvcf -V L291.gvcf -V L294.gvcf -V L295.gvcf -V L317.gvcf -V S167.gvcf -V S172.gvcf -V S188.gvcf -V S200.gvcf -V S224.gvcf -V S242.gvcf -V S256.gvcf -V S267.gvcf -V S273.gvcf -V S281.gvcf -V S283.gvcf -V S284.gvcf -O sample.gvcf

### genotype Gvcf
java -jar gatk-package-4.3.0.0-local.jar GenotypeGVCFs -R ~/OsNB1.0.fa -V sample.gvcf -O sample_genotype.gvcf

### choose only SNP varaints
java -jar gatk-package-4.3.0.0-local.jar SelectVariants -R ~/OsNB1.0.fa -V sample_genotype.gvcf -select-type SNP -O sample_snp.vcf

### filter out low quality SNP using GATK hard-filtering criteria
java -jar gatk-package-4.3.0.0-local.jar VariantFiltration -R ~/OsNB1.0.fa -V sample_snp.vcf -filter "QD < 2.0" --filter-name "QD2"     -filter "QUAL < 30.0" --filter-name "QUAL30"     -filter "SOR > 3.0" --filter-name "SOR3"     -filter "FS > 50.0" --filter-name "FS60"     -filter "MQ < 40.0" --filter-name "MQ40"     -filter "MQRankSum < -12.5" --filter-name "MQRankSum-4"     -filter "ReadPosRankSum < -4.0" --filter-name "ReadPosRankSum-4" -O sample_filtered.vcf

