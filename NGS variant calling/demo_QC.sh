#### Quality control ############
##write a for loop to conduct fastQC 
#!/bin/bash
run_fastqc() {
  if [ -z "$1" ]; then
        echo "please provide file name。"
        return 1
  fi
  echo "run FastQC: $1"
  ###如果沒有符合的檔案則變成空字串
  shopt -s nullglob
  files = ($1)
  if [ ${#files[@]} -eq 0 ]; then
  	echo "no file existed"
   	return 1
   fi
   ##執行fastqc
    for file in "${files[@]}"; 
    do  fastqc "$file"
    done
}

###檢查是否有安裝fastqc
whereis fastqc > check.txt
if ! grep -q "/" check.txt ; then
  read -p "Do you want to install FastQC？(y/n): " answer
  if [[ $answer == "y" || $answer == "Y" ]]; then
    sudo apt install fastqc
else
  echo -e "cannot conduct fastqc"
  exit 1
   fi
fi
###執行run_fastqc()
echo "Enter file pattern (e.g., L*.fq.gz):"
read file_pattern
run_fastqc "$file_pattern"

### trimmomatic (sample: L_294 for example)
trimmomatic PE ~/RILs/01.RawData/L_294/L_294_DKDN220005958-1A_HJTWTDSX3_L2_1.fq.gz ~/RILs/01.RawData/L_294/L_294_DKDN220005958-1A_HJTWTDSX3_L2_2.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:30 TRAILING:30 MINLEN:70

### snp quality check: counting ts/tv ratio using vcf file 
cd ~/vcf_file
for i in {1..12}; do
  input_file="final_sample_plus_parent.chr${i}.vcf"
  input_file_2="final_sample_plus_parent.chr${i}.vcf.gz"
  output_file="chr${i}.txt"
  bgzip $input_file
  bcftools index $input_file_2 
  bcftools stats "$input_file_2" > "$output_file"
done

