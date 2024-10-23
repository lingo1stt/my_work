#### Quality control ############
##write a for loop to conduct fastQC and trimmomatic at the same time
run_fastqc() {
  if [ -z "$1" ]; then
        echo "please provide file name。"
        return 1
  fi
  ###執行
  echo "run FastQC: $1"
    for file in L*.fq.gz
    do
	    if [ -e "$file"]; then
		    fastqc "$file"
	    else 
		    echo "NO files existed"
	    fi
    done
 
}


whereis fastqc > check.txt
if ! grep -q "/" check.txt ; then
  read -p "Do you want to install FastQC？(y/n): " answer
  if [[ $answer == "y" || $answer == "Y" ]]; then
    sudo apt install fastqc
    echo -e "Enter file name with regular expression (ex:*.fq.gz):"
    read answer_2
    run_fastqc "$answer_2"
	fi
else
  echo -e "cannot conduct fastqc, jump to next step \n"
fi


for file in L*.fq.gz
do
	if [ -e "$file"]; then
		echo "file existed"
		fastqc "$file"
	else 
		echo "NO files existed"
	fi
done

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

