$option = 0;
while ($option =~ /[0-8]/){
                if($option == 0){&menu();}
                if($option == 1){&input();}
                if($option == 2){& status();}
                if($option == 3){&remove();}
                if($option == 4){&blast();}
                if($option == 5){&reverse_complement();}
                if($option == 6){&amino();}
                if($option == 7){&GC_content();}
                if($option == 8){die();}
        }


sub menu{
        print "\nMain Menu\n";
        print "************************************\n";
        print "1) Input fasta file\n";
        print "2) Check status\n";
        print "3) Remove sequences\n";
        print "4) Local Blast\n";
        print "5) Reverse complement sequence\n";
        print "6) Convert sequences to amino acid\n";
        print "7) GC content\n";
	print "8) Exit\n";
        print "************************************\n";
        print "choose your option(1-8):";
        $option = <STDIN>;
        chomp $option;
        unless ($option =~ /^\d+$/ && $option >=1 && $option <=8){
	        print "Please select a number between 1~8!\n\n";
        	$option = 0;
        }
	return $option;
}
#假設input的fasta file以hash形式儲存:假設叫做%data = qw(SEQ1 seq1 SEQ2 seq2 SEQ3 seq3)

## 1) Input
sub input{
        ##現在有個問題是只要一進來上次的輸入數據就會被覆蓋
        $count = 0;
        while ($count+1){
                $num = $count + 1;
                print "Input sequence $num fasta file (F to finish): ";
                $seq[$count] = <STDIN>;
                chomp $seq[$count];
                if ($seq[$count] eq "F"){splice (@seq, -1); last;}
                unless ($seq[$count] ne ""){next;}
                $count++;
        }
        if ($count == 1){print "Warning: There is only one sequence inputted, [Blast] would unable to be use.\n"}
        
	for ($x = 0; $x < $count; $x++) {
		$name{$x} = $data{$x} = "";
                open FILE, $seq[$x];
		while(<FILE>) {
			chomp;
                        if (/>(\S+)/) {$name{$x} = $1; next;}
			$data{$x} .= $_;
		}
		close FILE;
                unless ($data{$x} ne ""| $data{$x} =~ /[ATCG]/){
                        print "There are some problem in [$seq[$x]], please check and try again.\n";
                        print "Warning: Only accept nucleotide sequence fasta file.\n\n";
                        $option = 1; 
                        return;
                }
	}
        $option = 0;
        return $input;
}

## 2) Status
sub status {
         unless (keys %data){
                print "Please input data first.\n";
                &input(); }
         print "Here are the sequences you have inputted:\n";
         for $key (sort {$a <=> $b} keys %data){
         my $num = $key +1;
        $length = length $data{$key};
        print " SEQ$num : $data{$key}\n";
        print "length : $length\n"; }
        $option = 0;
        return;}

## 3) Remove
sub remove {
        unless (keys %data){
                print "You don't have imputed any sequence\n";
                $option = 0;
                return;}
        print "Are you sure you want to remove the sequence? (yes or no)\n";
        $answer = <STDIN>;
        chomp $answer;
        if ($answer eq "yes"){
        for $key (keys %data){
                delete $data{$key};}}
        print "Your sequences status are as follow\n";
        if (keys %data){
                for $key (sort {$a <=> $b} keys %data){
                my $num = $key +1;
                print "SEQ$num :  $data{$key}\n";}}
                else {print "There's no sequence now!\n";}
        $option = 0;
        return;}



## 4) Blast
sub blast {
        unless (keys %data){
                print "Please input data first.\n";
                &input();
        }
        if ($count == 1){
                print "You cannot use Blast with only one sequence. Try other fuctions.\n";
                $option = 0;
                return;
        }
        print "\nHere are the sequences you have inputted:\n";
        for $key (sort {$a <=> $b} keys %data){
                my $num = $key +1;
                print "$num  $name{$key}\n";
        }
        print "Pick the query sequence(select by number):";
        chomp($input_Q = <STDIN>);
        unless ($input_Q =~ /^\d+$/ && $input_Q <= $count){print "No such sequence exist. Try again.\n"; &blast();}
        print "Pick the reference sequence(select by number):";
        chomp($input_R = <STDIN>);
        unless ($input_R =~ /^\d+$/ && $input_R <= $count){print "No such sequence exist. Try again.\n"; &blast();}

        while ($input_Q == $input_R){
                print "You sure you want to blast the sequence with itself? (Yes or No)";
                chomp($choice = <STDIN>);
                if ($choice eq "No"){print "Ok, try again.\n"; &blast();}
                elsif ($choice eq "Yes"){print "Ok, go ahead.\n"; last;}
                else {print "Please type 'Yes' or 'No'.\n"; next;}
                }

        $query = $ref = "";
        $query = $seq[$input_Q-1];
        $ref = $seq[$input_R-1];
        print "Query file: $query\n";
        print "Reference file: $ref\n";

        unless (-d './blast_output_dir') {
                system("mkdir blast_output_dir");
        }
        unless (-d './blast_db_dir') {
                system("mkdir blast_db_dir");
        }

        #/home/youylin/PERL/NCBI_BLAST/ncbi-blast-2.13.0/bin
        system ("/home/youylin/PERL/NCBI_BLAST/ncbi-blast-2.13.0/bin/makeblastdb -in $ref -dbtype nucl -out ./blast_db_dir/$ref");
        $output = "blast_SEQ$input_Q"."SEQ$input_R.txt";
        system ("/home/youylin/PERL/NCBI_BLAST/ncbi-blast-2.13.0/bin/blastn -db ./blast_db_dir/$ref -query $query -outfmt 6 -out ./blast_output_dir/$output");
        
        $option = 0;
        return;
}


## 5) Reverse compliment
sub reverse_complement {
        unless (keys %data){
                print "Please input data first.\n";
                &input();
        }
        for $key (sort {$a <=> $b} keys %data){
                $seq = $data{$key};
                print "$name{$key}\n";
                print "Sequence Input:     $seq\n";
                $seq =~ tr/ATCGBVDHKMRYSW/TAGCVBHDMKYRSW/; 
                $rcseq = reverse $seq;
                print "Reverse Compliment: $rcseq\n";
        }
        $option = 0;
        return;
}



## 6) amino acid
sub amino {
        unless (keys %data){
                print "Please input data first.\n";
                &input();
        }
        #%data = qw(seq1 ATTTATCCTTCTG seq2 ATTTGCATTATCATAGCT seq3 TTCCAGAATCAT);
        while(1){
                print"which ORF to convert (-3,-2,-1,1,2,3):\n";
                $input_2 = <STDIN>;
                chomp $input_2;
                if ( $input_2 == 1 || $input_2 == 2 || $input_2 ==3){
                        for $keys (sort {$a cmp $b}keys %data){
                                $seq = $data{$keys};
                                print"$name{$keys}:\n";

                                for($x = ($input_2 - 1); $x + 3 <= length$seq; $x +=3 ){
                                        $codon = substr($seq,$x,3);
                                        $codon =~ s/ATT|ATC|ATA/I/;
                                        $codon =~ s/CTT|CTC|CTA|CTG|TTA|TTG/L/;
                                        $codon =~ s/GTT|GTC|GTA|GTG/V/;
                                        $codon =~ s/TTT|TTC/F/;
                                        $codon =~ s/ATG/M/;
                                        $codon =~ s/TGT|TGC/C/;
                                        $codon =~ s/GCT|GCC|GCA|GACG/A/;
                                        $codon =~ s/GGT|GGC|GGA|GGG/G/;
                                        $codon =~ s/CCT|CCC|CCA|CCG/P/;
                                        $codon =~ s/ACT|ACC|ACA|ACG/T/;
                                        $codon =~ s/TCT|TCC|TCA|TCG|AGT|AGC/S/;
                                        $codon =~ s/TAT|TAC/Y/;
                                        $codon =~ s/TGG/W/;
                                        $codon =~ s/CAA|CAG/Q/;
                                        $codon =~ s/AAT|AAC/N/;
                                        $codon =~ s/CAT|CAC/H/;
                                        $codon =~ s/GAA|GAG/E/;
                                        $codon =~ s/GAT|GAC/D/;
                                        $codon =~ s/AAA|AAG/K/;
                                        $codon =~ s/CGT|CGC|CGA|CGG|AGA|AGG/R/;
                                        $codon =~ s/TAA|TAG|TGA/\*/;
                                        print"$codon";

                                }
                                print"\n";
                        }
                        last;
                }elsif($input_2 == -1 || $input_2 == -2 || $input_2 == -3){
                        for $keys (sort {$a cmp $b}keys %data){
                                $seq = $data{$keys};
                                $seq =~tr/ATCG/TAGC/;
                                $rev_seq = reverse $seq;
                                print"$keys:\n";

                                for($y = (abs($input_2) - 1); $y + 3 <= length $rev_seq; $y +=3 ){
                                        $codon = substr($rev_seq,$y,3);
                                        $codon =~ s/ATT|ATC|ATA/I/;
                                        $codon =~ s/CTT|CTC|CTA|CTG|TTA|TTG/L/;
                                        $codon =~ s/GTT|GTC|GTA|GTG/V/;
                                        $codon =~ s/TTT|TTC/F/;
                                        $codon =~ s/ATG/M/;
                                        $codon =~ s/TGT|TGC/C/;
                                        $codon =~ s/GCT|GCC|GCA|GACG/A/;
                                        $codon =~ s/GGT|GGC|GGA|GGG/G/;
                                        $codon =~ s/CCT|CCC|CCA|CCG/P/;
                                        $codon =~ s/ACT|ACC|ACA|ACG/T/;
                                        $codon =~ s/TCT|TCC|TCA|TCG|AGT|AGC/S/;
                                        $codon =~ s/TAT|TAC/Y/;
                                        $codon =~ s/TGG/W/;
                                        $codon =~ s/CAA|CAG/Q/;
                                        $codon =~ s/AAT|AAC/N/;
                                        $codon =~ s/CAT|CAC/H/;
                                        $codon =~ s/GAA|GAG/E/;
                                        $codon =~ s/GAT|GAC/D/;
                                        $codon =~ s/AAA|AAG/K/;
                                        $codon =~ s/CGT|CGC|CGA|CGG|AGA|AGG/R/;
                                        $codon =~ s/TAA|TAG|TGA/\*/;
                                        print"$codon";

                                }
                                print"\n";
                        }
                        last;
                }else {print"you idiot\n";}

        }
        $option = 0;
        return;
}


## 7) GC content
sub GC_content {
        unless (keys %data){
                print "Please input data first.\n";
                &input();
        }
        #%data = qw(seq1 ATCGG seq2 ATTTG seq3 TTCCA);
        $length = 0;
        for $keys (sort {$a cmp $b}keys %data){
                $num = $data{$keys} =~ tr/GC//;
                $length = length$data{$keys};
                $gc_cont{$keys} = ($num / $length) * 100;
                print"$name{$keys}\t $gc_cont{$keys}% \n";
        }
        $option = 0;
        return;
}
