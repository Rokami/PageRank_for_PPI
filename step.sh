#!bin/bash

#get disease_id file
cut -f 1 MimMiner_Exp_AC_T_TXCS_basedonACMESH_filt_RW.mat > disease_id.txt
#get phenotype_genotype disease interactions file
perl -ane 'next if /^#/;@t=split /\t/,$_;($s)=$t[0]=~/,\s*([0-9]{6,})/g;print "$s\t$t[1]\n";' morbidmap.txt  | perl -ane 'BEGIN{open IN,"disease_id.txt" or die $!;@t=<IN>;chomp @t;foreach(@t){$h{$_}="";};}chomp;@t=split /\t/,$_;$h{$t[0]}.="$t[1]," if exists $h{$t[0]};END{foreach(sort{$a<=>$b}keys %h){$h{$_}=~s/,$//;print "$_\t$h{$_}\n"}}' > phenotype_genotype.txt
#get protein_protein interactions file
perl -ane '@t=split /\t/,$_;print "$t[1]\t$t[4]\n"' BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt > protein_protein_interactions.txt
#get protein name match to HPRD ID list
perl -ane '@t=split /\t/,$_;$h{$t[0]}=$t[1];$h{$t[3]}=$t[4];END{foreach(sort{$a cmp $b}keys %h){print "$_\t$h{$_}\n"}}' BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt > all_protein_HPRD_ID.txt
#get relationship of OMIM disease ID with HPRD protein ID
perl -ane 'BEGIN{open IN,"all_protein_HPRD_ID.txt" or die $!;@t=<IN>;chomp @t;%h=map{(split /\t/)[0]=>(split /\t/)[1]}@t;}@t=split /\t/,$_;$t[1]=~s/\s+//g;@s=split /,/,$t[1];print "$t[0]\t";undef @d;foreach(@s){next if !exists $h{$_};push @d,$h{$_}};print join ",",@d;print "\n";' phenotype_genotype.txt > omim_HPRD.txt

#main program
time python pagerank.py > output.txt 
#or time python pagerank.py <disease_id> > <output_file_name>

#Number of HPRD gene symbol exsits in OMIM disease database: 1442
#Total gene number in HPRD network: 9617
#interaction number in HPRD network: 39240
#Total gene symbol in OMIM disease phenotype database: 6027