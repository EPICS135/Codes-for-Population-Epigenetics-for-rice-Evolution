#usr/bin/perl
#The program can be used for mapping BS-seq data with PE-reads and generate reports for QC, mapping and methylation level.

#usage:perl script.pl $ARGV[0]

#dependents
#bowtie2-2.2.9
#bismark
#fastq

#the pathway for index of genome for bismark
$path_to_bismark_index="../genome";# modification
$coverage_cutoff=3;# modification

#raw reads
$R1=$ARGV[0].".R1.fastq.gz";
$R2=$ARGV[0].".R2.fastq.gz";
#filtered reads
$CR1=$ARGV[0].".R1.clean.fastq.gz";
$CR2=$ARGV[0].".R2.clean.fastq.gz";
#unpaired reads
$unCR1=$ARGV[0].".R1.unpaired.fastq.gz";
$unCR2=$ARGV[0].".R2.unpaired.fastq.gz";
#QC with fastp
system("fastp --detect_adapter_for_pe -i $R1 -I $R2 -o $CR1 --unpaired1 $unCR1 -O $CR2 --unpaired2 $unCR2 -q 20 -w 2 -3 -M 30 -Y 20 -j $ARGV[0].json -h $ARGV[0].html");
#mapping with bismark
system("bismark -q -p 2 --score_min L,0,-0.6 -X 1000 $path_to_bismark_index -1 $CR1 -2 $CR2");
unlink $CR1;unlink $CR2;
#-N: mismatch£¬-p: thread
#remove PCR amplification
system("deduplicate_bismark -p --bam $ARGV[0].R1.clean_bismark_bt2_pe.bam");
#extractor£¬--multicore: thread
system("bismark_methylation_extractor -p --no_overlap --comprehensive --multicore 8 $ARGV[0].R1.clean_bismark_bt2_pe.deduplicated.bam");
unlink $ARGV[0].".R1.clean_bismark_bt2_pe.deduplicated.bam";
#merge results of extractor
$file1="CpG_context_".$ARGV[0].".R1.clean_bismark_bt2_pe.deduplicated.txt";
$file2="CHG_context_".$ARGV[0].".R1.clean_bismark_bt2_pe.deduplicated.txt";
$file3="CHH_context_".$ARGV[0].".R1.clean_bismark_bt2_pe.deduplicated.txt";
$extractor=$ARGV[0]."_all.txt";
system("cp $file1 $extractor");
system("grep -v Bismark $file2 >> $extractor");
system("grep -v Bismark $file3 >> $extractor");
unlink $file1;unlink $file2;unlink $file3;
#generate bedgraph
system("bismark2bedGraph -o $ARGV[0]_all_bedgraph.txt -CX --ample_memory $extractor");
unlink $extractor;
#generate methylation report
$Creport=$ARGV[0]."_all_genome_C_report.txt";
system("coverage2cytosine -o $Creport -genome_folder $path_to_bismark_index/ -CX $ARGV[0]_all_bedgraph.txt.gz.bismark.cov.gz");
system("mv $Creport.CX_report.txt $Creport");
#the summary of number of raw and clean reads
open(F1,"$ARGV[0].json") or die "$!";
$i=0;
while(my $line=<F1>){
	chomp $line;
	my @a=split/\:|\,/,$line;
	$i+=1;
	if($i==4){
		$total_reads=$a[1]/2;
	}else{
		if($i==15){
			$clean_reads=$a[1]/2;
		}
	}
}
close F1;
#the number of mapped reads
open(F1,"$ARGV[0].R1.clean_bismark_bt2_PE_report.txt") or die "$!";
$i=0;
while(my $line=<F1>){
	chomp $line;
	$line=~s/\s+|\t//g;
	my @a=split/\:/,$line;
	$i+=1;
	if($i==8){
		$map_read=$a[1];
	}else{
		if($i==9){
			$map_rate=$a[1];
		}
	}
}
close F1;
#the number of read after deduplication
open(F1,"$ARGV[0].R1.clean_bismark_bt2_pe.deduplication_report.txt") or die "$!";
$i=0;
while(my $line=<F1>){
	chomp $line;
	$line=~s/\s+//g;
	my @a=split/\:|\(/,$line;
	$i+=1;
	if($i==6){
		my @b=split/\s+/,$a[1];
		$deduplicated=$b[0];
	}
}
close F1;
#the reads with more than coverage_cutoff
open(F1,$Creport) or die "$!";
$Creport_select=$ARGV[0]."_all_genome_C_report_selected.txt";
open(O,">$Creport_select") or die "$!";
$mC_sequenced_in_chrC1=0;$C_sequenced_in_chrC1=0; 
$mC_sequenced_in_chrC2=0;$C_sequenced_in_chrC2=0;
$mC_sequenced_in_chrC3=0;$C_sequenced_in_chrC3=0;
while(my $line=<F1>){
	chomp $line;
	my @a=split/\t/,$line;
	#cytoplasmic genome for evaluating conversion error
	if (($arr[0] eq "ChrC") or ($arr[0] eq "ChrM")){ 
    if($arr[0] eq "ChrC"){$mC_sequenced_in_chrC2+=$arr[3];$C_sequenced_in_chrC2=$C_sequenced_in_chrC2+$arr[3]+$arr[4];}
	if($arr[0] eq "ChrM"){$mC_sequenced_in_chrC3+=$arr[3];$C_sequenced_in_chrC3=$C_sequenced_in_chrC3+$arr[3]+$arr[4];}
	}else{
		$Csite+=1;
		if($a[5] eq "CG"){
			$CGsite+=1;
		}else{
			if($a[5] eq "CHG"){
				$CHGsite+=1;
			}else{
				if($a[5] eq "CHH"){
					$CHHsite+=1;
				}
			}
		}

		if($a[3]+$a[4]>=$coverage_cutoff){#coverage
			print O $line,"\n";
			$cCsite+=1;$mC+=$a[3];$unmC+=$a[4];
			if($a[5] eq "CG"){
				$cCGsite+=1;$mCG+=$a[3];$unmCG+=$a[4];
			}else{
				if($a[5] eq "CHG"){
					$cCHGsite+=1;$mCHG+=$a[3];$unmCHG+=$a[4];
				}else{
					if($a[5] eq "CHH"){
						$cCHHsite+=1;$mCHH+=$a[3];$unmCHH+=$a[4];
					}
				}
			}
		}
	}
}
close F1;close O;
#the result for QC and mapping
open(O,">$ARGV[0]_qc_map_mCcover_mClevel_statistic_report.txt") or die "$!";
print O $ARGV[0],"\t","Total reads","\t",$total_reads,"\t",
	"Clean reads","\t",$clean_reads,"\t",
	"Clean reads%","\t",$clean_reads/$total_reads,"\t",
	"Unique map reads","\t",$map_read,"\t",
	"Unique map reads%","\t",$map_rate,"\t",
	"Reads after deduplication","\t",$deduplicated,"\t",

	#coverage of cytosine
	"total CHH sites","\t",$CHHsite,"\t",
	"total CHG sites","\t",$CHGsite,"\t",
	"total CG sites","\t",$CGsite,"\t",
	"total C sites","\t",$Csite,"\t",
	
	"cover CHH sites","\t",$cCHHsite,"\t",
	"cover CHG sites","\t",$cCHGsite,"\t",
	"cover CG sites","\t",$cCGsite,"\t",
	"cover C sites","\t",$cCsite,"\t",
	
	"cover CHH sites%","\t",$cCHHsite/$CHHsite,"\t",
	"cover CHG sites%","\t",$cCHGsite/$CHGsite,"\t",
	"cover CG sites%","\t",$cCGsite/$CGsite,"\t",
	"cover C sites%","\t",$cCsite/$Csite,"\t",
	#the methylation level
	"mCHH level","\t",$mCHH/($mCHH+$unmCHH),"\t",
	"mCHG level","\t",$mCHG/($mCHG+$unmCHG),"\t",
	"mCG level","\t",$mCG/($mCG+$unmCG),"\t",
	"mC level","\t",$mC/($mC+$unmC),"\t",
	#error
"error rate ChrC:","	",$mC_sequenced_in_chrC2/$C_sequenced_in_chrC2,"\t",
"error rate ChrM:","	",$mC_sequenced_in_chrC3/$C_sequenced_in_chrC3,"\n";

close O;

@name=split/\_/,$ARGV[0];
open(F1,$ARGV[0]_PE_pe_all_genome_C_report.txt) or die "$!";
open(OUT1,">$ARGV[0]_CG.txt") or die "$!";
open(OUT2,">$ARGV[0]_CHG.txt") or die "$!";
open(OUT3,">$ARGV[0]_CHH.txt") or die "$!";

while (my $line=<F1>){
	chomp $line;
	my @arr=split/\t/,$line;
	if($arr[5] eq "CG"){if(($arr[3]+$arr[4])>=3){$level=$arr[3]/($arr[3]+$arr[4]);print OUT1 $arr[0],"-",$arr[1],"\t",$level,"\n";}else{print OUT1 $arr[0],"-",$arr[1],"\t","NA","\n";}}
    if($arr[5] eq "CHG"){if(($arr[3]+$arr[4])>=3){$level=$arr[3]/($arr[3]+$arr[4]);print OUT2 $arr[0],"-",$arr[1],"\t",$level,"\n";}else{print OUT2 $arr[0],"-",$arr[1],"\t","NA","\n";}}
    if($arr[5] eq "CHH"){if(($arr[3]+$arr[4])>=3){$level=$arr[3]/($arr[3]+$arr[4]);print OUT3 $arr[0],"-",$arr[1],"\t",$level,"\n";}else{print OUT3 $arr[0],"-",$arr[1],"\t","NA","\n";}}
}
close F1;close OUT1;close OUT2;close OUT3;

system("sort -t "-" -k 1,1 -k 2n,1 $ARGV[0]_CG.txt >$ARGV[0]_CG-sort.txt");

open(F1,$ARGV[0]_CG-sort.txt) or die "$!";
open(OUT4,">$ARGV[0]_CG-sort-for_paste.txt") or die "$!";
print OUT4 $name[0],"\n",$name[0],"\n";
while (my $line=<F1>){
	chomp $line;
	my @arr=split/\t/,$line;
	print OUT4 $arr[1],"\n";
}
close F1;close OUT4;
