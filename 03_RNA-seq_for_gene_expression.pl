#usr/bin/perl
#The program can be used for mapping 3'RNA-seq data with single-reads and generate gene expression level.

#usage:perl script.pl $ARGV[0]

#dependents
#bowtie2-2.2.9
#hisat2
#fastq

#the pathway for index of genome
$path_to_histat2_index="../genome";# modification
$path_to_fa="../genome";# modification
$path_to_gff_index="../genome.gff";# modification

#raw reads
$R2=$ARGV[0]."_2.fq.gz";
#filtered reads
$CR2=$ARGV[0]."_2.clean.fastq.gz";
#QC with fastp
system("fastp --detect_adapter_for_pe -i $R2 -o $CR2 -q 20 -w 2 -3 -M 30 -Y 20 -j $ARGV[0].json -h $ARGV[0].html");
#mapping with bwa
system("hisat2 --dta -p 4 -x $path_to_bismark_index -U $CR2 -S $ARGV[0].sam");#-t: thread

#QC for mapping reads
#the number of mapped reads
open(F1,"$ARGV[0].sam") or die "$!";
open(UNIQUE,">$ARGV[0]_unique.sam") or die "$!";
$unique=0;$concordant=0;
while(my $line=<F1>){
	chomp $line;
	my @b=split/\t/,$line;
	if($b[0]=~"@"){print UNIQUE $line,"\n";}else{
	   my @arr=split/\t/,$line;
		if($arr[-1] eq "NH:i:1"){
			print UNIQUE $line,"\n";
			   $unique+=1;
			   }

	}
}
close F1;close UNIQUE;

##########################################################
#gene expression(RPKM)
open(F1,"$ARGV[0]_unique.sam") or die "$!";
%hash=();
%hash1=();
while(my $line=<F1>){
	chomp $line;
	my @arr=split/\t/,$line;
	
	if($arr[0]=~"@"){$arr[0]=1;}else{
    		
	$hash{$arr[2]."|".$arr[3]}=1;
	$a=$arr[3]+149;
    $hash1{$arr[2]."|".$a}=1;

	}
}

open(F1,$path_to_gff_index) or die "$!";
open(EXPR,">$ARGV[0]_EXPR.txt") or die "$!";
print EXPR "gene","\t",$ARGV[1],"\n";
while(my $line=<F1>){
	chomp $line;
	my @arr=split/\t/,$line;
	if($arr[2] eq "gene"){
	$j=0;
	$i = $arr[3];
	while($i <= $arr[4]){
	if((exists $hash{$arr[0]."|".$i}) or (exists $hash1{$arr[0]."|".$i})){
	  $j+=1;
         }
	  	  $i+=1;
         }
		$geneleng=abs($arr[3]-$arr[4]);
		$rpkm=$j/$geneleng*1000/$unique*1000000;
        print EXPR $arr[-1],"\t",$rpkm,"\n";#文件1
}
}
#######################################################################
system("samtools sort -@ 4 -o $ARGV[0]_sorted.bam $ARGV[0]_unique.sam");

#######################################################################

#IGV
system("samtools index $ARGV[0]_sorted.bam");
system("igvtools count $ARGV[0]_sorted.bam $ARGV[0].tdf $path_to_fa");

#the summary of number of raw and clean reads
open(F1,"$ARGV[0].json") or die "$!";
$i=0;
while(my $line=<F1>){
	chomp $line;
	my @a=split/\:|\,/,$line;
	$i+=1;
	if($i==4){
		$total_reads=$a[1];
	}else{
		if($i==15){
			$clean_reads=$a[1];
		}
	}
}
close F1;
##########################################################################
open(O,">$ARGV[0]_qc_map_statistic_report.txt") or die "$!";
print O $ARGV[0],"\t","Total reads","\t",$total_reads,"\t",
	"Clean reads","\t",$clean_reads,"\t",
	"Clean reads%","\t",$clean_reads/$total_reads,"\t",
	"map uniquely%","\t",$unique/$clean_reads,"\n";

close O;

