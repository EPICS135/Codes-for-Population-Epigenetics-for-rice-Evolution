#usr/bin/perl
#The program can be used for mapping DAP-seq with PE-reads and generate reports for QC, mapping and peaks-calling.

#usage:perl script.pl $ARGV[0]

#dependents
#fastq
#bwa
#samtools
#igvtools


#the path for index of genome for bwa
$path_to_bwa_genome_index="../genome";# modification
$picard="/MarkDuplicates.jar";# modification

#raw reads
$R1=$ARGV[0].".R1.fastq.gz";
$R2=$ARGV[0].".R2.fastq.gz";
#filtered reads
$CR1=$ARGV[0].".R1.clean.fastq.gz";
$CR2=$ARGV[0].".R2.clean.fastq.gz";

#QC with fastp
system("fastp --detect_adapter_for_pe -i $R1 -I $R2 -o $CR1 -O $CR2 -q 20 -w 2 -3 -M 30 -Y 20 -j $ARGV[0].json -h $ARGV[0].html");
#mapping with bwa
system("bwa mem -t 2 $path_to_bwa_genome_index $CR1 $CR2 >$ARGV[0].sam");#-t: thread
unlink $CR1;unlink $CR2;

#QC for mapping reads
system("samtools view -bhS -q 20 $ARGV[0].sam > $ARGV[0]_MAPQ20.bam");
system("samtools sort $ARGV[0]_MAPQ20.bam -o $ARGV[0]_sorted.bam");
system("java -Xmx40G -jar $picard REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true I=$ARGV[0]_sorted.bam M=$ARGV[0]_sorted.metrics O=$ARGV[0]_picard.bam");

#IGV
system("igvtools count $ARGV[0]_picard.bam $ARGV[0].tdf $path_to_bwa_genome_index");


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
open(F1,"$ARGV[0].sam") or die "$!";
$unique=0;$total=0;
while(my $line=<F1>){
	chomp $line;
	my @b=split/\t/,$line;
	if($b[0]=~"#"){$b[0]=1;}else{
	   my @arr=split/\t/,$line;
	      @first=split/:/,$arr[13];
			@second=split/:/,$arr[14];
			$total+=1;
		if($first[2] ne $second[2]){
			$unique+=1;
		}
	}
}

close F1;
open(O,">$ARGV[0]_qc_map_statistic_report.txt") or die "$!";
print O $ARGV[0],"\t","Total reads","\t",$total_reads,"\t",
	"Clean reads","\t",$clean_reads,"\t",
	"Clean reads%","\t",$clean_reads/$total_reads,"\t",
	"Unique map reads%","\t",$unique/$clean_reads,"\n";

close O;

