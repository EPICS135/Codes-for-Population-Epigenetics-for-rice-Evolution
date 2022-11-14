#usr/bin/perl
#The program can be used for mapping re-seq data with PE-reads and generate snp.

#usage:perl script.pl $ARGV[0]

#dependents
#bwa
#fastq
#gatk

#the pathway for index of genome
$path_to_bwa_index="../genome";# modification

#raw reads
$R1=$ARGV[0].".R1.fastq.gz";
$R2=$ARGV[0].".R2.fastq.gz";
#filtered reads
$CR1=$ARGV[0].".R1.clean.fastq.gz";
$CR2=$ARGV[0].".R2.clean.fastq.gz";
#QC with fastp
system("fastp --detect_adapter_for_pe -i $R1 -I $R2 -o $CR1 -O $CR2 -q 20 -w 2 -3 -M 30 -Y 20 -j $ARGV[0].json -h $ARGV[0].html");
#mapping with bismark
system("bwa mem -t 4 $path_to_bwa_index $CR1 $CR2 >$ARGV[0].sam");
unlink $CR1;unlink $CR2;
#
open (OUT,">$ARGV[0]_unique.sam") or die "$!";
open (F1,$ARGV[0]) or die "$!";
$j=0;$m=0;
while(my $line=<F1>){
	chomp $line;
	if($i=~ "#"){
		print OUT $line,"\n";
	}else{
			my @arr=split/\t/,$line;
			@first=split/:/,$arr[13];
			@second=split/:/,$arr[14];
			$j+=1;
		if($first[2] ne $second[2]){
			print OUT $line,"\n";
			$m+=1;
		}
	}
}
close F1;
#generate SNPs of each lines
system("samtools sort -@ 4 -o $ARGV[0].bam $ARGV[0]_unique.sam");
system("samtools rmdup $ARGV[0].bam $ARGV[0]_rmdup.bam");
system("samtools mpileup -B -C 50 -uf $path_to_bwa_index $ARGV[0]_rmdup.bam | bcftools call -c -v --output-type v > $ARGV[0].vcf");
system("gatk IndexFeatureFile -F $ARGV[0].vcf");
system("gatk BaseRecalibrator -R $path_to_bwa_index -I $ARGV[0]_unique_rmdup_add.bam -O $ARGV[0]_recal_data.grp -known-sites $ARGV[0].vcf --bqsr-baq-gap-open-penalty 30");
system("gatk ApplyBQSR -R $path_to_bwa_index -I $ARGV[0]_unique_rmdup_add.bam -bqsr $ARGV[0]_recal_data.grp -O $ARGV[0]_recal.bam");
system("gatk HaplotypeCaller -R $path_to_bwa_index -I $ARGV[0]_recal.bam -O $ARGV[0].g.vcf -ERC GVCF --genotyping-mode DISCOVERY --output-mode EMIT_VARIANTS_ONLY -stand-call-conf 30 -ploidy 2 -bamout $ARGV[0]_Haplotype.bam"); 

