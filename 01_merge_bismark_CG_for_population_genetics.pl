#usr/bin/perl
#The program can be used for merging cytosine methylation level of each accession.

#usage:perl script.pl $ARGV[0]

#the pathway for snps
$snps_file="../03_select.snp.vcf";# modification

system("paste *_CG-sort-for_paste.txt >$ARGV[0]-for_DMRs_and_population_genetics.txt");
#removed snps
open(F1,$snps_file) or die "$!";
open(OUT1,">$ARGV[0]_moveSNP.txt") or die "$!";
%hash=();
while(my $line=<F1>){
	chomp $line;
	my @arr=split/\t/,$line;
	$key=$arr[0]."-".$arr[1];
		$hash{$key}=$line;
	}

open(F1,$ARGV[0]) or die "$!";
$i=1;
while(my $line=<F1>){
	chomp $line;
	my @arr=split/\t/,$line;
	$key=$arr[0];
	if($i>2){
	if(exists $hash{$key}){$i+=1;}else{
		print OUT1 $line,"\n";}else{
         print OUT1 $line,"\n";
		 }
$i+=1;
}

}

close F1;

#conserved to methylation haplotype (meplotype)
open(F1,$ARGV[0]_moveSNP.txt) or die "$!";
open(CHR1,">$ARGV[0]_meplotype.vcf") or die "$!";
$ii=1;
while (my $line=<F1>){
	chomp $line;
	
	if($ii<=2){
	$ii+=1;
	}else{
my @arr=split/\t/,$line;
$i=1;
@c=();
$num=@arr;
while($i<=$num){
if($arr[$i] ne "NA"){if($arr[$i] <= 0.1){$a = "0/0";}else{if(($arr[$i] > 0.1)&&($arr[$i]<= 0.6)){$a = "0/1";}else{$a = "1/1";}}} else{$a = "./.";}
$i+=1;
push(@c,$a,"\t");

}
$ii+=1;
my @brr=split/\-/,$arr[0];
print CHR1 $brr[0],"\t",$brr[1],"\t","CG",$ii,"\t","T	C	.	.	PR	GT","\t",@c,"\n";
}
}
close F1;close CHR1;
#PCA
system("plink --vcf $ARGV[0]_meplotype.vcf --make-bed --out $ARGV[0]_meplotype");
system("plink --vcf $ARGV[0]_meplotype-plink.vcf --pca 10 --out $ARGV[0]_meplotype-plink-pca");
#Nucleotide diversity (p)was
system("vcftools --vcf $ARGV[0]_meplotype-plink.vcf --keep wild.txt --recode --recode-INFO-all --out wild");
system("vcftools --vcf wild.recode.vcf --out wild --window-pi 100000 --window-pi-step 10000");
system("paste wild.windowed.pi.txt aro.pi.txt aus.pi.txt geng.pi.txt indica.pi.txt weedyin.pi.txt weedyge.pi.txt >all_line.pi.txt");

