#usr/bin/perl -w
#The program can be used for peaks-calling.

#usage:perl script.pl $ARGV[0] $ARGV[1]

#dependents
#MACS2.1.2
#gem
#bedtools
#meme



#the path for index of genome for bwa
$path_to_bwa_genome_index="../genome.fa";# modification
$picard="/MarkDuplicates.jar";# modification
$genomesize="340000000";#rice
$chrom_rice_sizes="chrom_rice.sizes";
#mapped reads
$T=$ARGV[0].".T.picard.bam";
$C=$ARGV[1].".C.picard.bam";
#macs2
system("macs2 callpeak -t $T -c $C -f BAM -n NAME --down-sample --shift -100 --extsize 200 --nomodel -B --SPMR --keep-dup all --gsize $genomesize --outdir out_macs2");

#motifs of peaks from macs2
system("bedtools getfasta -fi $path_to_bwa_genome_index -bed $ARGV[0].T.picard_peaks.narrowPeak >$ARGV[0].T-narrowPeak.fasta");
system("meme-chip -oc . -time 30000 -ccut 100 -order 1 -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 3 -meme-searchsize 100000 -dreme-e 0.05 -centrimo-score 5.0 -centrimo-ethresh 10.0 $ARGV[0].T-narrowPeak.fasta -o $ARGV[0].T");


