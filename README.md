# Population-epigenetics-for-rice

1 Analysis of DNA mehylation of rice population

1.1 The program can be used for mapping BS-seq data with PE-reads and generate reports for QC, mapping and methylation level.

#dependents

#bowtie2-2.2.9

#bismark

#fastq

usage:perl 01_bismark_for_PE_reads.pl $ARGV[0]

1.2 Merged cytosines methylated state of each accession.

usage:perl 01_merge_bismark_CG_for_population_genetics.pl 01_CG

2 Analysis of SNPs

2.1 The program can be used for mapping re-seq data with PE-reads and generate snp.

#dependents

#bwa

#fastq

#gatk

usage:perl 02_reseq_generate_SNPs_for_PE_reads.pl $ARGV[0]

2.2 merged GenotypeGVCFs

gatk --java-options -Xmx60g CombineGVCFs -V line1.g.vcf -V line1.g.vcf -R ../genome -O 01_CombineGVCFs_out.vcf

gatk --java-options -Xmx60g GenotypeGVCFs -R ../genome -V 01_CombineGVCFs_out.vcf -O 01_GenotypeGVCFs_out.vcf

2.3 SplitVcfs

gatk SelectVariants -R ../genome -select-type SNP -V 01_GenotypeGVCFs_out.vcf -O 02_GenotypeGVCFs_out.snp.vcf

2.4 filt SNP

gatk VariantFiltration -R ../genome -V 02_GenotypeGVCFs_out.snp.vcf -O 02_VariantFiltration.snp.vcf --filter-expression "QD < 10.0 || MQ < 50.0 || FS > 30.0 || SOR > 3.0 || MQRankSum < -2.5 || ReadPosRankSum < -3.5 || DP<2.0 ||DP>100.0" --filter-name "Fail" -G-filter "DP<2 || isHet==1 || isHomRef==1" -G-filter-name "Failed" &

#and then selected filter SNPs

2.5 conversion using Plink

plink --vcf 02_VariantFiltration_select.snp.vcf --recode vcf-iid --out 03_select.snp --allow-extra-chr

2.6 conversion to phy

python vcf2phylip.py --input 03_select.snp.vcf

2.7 Phylogeny Reconstruction of Neighbor-joining

/DATA2/maize2/21_software/MegaX/megacc -a infer_NJ_nucleotide.mao -d 03_select.snp.meg -o 03_select.snp.recode

2.8 Nucleotide diversity (p)was

vcftools --vcf 03_select.snp.vcf --keep wild.txt --recode --recode-INFO-all --out wild

vcftools --vcf wild.recode.vcf --out wild --window-pi 100000 --window-pi-step 10000

paste wild.windowed.pi.txt aro.pi.txt aus.pi.txt geng.pi.txt indica.pi.txt weedyin.pi.txt weedyge.pi.txt >all_line.pi.txt

3 RNA-seq analysis

3.1 split each rnaseq data based on the indexs

python split_fastq_for_each_tissue.py $R1 $R2 index.txt

3.2The program can be used for mapping 3RNA-seq data with single-reads and generate gene expression level.

#dependents

#bowtie2-2.2.9

#hisat2

#fastq

usage:perl 03_RNA-seq_for_gene_expression.pl $ARGV[0]

4 DAP-seq analysis

4.1 The program can be used for mapping DAP-seq with PE-reads and generate reports for QC, mapping and IGV.

#dependents

#fastq

#bwa

#samtools

#igvtools

usage:perl 04_DAP-seq_with_PE_reads_mapping.pl $ARGV[0]

4.2 The program can be used for peaks-calling.

#dependents

#MACS2.1.2

#gem

#bedtools

#meme

usage:perl 04_DAP-seq_with_PE_reads_callpeaks_under_control.pl $ARGV[0] $ARGV[1]

5 ATAC-seq analysis

5.1 The program can be used for mapping ATAC-seq with PE-reads and generate reports for QC, mapping and IGV.

#dependents

#fastq

#bwa

#samtools

#igvtools

usage:perl 05_ATAC-seq_for_PE_reads_mapping.pl $ARGV[0]

5.2 The program can be used for peaks-calling.

#dependents

#MACS2.1.2

#bedtools

#meme

usage:perl 05_DAP-seq_with_PE_reads_callpeaks_under_nakedDNA.pl $ARGV[0] $ARGV[1]

