#stop pipeline if error comes up
set -e

echo "space usage on HCC..."
pquota -u khalil

mkdir Lenski_genomics
echo "organizing project into libraries..."
mkdir -p ~/Lenski_genomics/data
mkdir ~/Lenski_genomics/results
mkdir ~/Lenski_genomics/src
mkdir ~/Lenski_genomics/doc

echo "downloading data..."
mkdir -p ~/Lenski_genomics/data/untrimmed_fastq
cd ~/Lenski_genomics/data/untrimmed_fastq
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
curl -O ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz 

echo "quality control using FastQC..."
module load fastqc/0.11.7
fastqc ~/Lenski_genomics/data/untrimmed_fastq/*.fastq.gz
mkdir ~/Lenski_genomics/results/fastqc_untrimmed_reads
mv *.zip ~/Lenski_genomics/results/fastqc_untrimmed_reads/
mv *.html ~/Lenski_genomics/results/fastqc_untrimmed_reads/
#will need to download html files to local computer and run them in browser
echo "summarizing results using MultiQC..."
module load multiqc/1.4
cd ~/Lenski_genomics/results/fastqc_untrimmed_reads
multiqc .
mv multiqc_report.html  ~/Lenski_genomics/docs/
echo "documenting QC results from FastQC output..."
for filename in *.zip
do
unzip $filename
done
cat */summary.txt > ~/Lenski_genomics/docs/fastqc_summaries.txt
#run on local computer to analyze report:
#scp kiriakov@ :~/Lenski_genomics/docs/* ~/Downloads
#cd Downloads
#start chrome multiqc_report.html
#grep FAIL fastqc_summaries.txt

#Filtering low quality reads and trimming poor quality bases with Trimmomatic to reduce false positive rate due to sequencing errors

module load trimmomatic/0.36
cd ~/Lenski_genomics/data/untrimmed_fastq

echo "downloading Nextera adaptor sequences..."
curl -L -o ~/Lenski_genomics/data/untrimmed_fastq/NexteraPE-PE.fa https://repo.jbei.org/users/mwornow/repos/seqvalidation/raw/Trimmomatic-0.33/adapters/NexteraPE-PE.fa?at=40d6f90987fc825914667f80e38ce9ebc6a24734

#mv NexteraPE-PE.fa?at=40d6f90987fc825914667f80e38ce9ebc6a24734 NexteraPE-PE.fa

#trimmomatic PE SRR2589044_1.fastq.gz SRR2589044_2.fastq.gz SRR2589044_1.trim.fastq.gz SRR2589044_1orphan.trim.fastq.gz SRR2589044_2.trim.fastq.gz SRR2589044_2orphan.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
#trimmomatic PE SRR2584863_1.fastq.gz SRR2584863_2.fastq.gz SRR2584863_1.trim.fastq.gz SRR2584863_1orphan.trim.fastq.gz SRR2584863_2.trim.fastq.gz SRR2584863_2orphan.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
#trimmomatic PE SRR2584866_1.fastq.gz SRR2584866_2.fastq.gz SRR2584866_1.trim.fastq.gz SRR2584866_1orphan.trim.fastq.gz SRR2584866_2.trim.fastq.gz SRR2584866_2orphan.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15

echo "Running Trimmomatic..."
$ for infile in *_1.fastq.gz
do
base=$(basename ${infile} _1.fastq.gz)
trimmomatic PE ${infile} ${base}_2.fastq.gz ${base}_1.trim.fastq.gz ${base}_1orphan.trim.fastq.gz ${base}_2.trim.fastq.gz ${base}_2orphan.trim.fastq.gz SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
done

mkdir ~/Lenski_genomics/data/trimmed_fastq
mv *.trim* ../trimmed_fastq

#rerun QC to see if trimming improved quality
cd ~/Lenski_genomics/data/trimmed_fastq
fastq *fastq*

mkdir ~/Lenski_genomics/results/fastqc_trimmed_reads
mv *.zip ~/Lenski_genomics/results/fastqc_trimmed_reads/
mv *.html ~/Lenski_genomics/results/fastqc_trimmed_reads/
cd ~/Lenski_genomics/results/fastqc_trimmed_reads
multiqc .
mv multiqc_report.html multiqc_report2.html
mv multiqc_report2.html  ~/Lenski_genomics/docs/

echo "unzipping trimmed fastq files..."
cd ~/Lenski_genomics/data/trimmed_fastq
gunzip *.gz

echo "downloading reference genome..."
mkdir -p ~/Lenski_genomics/data/ref_genome
curl -L -o ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz
gunzip ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta.gz

echo "creating directories for the results..."
mkdir -p results/sam results/bam results/bcf results/vcf

echo "indexing reference genome..."
module load bwa/0.7.17
bwa index ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta

echo "aligning reads to reference genome..."
bwa mem ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta ~/Lenski_genomics/data/trimmed_fastq/SRR2589044_1.trim.fastq ~/Lenski_genomics/data/trimmed_fastq/SRR2589044_2.trim.fastq > ~/Lenski_genomics/results/sam/SRR2589044.aligned.sam
bwa mem ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta ~/Lenski_genomics/data/trimmed_fastq/SRR2584863_1.trim.fastq ~/Lenski_genomics/data/trimmed_fastq/SRR2584863_2.trim.fastq > ~/Lenski_genomics/results/sam/SRR2584863.aligned.sam
bwa mem ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta ~/Lenski_genomics/data/trimmed_fastq/SRR2584866_1.trim.fastq ~/Lenski_genomics/data/trimmed_fastq/SRR2584866_2.trim.fastq > ~/Lenski_genomics/results/sam/SRR2584866.aligned.sam

echo "converting SAM files into BAM files..."
module load samtools/1.8
samtools view -S -b ~/Lenski_genomics/results/sam/SRR2589044.aligned.sam > ~/Lenski_genomics/results/bam/SRR2589044.aligned.bam
samtools view -S -b ~/Lenski_genomics/results/sam/SRR2584863.aligned.sam > ~/Lenski_genomics/results/bam/SRR2584863.aligned.bam
samtools view -S -b ~/Lenski_genomics/results/sam/SRR2584866.aligned.sam > ~/Lenski_genomics/results/bam/SRR2584866.aligned.bam

echo "sorting BAM files..."
samtools sort -o ~/Lenski_genomics/results/bam/SRR2589044.aligned.sorted.bam ~/Lenski_genomics/results/bam/SRR2589044.aligned.bam 
samtools sort -o ~/Lenski_genomics/results/bam/SRR2584863.aligned.sorted.bam ~/Lenski_genomics/results/bam/SRR2584863.aligned.bam 
samtools sort -o ~/Lenski_genomics/results/bam/SRR2584866.aligned.sorted.bam ~/Lenski_genomics/results/bam/SRR2584866.aligned.bam 

echo "Getting information about BAM files..."
samtools flagstat ~/Lenski_genomics/results/bam/SRR2589044.aligned.sorted.bam > SRR2589044.aligned.sorted.bam.info.txt
samtools flagstat ~/Lenski_genomics/results/bam/SRR2584863.aligned.sorted.bam > SRR2584863.aligned.sorted.bam.info.txt
samtools flagstat ~/Lenski_genomics/results/bam/SRR2584866.aligned.sorted.bam > SRR2584866.aligned.sorted.bam.info.txt
#may check these with cat

#variant calling
module load bcftools/1.9
echo "Calculating read coverage ..."
bcftools mpileup -O b -o ~/Lenski_genomics/results/bcf/SRR2589044_raw.bcf -f ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta ~/Lenski_genomics/results/bam/SRR2589044.aligned.sorted.bam
bcftools mpileup -O b -o ~/Lenski_genomics/results/bcf/SRR2584863_raw.bcf -f ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta ~/Lenski_genomics/results/bam/SRR2584863.aligned.sorted.bam
bcftools mpileup -O b -o ~/Lenski_genomics/results/bcf/SRR2584866_raw.bcf -f ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta ~/Lenski_genomics/results/bam/SRR2584866.aligned.sorted.bam

echo "Calling SNPs..."
bcftools call --ploidy 1 -m -v -o ~/Lenski_genomics/results/bcf/SRR2589044_variants.vcf ~/Lenski_genomics/results/bcf/SRR2589044_raw.bcf
bcftools call --ploidy 1 -m -v -o ~/Lenski_genomics/results/bcf/SRR2584863_variants.vcf ~/Lenski_genomics/results/bcf/SRR2584863_raw.bcf
bcftools call --ploidy 1 -m -v -o ~/Lenski_genomics/results/bcf/SRR2584866_variants.vcf ~/Lenski_genomics/results/bcf/SRR2584866_raw.bcf 
echo "Filtering and reporting SNPs..."
vcfutils.pl varFilter ~/Lenski_genomics/results/bcf/SRR2589044_variants.vcf  > ~/Lenski_genomics/results/vcf/SRR2589044_final_variants.vcf
vcfutils.pl varFilter ~/Lenski_genomics/results/bcf/SRR2584863_variants.vcf  > ~/Lenski_genomics/results/vcf/SRR2584863_final_variants.vcf
vcfutils.pl varFilter ~/Lenski_genomics/results/bcf/SRR2584866_variants.vcf  > ~/Lenski_genomics/results/vcf/SRR2584866_final_variants.vcf

#quick assessment of number of variants
#-v is invert match in grep
grep -v "#" ~/Lenski_genomics/results/vcf/SRR2589044_final_variants.vcf | wc -l
grep -v "#" ~/Lenski_genomics/results/vcf/SRR2584863_final_variants.vcf | wc -l
grep -v "#" ~/Lenski_genomics/results/vcf/SRR2584866_final_variants.vcf | wc -l

#indexing alignment files
samtools index ~/Lenski_genomics/results/bam/SRR2589044.aligned.sorted.bam
samtools index ~/Lenski_genomics/results/bam/SRR2584863.aligned.sorted.bam
samtools index ~/Lenski_genomics/results/bam/SRR2584866.aligned.sorted.bam

#tview inspection
samtools tview ~/Lenski_genomics/results/bam/SRR2589044.aligned.sorted.bam ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta
samtools tview ~/Lenski_genomics/results/bam/SRR2584863.aligned.sorted.bam ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta
samtools tview ~/Lenski_genomics/results/bam/SRR2584866.aligned.sorted.bam ~/Lenski_genomics/data/ref_genome/ecoli_rel606.fasta

