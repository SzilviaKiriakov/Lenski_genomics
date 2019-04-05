#quality control using FastQC
module load fastqc/0.11.7
fastqc ~/Lenski_genomics/data/untrimmed_fastq/*.fastq.gz
mkdir ~/Lenski_genomics/results/fastqc_untrimmed_reads
mv *.zip ~/Lenski_genomics/results/fastqc_untrimmed_reads/
mv *.html ~/Lenski_genomics/results/fastqc_untrimmed_reads/
#will need to download html files to local computer and run them in browser
#summarizing results using MultiQC
module load multiqc/1.4
cd ~/Lenski_genomics/results/fastqc_untrimmed_reads
multiqc .
mv multiqc_report.html  ~/Lenski_genomics/docs/
#documenting QC results from FastQC output
for filename in *.zip
do
unzip $filename
done
cat */summary.txt > ~/Lenski_genomics/docs/fastqc_summaries.txt