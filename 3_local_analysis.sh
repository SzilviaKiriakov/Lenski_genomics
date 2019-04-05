#run on local computer to analyze report:
#scp kiriakov@ :~/Lenski_genomics/docs/* ~/Downloads
cd Downloads
start chrome multiqc_report.html
grep FAIL fastqc_summaries.txt