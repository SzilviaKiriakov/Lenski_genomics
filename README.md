# Lenski_genomics
Variant calling workflow from Lenski's long term evolution experiment NGS data.

Thid project downloads fastq files from Richard Lenski's long term laboratory evolution experiment, where E.coli bacteria were grown for 40000+ generations. Reference: https://www.pnas.org/content/105/23/7899
The data is paired-end NGS from 5000, 15000 and 50000 generation timepoints.
It organizes the project into directories, downloads sequencing data, and performs quality checks.
Then it trims, filters the reads, aligns them to the reference genome, and automates variant calling.
