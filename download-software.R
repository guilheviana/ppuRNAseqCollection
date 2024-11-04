#!/usr/bin/env Rscript

# installs sra toolkit
system("sudo apt install sra-toolkit")

# downloads and extracts trimmomatic into the "software" directory
system("mkdir software")
system("wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip")
system("mv Trimmomatic* software")

# downloads P. putida KT2440 cds from genome file from ncbi
system("mkdir -p software/rsem")
system("wget -O software/kt2440_ASM756v2_cds.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/565/GCF_000007565.2_ASM756v2/GCF_000007565.2_ASM756v2_cds_from_genomic.fna.gz && gzip -dk software/kt2440_ASM756v2_cds.fna.gz ")

# installs RSEM using conda
system('conda create -n rsem -c bioconda rsem')
# create RSEM reference file from KT2440 transcriptome file
system("rsem-prepare-reference --bowtie2 software/kt2440_ASM756v2_cds.fna software/rsem/kt2440")
