#!/usr/bin/env Rscript

system("mkdir -p software")

# installs esearch

# downloads and extracts trimmomatic into the "software" directory
system("wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && unzip Trimmomatic-0.39.zip")
system("mv Trimmomatic* software")
