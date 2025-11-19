# Install in the base environment
conda install motus

# OR, create a new environment
conda create -n motu-env motus
conda activate motu-env

# Download and install mOTUs
pip install motu-profiler
# Download the mOTUs database
motus downloadDB

# test that motus installed correctly
motus profile --test



# generating taxonomic profiles
# Using the mOTU profiler with test files contained in the installation/database directory
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test1.motus -n test1

#resulting profile reports the relative abundance for each mOTU (ref + meta + ext)
head -n 3 test1.motus
awk '{print $NF,$0}' test1.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# read counting
# The -c flag changes the output from relative abundance to number of assigned reads
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test1.motus -n test1 -c
head -n 3 test1.motus
awk '{print $NF,$0}' test1.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# taxonomic level
# The -k flag changes taxonomic level
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test1.motus -n test1 -k phylum
awk '{print $NF,$0}' test1.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# threading
multiple cores (-t flag) to accelerate the alignment process can be assigned
motus profile -f sample_R1.fq.gz -r sample_R2.fq.gz -t 8 -o test1.motus

# database selection
# the -e flag performs taxonomic profiling using only ref-mOTUs database
# motus profile -f sample_R1.fq.gz -r sample_R2.fq.gz -e -o test1.motus

# merging profiles
motus merge -i test1.motus,test2.motus -o test.motus
head -n 3 test.motus
awk '{print $NF,$0}' test.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# Adding bee samples to you profiles
motus merge -i test1.motus,test2.motus -a bee -o test.motus
# Adding bee and dog samples to you profiles
motus merge -i test1.motus,test2.motus -a bee,dog -o test.motus
# Adding all public profiles to your samples
motus merge -i test1.motus,test2.motus -a all -o test.motus

# view more options
motus profile

