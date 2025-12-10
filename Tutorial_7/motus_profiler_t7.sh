# 1. Generating taxonomic profiles

# Standard example
# Using the mOTU profiler with test files contained in the installation/database directory
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test1.motus -n test1
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test1a.motus -n test1a
# relative abundance for each mOTU
head -n 3 test1.motus
awk '{print $NF,$0}' test1.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# read counting: -c is number assigned reads
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o tp_read_count.motus -n test1 -c
head -n 3 tp_read_count.motus
awk '{print $NF,$0}' tp_read_count.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# taxonomic level: -k changes taxonomic level
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o tp_tax.motus -n test1 -k phylum
head -n 3 tp_tax.motus
awk '{print $NF,$0}' tp_tax.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# threading: -t assign multiple cores
motus profile -f sample_R1.fq.gz -r sample_R2.fq.gz -t 8 -o tp_thread.motus

# database selection: -e taxonomic profiling using ref-mOTUs database
motus profile -f sample_R1.fq.gz -r sample_R2.fq.gz -e -o tp_db_select.motus

# merging profiles
motus merge -i test1.motus,test2.motus -o tp_merge.motus
head -n 3 tp_merge.motus
awk '{print $NF,$0}' tp_merge.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# Adding bee samples to you profiles
motus merge -i test1.motus,test1a.motus -a bee -o tp_bee.motus
# Adding bee and dog samples to you profiles
motus merge -i test1.motus,test1a.motus -a bee,dog -o tp_beedog.motus
# Adding all public profiles to your samples
motus merge -i test1.motus,test1a.motus -a all -o tp_all.motus

motus profile


# 2. Generating metatranscriptomic profiles

# Standard example
# Using the mOTU profiler with test files contained in the installation/database directory
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test2.motus -n test2
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o test2a.motus -n test2a
head -n 3 test2.motus
awk '{print $NF,$0}' test2.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# read counting: -c changes the output from relative abundance to number of assigned reads:
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o mt_read_count.motus -n test2 -c
head -n 3 mt_read_count.motus
awk '{print $NF,$0}' mt_read_count.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# taxonomic level: -k changes taxonomic leve
motus profile -s mOTUs/db_mOTU/db_mOTU_test/test1_single.fastq -o mt_tax.motus -n test2 -k phylum
head -n 3 mt_tax.motus
awk '{print $NF,$0}' mt_tax.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# threading: -t accelerate the alignment process
motus profile -f sample_R1.fq.gz -r sample_R2.fq.gz -t 8 -o mt_thread.motus

# database selection: -e perform taxonomic profiling using only ref-mOTUs database
motus profile -f sample_R1.fq.gz -r sample_R2.fq.gz -e -o mt_db_select.motus

# merging profiles
motus merge -i test2.motus,test2a.motus -o mt_merge.motus
head -n 3 mt_merge.motus
awk '{print $NF,$0}' mt_merge.motus | sort -n -k 1 | cut -f 2- -d " " | tail -n 10

# Adding bee samples to you profiles
motus merge -i test1.motus,test2.motus -a bee -o mt_b.motus
# Adding bee and dog samples to you profiles
motus merge -i test1.motus,test2.motus -a bee,dog -o mt_bd.motus
# Adding all public profiles to your samples
motus merge -i test1.motus,test2.motus -a all -o mt_a.motus

motus profile


# 3. Generating single nucleotide variant (SNV) profiles using MGs

motus map_snv -s sample.fq.gz > sample.bam
motus map_snv -f sample_R1.fq.gz -r sample_R2.fq.gz > sample.bam

motus map_snv -f sample_R1.fq.gz -r sample_R2.fq.gz -l 100 -t 8> sample.bam

motus snv_call -d DIRECTORY -o OUTPUT_DIRECTORY

motus snv_call
