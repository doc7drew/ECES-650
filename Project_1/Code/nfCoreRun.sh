#!/bin/bash

nextflow run nf-core/ampliseq -profile conda --input_folder "project1DataRenamedSingleEnd/" --skip_cutadapt --extension "*_extendedFrags.fastq.gz" --single_end  --outdir output1