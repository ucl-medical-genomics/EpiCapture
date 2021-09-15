#!/usr/bin/bash

    BSGEN='/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/Bisulfite_Genome_bow2'
    FILEPATH='/home/regmani@ad.ucl.ac.uk/RDS_C2c/EpiCapture/Illumina_EPIC/fastq_files/merged/run-77-87'

#-------------------------------------------------
# Count number of reads and output to a file
#    for f in ./*.fastq.gz
#    do
#        echo $f >> filelist.txt
#        zgrep -c @D00 $f >> readcount.txt
#    done
#

#-------------------------------------------------
## Quality control:
    for f1 in ${FILEPATH}/*R1.fastq.gz
    do
        f2=${f1%%R1.fastq.gz}"R2.fastq.gz"
        fastqc $f1 $f2
    done

#-------------------------------------------------
## Quality score trimming and adapter removal (for RRBS use --rrbs option):
    for f1 in ${FILEPATH}/*1.fastq.gz
    do
        f2=${f1%%1.fastq.gz}"2.fastq.gz"
        trim_galore --paired --trim1 $f1 $f2
    done


#-------------------------------------------------
## Align to reference genome
    for f1 in ${FILEPATH}/*1_val_1.fq.gz
    do
        f2=${f1%%1_val_1.fq.gz}"2_val_2.fq.gz"
        bismark --bowtie2 ${BSGEN} -1 $f1 -2 $f2 --multicore 8
    done

#-------------------------------------------------
## Deduplication (not for RRBS!):
    for file in ${FILEPATH}/*val_1.fq.gz_bismark_bt2_pe.bam
    do
        deduplicate_bismark -p --bam $file
    done


#-------------------------------------------------
## Methylation extraction:
    
    for file in ${FILEPATH}/*_bismark_bt2_pe.deduplicated.bam
    do
        bismark_methylation_extractor -p  --ignore_r2 2 --no_overlap --comprehensive --bedGraph --gzip --parallel 4  $file
        ## A printout of the coverage at each cytosine can also be achieved by specifying the --bedGraph flag. This then needs to be passed through bedGraph2cytosine (also part of bismark) to produce the required methylation status and coverage report.
            # --ignore_r2 2 was introduced because of an observed bias in R2 introduced by end repair step during library perp

    done

#-------------------------------------------------
## Bismark HTML & summary report:
  bismark2report # generates a visual HTML report from the Bismark alignment, deduplication, methylation extraction and M-bias reports. It attempts to automatically detect all relevant files in the current working directory folder for you (which should work just fine here), but files may also be specified if desired

  bismark2summary ## Bismark summary report accepts Bismark BAM files as input. It will then try to  identify Bismark reports, and optionally deduplication reports or methylation extractor (splitting) reports automatically based the BAM file basename. It produces a tab delimited overview table (.txt) as well as a graphical HTML report.

    coverage2cytosine

#-------------------------------------------------
## MultiQC report  for single or multiple samples
    multiqc -f .

#--------------------
