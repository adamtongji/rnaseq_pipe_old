#!/usr/bin/env bash
STAR_PATH=$1
GENOME_DIR=$2
GENOME_FILE=$3
ANNOTATION=$4
READLENGTH=$5
SEQTYPE=$6
FILE1=$7
FILE2=$8
Outputdir=$9
if [ $SEQTYPE = "strand_specific" ];then
    $STAR_PATH --genomeDir $GENOME_DIR --runThreadN 8 --outSAMattributes Standard\
     --outSAMtype BAM Unsorted  --readFilesIn $FILE1 $FILE2 --outFileNamePrefix A
else
    $STAR_PATH --genomeDir $GENOME_DIR --runThreadN 8 --outSAMattributes Standard\
     --outSAMstrandField intronMotif --outSAMtype BAM Unsorted\
      --readFilesIn $FILE1 $FILE2  --outFileNamePrefix A
fi