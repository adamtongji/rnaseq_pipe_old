#!/usr/bin/env bash
STAR_PATH=$1
GENOME_DIR=$2
GENOME_FILE=$3
ANNOTATION=$4
READLENGTH=$5

# build index with default parameter
$STAR_PATH --runMode genomeGenerate --runThreadN 8 --genomeDir $GENOME_DIR\
 --genomeFastaFiles $GENOME_FILE --sjdbGTFfile $ANNOTATION --sjdbOverhang $READLENGTH-1