#!/usr/bin/env python
#coding:utf-8
import os
import sys
# import traceback

sh = os.system


def starindex(starpath, star_index, genome_file, annotation_file, readlength):
    sh("{0} --runMode genomeGenerate --runThreadN 8 --genomeDir {1}\
     --genomeFastaFiles {2} --sjdbGTFfile {3} --sjdbOverhang {4}"\
       .format(starpath, star_index,genome_file, annotation_file, str(int(readlength)-1)))
    print "STAR Index build successfully at {0}".format(star_index)
    return True


def star_run(starpath, star_index,outprefix, seqtype, *files):
    # support only pair end now

    if seqtype[0].lower()=="pairend" or seqtype[0].lower()=='pair_end':

        if files[0].endswith(".gz"):
            if seqtype[1]=="strand-specific":
                sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
                Standard --outSAMtype BAM Unsorted --readFilesIn {2} {4} --outFileNamePrefix {3}"\
               .format(starpath, star_index,files[0], outprefix,files[1]))
            else:
                sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
                Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --readFilesIn\
                 {2} {4} --outFileNamePrefix {3}".format(starpath, star_index,files[0], outprefix,files[1]))
        elif files[0].endswith(".fastq"):
            if seqtype[1]=="strand-specific":
                sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
                Standard --outSAMtype BAM Unsorted --readFilesIn {2} {4} --outFileNamePrefix {3}"\
               .format(starpath, star_index,files[0], outprefix,files[1]))
            else:
                sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
                Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --readFilesIn\
                 {2} {4} --outFileNamePrefix {3}".format(starpath, star_index,files[0], outprefix,files[1]))
        else:
            print "please give fastq or gzipped fastq files"
            sys.exit(1)
    elif seqtype[0].lower()=="singleend" or seqtype[0].lower()=='single_end':

        if files[0].endswith(".gz"):
            if seqtype[1]=="strand-specific":
                sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
                Standard --outSAMtype BAM Unsorted --readFilesIn {2} --outFileNamePrefix {3}"\
               .format(starpath, star_index,files[0], outprefix))
            else:
                sh("{0} --genomeDir {1} --runThreadN 8 --readFilesCommand zcat --outSAMattributes\
                Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --readFilesIn\
                 {2} --outFileNamePrefix {3}".format(starpath, star_index,files[0], outprefix))
        elif files[0].endswith(".fastq"):
            if seqtype[1]=="strand-specific":
                sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
                Standard --outSAMtype BAM Unsorted --readFilesIn {2} --outFileNamePrefix {3}"\
               .format(starpath, star_index,files[0], outprefix))
            else:
                sh("{0} --genomeDir {1} --runThreadN 8 --outSAMattributes\
                Standard --outSAMstrandField intronMotif --outSAMtype BAM Unsorted --readFilesIn\
                 {2} --outFileNamePrefix {3}".format(starpath, star_index,files[0], outprefix))
        else:
            print "please give fastq or gzipped fastq files"
            sys.exit(1)





