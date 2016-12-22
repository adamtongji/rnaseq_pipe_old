#!/usr/bin/env python
#coding:utf-8
"""
Licensed Materials - Property of Tongji University
(C) Copyright Tongji University LifeBio Dept. 2016, 2016 All Rights Reserved
--------------------------------------------------------------------------------------
File Name   : rainbow_rnaseq.py
Description : RNAseq pipeline

Author: Shaliu Fu
Change activity:
v1.0.1
Support single end data and mouse data
Generate final files format for customers

V1.0.2
Add support to circRNA, need to add samtools path for different version of samtools

Version: v1.0.2
"""
import os, sys,time
import collections
import functools
from lib.star import starindex, star_run
from lib.ciri import ciri_run, ciri_process
from lib.de_expr import expr
from multiprocessing import Pool

sh = os.system

def const():
    """
    this function defines inner constants, which will not change
    in whole scripts
        Usage:
                import const (class _const saved as const.py)
                const.ConstNames = ConstValues
        Notice:
                ConstNames should be all upperCase!
                this func(class) can save as another file for import

    """
    class _const:
        class ConstError(TypeError):
            pass

        class ConstCaseError(ConstError):
            pass

        def __setattr__(self, name, value):
            if self.__dict__.has_key(name):
                raise self.ConstError, "Can't change const.%s" % name
            if not name.isupper():
                raise self.ConstCaseError,\
                    'const name "%s" is not all uppercase!' % name
            self.__dict__[name] = value

    sys.modules[__name__] = _const()


def merge_mapping_rate(outputdir):
    rates=os.listdir("{0}/results/mapping_rate/".format(outputdir))
    temps = [i.rstrip().split('\t') for i in open("{0}/results/mapping_rate/{1}".format(outputdir, rates[0]))]
    prefix = rates[0].rstrip("Log.final.out")
    _total_mapping = [['',prefix]]+temps[5:]

    for _file in rates[1:]:
        temps = [i.rstrip().split('\t') for i in open("{0}/results/mapping_rate/{1}".format(outputdir, _file))]
        prefix = _file.rstrip("Log.final.out")
        for item in temps:
            if len(item)==1:
                item.append('')
        temp2 = [i[1] for i in temps]
        temp2 = temp2[5:]
        for i in range(len(_total_mapping)):
            if i == 0:
                _total_mapping[i].append(prefix)
            else:
                _total_mapping[i].append(temp2[i-1])
    _mapping_out = ['\t'.join(i) for i in _total_mapping]
    with open("{0}/results/mapping_rate/mapping_rate_summary.txt".format(outputdir),'w') as f:
        for line in _mapping_out:
            print >>f, line


def get_report(Outputdir):

    if not os.path.exists('{0}/final/'.format(Outputdir)):
        sh('mkdir {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
    else:
        sh('rm -rf {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
        sh('mkdir {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
    sh('mkdir {0}/final/phase1-AllExpGenes {0}/final/phase2-DiffExpGenes\
       {0}/final/phase3-GO_KEGG {0}/final/phase4-GSEA {0}/final/phase5-SignalNet'.format(Outputdir))
    sh('cp {0}/star/*.final.out {0}/results/mapping_rate/'.format(Outputdir))
    merge_mapping_rate(Outputdir)

    sh('cp {0}/results/treat_and_control_expression.txt {0}/final/phase1-AllExpGenes/expression_level_all.xls'\
       .format(Outputdir))
    sh('cp {0}/results/genes_up.txt {0}/final/phase2-DiffExpGenes/genes_up.xls'\
       .format(Outputdir))
    sh('cp {0}/results/genes_down.txt {0}/final/phase2-DiffExpGenes/genes_down.xls'\
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.txt {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.xls'\
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.pdf {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.pdf'\
       .format(Outputdir))
    sh('cp -r {0}/results/down/ {0}/final/phase3-GO_KEGG/down/'.format(Outputdir)\
       .format(Outputdir))
    sh('cp -r {0}/results/up/ {0}/final/phase3-GO_KEGG/up/'.format(Outputdir)\
       .format(Outputdir))
    sh('rm {0}/final/phase3-GO_KEGG/*/Rplots.pdf'.format(Outputdir))
    sh('cp -r {0}/results/gsea/bp*/ {0}/final/phase4-GSEA/BP'.format(Outputdir))
    sh('cp -r {0}/results/gsea/kegg*/ {0}/final/phase4-GSEA/KEGG'.format(Outputdir))
    sh('cp {0}/results/cytoscape/* {0}/final/phase5-SignalNet/'.format(Outputdir))


def get_ciri_report(Outputdir):
    if not os.path.exists('{0}/final/'.format(Outputdir)):
        sh('mkdir {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
    else:
        sh('rm -rf {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
        sh('mkdir {0}/final/ {0}/results/mapping_rate/'.format(Outputdir))
    sh('mkdir {0}/final/phase1-AllExpGenes {0}/final/phase2-DiffExpGenes\
       {0}/final/phase3-GO_KEGG {0}/final/phase4-GSEA {0}/final/phase5-SignalNet'.format(Outputdir))
    sh('cp {0}/star/*.final.out {0}/results/mapping_rate/'.format(Outputdir))
    # merge_mapping_rate(Outputdir)

    sh('cp {0}/results/treat_and_control_expression.txt {0}/final/phase1-AllExpGenes/expression_level_all.xls' \
       .format(Outputdir))
    sh('cp {0}/results/genes_up.txt {0}/final/phase2-DiffExpGenes/genes_up.xls' \
       .format(Outputdir))
    sh('cp {0}/results/genes_down.txt {0}/final/phase2-DiffExpGenes/genes_down.xls' \
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.txt {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.xls' \
       .format(Outputdir))
    sh('cp {0}/results/Treat_vs_control_diff.pdf {0}/final/phase2-DiffExpGenes/Treat_vs_control_diff.pdf' \
       .format(Outputdir))
    sh('cp -r {0}/results/down/ {0}/final/phase3-GO_KEGG/down/'.format(Outputdir) \
       .format(Outputdir))
    sh('cp -r {0}/results/up/ {0}/final/phase3-GO_KEGG/up/'.format(Outputdir) \
       .format(Outputdir))
    sh('rm {0}/final/phase3-GO_KEGG/*/Rplots.pdf'.format(Outputdir))
    sh('cp -r {0}/results/gsea/bp*/ {0}/final/phase4-GSEA/BP'.format(Outputdir))
    sh('cp -r {0}/results/gsea/kegg*/ {0}/final/phase4-GSEA/KEGG'.format(Outputdir))
    sh('cp {0}/results/cytoscape/* {0}/final/phase5-SignalNet/'.format(Outputdir))
    sh("mkdir {0}/final/phase6-circRNA".format(Outputdir))
    sh("cp {0}/results/circ_up_RBP.txt {0}/final/phase6-circRNA/up_circRNA_with_RBP.txt".format(Outputdir))
    sh("cp {0}/results/circ_down_RBP.txt {0}/final/phase6-circRNA/down_circRNA_with_RBP.txt".format(Outputdir))
    sh("cp {0}/results/annotation_table.txt {0}/final/phase6-circRNA/circRNA_summary.xls".format(Outputdir))
    sh("cp {0}/results/annotation_table_circbase.txt {0}/final/phase6-circRNA/circRNA_circbase.xls".format(Outputdir))
    sh("cp {0}/results/targetscan/circ_up_high.txt {0}/final/phase6-circRNA/down_circRNA_miRNA_top.xls".format(Outputdir))
    sh("cp {0}/results/targetscan/circ_down_high.txt {0}/final/phase6-circRNA/up_circRNA_miRNA_top.xls".format(Outputdir))
    sh("cp {0}/results/targetscan/circ_up_pass.txt {0}/final/phase6-circRNA/down_circRNA_miRNA_okay.xls".format(
        Outputdir))
    sh("cp {0}/results/targetscan/circ_down_pass.txt {0}/final/phase6-circRNA/up_circRNA_miRNA_okay.xls".format(
        Outputdir))


def configparser(config_file = './config.txt'):
    _file = [i.strip() for i in open(config_file)]
    configs = {}
    for _line in _file:
        if ":" in _line[:40]:
            _key = _line.split(":")[0]
            _value = _line.split(":")[1]
            configs[_key]=_value
    treatment = configs["Treatment"].split(",")
    control = configs["Control"].split(",")
    Starindex = configs["STARindex"]
    STAR_PATH = configs["STARpath"]
    Starindexdir = configs["STARindexdir"]
    Genomefile = configs["Genomefile"]
    Annotationfile = configs["Annotationfile"]
    Outputdir = configs["Outputdir"]
    Readlength = configs["Readlength"]
    Max_process = configs["Max_parrallel_process"]
    Seqtype = configs["Seqtype"].split(",")
    Pair_rep=configs["Paired_rep"]
    Run_type=configs["Runtype"]
    Genome=configs["Genome"].split(",")
    Input_check=configs["Inputcheck"]
    Htseq_PATH=configs["Htseq_path"]
    Pvalue=configs["pvalue"]
    if Max_process.lower()=='default':
        Max_process=len(control)+len(treatment)
    if Genome[1]=='default':
        if Genome[0].lower()=='hg19':
            Annotationfile = "/home/Public/igenome/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf"
            Genomefile = "/home/Public/igenome/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa"
        elif Genome[0].lower()=='hg38':
            Annotationfile = "/home/Public/igenome/Homo_sapiens/UCSC/hg38/Annotation/Genes/genes.gtf"
            Genomefile = "/home/Public/igenome/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa"
        elif Genome[0].lower()=='mm9':
            Annotationfile = "/home/Public/igenome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf"
            Genomefile = "/home/Public/igenome/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa"
        elif Genome[0].lower()=='mm10':
            Annotationfile = "/home/Public/igenome/Mus_musculus/UCSC/mm10/Annotation/Genes/genes.gtf"
            Genomefile = "/home/Public/igenome/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa"
        else:
            print "We do not have the species yet. Please take 'user' choice. "
    return treatment, control, Starindex, STAR_PATH, Starindexdir, Genomefile, Annotationfile, Outputdir,\
    Readlength, Max_process, Seqtype, Pair_rep,Run_type,Genome,Input_check, Htseq_PATH,Pvalue


def input_check(*files):
    _file0 = files[0]+files[1]
    if len(_file0)>len(set(_file0)):
            print "Same files in input.  Pleas check treat and control filenames!"
            sys.exit(1)

    for _file in files:
        for each in _file:
            if not os.path.isfile(each):
                print '{0} do not exist!'.format(each)
                sys.exit(1)


def star_main(star_path,star_index, starindexdir, genomefile,annotationfile,readlength,\
              treatments,controls,outputdir, seqtype):

    print "outputdir",outputdir
    if not os.path.exists(outputdir):
        sh('mkdir {0}'.format(outputdir))
    if not os.path.exists('{0}/star/'.format(outputdir)):
        sh('mkdir {0}/star/'.format(outputdir))
    if star_index.lower() =='false':
        if not os.path.exists(starindexdir):
            sh('mkdir {0}'.format(starindexdir))
        starindex(star_path, starindexdir, genomefile, annotationfile, readlength)
    else:
        if not os.path.exists(starindexdir):
            print "cannot find your star index! exit!"
            sys.exit(1)

    if seqtype[0].lower()=='pairend' or seqtype[0].lower()=='pair_end':
        for index in range(0,len(treatments),2):
            treatment = treatments[index:index+2]
            outprefix = outputdir+"/star/treat{0}".format(str((index/2)+1))
            star_run(star_path, starindexdir,outprefix, seqtype,treatment[0],treatment[1])

        for index in range(0,len(controls),2):
            control = controls[index:index+2]
            outprefix = outputdir+"/star/control{0}".format(str((index/2)+1))
            star_run(star_path, starindexdir,outprefix, seqtype,control[0],control[1])
    elif seqtype[0].lower()=='singleend' or seqtype[0].lower()=='single_end':
        for index in range(len(treatments)):
            treatment = treatments[index]
            outprefix = outputdir+"/star/treat{0}".format(str(index+1))
            star_run(star_path, starindexdir,outprefix, seqtype,treatment)

        for index in range(len(controls)):
            control = controls[index]
            outprefix = outputdir+"/star/control{0}".format(str(index+1))
            star_run(star_path, starindexdir,outprefix, seqtype,control)
    else:
        print "Please give your sequencing type!"
        sys.exit(1)


def htseq_main(outputdir, Annotationfile, group, sampletype, stranded, htseq_path):
    if not os.path.exists('{0}/expr'.format(outputdir)):
        sh('mkdir {0}/expr'.format(outputdir))
    sh("samtools sort -n {0}/star/{2}{1}Aligned.out.bam -o {0}/star/{2}{1}.sorted.bam  "\
       .format(outputdir,group,sampletype))
    # sh("samtools index {0}/star/{2}{1}.sorted.bam".format(outputdir,group, sampletype))
    if stranded.lower()=='strand-specific':
        sh("{4} -i gene_id {0}/star/{3}{2}.sorted.bam {1} -r name -s yes -f bam -q \
        >{0}/expr/{3}{2}_transcript.txt".format(outputdir, Annotationfile, group, sampletype,htseq_path))
    elif stranded.lower()=='non-strand-specific':
        sh("htseq-count -i gene_id {0}/star/{3}{2}.sorted.bam {1} -r name -s no -f bam -q \
        >{0}/expr/{3}{2}_transcript.txt".format(outputdir, Annotationfile, group, sampletype,htseq_path))
    else:
        print "We are running with not strand-specific library!!"
        sh("htseq-count -i gene_id {0}/star/{3}{2}.sorted.bam {1} -r name -s no -f bam -q \
        >{0}/expr/{3}{2}_transcript.txt".format(outputdir, Annotationfile, group, sampletype,htseq_path))


def ciri_main(outputdir, genomefile,annotationfile,treatments, controls):
    print "outputdir:",outputdir
    if not os.path.exists(outputdir):
        sh('mkdir {0}'.format(outputdir))
    if not os.path.exists('{0}/ciri_out/'.format(outputdir)):
        sh('mkdir {0}/ciri_out/'.format(outputdir))
    pool=Pool(processes=4)
    for index in range(0,len(treatments),2):
        treatment = treatments[index:index+2]
        outprefix = outputdir+"/ciri_out/treat{0}".format(str((index/2)+1))
        ciri_run(outprefix, genomefile,annotationfile,treatment[0],treatment[1])

    for index in range(0,len(controls),2):
        control = controls[index:index+2]
        outprefix = outputdir+"/ciri_out/control{0}".format(str((index/2)+1))
        ciri_run(outprefix, genomefile,annotationfile,control[0],control[1])


if __name__=='__main__':
    const.TITLE = 'RNAseq pipeline'
    const.AUTHOR = os.popen("whoami").read().rstrip()
    const.TEMPLATEDATE = 'April 7, 2016'
    const.CODEFUNCTION = 'This program is for rnaseq.'
    begin_time = time.clock()
    print "\033[1;31;38m"
    print "#" * 50
    print 'User ' + const.AUTHOR + ':'
    print const.TITLE + ' runs at ' + \
        time.strftime("%Y-%m-%d %X", time.localtime()) + '...'
    print 'Pipeline made: ' + const.TEMPLATEDATE
    print const.CODEFUNCTION
    print "#" * 50
    print "\033[0m"


    # Parsing argvs and configs
    if len(sys.argv)==2:
        config_file = sys.argv[1]
    else:
        config_file = False

    if config_file:
        treatments, controls, Starindex, STAR_PATH, Starindexdir, Genomefile, Annotationfile, Outputdir,\
        Readlength, Max_process, Seqtype, Pair_rep,Run_type,Genome,Input_check, Htseq_PATH, P_value = configparser(config_file)
    else:
        try:
            f = open('config.txt')
        except IOError:
            print "Please give config files"
            sys.exit(1)

        treatments, controls, Starindex, STAR_PATH, Starindexdir, Genomefile, Annotationfile, Outputdir,\
        Readlength, Max_process, Seqtype, Pair_rep,Run_type,Genome,Input_check, Htseq_PATH, P_value = configparser()

    if Input_check.lower()=="false":
        print "Skipping checking of input files"
    else:
        input_check(treatments, controls)

    # Mapping


    if Run_type.lower()=='full' or Run_type.lower()=='mapping':
        star_main(STAR_PATH,Starindex, Starindexdir, Genomefile, Annotationfile, Readlength,treatments
                  , controls, Outputdir, Seqtype)
        pool = Pool(processes=int(Max_process))
        result=[]
        if Seqtype[0].lower()=='pairend' or Seqtype[0].lower()=='pair_end':
            for group in range(1,(len(treatments)/2)+1):
                result.append(pool.apply_async(htseq_main,(Outputdir, Annotationfile,group,"treat",Seqtype[1],Htseq_PATH)))
            for group in range(1,(len(controls)/2)+1):
                result.append(pool.apply_async(htseq_main,(Outputdir, Annotationfile,group,"control",Seqtype[1],Htseq_PATH)))
        elif Seqtype[0].lower()=='singleend' or Seqtype[0].lower()=='single_end':
            for group in range(1,len(treatments)+1):
                result.append(pool.apply_async(htseq_main,(Outputdir, Annotationfile,group,"treat",Seqtype[1],Htseq_PATH)))
            for group in range(1,len(controls)+1):
                result.append(pool.apply_async(htseq_main,(Outputdir, Annotationfile,group,"control",Seqtype[1],Htseq_PATH)))
        else:
            print("Please give your sequencing type!")
            sys.exit(1)

        pool.close()
        pool.join()

    if Run_type == "circRNA_full" or Run_type == "circRNA_mapping":
        print "circRNA mode is only for testing pair-end data! only fastq input is allowed!"
        ciri_main(Outputdir,Genomefile, Annotationfile,treatments, controls)
    if Run_type == "circRNA_full" or Run_type == "circRNA_process":
        nrep = len(treatments)/2
        ciri_process(Outputdir,nrep,Pair_rep,P_value, Genomefile)

    # DE analysis
    if Run_type.lower()=='full' or Run_type.lower()=='de_full' or Run_type.lower()=='de_only':
        expr(Outputdir, controls, Pair_rep,Genome[0], P_value, Run_type, Seqtype[0])
        pass


    # collect report only after finish mouse or human project
    if Run_type.lower()=='full' or Run_type.lower()=='de_full':
        get_report(Outputdir)
    if Run_type == "circRNA_full" or Run_type == "circRNA_process":
        # get_report(Outputdir)
        get_ciri_report(Outputdir)


    end_time = time.clock()
    lapsed_time = end_time - begin_time

    print "\033[1;31;38m"
    print "#" * 50
    print const.TITLE + ' ends at ' + \
        time.strftime("%Y-%m-%d %X", time.localtime()) + '...'
    print "Totally \t%.03f seconds lapsed" % lapsed_time
    print "#" * 50
    print "\033[0m"