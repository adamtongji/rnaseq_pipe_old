#!/usr/bin/env python
#coding:utf-8
import os
import sys

sh = os.system


def further(outputdir, genome):
    if not os.path.exists('{0}/results/gsea/'.format(outputdir)):
        sh('mkdir {0}/results/gsea/'.format(outputdir))
    else:
        sh('rm -rf {0}/results/gsea/'.format(outputdir))
        sh('mkdir {0}/results/gsea/'.format(outputdir))

    gsea_file = [i.rstrip().split('\t') for i in open('{0}/results/Treat_vs_control_diff.txt'.format(outputdir))]
    gsea_file = [i[:1]+i[4:5] for i in gsea_file]
    gsea_file = gsea_file[1:]
    gsea_file2 = []
    for _line in gsea_file:
        if _line[1]=='NA':
            pass
        else:
            _line[0]=_line[0].upper()
            gsea_file2.append(_line)
    gsea_file = ['\t'.join(i) for i in gsea_file2]
    with open('{0}/results/gsea_input.txt'.format(outputdir),'w') as f:
        for line in gsea_file:
            print >>f, line

    sh("sort -k 2gr {0}/results/gsea_input.txt >{0}/results/gsea/gsea_input.rnk".format(outputdir))
    sh("export LANG=en_US.UTF-8;java -cp /usr/local/src/gsea2-2.2.2.jar -Xmx2g xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c5.bp.v5.1.symbols.gmt\
     -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk {0}/results/gsea/gsea_input.rnk -scoring_scheme weighted\
      -rpt_label bp -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500\
       -set_min 15 -zip_report false -out {0}/results/gsea/ -gui false".format(outputdir))
    sh("export LANG=en_US.UTF-8; java -cp /usr/local/src/gsea2-2.2.2.jar -Xmx2g xtools.gsea.GseaPreranked -gmx gseaftp.broadinstitute.org://pub/gsea/gene_sets/c2.cp.kegg.v5.1.symbols.gmt\
     -collapse false -mode Max_probe -norm meandiv -nperm 1000 -rnk {0}/results/gsea/gsea_input.rnk -scoring_scheme weighted\
      -rpt_label kegg -include_only_symbols true -make_sets true -plot_top_x 20 -rnd_seed timestamp -set_max 500\
     -set_min 15 -zip_report false -out {0}/results/gsea/ -gui false".format(outputdir))

    # step 4 cytoscape
    if not os.path.exists('{0}/results/cytoscape/'.format(outputdir)):
        sh('mkdir {0}/results/cytoscape'.format(outputdir))
    if genome.lower()=='hg19' or genome.lower()=='hg38':
        sh('/home/Public/software/rnaseq2report/lib/DEG2network.py -p 0.05 -n 5 -k\
         /home/Public/software/rnaseq2report/lib/merged_KEGG.txt -i {0}/results/Treat_vs_control_diff.txt -d \
           {0}/results/cytoscape'.format(outputdir))
    elif genome.lower()=='mm9' or genome.lower()=='mm10':
        sh('/home/Public/software/rnaseq2report/lib/DEG2network_mouse.py -p 0.05 -n 5 -k\
         /home/Public/software/rnaseq2report/lib/mouse_merged_KEGG.txt -i {0}/results/Treat_vs_control_diff.txt -d \
           {0}/results/cytoscape'.format(outputdir))


def expr(outputdir, control, pairtype,genome, pvalue, runtype,seqtype):
     # step2 DEseq2
    control_names= [i.split('.')[0] for i in control]
    if not os.path.exists('{0}/results'.format(outputdir)):
        sh("mkdir {0}/results".format(outputdir))
    if seqtype.lower()=='pairend' or seqtype.lower()=='pair_end':
        sh("Rscript /home/Public/software/rnaseq2report/lib/deseq2.r {0}/expr/ {1} {2} {3}".format(outputdir, len(control_names)/2, pairtype, pvalue))
    elif seqtype.lower()=='singleend' or seqtype.lower()=='single_end':
        sh("Rscript /home/Public/software/rnaseq2report/lib/deseq2.r {0}/expr/ {1} {2} {3}"\
           .format(outputdir, len(control_names), pairtype, pvalue))
    else:
        print("Please give your sequencing type!")
        sys.exit(1)

    sh('mkdir {0}/results/up;mkdir {0}/results/down'.format(outputdir))
    #
    if genome.lower()=='hg19' or genome.lower()=='hg38':
        sh("Rscript /home/Public/software/rnaseq2report/lib/pathway.r {0}/results/".format(outputdir))
    if genome.lower()=='mm9' or genome.lower()=='mm10':
        sh("Rscript /home/Public/software/rnaseq2report/lib/mouse_pathway.r {0}/results".format(outputdir))

    if runtype=='de_full' or runtype=='full':
        further(outputdir,genome)
