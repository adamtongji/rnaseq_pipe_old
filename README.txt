RNAseq pipeline用法（任意路径）：
rainbowRNA your_config.txt

跑错了停止不了：
pkill -u abc   (abc为你的用户名）

关于Configs：
必定修改项：
Treatment， Control, Outputdir
必须检查项：
Seqtype ,Pair_rep, Read_length, Genome

单独只跑star和htseq，将Runtype改为mapping。
STAR和htseq跑完并确定正确后，修改重跑后续的，将差异表达分析重跑，把Runtype改为de_full。如果不是人类或者小鼠，则改为de_only。
Inputcheck是检查输入fastq文件是否存在，是否重复，建议设定为True。


输出结果：
直接的输出结果在results文件夹里，其中mapping率在results/mapping_rate里，用excel打开即可。
给客户的结果在final里的5个phase文件夹，里面图片可以用来填写报告模板。
cytoscape的结果请手动生成。
所有输出表格都已改成xls。源txt文件都在results文件夹都能够找到。
重点检查：mapping rate是否大于80%， heatmap是否大小和聚类合理， 差异基因数量是否合理， 结果中GO和KEGG结果是否和预期相符。
生成的中文几月几日文件夹为gsea自动生成空文件夹，记录你今天干了什么坏事，请手动删除。

