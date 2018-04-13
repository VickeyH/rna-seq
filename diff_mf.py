#!/usr/bin/env python
#coding:utf-8

import os,argparse,re
#from pathlib2 import Path
from os.path import join,abspath

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is used to do DEG analysis")
    parser.add_argument("-c","--count",type=str,help="The read count file",required=True)
    parser.add_argument("-sg","--sample_group",type=str,help="The all samples and groups rep. like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
    parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB .... if paired group, the vs must be like 'CasegroupA-VS-ControlgroupB,paired' and deg method will be limma", required = True,nargs="+")
    parser.add_argument("-o","--output_dir",type=str,help="The DEG output directory",required = True)
    parser.add_argument("-anno","--anno",type=str,help="the track_ref anno file",required = True)
    parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def write_deg_makeflow(makeflow,count,sg,vslist,anno,outputdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=do_deg\n")
    fo.write(" ".join([join(outputdir,vs,vs+".txt") for vs in [v.split(",")[0] for v in vslist]]) + " " + join(outputdir,"gene_count_matrix.normalized.txt : ") + count + "\n")
    fo.write('\tpython /lustre/work/yongdeng/software/protokaryon/flow/get_deg.py -c %s -sg %s -o %s -vs %s &> %s\n\n'%(count," ".join(sg),outputdir," ".join(vslist),join(outputdir,"get_deg.log")))
    vslist = [i.split(",")[0] for i in vslist]
    fo.write("CATEGORY=parse_deg\n")
    fo.write(" ".join([join(outputdir,vs,vs + i) for vs in vslist for i in [".diff.stat",".plot_valcano.file"]]) + " " + " ".join([join(outputdir,vs,i) for vs in vslist for i in ["deg.down.txt","deg.num.txt","deg.up_down.txt","deg.up.txt"]]) + " " + join(outputdir,"deg.num.txt : ") + " ".join([join(outputdir,vs,vs+".txt") for vs in vslist]) + "\n" )
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/parse_deg.py -i %s -c %s &> %s\n\n"%(" ".join([join(outputdir,vs,vs+".txt") for vs in vslist]),count,join(outputdir,"parse_deg.log")))
    fo.write("CATEGORY=anno_deg\n")
    fo.write(" ".join([join(outputdir,vs,vs +".anno.txt") for vs in vslist]) + " " + " ".join([join(outputdir,vs,"deg.down.anno.txt") for vs in vslist]) + " " +  " ".join([join(outputdir,vs,"deg.up.anno.txt") for vs in vslist])  + " : " +  " ".join([join(outputdir,vs,vs +".txt") for vs in vslist]) + " " + " ".join([join(outputdir,vs,"deg.down.txt") for vs in vslist]) + " " +  " ".join([join(outputdir,vs,"deg.up.txt") for vs in vslist]) + " " + anno + " " + join(outputdir,"gene_count_matrix.normalized.txt") + " \n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/anno_deg.py -i %s -sg %s -nc %s -anno %s &> %s\n\n"%(" ".join([join(outputdir,vs,vs +".txt") for vs in vslist]) + " " + " ".join([join(outputdir,vs,"deg.down.txt") for vs in vslist]) + " " +  " ".join([join(outputdir,vs,"deg.up.txt") for vs in vslist])," ".join(sg),join(outputdir,"gene_count_matrix.normalized.txt"),anno,join(outputdir,"anno_deg.log")))
    fo.write("CATEGORY=get_deg_list\n")
    fo.write(" ".join([join(outputdir,vs,"deg.up_down.list") for vs in vslist]) + " : " + " ".join([join(outputdir,vs,"deg.up_down.txt") for vs in vslist]) + " " + anno + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/get_deg_list.py %s %s &> %s\n\n"%(" ".join([join(outputdir,vs,"deg.up_down.txt") for vs in vslist]),anno,join(outputdir,"get_deg_list.log")))
    fo.write("CATEGORY=deg_stat\n")
    fo.write(join(outputdir,"deg.num.png ") + join(outputdir,"deg.num.pdf : ") + join(outputdir,"deg.num.txt\n"))
    fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/flow/deg_stat.r %s %s %s &> %s\n\n"%(join(outputdir,"deg.num.txt"),join(outputdir,"deg.num.png"),join(outputdir,"deg.num.pdf"),join(outputdir,"deg.num.log")))

    for vs in vslist:
        fo.write("CATEGORY=venn_%s\n"%vs)
        fo.write(join(outputdir,vs,vs + ".venn.png ") + join(outputdir,vs,vs + ".venn.pdf : ") + join(join(outputdir,vs,vs+ ".diff.stat\n")))
        fo.write("\tperl /lustre/work/zhonghuali/software/rna.ref/bin/cuffdiff/venn.expression.pl %s %s %s &> %s\n\n"%(join(outputdir,vs,vs+ ".diff.stat"),join(outputdir,vs,vs + ".venn.png"),join(outputdir,vs,vs + ".venn.pdf"),join(outputdir,vs,vs + ".venn.log")))
        fo.write("CATEGORY=valcano_%s\n"%vs)
        fo.write(join(outputdir,vs,vs+".valcano.png : ") + join(outputdir,vs,vs+".plot_valcano.file\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/flow/Scatter_volcano.r %s %s &> %s\n\n"%(join(outputdir,vs,vs+".plot_valcano.file"),join(outputdir,vs,vs+".valcano.png"),join(outputdir,vs,"valcano.log")))
    
    fo.write("CATEGORY=venn_stat_draw\n")   
    fo.write(join(outputdir,"venn_DEG","stat.txt ") + join(outputdir,"venn_DEG","list.txt : ") + " ".join([join(outputdir,vs,"deg.up_down.txt") for vs in vslist]) + "\n")
    fo.write("\tperl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/venn/venn_stat_draw.pl -DC -sample_file %s -sample_file_head -od %s &> %s\n\n"%(" ".join([vs.replace("-","_") + "," + join(outputdir,vs,"deg.up_down.txt") for vs in vslist]),join(outputdir,"venn_DEG"),join(outputdir,"venn_DEG.log")))

    fo.write("CATEGORY=venn_stat_draw\n")   
    fo.write(join(outputdir,"venn_UP","stat.txt ") + join(outputdir,"venn_UP","list.txt : ") + " ".join([join(outputdir,vs,"deg.up.txt") for vs in vslist]) + "\n")
    fo.write("\tperl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/venn/venn_stat_draw.pl -DC -sample_file %s -sample_file_head -od %s &> %s\n\n"%(" ".join([vs.replace("-","_") + "," + join(outputdir,vs,"deg.up.txt") for vs in vslist]),join(outputdir,"venn_UP"),join(outputdir,"venn_UP.log")))

    fo.write("CATEGORY=venn_stat_draw\n")   
    fo.write(join(outputdir,"venn_DOWN","stat.txt ") + join(outputdir,"venn_DOWN","list.txt : ") + " ".join([join(outputdir,vs,"deg.down.txt") for vs in vslist]) + "\n")
    fo.write("\tperl /lustre/database/rnaseq/program/workflow/lncRNA/stringtie/venn/venn_stat_draw.pl -DC -sample_file %s -sample_file_head -od %s &> %s\n\n"%(" ".join([vs.replace("-","_") + "," + join(outputdir,vs,"deg.down.txt") for vs in vslist]),join(outputdir,"venn_DOWN"),join(outputdir,"venn_DOWN.log")))
    fo.close()
 
def main():
    args = parseArg()
    count = abspath(args.count)
    anno = abspath(args.anno)
    out = abspath(args.output_dir)
    mkdir(out,[i.split(",")[0] for i in args.vs])
    write_deg_makeflow(args.makeflow,count,args.sample_group,args.vs,anno,out)

if __name__ == "__main__":
    main()

