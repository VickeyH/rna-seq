#!/usr/bin/env python
#coding:utf-8
### 比对步骤
## cd /mnt/icfs/work/yongdeng/Project/161565_test && python /lustre/work/yongdeng/software/protokaryon/flow/coexpression.py -eq analysis/assemble/gene.expression.txt -vs IB-VS-NT IR-VS-NT -deg analysis/deg_hisat -config config_illuqc.txt -o analysis/coexpression -mf flow/coexpression.mf
import argparse,os,sys,glob
import RNA

parser = argparse.ArgumentParser(description="This Script is used to align cleandata to reference genome by bowtie2 program.")
parser.add_argument("-deg","--deg_dir",type=str,help="the deg output dir, 'vs' directory must under this dir",required = True)
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB",nargs="*",required = True)
parser.add_argument("-o","--output_dir",type=str,help="The output coexpression directory",required = True)
parser.add_argument("-db","--db_file",type=str,help="The pro dbfile, eg.: pro_pro_interaction.norep.db",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
            
def write_propro_makeflow(makeflow,vslist,degdir,outputdir,prodb):
    fo = open(makeflow,"a+")
    deg_pro_pro_interaction = RNA.get_bin_abspath("deg_pro-pro_interaction.pl")
    coexp_propro_forCytoscape = RNA.get_bin_abspath("filter.coexp.propro.forCytoscape.pl")
    fo.write("CATEGORY=get_diff_info\n")
    fo.write(" ".join([os.path.join(degdir,vs,vs+".diff_info") for vs in vslist]) + " : " + " ".join([os.path.join(degdir,vs,i) for vs in vslist for i in ["deg.up.anno.txt","deg.down.anno.txt"]]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/get_diff_info.py -vsdir %s &> %s\n\n"%(" ".join([os.path.join(degdir,vs) for vs in vslist]),os.path.join(degdir,"diff_info.log")))
    fo.write("CATEGORY=pro_pro_interaction\n")
    fo.write(" ".join([os.path.join(outputdir,vs,"deg.pro-pro.txt") for vs in vslist]) + " : " + " ".join([os.path.join(degdir,vs,vs+".diff_info") for vs in vslist]) + "\n")
    fo.write("\tperl %s -deg %s -pp %s -od %s &> %s\n\n"%(deg_pro_pro_interaction,degdir,prodb,outputdir,os.path.join(outputdir,"deg.propro.log")))
    fo.write("CATEGORY=pro_pro_interaction\n")
    fo.write(" ".join([os.path.join(outputdir,vs,"filter.propro.txt") for vs in vslist]) + " : " + " ".join([os.path.join(outputdir,vs,"deg.pro-pro.txt") for vs in vslist]) + "\n")
    fo.write("\tperl %s -input %s,deg.pro-pro.txt -name filter.propro.txt -type propro &> %s\n\n"%(coexp_propro_forCytoscape,outputdir,os.path.join(outputdir,"filter.log")))
    fo.close()
    
def main():
    vslist = args.vs
    degdir = os.path.abspath(args.deg_dir)
    outputdir = os.path.abspath(args.output_dir)
    prodb = os.path.abspath(args.db_file)
    mkdir(outputdir,vslist)
    write_propro_makeflow(args.makeflow,vslist,degdir,outputdir,prodb)
    
if __name__ == "__main__":
    main()
