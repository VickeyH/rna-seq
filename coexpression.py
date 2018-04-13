#!/usr/bin/env python
#coding:utf-8
### 比对步骤
## cd /lustre/work/yongdeng/Project && python ~/align_clean.py -i 160761_test/QC/ -g 160761_test/genome/ -mf 160761_test/align_clean.mf -o 160761_test/align
import argparse,os,sys,glob
import RNA

parser = argparse.ArgumentParser(description="This Script is used to align cleandata to reference genome by bowtie2 program.")
parser.add_argument("-eq","--expression",type=str,help="The input expression table of all sampls",required = True)
parser.add_argument("-deg","--deg_dir",type=str,help="the deg output dir, 'vs' directory must under this dir",required = True)
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB",nargs="*",required = True)
parser.add_argument("-sg","--sample_group",type=str,help="The sample info of this project, like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
parser.add_argument("-o","--output_dir",type=str,help="The output coexpression directory",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    if sampls:
        for n in sampls:
            d = os.path.join(dir,n)
            if not os.path.exists(d):
                os.makedirs(d)
            else:
                continue
    else:
        if not os.path.isdir(dir):
            os.makedirs(dir)
            
def get_group_dict(sampleinfo):
    groupDict = {}
    for i in sampleinfo:
        k = i.split("=")[0]
        v = i.split("=")[-1].split(",")
        groupDict[k] = v
    return groupDict        
            
def write_coexpression_makeflow(makeflow,vslist,degdir,config,outputdir,expressionmatrix,vs):
    paser_for_coexpression = RNA.get_bin_abspath("paser_for_coexpression.py")
    CoExpression_denovo = RNA.get_bin_abspath("CoExpression_denovo.pl")
    coexp_propro_forCytoscape = RNA.get_bin_abspath("filter.coexp.propro.forCytoscape.pl")
    fo =open(makeflow,"a+")
    fo.write("CATEGORY=pre_coexpression\n")
    fo.write((" ".join([os.path.join(outputdir,v,"deg.exp.list") for v in vslist]) + " " +  os.path.join(outputdir,"deg.exp.list : ") + expressionmatrix + " " + " ".join([os.path.join(degdir,v,"deg.up_down.txt") for v in vs]) + "\n").lstrip())
    fo.write("\tpython %s -deg %s -sg %s -eq %s -o %s &> %s\n\n"%(paser_for_coexpression,degdir,config,expressionmatrix,outputdir,os.path.join(outputdir,"get_deg_exp_list.log")))
    for vs in vslist:
        group1,group2 = vs.split("-VS-")
        if len(groupDict[group1]) + len(groupDict[group2]) <= 2:
            os.system('rm -fr %s'%(os.path.join(outputdir,vs)))
            continue
        fo.write("CATEGORY=coexpression_%s\n"%vs)
        fo.write(os.path.join(outputdir,vs,"Co-Expression.xls : ") + os.path.join(outputdir,vs,"deg.exp.list\n"))
        fo.write("\tperl %s -e %s -o %s &> %s\n\n"%(CoExpression_denovo,os.path.join(outputdir,vs,"deg.exp.list"),os.path.join(outputdir,vs),os.path.join(outputdir,vs,"Coexpression.log")))
        fo.write("CATEGORY=coexpression_%s\n"%vs)
        fo.write(os.path.join(outputdir,vs,"filter.coexpression.txt : ") + os.path.join(outputdir,vs,"Co-Expression.xls\n"))
        fo.write("\tperl %s -input %s -name %s -type coexpression &> %s\n\n"%(coexp_propro_forCytoscape,os.path.join(outputdir,vs,"Co-Expression.xls"),"filter.coexpression.txt",os.path.join(outputdir,vs,"filter.coexp.log")))
    fo.write("CATEGORY=coexpression_alldeg\n")
    fo.write(os.path.join(outputdir,"Co-Expression.xls : ") + os.path.join(outputdir,"deg.exp.list\n"))
    fo.write("\tperl %s -e %s -o %s &> %s\n\n"%(CoExpression_denovo,os.path.join(outputdir,"deg.exp.list"),os.path.join(outputdir),os.path.join(outputdir,"Coexpression.log")))
    fo.write("CATEGORY=coexpression_alldeg\n")
    fo.write(os.path.join(outputdir,"filter.coexpression.txt : ") + os.path.join(outputdir,"Co-Expression.xls\n"))
    fo.write("\tperl %s -input %s -name %s -type coexpression &> %s\n\n"%(coexp_propro_forCytoscape,os.path.join(outputdir,"Co-Expression.xls"),"filter.coexpression.txt",os.path.join(outputdir,"filter.coexp.log")))
    fo.close()
    
def main():
    global groupDict
    groupDict = get_group_dict(args.sample_group)
    vslist = []
    for i in args.vs:
        group1,group2 = i.split("-VS-")
        if group1 in groupDict.keys() and group2 in groupDict.keys() and len(groupDict[group1]) + len(groupDict[group2]) >= 3:
            vslist.append(i)
    outputdir = os.path.abspath(args.output_dir)
    degdir = os.path.abspath(args.deg_dir)
    config = " ".join(args.sample_group)
    expressionmatrix = os.path.abspath(args.expression)
    mkdir(outputdir,vslist)
    write_coexpression_makeflow(args.makeflow,vslist,degdir,config,outputdir,expressionmatrix,args.vs)
 
if __name__ == "__main__":
    main()
    
