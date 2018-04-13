#!/usr/bin/env python
#coding:utf-8
import argparse,os,sys,re
import RNA

parser = argparse.ArgumentParser(description="This Script is used to align cleandata to reference genome by bowtie2 program.")
parser.add_argument("-eq","--expression",type=str,help="The input expression table of all sampls",required = True)
parser.add_argument("-deg","--deg_dir",type=str,help="the deg output dir, 'vs' directory must under this dir",required = True)
parser.add_argument("-sg","--sample_group",type=str,help="The sample info of this project, like groupA=sample1,sample2,... groupB=sample3,sample4,...",required = True,nargs="+")
parser.add_argument("-anno","--anno",type=str,help="The track-ref anno file")
parser.add_argument("-o","--output_dir",type=str,help="The output coexpression directory",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
            
def get_group_dict(sampleinfo):
    groupDict = {}
    for i in sampleinfo:
        k = i.split("=")[0]
        v = i.split("=")[-1].split(",")
        groupDict[k] = v
    return groupDict 

def getanno(anno):
    d = {}
    if anno == None:return d
    with open(anno) as fi:
        for line in fi:
         k = line.split("\t")[0]
         v = line.split("\t")[2]
         d[k] = k + "(" + v + ")"
    return d

def get_expression_dict(expression):
    expressionDict = {}
    with open(expression) as ex:
        head = ex.next()
        sample = head.strip().split("\t")[1:]
        for line in ex:
            k = line.split("\t")[0]
            v = dict(zip(sample,line.strip().split("\t")[1:]))
            expressionDict[k] = v
    return expressionDict
    
sampleinfo = args.sample_group
groupDict = get_group_dict(sampleinfo)

vs_l =[i for i in next(os.walk(args.deg_dir))[1] if "-VS-" in i]
vs_list = []
for i in vs_l:
    group1,group2 = i.split("-VS-")
    if group1 in groupDict.keys() and group2 in groupDict.keys() and len(groupDict[group1]) + len(groupDict[group2]) >= 3:
        vs_list.append(i)

mkdir(args.output_dir,vs_list)
expressionDict = get_expression_dict(args.expression)
anno = getanno(args.anno)
if vs_list:
    for vs in vs_list:
        group1,group2 = vs.split("-VS-")
        degfile = os.path.join(args.deg_dir,vs,"deg.up_down.txt")
        outfile = os.path.join(args.output_dir,vs,"deg.exp.list")
        with open(degfile) as df,open(outfile,"w") as fo:
            header = df.next()
            fo.write("ID\t" + "\t".join(sorted(groupDict[group1])) + "\t" + "\t".join(sorted(groupDict[group2])) + "\n")
            for line in df:
                #if line.strip().endswith("no"):
                #    continue
                #else:
                    geneid = line.split("\t")[0]
                    fo.write(anno.get(geneid,geneid) + "\t" + "\t".join([expressionDict[geneid][i] for i in sorted(groupDict[group1])+sorted(groupDict[group2])]) + "\n")
                
    all_group = groupDict.keys()
    all_sample = sorted(set(("\t".join(["\t".join(sorted(groupDict[i])) for i in all_group])).split("\t")))
    all_deg = set(os.popen("cut -f1 %s"%os.path.join(args.output_dir,"*","deg.exp.list")).read().strip("\n").split("\n"))
    with open(os.path.join(args.output_dir,"deg.exp.list"),"w") as all_exp:
        all_exp.write("ID\t" + "\t".join(all_sample) + "\n" )
        for geneid in all_deg:
            if geneid in expressionDict:
                all_exp.write(anno.get(geneid,geneid) + "\t" + "\t".join([expressionDict[geneid][i] for i in all_sample]) + "\n")
else:
    all_group = groupDict.keys()
    all_sample = sorted(set(("\t".join(["\t".join(sorted(groupDict[i])) for i in all_group])).split("\t")))
    all_deg = set(os.popen('cut -f1 %s'%os.path.join(args.deg_dir,"*","deg.up_down.txt")).read().strip("\n").split("\n"))
    with open(os.path.join(args.output_dir,"deg.exp.list"),"w") as all_exp:
        all_exp.write("ID\t" + "\t".join(all_sample) + "\n" )
        for geneid in all_deg:
            if geneid in expressionDict:
                all_exp.write(anno.get(geneid,geneid) + "\t" + "\t".join([expressionDict[geneid][i] for i in all_sample]) + "\n")
    
    
