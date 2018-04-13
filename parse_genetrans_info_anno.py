#!/usr/bin/env python
#coding:utf-8

import sys,os,argparse,re
from collections import OrderedDict,defaultdict

def parseArg():
    parser = argparse.ArgumentParser(description="This Script is get all gene anno table from total gtf file and anno file")
    parser.add_argument("-gtf","--gtf",type=str,help="The input total gtf file",required = True)
    parser.add_argument("-anno","--anno",type=str,help="The input anno file, second clonum will be used for backlist and kobas id",required = True)
    parser.add_argument("-g","--output_gene",type=str,help="The output all genes table")
    parser.add_argument("-t","--output_trans",type=str,help="The output all transcripts table")
    parser.add_argument("-b","--back",type=str,help="The output back list")
    parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
    return parser.parse_args()
    
def get_gene_table(gtf,anno):
    backlist=[]
    annogenetable = OrderedDict()
    genelist = OrderedDict()       # trackgene_list
    tran_gene_dict = OrderedDict()
    gene_dict = OrderedDict()    ## trackgene:refgene
    gene_info = defaultdict(dict) ## trackgene:trackgene_pos
    trans_info = defaultdict(dict)
    with open(gtf) as fi:
        for line in fi:
            if line.strip() and not line.startswith("#"):
                line = line.strip().split("\t")
                description = line[-1]
                geneid = re.search('gene_id "(.+?)";',description).group(1)
                genelist[geneid] = ""
                gene_info[geneid].setdefault("pos",[]).extend([int(line[3]),int(line[4])])
                gene_info[geneid].setdefault("chr",set()).add(line[0])
                gene_info[geneid].setdefault("strand",set()).add(line[6])
                #if re.search('gene_id ".+?";',description):
                if re.search('ref_gene_id ".+?";',description):
                    gene_dict.setdefault(geneid,set()).add(re.search('ref_gene_id "(.+?)";',description).group(1))
                if re.search('transcript_id ".+?";',description):
                    transid = re.search('transcript_id "(.+?)";',description).group(1)
                    trans_info[transid].setdefault("pos",[]).extend([int(line[3]),int(line[4])])
                    trans_info[transid].setdefault("chr",set()).add(line[0])
                    trans_info[transid].setdefault("strand",set()).add(line[6])
                    tran_gene_dict.setdefault(transid,geneid)

    anno_dict = {}
    with open(anno) as fa:
        for line in fa:
            if line.strip() and not line.startswith("#"):
                k = line.split("\t")[0]
                anno_dict[k] = line.strip().split("\t")
    if len(gene_dict) == 0:                                     ### 表达谱gtf为refgtf的时候
        gene_dict = {k:[k] for k in genelist}
    for k in genelist:
        Gene_info = list(gene_info[k]["chr"])[0] + ":" + list(gene_info[k]["strand"])[0] + ":" + str(min(gene_info[k]["pos"])) + "-" + str(max(gene_info[k]["pos"]))
        if gene_dict.has_key(k):
            Ref_gene = "|".join(list(gene_dict[k]))
            Ref_gene_info = "|".join([anno_dict[g][3] + ":" + list(gene_info[k]["strand"])[0] + ":" + anno_dict[g][4] + "-" + anno_dict[g][5] for g in gene_dict[k]])
            Entrez = "|".join([anno_dict[g][1] for g in gene_dict[k]])
            Uniprot = "|".join([anno_dict[g][2] for g in gene_dict[k]])
            Symbol = "|".join([anno_dict[g][-2] for g in gene_dict[k]])
            Description = "|".join([anno_dict[g][-1] for g in gene_dict[k]])
        else:
            Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description = list("-")*6
        annogenetable[k] = [Gene_info,Ref_gene,Ref_gene_info,Entrez,Uniprot,Symbol,Description]    
    for i in gene_dict:
        for j in gene_dict[i]:
            if anno_dict.has_key(j) and anno_dict[j][1] != "-":
                backlist.append(anno_dict[j][1])
    return annogenetable,tran_gene_dict,trans_info,backlist

def main():
    args = parseArg()
    annoGeneTable,Tran_Gene_Dict,trans_info,backlist=get_gene_table(args.gtf,args.anno)
    if args.output_gene:
        with open(args.output_gene,"w") as fg:
            fg.write("ID\tGene_info\tRef_gene\tRef_gene_info\tEntrez\tUniprot\tSymbol\tDescription\n")
            for k in annoGeneTable:
                fg.write(k + "\t" + "\t".join(annoGeneTable[k]) + "\n")
    if args.output_trans:
        with open(args.output_trans,"w") as ft:
            ft.write("ID\tTrans_info\tTrack_gene\tTrack_gene_pos\tRef_gene\tRef_gene_info\tEntrez\tUniprot\tSymbol\tDescription\n")
            for t in Tran_Gene_Dict:
                tinfo = list(trans_info[t]["chr"])[0] + ":" + list(trans_info[t]["strand"])[0] + ":" + str(min(trans_info[t]["pos"])) + "-" + str(max(trans_info[t]["pos"]))
                ft.write(t + "\t" + tinfo +"\t" + Tran_Gene_Dict[t] + "\t" + "\t".join(annoGeneTable[Tran_Gene_Dict[t]]) + "\n")
    if args.back:
        with open(args.back,"w") as fb:
            for b in backlist:
                fb.write(str(b) + "\n")
            
                
if __name__ == "__main__":
    main()
        
        
