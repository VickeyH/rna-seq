#!/usr/bin/env python
#coding:utf-8
## 当蛋白互作的蛋白id和gtf文件中的蛋白id不对应的时候，用此脚本对应上注释文件的第一列id和蛋白库中的id
## get gene_protein.txt
## python scripts.py entrez_ref.txt protein_entrez.txt
import sys
e2r={}
for i in open(sys.argv[1]):
    if i.strip().split("\t")[-1] != "-":
        e2r[i.strip().split("\t")[-1]] = i.strip().split("\t")[0]
f=open(sys.argv[2])
p2r = {}
for line in f:           
    if bool(line.strip("\n").split("\t")[-1]) and e2r.has_key(line.strip("\n").split("\t")[-1]) and line.strip("\n").split("\t")[0]:
        p2r.setdefault(e2r[line.strip("\n").split("\t")[-1]],[]).append(line.strip("\n").split("\t")[0])
fo = open("gene_protein.txt","w")
for i in p2r:
    fo.write(i + "\t" + ",".join(p2r[i]) + "\n")
