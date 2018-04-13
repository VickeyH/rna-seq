#!/usr/bin/env python
#coding:utf-8
## python scripys.py sample1.count.txt sample2.count.txt sample3.count.txt.. gene_count.txt

import sys,os
from collections import OrderedDict
info = os.popen("grep -v '^__' %s"%" ".join(sys.argv[1:-1])).read().strip().split("\n")
info_dict = OrderedDict()
for line in info:
    k = os.path.basename(line.split(":")[0])
    info_dict.setdefault(k,[]).extend(line.split(":")[-1].split())
gene = info_dict.values()[0][::2]
count = [gene]
for l in info_dict.values():
    count.append(l[1::2])
count_line = []
for i in range(len(gene)):
    a = []
    for n in count:
        a.append(n[i])
    count_line.append(a)
   
with open(sys.argv[-1],"w") as fo:
    for i in count_line:
        fo.write("\t".join(i) + "\n")


