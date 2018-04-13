#!/usr/bin/env python
#coding:utf-8

### pick_deg_by_Fc.py fc group1-VS-group2.txt group1-VS-group2.new.txt

from math import log
import sys
if len(sys.argv) != 4:
    print "USAGE: python pick_deg_by_Fc.py fc group1-VS-group2.txt group1-VS-group2.new.txt"
    sys.exit(-1)
with open(sys.argv[2]) as fi, open(sys.argv[3],"w") as fo:
    head = fi.next()
    fo.write(head.replace("Log2FC","Log%sFC"%sys.argv[1]))
    for line in fi:
        line_list = line.split("\t")
        if line_list[-2] == "NA":
            fo.write("\t".join(line_list[:4]) + "\t"  + str(float(line_list[4])/log(float(sys.argv[1]),2)) + "\t" + "\t".join(line_list[5:-1]) + "\tno\n")
            continue
        if float(line_list[-2]) <= 0.05 and float(line_list[4]) >= log(float(sys.argv[1]),2):  ## 上调
            fo.write("\t".join(line_list[:4]) + "\t"  + str(float(line_list[4])/log(float(sys.argv[1]),2)) + "\t" + "\t".join(line_list[5:-1]) + "\tup\n")
        elif float(line_list[-2]) <= 0.05 and float(line_list[4]) <= -log(float(sys.argv[1]),2):
            fo.write("\t".join(line_list[:4]) + "\t"  + str(float(line_list[4])/log(float(sys.argv[1]),2)) + "\t" + "\t".join(line_list[5:-1]) + "\tdown\n")
        else:
            fo.write("\t".join(line_list[:4]) + "\t"  + str(float(line_list[4])/log(float(sys.argv[1]),2)) + "\t" + "\t".join(line_list[5:-1]) + "\tno\n")
           
            