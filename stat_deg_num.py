#!/usr/bin/env python
#coding:utf-8
## python scripts.py deg.num.txt1 deg.num.txt2.... deg.state.txt
import os,sys
if len(sys.argv) == 1:
    print "python stat_deg_num.py deg.num.txt1 deg.num.txt2.... deg.state.txt"
    sys.exit(0)
for f in sys.argv[1:-1]:
        if not os.path.exists(f):
            print "%s file not exists!" %f
            sys.exit()
state_txt = os.popen("cat %s"%" ".join(sys.argv[1:-1])).read().strip().split("\n")
fo = open(sys.argv[-1],"w")
for stat in state_txt:
    if stat.strip() and not stat.endswith("Down"):
        fo.write(stat.split("\t")[0] + "\tup\t" + stat.split("\t")[-2] + "\n")
        fo.write(stat.split("\t")[0] + "\tdown\t" + stat.split("\t")[-1] + "\n")
fo.close()
