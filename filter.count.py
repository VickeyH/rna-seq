#!/usr/bin/env python

import sys,os

if len(sys.argv) != 3 or sys.argv[-1] == "-h" or sys.argv[-1] == '--help':
    print "USAGE: python scripts.py fpkm outputdir"
    sys.exit(0)

f = open(sys.argv[1])
if not os.path.isdir(sys.argv[-1]):
    os.makedirs(sys.argv[-1])
fcor = open(os.path.abspath(sys.argv[-1])+ "/gene.corfpkm.txt","w")
fany = open(os.path.abspath(sys.argv[-1])+ "/gene.anyfpkm.txt","w")
f0 = open(os.path.abspath(sys.argv[-1])+ "/gene.all0.txt","w")
header = f.next()
#header = header.replace("\tlength","")
fcor.write(header);fany.write(header);f0.write(header)
for line in f:
    l = line.strip().split("\t")
    if all(map(float,l[1:])):
        fcor.write("\t".join(l[:]) + "\n")
    if any(map(float,l[1:])):
        fany.write("\t".join(l[:]) + "\n")
    if not any(map(float,l[1:])):
        f0.write("\t".join(l[:]) + "\n")
fany.close()
fcor.close()
f0.close()
