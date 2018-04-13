#!/usr/bin/env python

import sys,os
if len(sys.argv) == 1 or sys.argv[-1] == "--help" or sys.argv[-1] == "-h" :
    print "python scripts.py deg.txt ... deg.anno.txt"
    sys.exit(1)
f2 = open(sys.argv[-1])

annodict = {}
for line in f2:
    k = line.split("\t")[0]
    v = line.split("\t")[4].split("|")
    for i in v:
        if i != "-":
            annodict.setdefault(k,[]).append(i)
f2.close()            

for degfile in sys.argv[1:-1]:
    if os.path.basename(degfile) == "deg.up_down.list":
        print degfile + " will be overwrite, skip"
        continue
    with open(degfile) as fi, open(os.path.dirname(degfile) + "/deg.up_down.list","w") as fo:
        header = fi.next()
        fo.write(header)
        for line in fi:
            geneid = line.split("\t")[0]
            if annodict.has_key(geneid):
                for e in annodict[geneid]:
                    fo.write(line.replace(geneid,e))
                    
        
        
        
    
