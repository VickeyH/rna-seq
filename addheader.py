#!/usr/bin/env python
import sys
if len(sys.argv) != 2:
    print "python scripts.py ASprofile.as.info"
    sys.exit(1)
f = open(sys.argv[-1])
txt = f.readlines()
f.close()
f=open(sys.argv[-1],"w")
f.write("#event_id\tevent_type\tgene_id\ttrans_id\n")
f.writelines(txt)
f.close()
