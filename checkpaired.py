#!/usr/bin/env python
#coding:utf-8
import os,sys,re
from commands import getstatusoutput

if len(sys.argv)!=3:
    print "USAGE: python script.py r1.fq r2.fq"
    sys.exit(1)

if not os.path.exists(sys.argv[1]) or not os.path.exists(sys.argv[2]):
    print "Error: file do not exist"
    sys.exit(1)

s,out = getstatusoutput('python /lustre/work/yongdeng/software/protokaryon/flow/sms3.py -fq %s %s'%(sys.argv[1],sys.argv[2]))
if re.search('\nsequence number:.+?\d+?\n',out):
    r = re.findall('\nsequence number:.+?(\d+?)\n',out)
    if len(r) != 2:
        print "Error: paired sequence not match [%s,%s], please check files [ '%s', '%s' ]"%(r[0],r[1],os.path.abspath(sys.argv[1]),os.path.abspath(sys.argv[2]))
        sys.exit(1)
    if r[0] != r[1]:
        print "Error: paired sequence not match [%s,%s], please check files [ '%s', '%s' ]"%(r[0],r[1],os.path.abspath(sys.argv[1]),os.path.abspath(sys.argv[2]))
        sys.exit(1)
    print "Paired match, looks fine."
    sys.exit(0)
else:
    print "Error: can not read sequence file, please check files [ '%s', '%s' ]"%(os.path.abspath(sys.argv[1]),os.path.abspath(sys.argv[2]))
    sys.exit(1)

