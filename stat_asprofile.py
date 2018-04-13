#!/usr/bin/env python
import os
import re
import sys

usge = """
usge :
        python xx.py -dir inputdir -o out.file
"""
def argvlist(argv):
        argvlist = []
        for arg in argv:
                argvlist.append(arg)
        return argvlist
def contdict(filepath):
        filedict = {}
        fileopen = open(filepath).readlines()
        for item in fileopen:
                key = item.split()[0].strip()
                filedict[key] = item.split()[1].strip()
        return filedict
def samplename(dir):
        samplelist = os.listdir(dir)
        return samplelist

def isdir(listdir,path):
        newlistdir = []
        for dir in listdir:
                dirpath = path+"/"+dir
                if os.path.isdir(dirpath):
                        newlistdir.append(dir)
        return newlistdir

def main():
	mainlist = argvlist(sys.argv)[1:]
	if len(mainlist)<1:
                print usge
                sys.exit()
        elif mainlist[0] == "-h" or mainlist[0] == "--help":
                print usge
                sys.exit()
        elif len(mainlist)!=4:
                print usge
                sys.exit()
        else:
		samplelist = samplename(mainlist[1])
		samplelist = isdir(samplelist,mainlist[1])
		if "merge" in samplelist:
                        samplelist.remove("merge")
                else:
                        pass
		samplelist = sorted(samplelist)
		filesave = open(mainlist[3],"w+")
		filesave.write("Samples\tTSS\tTTS\tAE\tSKIP\tIR\n")
		for samplenames in samplelist:
			samplenamestat = mainlist[1]+"/"+samplenames+"/ASprofile.as.types_stat"
			sampledict = contdict(samplenamestat)
			filesave.write(samplenames+"\t"+sampledict.get("TSS","0")+"\t"+sampledict.get("TTS","0")+"\t"+sampledict.get("AE","0")+"\t"+sampledict.get("SKIP","0")+"\t"+sampledict.get("IR","0")+"\n")
		cuffmergepath = mainlist[1]+"/merge/ASprofile.as.types_stat"
		cuffmergedict = contdict(cuffmergepath)
		filesave.write("merge"+"\t"+cuffmergedict.get("TSS","0")+"\t"+cuffmergedict.get("TTS","0")+"\t"+cuffmergedict.get("AE","0")+"\t"+cuffmergedict.get("SKIP","0")+"\t"+cuffmergedict.get("IR","0")+"\n")
		filesave.close()
			
if __name__ == '__main__':
        main()
