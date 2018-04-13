#!/usr/bin/env python
import argparse,os,re

parser = argparse.ArgumentParser(description="This Script is used to combine all raw fq file to one fq file per sample.")
parser.add_argument("-i","--input",type=str,help="The input rna_seq_metrics file generated",required = True)
parser.add_argument("-u","--uniform",type=str,help="The output uniform distribution file")
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-e","--element",type=str,help="The output element distribution file")
args = parser.parse_args()

def main():
    if not os.path.exists(args.input):
        print "%s file not exists!"%args.input
        sys.exit(0)
    ele_dict = {}
    p,c = [],[]
    with open(args.input) as m:
        for line in m:
            if line.startswith("PF_BASES"):
                k = line.strip().split()
                v = m.next().strip("\n").split("\t")
            if line.startswith("normalized_position"):
                h = line
                break
        ele_dict = dict(zip(k,v))
        for line in m:
            if line.strip():
                p.append(line.split()[0])
                c.append(line.strip().split()[1])
    uniform = zip(p,c)
    if args.element:
        with open(args.element,"w") as ele_file:
            ele_file.write("CODING_BASES\t" + ele_dict["CODING_BASES"] + "\n")
            ele_file.write("UTR_BASES\t" + ele_dict["UTR_BASES"] + "\n")
            ele_file.write("INTRONIC_BASES\t" + ele_dict["INTRONIC_BASES"] + "\n")
            ele_file.write("INTERGENIC_BASES\t" + ele_dict["INTERGENIC_BASES"] + "\n")
    if args.uniform:
        with open(args.uniform,"w") as uni_file:
            uni_file.write(h)
            for i in uniform:
                uni_file.write("\t".join(i) + "\n")

if __name__ == "__main__":
    main()
            