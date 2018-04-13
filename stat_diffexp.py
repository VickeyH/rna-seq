#!/usr/bin/env python
#coding:utf-8
## python stat_diffexp.py group1-VS-group2.txt group3-VS-group4.txt.... 
import sys,os
def stat_diffexp(diffinfo):
    with open(diffinfo) as df:
        header = df.next()
        group1 = header.split("\t")[2]
        group2 = header.split("\t")[3]
        total_gene=0
        Expressed_Genes = 0 
        Expressed_g1 = [group1] 
        Expressed_g2 = [group2]
        Up_Diff_Expressed_Genes = 0
        Down_Diff_Expressed_Gene = 0
        for line in df:
            total_gene += 1
            line_list = line.strip().split("\t")
            if float(line_list[2]) + float(line_list[3]) > 0:
                Expressed_Genes += 1
            if float(line_list[2]) > 0:
                Expressed_g1.append(line_list[0])
            if float(line_list[3]) > 0:
                Expressed_g2.append(line_list[0])
            if line_list[-1] == "up":
                Up_Diff_Expressed_Genes += 1
            if line_list[-1] == "down":
                Down_Diff_Expressed_Gene += 1
        Total_Diff_Expressed_Genes = Up_Diff_Expressed_Genes + Down_Diff_Expressed_Gene
    return total_gene,Expressed_Genes,Expressed_g1,Expressed_g2,Total_Diff_Expressed_Genes,Up_Diff_Expressed_Genes,Down_Diff_Expressed_Gene

if __name__ == "__main__":
    for filename in sys.argv[1:]:
        pre = filename[:filename.rfind(".txt")] if filename.rfind(".txt") > 0 else filename
        with open(pre+".diff_stat","w") as ds:
            total_gene,Expressed_Genes,Expressed_g1,Expressed_g2,Total_Diff_Expressed_Genes,Up_Diff_Expressed_Genes,Down_Diff_Expressed_Gene = stat_diffexp(filename)
            eg1_num = len(Expressed_g1[1:])
            eg2_num = len(Expressed_g2[1:])
            eg12_num = len(set.intersection(set(Expressed_g1[1:]),set(Expressed_g2[1:])))
            eg1_no2_num = eg1_num - eg12_num
            eg2_no1_num = eg2_num - eg12_num
            ds.write('Discription\tNumber\tRatio(%)\n')
            ds.write('Total Genes\t%d\t%.2f\n'%(total_gene,100.0))
            ds.write('Expressed Genes\t%d\t%.2f\n'%(Expressed_Genes,float(Expressed_Genes)/total_gene*100))
            ds.write('Expressed In %s\t%d\t%.2f\n'%(Expressed_g1[0],eg1_num,float(eg1_num)/total_gene*100))
            ds.write('Expressed In %s\t%d\t%.2f\n'%(Expressed_g2[0],eg2_num,float(eg2_num)/total_gene*100))
            ds.write('Expressed In Both\t%d\t%.2f\n'%(eg12_num,float(eg12_num)/total_gene*100))
            ds.write('Expressed Only In %s\t%d\t%.2f\n'%(Expressed_g1[0],eg1_no2_num,float(eg1_no2_num)/total_gene*100))
            ds.write('Expressed Only In %s\t%d\t%.2f\n'%(Expressed_g2[0],eg2_no1_num,float(eg2_no1_num)/total_gene*100))
            ds.write('Total Diff Expressed Genes\t%d\t%.2f\n'%(Total_Diff_Expressed_Genes,float(Total_Diff_Expressed_Genes)/total_gene*100))
            ds.write('Up Diff Expressed Genes\t%d\t%.2f\n'%(Up_Diff_Expressed_Genes,float(Up_Diff_Expressed_Genes)/total_gene*100))
            ds.write('Down Diff Expressed Genes\t%d\t%.2f\n'%(Down_Diff_Expressed_Gene,float(Down_Diff_Expressed_Gene)/total_gene*100))
            