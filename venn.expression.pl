#! perl -w 
use File::Basename;
#this program is for valuate the expression in each sample and draw the Venn graph
if(!$ARGV[0]){print "Usage: perl $0 cuffdiff/gene_exp.diff.diff_stat venn.png venn.pdf\n";exit};
my $info=$ARGV[0];
my $dire=dirname $info;
#my $png=$ARGV[1];
#my $pdf=$ARGV[2];
my ($sample1,$sample2)=("","");
my ($num1,$num2,$num3)=(0,0,0);
open IN,$info;
my $i=0;
while(<IN>){
	chomp;
	if(/^Expressed\s+In\s+Both\t(\d+)/){$num3=$+;next};
	if($i==2){last}
	elsif($i==0 && $_=~/^Expressed\s+In\s+(\S+)\t(\d+)/){$sample1=$1;$num1=$2;$i++;next}
	elsif($i==1 && $_=~/^Expressed\s+In\s+(\S+)\t(\d+)/){$sample2=$1;$num2=$2;$i++;next}
	}
print "$num1 $num2 $num3 $sample1 $sample2\n";
#`Rscript /lustre/work/zhonghuali/software/rna.ref/bin/bin/R/Venn.r $num1 $num2 $num3 $sample1 $sample2 $ARGV[1] $ARGV[2] >$dire/venn.out 2> $dire/venn.err`;
`Rscript /lustre/work/yongdeng/software/protokaryon/flow/Venn.r $num1 $num2 $num3 $sample1 $sample2 $ARGV[1] $ARGV[2] >$dire/venn.out 2> $dire/venn.err`;
