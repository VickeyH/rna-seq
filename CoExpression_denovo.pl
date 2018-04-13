#!/usr/bin/perl

# Script Name: CoExpression.pl
# Version: v1.0
# Usage: just run this script.
# Author: yujingzhang@capitalbio.com

# Version: v1.01
# Date: 2014-06-26

use FileHandle;
use File::Basename qw(basename dirname);
use Cwd;
use warnings;
use strict;
use File::Spec;
use Data::Dumper;
use Getopt::Long;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use CapitalBioTools;
use PlotInter;
use List::Util qw/max/;
use Cwd;
sub getImage;

########## Configuration & Parameters ##########

#my $pgmPath="/home/mammoth/zhangyujing/Coexpression/Programme/";
my $pgmPath="$RealBin/"; 
my $nodenum=10000;
my $edgenum=10000;
my $filterout="filter.file";
my $config  = "${pgmPath}kointerfigconfig.txt";
my $configOlp="${pgmPath}kointerfigconfigOlp.txt";
my $configPro="${pgmPath}kointerfigconfigPro.txt";
print "$config\n$configOlp\n";

my ($help,$exp_sgl,$gene_filter,$sample_filter,$rcut,$pcut,$out,$name1,$name2);

GetOptions(
			
                        "h|help" => \$help,
                        "e|exp_signal=s" => \$exp_sgl,
						"r|r_cutoff=s" => \$rcut,
						"p|p_cutoff=s" => \$pcut,
						"o|out=s" => \$out,

                       
);

my $parametersInvalid = 0;

if ( defined($help) ) {
	prtUsage();
	exit;
}
if(!defined $exp_sgl){
	CapitalBioTools::logError("NO input expression_signal file!");
	$parametersInvalid = 1;
}

$gene_filter="abcdefg---------a";
if(!defined $rcut){
     $rcut=0.99;
}
if(!defined $pcut){
     $pcut=0.05;
}
$sample_filter="";
if(!defined $out){
	$out=getcwd();
}
if(!defined $name1){
	$name1="case";
}
if(!defined $name2){
	$name2="control";
}
if ($parametersInvalid) {
	prtUsage();
	exit;
}
########## Log Configurations & Parameters ##########
&MKDIR($out);
$out=ABSOLUTE_DIR($out);
$exp_sgl=ABSOLUTE_DIR($exp_sgl);
chdir($out);
########## Main ##########

CapitalBioTools::logInfo("Program started");
#构建共表达网络

open(RPROG, "|R --vanilla --slave ") || failmessage($!);
print RPROG <<"CODE";
source("${pgmPath}CoExpression.R")
CoExpression("$rcut","$pcut" ,"$exp_sgl" ,"$gene_filter","$sample_filter","$out","$name1","$name2")
q()
CODE
close RPROG;


sub MKDIR
{ # &MKDIR($out_dir);
	my ($dir)=@_;
	rmdir($dir) if(-d $dir);
	mkdir($dir) if(!-d $dir);
}
########## Subs ##########
sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $dir=dirname($in);
		my $file=basename($in);
		chdir $dir;$dir=`pwd`;chomp $dir;
		$return="$dir/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and dir in [sub ABSOLUTE_DIR]\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}

# 网络作图
sub getImage(@){
	my ($path,$inputFile,$outputPng,$outputFile,$P,$R,$config) =@_;
	my %degP = %$P;
	my %degR;
	if ($R ne "") {
		%degR = %$R;
	}
	
	#读取共表达网络
	open(NET,$inputFile)||die "无法打开文件$inputFile\n$!";
	my (%edgeinfo,%degq,%itainfo);
	my ($sc,$tg,$scgs,$tggs,$cor);
	while (<NET>) {
		chomp;
		tr/\r'//d;
		my @array = split /\t/;
		if (/^Source/i){
			my $i=0;
			foreach my $n (@array) {
				$sc = $i if $n =~ /^Source$/i;
				$tg = $i if $n =~ /^Target$/i;
				$scgs = $i if $n =~ /^Source.GeneSymbol$/i;
				$tggs = $i if $n =~ /^Target.GeneSymbol$/i;
				$cor = $i if $n =~ /^Correlation$/i;
				$i++;
			}
			next;
		};

		$itainfo{$array[$sc]}{$array[$tg]}=join("\t",@array[2..$#array]);
		$edgeinfo{$array[$sc]}{$array[$tg]}=$array[$cor];
		$degq{$array[$sc]}{"q"}=1/$degP{$array[$sc]};
		$degq{$array[$tg]}{"q"}=1/$degP{$array[$tg]};
		if (defined $scgs) {
			$degq{$array[$sc]}{"des"}=$array[$scgs];
			$degq{$array[$tg]}{"des"}=$array[$tggs];
		}
		if ($R ne "") {
			$degq{$array[$sc]}{"h"}=$degR{$array[$sc]};
			$degq{$array[$tg]}{"h"}=$degR{$array[$tg]};
		}
	}
	close NET;

	#绘图
	my @png=qw(png);
	push my @pngout,File::Spec->catfile($path,$outputPng);
	my(%newinter2)=FilterInteraction($nodenum,$edgenum,$filterout,\%edgeinfo,\%degq,$config);
	PlotInter::PlotInteraction(\@png,\@pngout,$config,\%newinter2,\%degq);
	my $output=File::Spec->catfile($path,$outputFile);
	my %newintsym;
	open(OUT,">$output");
	if ($R ne "") {
		print OUT "Source\tTarget\tCorrelation\tP.value\tSource.GeneSymbol\tTarget.GeneSymbol\n";
	}else{
		print OUT "Source\tTarget\tCorrelation\tP.value\n";
	}
	foreach my $A (keys %newinter2) {
		foreach my $B (keys %{$newinter2{$A}}) {
			if (defined $degq{$A}{"des"}) {
				$newintsym{$degq{$A}{"des"}}{$degq{$B}{"des"}}=1;
			}
			print OUT "$A\t$B\t$itainfo{$A}{$B}\n";
		}
	}
	close OUT;
	%newintsym;
}

sub getDegree(@){
	my $inputFile = shift @_;
	my %nodeDegree;
	open(IN,$inputFile)||die $!;;
	while (<IN>) {
		chomp;
		tr/\r'//d;
		my @array = split /\t/;
		next if /^node/;
		$nodeDegree{$array[1]}=$array[2];
	}
	close IN;
	%nodeDegree;
}

#perl程序参数使用说明
sub prtUsage {
    print "\n\n$0: The program is to plot Co-Expression network\n\n";
	print "\n$0 <require> \n\n";
	print "  -e | -exp_signal (all expression signal file)\n";
	print "  -r | -r_cutoff correlation coefficient \n";
	print "  -p | -p_cutoff p value\n";
	print "  -o | output dir\n";
	print "### [Optional]\n";
	print "  -h | -help\n";
	print "\n";
}

