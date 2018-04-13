#!/usr/bin/perl
#####程序所需加载的包################################################################################################
use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Data::Dumper;
use File::Basename qw(basename dirname);
use Cwd;
use FindBin qw($Bin $Script);

chomp $Bin;
#####时间和版本#######################################################################################################
my $BEGIN_TIME=time();
my $version="1.1";
########################################
####使用程序需要输入的参数############################################################################################
my ($deg,$pro_pro,$od);
GetOptions(		
	"help|?" =>\&USAGE,
	"deg:s"=>\$deg, 
	"pp:s"=>\$pro_pro, 
	"od:s"=>\$od,
) or &USAGE;

&USAGE unless ($deg and $pro_pro and $od);

$od=ABSOLUTE_DIR($od);
$pro_pro=ABSOLUTE_DIR($pro_pro);
$deg=ABSOLUTE_DIR($deg);

my @degfile;
my @deg;
if(-d $deg){
	@degfile=glob "$deg/*/*.diff_info";
}
if(-f $deg){
	push @degfile,$deg;
	push @deg,$deg;
}

for my $ele(@degfile){
	my $elename=basename $ele;
	if($elename=~/-VS-/){
		push @deg,$ele;
	}
}

my %degfile_list;
foreach my $degfile (@deg) {
	#print "$degfile\n";
	my $vsname;
	if(-f $deg){
		$vsname="inter";
	}
	else{
		$vsname=(split /\//,dirname $degfile)[-1];
	}
	open (IN,$degfile) or die $!;
	#my $head=<IN>;
	while (<IN>) {
		chomp;
		my $id=(split /\t/,$_)[0];
		if($id=~/^-/){next;}
		else{
			$degfile_list{$vsname}{$id}=1;
		}
	}
close IN;
}


open (PP,$pro_pro) or die $!;
my $head=<PP>;
chomp $head;
foreach my $vsname1 (sort keys %degfile_list) {
	`mkdir $od/$vsname1` if(!-d "$od/$vsname1");
	open OUT,">","$od/$vsname1/deg.pro-pro.txt";
	print OUT "source_id\ttarget_id\tcombine_score\tadj_score\n";
	close OUT;
}
while (<PP>) {
	chomp;
	my ($gene11,$gene12,$score)=(split /\t/,$_)[0,1,2];
#	my $gene1=uc($gene11);
#	my $gene2=uc($gene12);
	my $gene1=$gene11;
	my $gene2=$gene12;
	my $score1=$score/200.0;
#print "$gene1\t$gene2\n ";
	foreach my $vsname (sort keys %degfile_list) {
		open OUT,">>","$od/$vsname/deg.pro-pro.txt";
		if ($degfile_list{$vsname}{$gene1} && $degfile_list{$vsname}{$gene2} ) {
			print OUT"$gene1\t$gene2\t$score\t$score1\n";
		}else{
			#print "$vsname\t$gene11\t$gene12\n";
		}
		close OUT;
	}
}
close PP;


#####输入说明子程序####################################################################################################
sub USAGE
{
my $usage=<<"USAGE";
Program: combine fpkm deg
Version: $version
Contact: lv ran <lvran18931993760\@163.com>
Description:

Usage:

	"help|?" =>\&USAGE,
	"deg:s"=>\$deg,	deg dir or deg file 
	"pp:s"=>\$pro_pro, protein interaction db
	"od:s"=>\$od,	output directory
	
#############################################
USAGE
	print $usage;
	exit;
}

#####获取文件绝对路径####################################################################################################
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
