#! perl -w
use File::Basename;
use Getopt::Long;

my ($species,$od,$name,$snpdir,$fa,$help,$aligndir)=("","","","","","","","");
GetOptions(
	"species=s"	=>\$species,
	"samples=s"	=>\$samples,
	"od=s"	=>\$od,
	"snpdir=s"	=>\$snpdir,
	"fa=s"	=>\$fa,
	"name=s"	=>\$name,
	"h|help!"	=>\$help,
	"aligndir=s"	=>\$aligndir,
);

sub usage {
	print<<DO;
	This program is for calling SNP and InDel
	Usage:perl $0 -species [species] -samples [samples] -od [makeflow output dir] -snpdir [SNP output dir] -fa [fa genome] -aligndir [alignment results dir] -name [makeflow name]
	
	-species	species name same with snpEff annotation, required; Recommand info:/lustre/work/minzhao/software/snpEff
	-samples	samples info, sep by ","
	-od	makeflow output dir, required;
	-snpdir	SNP output directroy, required;
	-fa	fa reference sequence file, required;
	-aligndir	alignment output directory, required;
	-name	makeflow name, required;
	-h|help	help option
DO
exit;
}

my $ooops=0;
if(!$od){
	print "Please give a makeflow output directory\n";
	$ooops=1;
}
if(!$snpdir){
	print "please give a SNP output directory\n";
	$ooops=1;
}
if(!$fa){
	print "please give the reference genome sequences fa file\n";
	$ooops=1;
}
if(!$aligndir or !-d $aligndir){
	print "please give the alignment results directory\n";
	$ooops=1;
}
if($ooops==1){
	&usage;
	exit;
}

my $samtools="/lustre/software/target/samtools-1.2/bin/samtools";

`mkdir $od` if(!-d "$od");
`mkdir $snpdir` if(!-d "$snpdir");

$od=&ABSOLUTE_DIR($od);
$snpdir=&ABSOLUTE_DIR($snpdir);
my $Bin="/lustre/work/zhonghuali/software/rna.ref/bin/snp";

my @samples;
if(-f $samples){
	$samples=&ABSOLUTE_DIR($samples);
	open SAM,"<",$samples;
	while(<SAM>){
		chomp;
		next if(/^#|^\s+/);
		if(/^rep=(\S+)/){
			my $sample=$+;
			my @rep=split /;/,$sample;
			for my $ele(@rep){
				if($ele=~/=(\S+?)\)/){
					my $samplesss=$+;

					my @samplessss=split /,/,$samplesss;
					for my $samele(@samplessss){
						push @samples,$samele;
					}
				}
				else{
					push @samples,$ele;
				}
			}
			last;	
		}
		elsif($_=~/;/ and $_!~/-VS-/){
			my @rep=split /;/;
			for my $repele(@rep){
				push @samples, $repele;
			}
			last;
		}	
	}
}
else{
	my @rep=split /,/,$samples;
	for my $repele(@rep){
		push @samples, $repele;
	}
}
my %hash;
@samples=grep {++$hash{$_}<2} @samples;

$aligndir=&ABSOLUTE_DIR($aligndir);

my @snpstat;
open FILE,">$snpdir/bam.file.list";
my @dir;

if($name){
	open MAKEFLOW,">>", "$od/$name";
}
else{
	open MAKEFLOW , ">","$od/snp.Makeflow";
}

for my $key( sort @samples){

	push @snpstat,"$snpdir/$key.bam";
	print FILE "$snpdir/$key.bam\n";

	print MAKEFLOW "CATEGORY=SNP\n";
	print MAKEFLOW "$snpdir/$key.bam : $aligndir/$key/accepted_hits.bam\n";
	print MAKEFLOW "\tln -s $aligndir/$key/accepted_hits.bam $snpdir/$key.bam 2> $snpdir/$key.snp.err\n\n";
	}
close FILE;

my $snpstat=join " ",@snpstat;

print MAKEFLOW "CATEGORY=SNP\n";
print MAKEFLOW "$snpdir/total.bam : $snpstat\n\@BATCH_OPTIONS = -l h_vmem=30G\n";
print MAKEFLOW "\t$samtools mpileup -f $fa -q 1 -b $snpdir/bam.file.list -o $snpdir/total.bam > $snpdir/total.bam.out 2> $snpdir/total.bam.err\n\n";

print MAKEFLOW "CATEGORY=SNP\n";
print MAKEFLOW "$snpdir/snp.vcf : $snpdir/total.bam\n\@BATCH_OPTIONS = -l h_vmem=30G\n";
print MAKEFLOW "\tjava -Xmx20g -jar /lustre/software/target/varscan-2.3.7/VarScan.v2.3.7.jar mpileup2snp $snpdir/total.bam --min-reads2 3 --min-avg-qual 20 --min-var-freq 0.2 --p-value 0.05 --output-vcf 1 > $snpdir/snp.vcf 2> $snpdir/snp.vcf.err\n\n";

print MAKEFLOW "CATEGORY=SNP\n";
print MAKEFLOW "$snpdir/snp.readcount.txt : $snpdir/snp.vcf $snpdir/bam.file.list\n";
print MAKEFLOW "\tperl $Bin/get.SNP.InDel.readcounts.pl -snp $snpdir/snp.vcf -out $snpdir/snp.readcount.txt -bamlist $snpdir/bam.file.list > $snpdir/snp.readcount.out 2> $snpdir/snp.readcount.err\n\n";

print MAKEFLOW "CATEGORY=SNP\n";
print MAKEFLOW "$snpdir/indel.vcf : $snpdir/total.bam\n\@BATCH_OPTIONS = -l h_vmem=30G\n";
print MAKEFLOW "\tjava -Xmx20g -jar /lustre/software/target/varscan-2.3.7/VarScan.v2.3.7.jar mpileup2indel $snpdir/total.bam --min-reads2 5 --min-avg-qual 20 --min-var-freq 0.2 --p-value 0.05 --output-vcf 1 > $snpdir/indel.vcf 2> $snpdir/indel.vcf.err\n\n";

print MAKEFLOW "CATEGORY=SNP\n";
print MAKEFLOW "$snpdir/indel.readcount.txt : $snpdir/indel.vcf $snpdir/bam.file.list\n";
print MAKEFLOW "\tperl $Bin/get.SNP.InDel.readcounts.pl -snp $snpdir/indel.vcf -out $snpdir/indel.readcount.txt -bamlist $snpdir/bam.file.list > $snpdir/indel.readcount.out 2> $snpdir/indel.readcount.err\n\n";

print MAKEFLOW "CATEGORY=SNP\n";
print MAKEFLOW "$snpdir/rm.bam.out : $snpdir/snp.vcf $snpdir/indel.vcf $snpdir/total.bam\n";
print MAKEFLOW "\trm -rf $snpdir/total.bam > $snpdir/rm.bam.out 2> $snpdir/rm.bam.err\n\n";


print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/snpeff.vcf $snpdir/snpeff_summary.csv: /lustre/work/minzhao/software/snpEff/snpEff.config $snpdir/snp.vcf\n\@BATCH_OPTIONS = -l h_vmem=30G\n";
print MAKEFLOW "\tjava -Xmx20g -jar /lustre/work/minzhao/software/snpEff/snpEff.jar eff -i vcf -c /lustre/work/minzhao/software/snpEff/snpEff.config -csvStats $snpdir/snpeff_summary.csv -s $snpdir/snpeff_summary.html -ud 3000 $species $snpdir/snp.vcf > $snpdir/snpeff.vcf 2> $snpdir/snpeff.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/indeleff.vcf $snpdir/indeleff_summary.csv: /lustre/work/minzhao/software/snpEff/snpEff.config $snpdir/indel.vcf\n\@BATCH_OPTIONS = -l h_vmem=30G\n";
print MAKEFLOW "\tjava  -Xmx20g -jar /lustre/work/minzhao/software/snpEff/snpEff.jar eff -i vcf -c /lustre/work/minzhao/software/snpEff/snpEff.config -csvStats $snpdir/indeleff_summary.csv -s $snpdir/indeleff_summary.html -ud 3000 $species $snpdir/indel.vcf > $snpdir/indeleff.vcf 2> $snpdir/indeleff.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/snpeff.filter.txt $snpdir/snpeff.high.filter.txt : $snpdir/snpeff.vcf $snpdir/bam.file.list\n";
print MAKEFLOW "\tperl $Bin/snp.indel.filter.pl -vcf $snpdir/snpeff.vcf -od $snpdir -out snpeff.filter.txt -bamlist $snpdir/bam.file.list > $snpdir/snpeff.high.filter.txt 2> $snpdir/snpfilter.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/indeleff.filter.txt $snpdir/indeleff.high.filter.txt : $snpdir/indeleff.vcf $snpdir/bam.file.list\n";
print MAKEFLOW "\tperl $Bin/snp.indel.filter.pl -vcf $snpdir/indeleff.vcf -od $snpdir -out indeleff.filter.txt -bamlist $snpdir/bam.file.list > $snpdir/indeleff.high.filter.txt 2> $snpdir/indelfilter.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/snpeff.csv $snpdir/indeleff.csv : $snpdir/snpeff.vcf $snpdir/indeleff.vcf\n";
print MAKEFLOW "\tperl $Bin/vcf2csv.pl -v $snpdir >$snpdir/csv.out 2> $snpdir/csv.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/indeleff.state.txt $snpdir/snpeff.state.txt : $snpdir/snpeff_summary.csv $snpdir/indeleff_summary.csv\n";
print MAKEFLOW "\tperl $Bin/call.summary.info.pl -d $snpdir >$snpdir/stat.out 2> $snpdir/stat.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/snpEFF.png $snpdir/snpEFF.pdf : $snpdir/snpeff.state.txt\n";
print MAKEFLOW "\tRscript  $Bin/snp.histgram.r -f $snpdir/snpeff.state.txt -o $snpdir/snpEFF -x \\\'\\\"genomic regions\\\"\\\' -y \\\'\\\"SNP Number\\\"\\\' -l 6 -a 4 > $snpdir/snp.png.out 2> $snpdir/snp.png.err\n\n";

print MAKEFLOW "CATEGORY=SNPEFF\n";
print MAKEFLOW "$snpdir/indelEFF.png $snpdir/indelEFF.pdf : $snpdir/indeleff.state.txt\n";
print MAKEFLOW "\tRscript  $Bin/snp.histgram.r -f $snpdir/indeleff.state.txt -o $snpdir/indelEFF -x \\\'\\\"genomic regions\\\"\\\' -y \\\'\\\"InDel Number\\\"\\\' -l 6 -a 4 > $snpdir/indel.png.out 2> $snpdir/indel.png.err\n\n";

close MAKEFLOW;

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
