#! perl -w 
use Getopt::Long;
use File::Basename;

my  ($input,$help,$name,$type,$expression)=("","","","","");

GetOptions(
	"input=s"	=>\$input,
	"name=s"	=>\$name,
	"type=s"	=>\$type,
	"expression!"	=>\$expression,
	"h|help!"	=>\$help,
);

sub usage{
	print<<DO;
	
	This program is to the coexpression or protein-protein interaction info
	version 2.0 : first filter according to correlation and p-val, then order them by degree
	
	perl $0 -input [coexpression file] -name [output file name] -type [coexpression or pro-pro] -expression
	
	-input	infomation, file or directory,file eg:dir/filename; dir eg: dir,filenme which can get sample info by dir/*/filename
	-name	output file name
	-type	pro-pro/coexpression; pro-pro should be 3 columns; but coexpression 4
	-expression	if the table delete the tracking_id with brackets
	-h|help	help option
	
DO
}

my $ooops=0;
if(!$input){
	print "NO input file(coexpression/pro-pro)\n";
	$ooops=1;
}
if(!$name){
	print "NO output file name\n";
	$ooops=1;
}
if(!$type){
	print "NO input data type\n";
	$ooops=1;
}
if($ooops==1){
	&usage;
	exit
}

my @file;
if(-f $input){
	$input=&ABSOLUTE_DIR($input);
	push @file,$input; 
}
elsif($input=~/,/){
	my @infotmp=split /,/,$input;
	@file=glob "$infotmp[0]/*/$infotmp[1]";
}
my %repdele;
@file=grep {++$repdele{$_}<2} @file;
for my $file(@file){
	chomp $file;
	my $od=dirname $file;
	open IN,"<",$file;
	open OUT,">","$od/$name";
	print OUT "Source\tTarget\tCorrelation\tP.value\tadj_cor\n";
	my (%inf,%num,%p,%cor,%score,%done);
	my @total;
	my ($n,$m)=(0,0);
	if($type=~/coexp/){
		while(<IN>){
			chomp;
			next if(/^\s+|^#|^Source/);
			my($gene1,$gene2,$cor,$p)=split /\t/;
			if($expression){
					if($gene1=~/^(\S+?)\(/){
						$gene1=$+;
					}
					if($gene2=~/^(\S+?)\(/){
						$gene2=$+;
					}
			}
			if($m<1000){
				my $adj_cor;
				if($cor>0){
						$adj_cor=1;
				}
				else{$adj_cor=-1;}
				push @total,"$_\t$adj_cor";
			}
			#else{@total=();}
			$m++;
			if(abs($cor)>=0.99 and $p<=0.001){
				($gene1,$gene2)=sort ($gene1,$gene2);
				$num{$gene1}++;$num{$gene2}++;
				$inf{$gene1}{"$gene1 $gene2"}=1;
				$inf{$gene2}{"$gene1 $gene2"}=1 ;
				$p{"$gene1 $gene2"}=$p;
				$cor{"$gene1 $gene2"}=$cor;
				$n++;
			}
		}
		close IN;
		if($m<1000 or ($m>1000 and $n<100)){
			for my $ele(@total){
				print OUT "$ele\n";
			}
		}
		elsif($m>1000 and $n<1000){
			for my $ele(keys %p){
				my $eleid=$ele;
				$eleid=~s/ /\t/g;
				if($cor{$ele}>0){
						$adj_cor=1;
				}
				else{$adj_cor=-1;}
				print OUT "$eleid\t$cor{$ele}\t$p{$ele}\t$adj_cor\n";
			}
		}
		elsif($n>1000){
			my $i=0;
			for my $id1(sort {$num{$b} <=> $num{$a}} keys %num){
				print "$id1\t$num{$id1}\n";
				for my $idx (sort {$p{$a} <=> $p{$b}}(sort {abs($cor{$b}) <=> abs($cor{$a})} keys $inf{$id1})){
					if($done{$idx}){ next };
					my $id=$idx;
					$idx=~s/ /\t/g;
					my $adj_cor;
					if($i<1000){
						if($cor{$id}>0){
							$adj_cor=1;
						}
						else{$adj_cor=-1;}
						if(abs($cor{$id})<0.99 or $p{$id}>0.001){next};
						print OUT "$idx\t$cor{$id}\t$p{$id}\t$adj_cor\n";
						$done{$id}=1;
						$i++;
					}
				}
			}
		}
	}
	if($type=~/pro/){
		while(<IN>){
			chomp;
			next if(/^\s+|^#|^source/);
			my($gene1,$gene2,$score)=split /\t/;
			if($m<1000){
				my $adj_score=$score/200.0;
				push @total,"$_\t$adj_score";
			}
		#	else{@total=()}
			$m++;
			if($score>=400){
				($gene1,$gene2)=sort ($gene1,$gene2);
				$num{$gene1}++;$num{$gene2}++;
				$inf{$gene1}{"$gene1 $gene2"}=1;
				$inf{$gene2}{"$gene1 $gene2"}=1;
				$score{"$gene1 $gene2"}=$score;
				$n++;
			}
		}
		close IN;
		if($m<1000 or ($m>1000 and $n<100)){
			for my $ele(@total){
				print OUT "$ele\n";
			}
		}
		elsif($m>1000 and $n<1000){
			for my $ele(keys %score){
				my $eleid=$ele;
				$eleid=~s/ /\t/g;
				my $adj_score=$score{$ele}/200.0;
				print OUT "$eleid\t$score{$ele}\t$adj_score\n";
			}
		}
		elsif($n>1000){
			my $i=0;
			for my $id1(sort {$num{$b} <=> $num{$a}} keys %num){
				print "$id1\t$num{$id1}\n";
				my @linkid;
				for my $idx (sort {$score{$b} <=> $score{$a}} keys $inf{$id1}){
					if($done{$idx}){ next };
					my $id=$idx;
					$idx=~s/ /\t/g;
					my $adj_score=$score{$id}/200.0;
					if($adj_score<2.0){next};
					if($i<1000){
						print OUT "$idx\t$score{$id}\t$adj_score\n";
						$done{$id}=1;
						$i++;
					}
				}
			}
		}
	}
}

#.............................................................................

#.............................................................................

sub ABSOLUTE_DIR
{ #$pavfile=&ABSOLUTE_DIR($pavfile);
	my $cur_dir=`pwd`;chomp($cur_dir);
	my ($in)=@_;
	my $return="";
	
	if(-f $in)
	{
		my $od=dirname($in);
		my $file=basename($in);
		chdir $od;$od=`pwd`;chomp $od;
		$return="$od/$file";
	}
	elsif(-d $in)
	{
		chdir $in;$return=`pwd`;chomp $return;
	}
	else
	{
		warn "Warning just for file and od in [sub ABSOLUTE_DIR] $in\n";
		exit;
	}
	
	chdir $cur_dir;
	return $return;
}
