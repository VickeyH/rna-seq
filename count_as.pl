#!/usr/bin/perl -w
use strict;
my $file=shift;   #### *.ASprofile.as
my %hash;
my %iden_tss;
my %iden_tts;
open IN,$file;
open OUT1,">$file.info\n" or die "$!\n";
while(<IN>){
	chomp;
	next if(/^#/);
	my @c=split(/\t/);
	if($c[1]=~/TSS|TTS/){
	#	next if($c[8]=~/^ENS/);
		@{$iden_tss{$c[8]}}=($c[0],$c[1],$c[2]) if($c[1]=~/TSS/);
		@{$iden_tts{$c[8]}}=($c[0],$c[1],$c[2]) if($c[1]=~/TTS/);
#		$hash{$c[1]}++;
#		print OUT1"$c[0]\t$c[1]\t$c[2]\t$c[8]\n";
	}
	if($c[1]=~/AE|XAE/){
		my $line=<IN>;
		chomp $line;
		my @d=split("\t",$line);
	#	next if($c[8]=~/^ENS/ && $d[8]=~/^ENS/);
		$hash{'AE'}++;
		print OUT1"$c[0]\tAE\($c[1]\)\t$c[2]\t$c[8],$d[8]\n";	
	}if($c[1]=~/SKIP_ON$|IR_ON$/){
		my $line=<IN>;
		chomp $line;
		my @d=split("\t",$line);
	#	next if($c[8]=~/^ENS/ && $d[8]=~/^ENS/);
		my $sub_type=$1 if($c[1]=~/(\S+)_/);
		my $type=($c[1]=~/SKIP/)?"SKIP":"IR";
		print OUT1"$c[0];$d[0]\t$type\($sub_type\)\t$c[2]\t$c[8],$d[8]\n";
		$hash{$type}++;
	}
}
close IN;

foreach my $key(sort keys %iden_tss){
	if(!exists $iden_tts{$key}){
		$hash{'TSS'}++;
		print OUT1"$iden_tss{$key}[0]\t$iden_tss{$key}[1]\t$iden_tss{$key}[2]\t$key\n";
	}
}

foreach my $key(sort keys %iden_tts){
	if(!exists $iden_tss{$key}){
		$hash{'TTS'}++;
		print OUT1"$iden_tts{$key}[0]\t$iden_tts{$key}[1]\t$iden_tts{$key}[2]\t$key\n";
	}
}
close OUT1;
open OUT,">$file.types_stat";
foreach my $key(sort keys %hash){
	print OUT"$key\t$hash{$key}\n";
}
close OUT;
