#!/usr/bin/perl -w
use strict;
if($#ARGV<1){ print "Usage: perl $0 <fa> <species> <output>\n";exit}
my $file=shift;
my $org=shift;  ## eg. org=Homo_Sapiens
my $out=shift;
open IN,$file;
open OUT,">$out";
$/=">";<IN>;
while(<IN>){
	chomp;
	my @array=split("\n",$_);
	my $head=shift @array;
	my $id=$1 if($head=~/^(\S+)/);
	my $seq=join ("\n",@array);
    $seq=~s/\s//g;
    my $N_num=($seq =~ tr/Nn/Nn/);
    my $len=length $seq;
    my $nonNlen=$len-$N_num;
    print OUT">$id /len=$len /nonNlen=$nonNlen /org=$org\n";
}
close IN;
