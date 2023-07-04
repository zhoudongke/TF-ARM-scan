#!/bin/perl -w
use strict;
die "usage: <input raw basic patches file> <IDRs info file>\n" unless @ARGV == 2;

my %idr;
open IDR,"$ARGV[1]" or die "file $ARGV[1] not found!\n";
#open IDR,"$idr" or die "IDRs file not found!\n";
while (<IDR>){
	chomp;
	my @a=split /,/;
	my $tf_name=(split /\s+/,$a[1])[1];
	my $n=0;
	foreach (@a[3..$#a]){
		$n++;
		$idr{$tf_name}{$n}=1 if $_ > 0.2;
	}
}
close IDR;


open IN,"$ARGV[0]" or die "file $ARGV[0] not found!\n";
open OUT,">$ARGV[0]\.filtered.by.IDRs" or die "permission denied!\n";
my $switch = 0;
my $head;
my $tf = '';

READ:
while (<IN>){
    if (/^>/){
        $head = $_;
        s/^>//;
        $tf = (split /\s+/)[0];
        $switch = 1;
        next;
    }
    chomp;
    my $bp = $_;
    my @a=split /\:/;
    my $order=0;
    foreach ($a[0]..$a[1]){
        unless (exists $idr{$tf}{$_}){
            $order++;
        }
    }
    my $ratio = $order / $a[-1] if $a[-1] != 0;
    if ($ratio > 0){
        next READ;
    }
    if ($switch == 1){
        print OUT "$head$bp\n";
        $switch = 0;
    }else{
        print OUT "$bp\n";
    }
}
close OUT;
close IN;