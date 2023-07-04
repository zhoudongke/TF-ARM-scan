#!/bin/perl -w
use strict;
die "usage: <input raw .tsv file from InterPro analysis>\n" unless @ARGV == 1;

open DBDENTRY,"/home/dongke/TF_RNA/prediction.of.TF.DNA.binding.domain.by.InterPro/entry.id.TF.with.only.1.list" or die "DBD entryID list file not found!\n";
my %entry;
while (<DBDENTRY>){
    chomp;
    $entry{$_}=1;
}
close DBDENTRY;

open IN,"$ARGV[0]" or die "raw .tsv file not found!\n";
my %dbd;
while (<IN>){
	chomp;
	my @a=split /\t/;
    next unless exists $entry{$a[-2]};
	if (defined $dbd{$a[0]}){
		$dbd{$a[0]}{'start'} = $a[6] if $dbd{$a[0]}{'start'} > $a[6];
		$dbd{$a[0]}{'end'} = $a[7] if $dbd{$a[0]}{'end'} < $a[7];
	}else{
		$dbd{$a[0]}{'start'} = $a[6];
		$dbd{$a[0]}{'end'} = $a[7];
	}
	$dbd{$a[0]}{'TF_len'}=$a[2];
	$dbd{$a[0]}{'describe'}=$a[-1];
    $dbd{$a[0]}{'entryid'}=$a[-2];
}
close IN;

open OUT,">$ARGV[0]\.DBD" or die "permission denied!\n";
foreach my $tf (keys %dbd){
    print OUT "$tf\t$dbd{$tf}{'entryid'}\t$dbd{$tf}{'start'}\t$dbd{$tf}{'end'}\t$dbd{$tf}{'TF_len'}\t$dbd{$tf}{'describe'}\n";
}
close OUT;
