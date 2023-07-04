#!/bin/perl -w
use strict;
die "usage: <input domain region info file(6 columns)> <IDRs file>\n" unless @ARGV == 2;

my %idr;
open IDR,"$ARGV[1]" or die "$ARGV[1] not found!\n";
#open IDR,"$idr" or die "IDRs file not found!\n";
while (<IDR>){
	chomp;
	my @a=split /,/;
	my $tf_name=(split /\s+/,$a[1])[1];
	my $n=0;
	foreach (@a[3..$#a]){
		$n++;
		$idr{$tf_name}{$n}=$_;
	}
}
close IDR;

open IN,"$ARGV[0]" or die "file $ARGV[0] not found!\n";
open OUT,">$ARGV[0]\.adjusted.by.IDRs.new.region" or die "permission denied!\n";
open OUT2,">$ARGV[0]\.adjusted.by.IDRs.sta" or die "sta permission denied!\n";
while (<IN>){
    #chomp;
    my @a=split /\t/;
    my ($m,$n,$len_up,$len_down)=&Adjust($a[0],$a[2],$a[3],\%idr);
    $a[2] = $m ;
    $a[3] = $n ;
    my $line = join ("\t",@a);
    print OUT $line;
    print OUT2 "$a[0]\t$m\t$n\t$len_up\t$len_down\n";
}

sub Adjust{
    my ($tf,$start,$end,$idr)=@_;
    my $up_index = $start;
    #print "$idr->{$tf}{$up_index}";
    while ($idr->{$tf}{$up_index} > 0.5 && $up_index < $end){
        $up_index ++;
    }
    my $new_start = $up_index;
    my $len_up_adjust = $new_start - $start;
    my $down_index = $end;
    while ($idr->{$tf}{$down_index} > 0.5 && $down_index > $new_start){
        $down_index += -1;
    }
    my $new_end = $down_index;
    my $len_down_adjust = $end - $new_end;
    my @new= ($new_start,$new_end,$len_up_adjust,$len_down_adjust);
    return @new;
}