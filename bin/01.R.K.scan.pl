#!/bin/perl -w
use strict;
die "usage:<input .fa file>\n" unless @ARGV == 1;
open IN,"$ARGV[0]" or die "file $ARGV[0] not found!\n";
open OUT,">RK.raw.basic.patches.info" or die "permission denied!\n";
local $/='>';
<IN>;
while (<IN>){
	s/\r?\n>?$//;
	my ($head,$seq)=split /\r?\n/,$_,2;
	my $switch = 1;
	$seq=~s/\r?\n//g;
	my $start=0;
    my $end = 5;
    my $bp= 0 ;
    my $bpseq='';
	my $len=length $seq;
	while ($end < $len){
		my $cut=substr($seq,$start,$end-$start+1);
		my $score= &RKratio($cut);
        if ($score > 0.5){
            $bp=1;
            $end+=1;
            $bpseq=$cut;
        }elsif($bp==1){
            my $bplen=length $bpseq;
            if ($switch == 1){
                print OUT ">$head\n$start:$end:$bpseq:$bplen\n";
                $switch = 0;
            }else{
                print OUT "$start:$end:$bpseq:$bplen\n";
            }
            #$start += 1 ;
            $start = $end + 1 ;
            $end = $start + 4;
            $bp = 0;
        }else{
            $start += 1;
            $end = $start + 4;
        }
	}
}
$/="\n";
close IN;

sub RKratio{
    my ($seq)=@_;
    my $c=()=$seq=~/[RK]/ig;
    my $len = length $seq;
    my $ratio =sprintf ("%.1f" ,$c / $len) if $len != 0;
    return $ratio;
}
