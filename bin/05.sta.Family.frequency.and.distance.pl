#!/bin/perl -w
use strict;
die "usage: <input filtered basic patches file> <DBD info file> <family sta file>\n" unless @ARGV == 3;

### read DBD info ###
#my $dbd='/home/dongke/TF_RNA/prediction.of.TF.DNA.binding.domain/InterPro/Ath_pep.fas.tsv.ProSiteProfiles';
my $dbd= $ARGV[1];
open DBD,"$dbd" or die "file $dbd not found!\n";
my %dbd;
while (<DBD>){
	chomp;
	my @a=split /\t/;
	$dbd{$a[0]}{'start'} = $a[2];
	$dbd{$a[0]}{'end'} = $a[3];
	$dbd{$a[0]}{'TF_len'}=$a[4];
	$dbd{$a[0]}{'family'}=$a[5];
}
my $tf_counts_with_dbd = () = keys %dbd;
close DBD;

### read RK score sta file ###
open IN,"$ARGV[0]" or die "input filtered basic patches file not found!\n";
open OUT,">$ARGV[0]\.distance.to.dbd.lt.20.family.sta" or die "permission denied!\n";
open OUT2,">$ARGV[0]\.distance.to.dbd.lt.20.with.region.info" or die "permission denied!\n";
my $name;
my $switch;
my $tf_with_dbd;
my $total_gene;
my $win_o;
my %pro;
my $ordered_win;
my $disorded_lt20;

my %tf_with_bb;

my %tf_with_bp_round_20aa;
my $bp_sum_count;
my %dis_pro;
my $head;
my $family;
READ:
while (<IN>){
	if (/^>/){
        $head = $_;
		s/^>//;
        chomp;
        if (/\|/){
            $family = (split /\|/,$_)[1];
        }else{
            $family = (split /\t/,$_)[-1];
        }
        #print "$family\n";
		$name=(split /\s+/)[0];
		$total_gene++;
		$tf_with_dbd++ if exists $dbd{$name};
		#print OUT2 ">$_\n";
		$switch=1;
		next;
	}
	$bp_sum_count++;
	chomp;
	my $info=$_;
	my @bp=split /:/,$info;
    my $start=$bp[0];
    my $end=$bp[1];
	my $dis;
    my $weather_lt20=0;
	if (exists $dbd{$name}){
		my $m=$dbd{$name}{'start'};
		my $n=$dbd{$name}{'end'};
        my $len=$dbd{$name}{'TF_len'};
        foreach (1..$m){
            $dis_pro{$name}{$m-$_}=0;
        }
        foreach ($n..$len){
            $dis_pro{$name}{$_-$n}=0;
        }
		my $region="$m\-$n";
		next if $n-$m <= 10;
        foreach ($start..$end){
            my $loc = $_;
            if ($loc>$m && $loc<$n){
                next;
            }elsif( $loc<=$m){
                $dis=$m-$loc;
                $weather_lt20=1 if $dis <=20;
                $dis_pro{$name}{$dis}=1;
            }elsif( $loc>=$n){
                $dis=$loc-$n;
                $weather_lt20=1 if $dis <=20;
                $dis_pro{$name}{$dis}=1;
            }
        }
        if ($weather_lt20==1){
            $tf_with_bp_round_20aa{$name}=$family;
            if ($switch == 1){
                print OUT2 $head;
                print OUT2 "$info\t$region\n";
                $switch = 0;
            }else{
                print OUT2 "$info\t$region\n";
            }
        }
    }
}

close IN;


### print family sta ###
#my
#$family_file='/home/dongke/TF_RNA/fasta.and.gene.info.of.whole.proteome.and.TF/TF.family.sta';
my $family_file=$ARGV[2];
open FAM,"$family_file" or die "file $family_file not found!\n";
my %fam;
while (<FAM>){
    chomp;
    my @a=split;
    $fam{$a[0]}=$a[1];
}
close FAM;
my %family;
foreach (keys %tf_with_bp_round_20aa){
    $family{$tf_with_bp_round_20aa{$_}}++;
}
foreach (keys %family){
    my $percentage_of_family=sprintf ("%.2f",$family{$_} / $fam{$_}) if exists $fam{$_};
	print OUT "$_\t$family{$_}\t$percentage_of_family\n";
}



### print distance sta ###
open PRO,">$ARGV[0]\.distance.sta.probability" or die "permission denied!\n";
my $tf_n;
my $pro_sum;
foreach my $dis (0..50){
    foreach my $tf (keys %dis_pro){
        if (defined $dis_pro{$tf}{$dis}){
            $pro_sum+=$dis_pro{$tf}{$dis};
            $tf_n ++;
        }
    }
    my $probability = sprintf ("%.2f",$pro_sum / $tf_n) if $tf_n != 0;
    print PRO "$dis\t$probability\n";
}
close PRO;


### print TFs with basic batches IDs and their family info ###
#open GENE,">$ARGV[0]\.TFs.with.basic.batches.ID.family.sta";
#my @tf_with_bb = keys %tf_with_bb;
#my $tf_with_bb = @tf_with_bb;
#foreach (@tf_with_bb){
#    my $tfs=$_;
#    print GENE "$tfs\t$tf_with_bb{$tfs}{'bbcount'}\t$tf_with_bb{$tfs}{'family'}\n";
#}
#close GENE;


### print windows count and TF counts info to stdout ###
#open LOG,">$ARGV[0]\.out.log" or die "permission denied!\n";
#print LOG "Total TFs with DBD counts $tf_counts_with_dbd\n";
#print LOG "$tf_with_dbd in whole input $total_gene TFs found with DBD info!\n";
#print LOG "we count total $win_a windows and $win_o windows found outside DBD,including $ordered_win windows in ordered region. Finaly, $tf_with_bb TFs were detected by $disorded_lt20 basic batches windows which are in IDRs and at the distance to DBD less than 20 aa.\n";

=pod
sub Min{
	my @a=@_;
	my $min = shift @a;
	foreach (@a){
		$min = $_ if $_ < $min;
	}
	return $min;
}

sub Max{
	my @a=@_;
        my $max = shift @a;
        foreach (@a){
                $max = $_ if $_ > $max;
        }
        return $max;
}
=cut