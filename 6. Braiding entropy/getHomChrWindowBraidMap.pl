#!/usr/bin/perl -w
use strict;
use Statistics::RankCorrelation;

if ($#ARGV==-1) {
 die(" Usage: [psynt.tab] [pairwise alignment] [window size] [minlines]\n");
}

my $win=$ARGV[2];
my $minlines=$ARGV[3];

my %hom=();
open(I,"<$ARGV[0]");
while (<I>) {
	chomp;
	my @tmp = split /\t/;
	if ($tmp[0]=~/(Lachesis_group\d+)\_\_/) { $tmp[0]=$1 }
        if ($tmp[1]=~/(Lachesis_group\d+)\_\_/) { $tmp[1]=$1 }
	$hom{$tmp[0]}{$tmp[1]}=1;
}
close I;

open(I,"<$ARGV[1]");
my %data=();
while (<I>) {
	chomp;
	my @tmp = split /\t/;
	if (exists $hom{$tmp[0]}{$tmp[3]}) {
		push @{$data{$tmp[0]}{$tmp[3]}}, [ ( @tmp ) ];
	}
}
close I;

my $ID=0;
for my $x (keys %data) {
 for my $y (keys %{$data{$x}}) {
  my @homs=@{$data{$x}{$y}};
  my @shom=sort { ${$a}[1] <=> ${$b}[1] } @homs;
	  my $pos=-1;
  while ($pos<$#shom) {
	  $pos++;
   my $stx=$shom[$pos][1];
   my $dist=0;
   my $posend=$pos+1;
   my @poswin=();
   push @poswin, $pos;
   while (($dist<=$win)&&($posend<=$#shom)) {
	   $dist=$shom[$posend][1]-$stx;
	   if ($dist<=$win) { push @poswin, $posend }
	   $posend++;
   }
   if (($#poswin+1)>=$minlines) {
    #hom cluster search
    my %hit=();
    for my $hompos (@poswin) {
	    for my $neighbor (@poswin) {
	     if (abs($shom[$hompos][4]-$shom[$neighbor][4])<=$win) { push @{$hit{$hompos}}, $neighbor }
	    }
    }
    my @shit=reverse sort {$#{$hit{$a}} <=> $#{$hit{$b}}} keys %hit;
    if (($#{$hit{$shit[0]}}+1)>=$minlines) {
     #braid map 
     $pos=$posend-1;
     $ID++;
      my @bdata=();
      for my $e (@{$hit{$shit[0]}}) { push @bdata, [ @{$shom[$e]} ] } 
      my @pos1=sort {$bdata[$a][1] <=> $bdata[$b][1]} 0..$#bdata;
      my @pos2=sort {$bdata[$a][4] <=> $bdata[$b][4]} 0..$#bdata;
	for my $xx (0..$#pos1) {
 	 print "$ID\t$x\t$bdata[$pos1[0]][1]\t$bdata[$pos1[$#pos1]][1]\t$y\t$bdata[$pos2[0]][4]\t$bdata[$pos2[$#pos2]][4]\t$x:$bdata[$pos1[0]][1]..$bdata[$pos1[$#pos1]][1]\t$y:$bdata[$pos2[0]][4]-$bdata[$pos2[$#pos2]][4]\t".($xx+1)."\t".$pos1[$xx]."\t$pos2[$xx]\n";
	}
      my $cor = Statistics::RankCorrelation->new( \@pos1, \@pos2, sorted => 1 );
      print STDERR "$x\t$y\t$ID\t".$cor->spearman."\t".$cor->kendall."\t".$cor->csim."\t".$cor->size."\n";
    }
   }
  }
 }
}
