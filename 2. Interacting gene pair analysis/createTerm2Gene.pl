#!/usr/bin/perl -w
use strict;

open(I,"<$ARGV[0]");
my %s=();
open(O,">gid.go");
print O "GID\tGO\tEVIDENCE\n";
open(O2,">gname.go");
print O2 "GID\tSYMBOL\tGENENAME\n";
my %gid=();
my $c=0;
while (<I>) {
	chomp;
	my @tmp = split /\t/;
	if (not exists $gid{$tmp[0]}) { $c++; $gid{$tmp[0]}=$c; print O2 "$c\t$tmp[0]\t$tmp[0]\n" }
	while (/(GO:\d+)/g) {
		if (exists $s{"$1-$tmp[0]"}) { next }
		print "$1\t$tmp[0]\n";
		print O "$gid{$tmp[0]}\t$1\tIEA\n";
		$s{"$1-$tmp[0]"}=1;
	}
}
close I;

close O;
close O2;
