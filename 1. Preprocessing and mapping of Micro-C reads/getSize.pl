#!/usr/bin/perl -w
use strict;

open(I,"<$ARGV[0]");
my %len=();
my $n="";
while (<I>) {
 chomp;
 if (/^>([^\s]*)/) { $n=$1 } else { $len{$n}+=length $_ }
}
close I;

for my $x (keys %len) {
 print "$x\t$len{$x}\n";
}
