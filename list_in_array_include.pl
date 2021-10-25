#!/usr/bin/perl

use warnings;
use strict;

my $sub = shift or die $!;
my $main = shift or die $!;


my %hash;
open IN, "$main" or die $!;
while(<IN>){
	chomp;
	my @tmp =  split /\t/;
	$hash{$tmp[0]} = $_;
}
close IN;

open IN, "$sub" or die $!;
while(<IN>){
	chomp;
	my @tmp = split /\t/;
	if(defined $hash{$tmp[0]}){
		print "$hash{$tmp[0]}\n";
	}else{
		next;
	}
}
close IN;
