#!/usr/bin/env perl

use strict;
use warnings;

my $header = <STDIN>;
chomp $header;
my @h = split("\t", $header);

my @keep_i = (0,1);
my @n = ();

my @keep = ('Quant\. meas\. - Date_(\d+)','Quant\. meas\. - Volume_(\d+)');
	  #T:Quant. meas. - Date_0
for (my $i=1; $i < scalar(@h); $i++) {
	foreach my $regex (@keep) {
		if ($h[$i] =~ /$regex/) {
			push(@keep_i, $i);
			push(@n, $1);
			#print $h[$i] . "\t" . $1 . "\n";
		}
	}
}

# Consistency check of having found same n. of Dates and Volumes
for (my $i = 0; $i < scalar(@n) -1; $i += 2) {
	die "something nasty! $n[$i] $n[$i+1]\n" if $n[$i] != $n[$i+1];
}

while(<STDIN>) {
	chomp;
	my @line = split("\t", $_);
	my @res = ();
	foreach my $col (@keep_i) {
		if ($line[$col] ne "") {
			push(@res, $line[$col]);
		}
	}
	my $gen = shift @res;
	my $exp = shift @res;
	for (my $i = 0; $i < scalar(@res) -1; $i += 2) {
		print $gen . "\t" . $exp . "\t" .  $res[$i] . "\t" . $res[$i+1] . "\n";
	}
}
