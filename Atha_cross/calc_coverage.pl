#!/usr/bin/perl
use strict;
use warnings;

my ($sum, $count) = (0, 0);
while (<STDIN>) {
    my @f = split;
    next unless @f >= 3 && $f[2] =~ /^\d+$/;
    $sum += $f[2];
    $count++;
}
printf "%.2f\n", $count > 0 ? $sum/$count : 0;
