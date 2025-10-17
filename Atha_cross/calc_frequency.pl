#!/usr/bin/perl
use strict;
use warnings;

while (<STDIN>) {
    chomp;
    my @f = split(/\t/);
    next unless @f >= 6;

    my ($chrom, $pos, $ref, $alt_str, $ad_str, $qual) = @f[0..5];
    $qual = 0 unless (defined $qual && $qual =~ /^[-+]?\d+(\.\d+)?$/);

    my @alts = split(/,/, $alt_str // "");
    my @ad = split(/,/, $ad_str // "");

    next unless @ad == @alts + 1;   # 确保 AD 长度正确

    my $ref_cnt = $ad[0] // 0;
    my $total = $ref_cnt + 0;
    $total += $_ for @ad[1..$#ad];
    next if $total == 0;

    for my $i (0..$#alts) {
        my $alt_cnt = $ad[$i+1] // 0;
        my $freq = $alt_cnt / $total;
        next if $freq < 0.01;   # 过滤 <1%
        printf "%s\t%s\t%s\t%s\t%.6f\t%d\n", $chrom, $pos, $ref, $alts[$i], $freq, $total;
    }
}
