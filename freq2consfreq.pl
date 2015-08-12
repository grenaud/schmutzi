#!/usr/bin/perl

#use bignum;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use File::Copy;
use Time::HiRes qw/ time sleep /;


open(FILE,$ARGV[0]) or die "cannot open ".$ARGV[0];
while(my $line = <FILE>){
  chomp($line);
  my @arrtmp=split("\t",$line);
  my $maxF=max( $arrtmp[1],$arrtmp[2],$arrtmp[3],$arrtmp[4] );
  if(0          == $maxF){  print $arrtmp[0]."\t0.0\t0.0\t0.0\t0.0\n"; next; }
  if($arrtmp[1] == $maxF){  print $arrtmp[0]."\t1.0\t0.0\t0.0\t0.0\n"; next; }
  if($arrtmp[2] == $maxF){  print $arrtmp[0]."\t0.0\t1.0\t0.0\t0.0\n"; next; }
  if($arrtmp[3] == $maxF){  print $arrtmp[0]."\t0.0\t0.0\t1.0\t0.0\n"; next; }
  if($arrtmp[4] == $maxF){  print $arrtmp[0]."\t0.0\t0.0\t0.0\t1.0\n"; next; }

  #print $line;
}
close(FILE);
