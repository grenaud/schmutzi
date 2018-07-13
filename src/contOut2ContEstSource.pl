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

my $mock =0;
my $lengthDeam=2;

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}




sub computeIntervalCont{
  my ($mtcontOutputLog,$source) = @_;

  #print "reading the contamination log ".$mtcontOutputLog."\n";

  if ($mock != 1) {
    open(FILE,$mtcontOutputLog) or die "cannot open ".$mtcontOutputLog."";
    my $maxL=10;
    my $maxLi=0;

    my @arrayOfValues;
    while (my $line = <FILE>) {
      chomp($line);
      #print $line;

      my @array = split("\t",$line);
      if($array[0] eq $source){
      if ($maxL == 10) {
	$maxL = $array[2];
      }

      if ($maxL < $array[2]) {
	$maxL = $array[2];
	$maxLi = ($#arrayOfValues+1);
      }
      my $hashVal =  { 'cont'   => $array[1],
		       'logL'   => $array[2]};

      push(@arrayOfValues,$hashVal);
    }
    }
    close(FILE);


    my $sum=0;
    foreach my $hashVal (@arrayOfValues) {

      $hashVal->{'logLS'}=$hashVal->{'logL'}-$maxL;
      # print "".($hashVal->{'logLS'})."\t".(10**($hashVal->{'logLS'}))."\t".$sum."\n";
      $sum+=(10**($hashVal->{'logLS'}));
    }
    #print $sum;
    my $targetSum = 0.95 * $sum;


    my $il=$maxLi;
    my $ih=$maxLi;

    while (1) {
      $il=max($il-1,              0);
      $ih=min($ih+1, $#arrayOfValues );
      #print "i ".$il."\t".$ih."\n";
      my $subsum = 0;
      for (my $i=$il;$i<=$ih;$i++) {
	$subsum += (10**($arrayOfValues[$i]->{'logLS'}));
      }
      #print $targetSum."\t".$subsum."\t".$il."\t".$ih."\n";

      if ($subsum>$targetSum) {
	last;
      }

    }

    return ($arrayOfValues[$maxLi]->{'cont'},$arrayOfValues[$il]->{'cont'},$arrayOfValues[$ih]->{'cont'});

  }else{
    #mock, nothing to do
    return (-1,-1,-1);
  }

}



sub readMTCONTRetainMostlikely{
  my ($mtcontOutputLog) = @_;

#  print "reading the log file ".$mtcontOutputLog." and print most likely source to ".$outputToWrite."\n";


  my $sourceALL   = -1;
  my $contMaxSALL = -1;
  my $llikMaxSALL = -1;
  my $firstRecord = 1;
  if ($mock != 1) {

    open(FILEoutlogmtcont, $mtcontOutputLog) or die "cannot open ".$mtcontOutputLog."\n";

    while(my $line = <FILEoutlogmtcont>){
      my @arrayTemp = split("\t",$line);
      if($#arrayTemp != 2){
	die "Error line $line in $mtcontOutputLog does not have 3 fields";
      }


      if( $firstRecord == 1 ){#first record
	$sourceALL   = $arrayTemp[0];
	$contMaxSALL = $arrayTemp[1];
	$llikMaxSALL = $arrayTemp[2];
	$firstRecord = 0;
      }else{

	if($llikMaxSALL < $arrayTemp[2]){
	  $sourceALL   = $arrayTemp[0];
	  $contMaxSALL = $arrayTemp[1];
	  $llikMaxSALL = $arrayTemp[2];
	}

      }

    }

    close(FILEoutlogmtcont);

    if ( $contMaxSALL == -1) {
      die "Cannot parse the contamination output file = ".$mtcontOutputLog."\n";
    }

    return $sourceALL;
  } else {
    #nothing to do
  }

}






my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);




my $contFile = $ARGV[ $#ARGV  ];
my $source = readMTCONTRetainMostlikely($contFile);
my @arrayCont = computeIntervalCont($contFile,$source);
my @artemp = split("/",$source);
$source = $artemp[ $#artemp ];
print $source."\t".$arrayCont[0]."\t".$arrayCont[1]."\t".$arrayCont[2]."\n";


