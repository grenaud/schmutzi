#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;
use POSIX;
use Cwd 'abs_path';

my $dir = getcwd;
my $abs_path = abs_path($dir);

#ID ID of the sample
#deam5p rate of deamination at the 5' end
#deam3p rate of deamination at the 5' end
#contDeam contamination estimate obtained via deamination patterns
#contDeamL contamination estimate obtained via deamination patterns lower bbound
#contDeamH contamination estimate obtained via deamination patterns upper bound
#contWPrd contamination estimate obtained by predicting the contaminant+DB of putative contaminant
#contWPrdL contamination estimate obtained by predicting the contaminant+DB of putative contaminant lower bound
#contWPrdH contamination estimate obtained by predicting the contaminant+DB of putative contaminant upper bound
#contNPrd contamination estimate obtained by DB only
#contNPrdL contamination estimate obtained by DB only lower bound
#contNPrdH contamination estimate obtained by DB only upper bou
#avgCov average coverage
#hpGrq10w haplogrep for endogenous consensus obtained by predicting the contaminant and QC>10
#hpGrq10wq haplogrep quality for endogenous consensus obtained by predicting the contaminant and QC>10
#hpGrq30w haplogrep for endogenous consensus obtained by predicting the contaminant and QC>30
#hpGrq30wq haplogrep quality for endogenous consensus obtained by predicting the contaminant and QC>30
#hpGrq50w haplogrep for endogenous consensus obtained by predicting the contaminant and QC>50
#hpGrq50wq haplogrep quality for endogenous consensus obtained by predicting the contaminant and QC>50
#hpGrq10n haplogrep for endogenous consensus obtained by DB alone and QC>10
#hpGrq10nq haplogrep quality for endogenous consensus obtained by DB alone and QC>10
#hpGrq30n haplogrep for endogenous consensus obtained by DB alone and QC>30
#hpGrq30nq haplogrep quality for endogenous consensus obtained by DB alone and QC>30
#hpGrq50n haplogrep for endogenous consensus obtained by DB alone and QC>50
#hpGrq50nq haplogrep quality for endogenous consensus obtained by DB alone and QC>50

if($#ARGV==-1){
print "#ID\tdeam5p\tdeam3p\tcontDeam\tcontDeamL\tcontDeamH\tcontWPrd\tcontWPrdL\tcontWPrdH\tcontNPrd\tcontNPrdL\tcontNPrdH\tavgCov\thpGrq10w\thpGrq10wq\thpGrq30w\thpGrq30wq\thpGrq50w\thpGrq50wq\thpGrq10n\thpGrq10nq\thpGrq30n\thpGrq30nq\thpGrq50n\thpGrq50nq\n";
exit;
}

my $filebam    = $ARGV[0];
my @filea      = split("/",$filebam);
#my $filebamname=$filea[$#filea];
my $filebamname=$filebam;


my $name="output";
#print $filebam."\n";
#print $filebamname."\n";


if($#ARGV>0){
  $name=$ARGV[1];
}

#die;


my $line;
print $name."\t";
my $filename5p=$filebamname.".endo.5p.prof";
#die $filename5p;
if (! -e $filename5p) {
  $filename5p="contdeam/".$filebamname.".endo.5p.prof";
}

if (-e $filename5p) {
  open(FILE,$filename5p) or die "Cannot open file $filename5p\n";
  $line=<FILE>;			#header
  while ($line=<FILE>) {
    my @array=split("\t",$line);
    print "".100*sprintf("%.3f",$array[5])."\t";
    last;
  }
  close(FILE);
} else {
  print "NA\t";
}


my $filename3p=$filebamname.".endo.3p.prof";
if (! -e $filename3p) {
  $filename3p="contdeam/".$filebamname.".endo.3p.prof";
}


if (-e $filename3p) {
  open(FILE,$filename3p) or die "Cannot open file $filename3p\n";
  $line=<FILE>;			#header
  while ($line=<FILE>) {
    my @array=split("\t",$line);
    print 100*sprintf("%.3f",$array[6])."\t";
    last;
  }
  close(FILE);
} else {
  print "NA\t";
}



my $filenamedcont=$filebamname.".cont.est";
if (! -e $filenamedcont) {
  $filenamedcont="contdeam/".$filebamname.".cont.est";
}


if (-e $filenamedcont) {
  open(FILE,$filenamedcont) or die "Cannot open file $filenamedcont\n";
  while ($line=<FILE>) {
    my @array=split("\t",$line);
    print 100*sprintf("%.3f",$array[0])."\t".100*sprintf("%.3f",$array[1])."\t".100*sprintf("%.3f",$array[2])."\t";
    last;
  }
  close(FILE);
} else {
  print "NA\tNA\tNA\t";
}


my $filenamescontw=$filebamname."_wtpred_final.cont.est";
if (! -e $filenamescontw) {
  $filenamescontw="wtpred/".$filebamname."_wtpred_final.cont.est";
}


if(-e $filenamescontw){
  open(FILE,$filenamescontw) or die "Cannot open file $filenamescontw\n";
  while($line=<FILE>){
    my @array=split("\t",$line);
    print 100*sprintf("%.3f",$array[0])."\t".100*sprintf("%.3f",$array[1])."\t".100*sprintf("%.3f",$array[2])."\t";
    last;
  }
  close(FILE);
}else{
  print "NA\tNA\tNA\t";
}


my $filenamescontn=$filebamname."_nopred_final.cont.est";
if (! -e $filenamescontn) {
  $filenamescontn="nopred/".$filebamname."_nopred_final.cont.est";
}

if(-e $filenamescontn){
  open(FILE,$filenamescontn) or die "Cannot open file $filenamescontn\n";
  while($line=<FILE>){
    my @array=split("\t",$line);
    print 100*sprintf("%.3f",$array[0])."\t".100*sprintf("%.3f",$array[1])."\t".100*sprintf("%.3f",$array[2])."\t";
    last;
  }
  close(FILE);
}else{
  print "NA\tNA\tNA\t";
}



my $filenamescov=$filebamname."_nopred_final_endo.cov";

if (! -e $filenamescov) {
  $filenamescov="nopred/".$filebamname."_nopred_final_endo.cov";
}

if(-e $filenamescov){
  open(FILE,$filenamescov) or die "Cannot open file $filenamescov\n";
  while($line=<FILE>){
#    warn "#".$line."#\n";
    my @array=split("\t",$line);
    print sprintf("%.3f",$array[0])."\t";
    last;
  }
  close(FILE);
}else{
  print "NA\t";
}





my $filenamesq10hsd=$filebamname."_wtpred_final_endo.q10.hsd";

if (! -e $filenamesq10hsd) {
  $filenamesq10hsd="wtpred/".$filebamname."_wtpred_final_endo.q10.hsd";
}

if(-e $filenamesq10hsd){
  open(FILE,$filenamesq10hsd) or die "Cannot open file $filenamesq10hsd\n";
  $line=<FILE>;
  my $print=0;
  while($line=<FILE>){
    my @array=split("\t",$line);
    #print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    if($#array>3){
      print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    }else{
      print "NA\tNA\t";
    }
    $print=1;
    last;
  }
  close(FILE);
  if($print == 0){
    print "NA\tNA\t";
  }
}else{
  print "NA\tNA\t";
}

my $filenamesq30hsd=$filebamname."_wtpred_final_endo.q30.hsd";

if (! -e $filenamesq30hsd) {
  $filenamesq30hsd="wtpred/".$filebamname."_wtpred_final_endo.q30.hsd";
}


if(-e $filenamesq30hsd){
  open(FILE,$filenamesq30hsd) or die "Cannot open file $filenamesq30hsd\n";
  $line=<FILE>;
  my $print=0;
  while($line=<FILE>){
    my @array=split("\t",$line);
    #print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    if($#array>3){
      print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    }else{
      print "NA\tNA\t";
    }
    $print=1;
    last;
  }
  close(FILE);
  if($print == 0){
    print "NA\tNA\t";
  }

}else{
  print "NA\tNA\t";
}

my $filenamesq50hsd=$filebamname."_wtpred_final_endo.q50.hsd";

if (! -e $filenamesq50hsd) {
  $filenamesq50hsd="wtpred/".$filebamname."_wtpred_final_endo.q50.hsd";
}

if(-e $filenamesq50hsd){
  open(FILE,$filenamesq50hsd) or die "Cannot open file $filenamesq50hsd\n";
  $line=<FILE>;
  my $print=0;
  while($line=<FILE>){
    my @array=split("\t",$line);
    if($#array>3){
      print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    }else{
      print "NA\tNA\t";
    }
    $print=1;
    last;
  }
  close(FILE);
  if($print == 0){
    print "NA\tNA\t";
  }

}else{
  print "NA\tNA\t";
}




my $filenamesq10hsdn=$filebamname."_nopred_final_endo.q10.hsd";


if (! -e $filenamesq10hsdn) {
  $filenamesq10hsdn="nopred/".$filebamname."_nopred_final_endo.q10.hsd";
}

if(-e $filenamesq10hsdn){
  open(FILE,$filenamesq10hsdn) or die "Cannot open file $filenamesq10hsdn\n";
  $line=<FILE>;
  my $print=0;
  while($line=<FILE>){
    my @array=split("\t",$line);
    #print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    if($#array>3){
      print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    }else{
      print "NA\tNA\t";
    }
    $print=1;
    last;
  }
  close(FILE);
  if($print == 0){
    print "NA\tNA\t";
  }

}else{
  print "NA\tNA\t";
}

my $filenamesq30hsdn=$filebamname."_nopred_final_endo.q30.hsd";
if (! -e $filenamesq30hsdn) {
  $filenamesq30hsdn="nopred/".$filebamname."_nopred_final_endo.q30.hsd";
}


if(-e $filenamesq30hsdn){
  open(FILE,$filenamesq30hsdn) or die "Cannot open file $filenamesq30hsdn\n";
  $line=<FILE>;
  my $print=0;
  while($line=<FILE>){
    my @array=split("\t",$line);
    #print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    if($#array>3){
      print $array[2]."\t".100*sprintf("%.3f",$array[3])."\t";
    }else{
      print "NA\tNA\t";
    }
    $print=1;
    last;
  }
  close(FILE);
  if($print == 0){
    print "NA\tNA\t";
  }

}else{
  print "NA\tNA\t";
}

my $filenamesq50hsdn=$filebamname."_nopred_final_endo.q50.hsd";
if (! -e $filenamesq50hsdn) {
  $filenamesq50hsdn="nopred/".$filebamname."_nopred_final_endo.q50.hsd";
}


if(-e $filenamesq50hsdn){
  open(FILE,$filenamesq50hsdn) or die "Cannot open file $filenamesq50hsdn\n";
  $line=<FILE>;
  my $print=0;
  while($line=<FILE>){
    my @array=split("\t",$line);
    if($#array>3){
      print $array[2]."\t".100*sprintf("%.3f",$array[3])."";
    }else{
      print "NA\tNA";
    }
    $print=1;
    last;
  }
  close(FILE);
  if($print == 0){
    print "NA\tNA\t";
  }

}else{
  print "NA\tNA";
}

print "\n";
