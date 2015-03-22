#!/usr/bin/perl

#use bignum;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;

my $mock =0;

sub runcmd{
  my ($cmdtorun) = @_;

  print "running cmd ". $cmdtorun."\n";
  if($mock != 1){
    my @argstorun = ( "bash", "-c", $cmdtorun );

    if(system(@argstorun) != 0){
      die "system  cmd $cmdtorun failed: $?"
    }else{
    }
  }

}

sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $$exeFile)) {
    if (!( -e (($$exeFile).".exe")) ) {
      $$exeFile=(($$exeFile).".exe");
    }else{
      die "Executable ".$$exeFile." does not exist\n";
    }
  }
}

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $insertSize  = $pathdir."/insertSize";
my $approxDist  = $pathdir."/approxDist.R";
my $bam2prof    = $pathdir."/bam2prof";
my $contDeam    = $pathdir."/contDeam";
my $contDeamR   = $pathdir."/posteriorDeam.R";
my $splitEndo   = $pathdir."/splitEndoVsCont/poshap2splitbam";

fileExists(\$insertSize);
fileExists(\$approxDist);
fileExists(\$bam2prof);
fileExists(\$contDeam);
fileExists(\$contDeamR);
fileExists(\$splitEndo);

#
# protocol s or d
# output prefix

my $lengthDeam=2;

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\nThis script is a wrapper to call the core program and allows for the estimation of endogenous\ndeamination andsubsequent contamination estimate. By default it will condition on seeing a \ndeaminated base on the 5' end to measure endogenous deamination rates on the 3' end and\nvice-versa. Users have also the possibility of conditioning on diagnostic positions to \nmeasure endogenous deamination.\n\nusage:\t".$0." <options> input.bam\n\n".
"Options:\n".
"\n\t--library (single|double)\tType of library used".
"\n\t--mock\t\t\t\tDo nothing, just print the commands used\n".

"\t--lengthDeam (bp)\t\tOnly consider this about of bases to be deaminated on each end (Default : $lengthDeam)\n\n".

"Output Options:\n".
"\t--out (output prefix)\t\tAll output files will share this prefix\n".
"\t--title (title)\t\t\tTitle for the graph of the posterior distribution\n".
"\t--cont (cont)\t\t\tIf you have prior knowledge about the contamination\n\t\t\t\t\trate, enter it here [0-1]\n".
"\nSplit into endogenous/contaminant options:\n".
"\t--split (file)\t\t\tSplit endogenous/contaminant according to diagnostic positions\n".
"\t\t\t\t\tThe file must have the following format:\n".
"\t\t\t\t\t\t[coord]tab[nucleotide]tab[endo or cont]\n".
"\t\t\t\t\tWhere the coordinate is on the reference genome\n".
"\t\t\t\t\tex:\t385\tA\tendo\n".
"\t--splitdeam\t\t\tEstimate the deamination rates using the reads split into endogenous/contaminant\n".
"\n\t--uselength\t\t\tUse length of the molecules as well\n\n".
  "\nOther options:\n".
"\t--ref (reference genome)\tThe fasta file used for alignment\n".

#"\t--help|-?".
"\n\n";
  exit;
}

my $help;
my $library        = "none";
my $outputPrefix   = "outputdeam";
my $inbam          = "none";
my $referenceFasta = "";
my $contPriorKnow  = -1;
my $textGraph      = "Posterior probability for contamination\\nusing deamination patterns";
my $splitPos       = "";
my $useLength      = 0;
my $splitDeam      = 0;

usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'split=s' => \$splitPos,'library=s' => \$library,'ref=s' => \$referenceFasta,'title=s' => \$textGraph,'cont=f' => \$contPriorKnow, 'out=s' => \$outputPrefix,'mock' => \$mock, 'uselength' => \$useLength,'lengthDeam=i' => \$lengthDeam,'splitDeam' => \$splitDeam )
          or defined $help );


#die $contPriorKnow;
#print $data."\n";
#run deam

if($library ne "single" &&
   $library ne "double" ){
  die "Please enter the type of library as either double or single\n";
}

#if($referenceFasta eq "none" ){
#  die "Please enter the fasta reference used for mapping\n";
#}

$inbam = $ARGV[ $#ARGV ];

my $scaleLocSpecified=0;
my $locE=1;
my $scaE=1;
my $locC=1;
my $scaC=1;


my $cmdBam2Prof = $bam2prof." -length  $lengthDeam -endo -".$library." -5p ".$outputPrefix.".endo.5p.prof  -3p ".$outputPrefix.".endo.3p.prof $inbam";

runcmd($cmdBam2Prof);

if($splitPos ne ""){

  print "Examining file ".$splitPos."....\n";
  open(FILEdiag,$splitPos) or die "cannot open ".$splitPos;

  while (my $line = <FILEdiag>) {
    chomp($line);

    my @array = split("\t",$line);
    if ($#array != 2) {
      die "Line ".$line." does not have 3 fields\n";
    }

    if ($array[2] ne "endo"  &&
	$array[2] ne "cont"  ) {
      die "The third field in ".$line." is not either \"endo\" or \"cont\"\n";
    }

  }
  close(FILEdiag);
  print ".... fine\n";

  my $cmdBamSplit = $splitEndo."  ".$splitPos." $inbam ".$outputPrefix." > ".$outputPrefix."_split.log 2> /dev/null ";
  runcmd($cmdBamSplit);
  #evaluate deamination

  if($splitDeam){
    my $cmdBam2ProfEndo = $bam2prof." -length $lengthDeam -".$library." -5p ".$outputPrefix.".endo.5p.prof  -3p ".$outputPrefix.".endo.3p.prof ".$outputPrefix."_endo.bam";
    runcmd($cmdBam2ProfEndo);

    my $cmdBam2ProfCont = $bam2prof." -length $lengthDeam -".$library." -5p ".$outputPrefix.".cont.5p.prof  -3p ".$outputPrefix.".cont.3p.prof ".$outputPrefix."_cont.bam";
    runcmd($cmdBam2ProfCont);
  }

  if ($useLength) {
    #evaluate size
    #my $inserSize  = $pathdir."/insertSize";
    my $cmdBam2InsertsizeEndo = $insertSize."  ".$outputPrefix."_endo.bam |gzip  > ".$outputPrefix."_endo.size.gz";
    runcmd($cmdBam2InsertsizeEndo);
    my $cmdBam2InsertsizeCont = $insertSize."  ".$outputPrefix."_cont.bam |gzip  > ".$outputPrefix."_cont.size.gz";
    runcmd($cmdBam2InsertsizeCont);

    #Log normal fit
    my $cmdinsertsize2LognormEndo = $approxDist."  ".$outputPrefix."_endo.size.gz >  ".$outputPrefix."_endo.size.param";
    runcmd($cmdinsertsize2LognormEndo);
    my $cmdinsertsize2LognormCont = $approxDist."  ".$outputPrefix."_cont.size.gz >  ".$outputPrefix."_cont.size.param";
    runcmd($cmdinsertsize2LognormCont);


    if ($mock != 1) {


      open(FILEparamEndo,$outputPrefix."_endo.size.param") or die "cannot open ".$outputPrefix."_endo.size.param";
      my @linesEndo = <FILEparamEndo>;
      close(FILEparamEndo);

      if ($#linesEndo != 2) {
	die "ERROR: Parameter file from the log normal distribution does not have 3 lines, it has ".$#linesEndo." lines check the size distribution of the endogenous and contaminant molecules";
      }

      my $lineParamEndo= $linesEndo[1];
      chomp($lineParamEndo);
      $lineParamEndo =~ s/^\s+//;
      $lineParamEndo =~ s/\s+$//;

      my @arrayEndo = split(/\s+/,$lineParamEndo);

      if (!looks_like_number( $arrayEndo[0] )) {
	die "ERROR: Parameter file from the log normal distribution contains a line that does not have the expected numerical parameters: ".$lineParamEndo."\n";
      }
      if (!looks_like_number( $arrayEndo[1] )) {
	die "ERROR: Parameter file from the log normal distribution contains a line that does not have the expected numerical parameters: ".$lineParamEndo."\n";
      }

      $locE=$arrayEndo[0];
      $scaE=$arrayEndo[1];



      open(FILEparamCont,$outputPrefix."_cont.size.param") or die "cannot open ".$outputPrefix."_cont.size.param";
      my @linesCont = <FILEparamCont>;
      close(FILEparamCont);

      if ($#linesCont != 2) {
	die "ERROR: Parameter file from the log normal distribution does not have 3 lines, it has ".$#linesCont." lines check the size distribution of the contgenous and contaminant molecules";
      }

      my $lineParamCont= $linesCont[1];
      chomp($lineParamCont);
      $lineParamCont =~ s/^\s+//;
      $lineParamCont =~ s/\s+$//;

      my @arrayCont = split(/\s+/,$lineParamCont);

      if (!looks_like_number( $arrayCont[0] )) {
	die "ERROR: Parameter file from the log normal distribution contains a line that does not have the expected numerical parameters: ".$lineParamCont."\n";
      }
      if (!looks_like_number( $arrayCont[1] )) {
	die "ERROR: Parameter file from the log normal distribution contains a line that does not have the expected numerical parameters: ".$lineParamCont."\n";
      }

      $locC=$arrayCont[0];
      $scaC=$arrayCont[1];

      $scaleLocSpecified=1;
    }
  }

}

my $cmdcontdeam = $contDeam." ";
if($scaleLocSpecified){
  $cmdcontdeam.=" --loce ".$locE." --scalee ".$scaE." --locc ".$locC." --scalec ".$scaC." ";
}

$cmdcontdeam.=" -deamread -deam5p ".$outputPrefix.".endo.5p.prof  -deam3p ".$outputPrefix.".endo.3p.prof  -log  ".$outputPrefix.".cont.deam ";
if($referenceFasta ne ""){
  $cmdcontdeam.=" -r ".$referenceFasta." ";
}
$cmdcontdeam.="  $inbam";

runcmd($cmdcontdeam);


#print "Wrapper finished succesfully output available as ".$outputPrefix.".cont.deam\n";
#
open(FILE,$outputPrefix.".cont.deam") or die "cannot open ".$outputPrefix.".cont.deam";
my $maxL=10;
my $maxLi=0;

my @arrayOfValues;
while(my $line = <FILE>){
  chomp($line);
  #print $line;
  my @array = split("\t",$line);
  if($maxL == 10){
    $maxL = $array[2];
  }

  if($maxL < $array[2]){
    $maxL = $array[2];
    $maxLi = ($#arrayOfValues+1);
  }
  my $hashVal =  { 'cont'   => $array[1],
		   'logL'   => $array[2]};

  push(@arrayOfValues,$hashVal);
}
close(FILE);


my $sum=0;
foreach my $hashVal (@arrayOfValues){

  $hashVal->{'logLS'}=$hashVal->{'logL'}-$maxL;
 # print "".($hashVal->{'logLS'})."\t".(10**($hashVal->{'logLS'}))."\t".$sum."\n";
  $sum+=(10**($hashVal->{'logLS'}));
}
#print $sum;
my $targetSum = 0.95 * $sum;


my $il=$maxLi;
my $ih=$maxLi;

while(1){
  $il=max($il-1,              0);
  $ih=min($ih+1, $#arrayOfValues+1 );
  #print "i ".$il."\t".$ih."\n";
  my $subsum = 0;
  for(my $i=$il;$i<=$ih;$i++){
    $subsum += (10**($arrayOfValues[$i]->{'logLS'}));
  }
  #print $targetSum."\t".$subsum."\t".$il."\t".$ih."\n";

  if($subsum>$targetSum){
    last;
  }

}

if ($mock != 1) {
  open(FILEOUT,">".$outputPrefix.".cont.est") or die "cannot write to ".$outputPrefix.".cont.est";
  print FILEOUT $arrayOfValues[$maxLi]->{'cont'}."\t".$arrayOfValues[$il]->{'cont'}."\t".$arrayOfValues[$ih]->{'cont'}."\n";
  close(FILEOUT);

  open(FILECONFIGOUT,">".$outputPrefix.".deam.config") or die "cannot write to ".$outputPrefix.".deam.config";
  print FILECONFIGOUT "library\t$library\n";
  print FILECONFIGOUT "outputPrefix\t$outputPrefix\n";
  #print FILECONFIGOUT "inbam\t$inbam\n";
  print FILECONFIGOUT "referenceFasta\t$referenceFasta\n";
  print FILECONFIGOUT "contPriorKnow\t$contPriorKnow\n";
  print FILECONFIGOUT "textGraph\t$textGraph\n";
  print FILECONFIGOUT "splitPos\t$splitPos\n";
  print FILECONFIGOUT "useLength\t$useLength\n";
  print FILECONFIGOUT "splitDeam\t$splitDeam\n";
  close(FILECONFIGOUT);

}


my $cmdPlot = $contDeamR." ".$outputPrefix.".cont.deam ".$outputPrefix.".cont.pdf  \"$textGraph\" ";
if($contPriorKnow != -1){
  $cmdPlot =  $cmdPlot." ".$contPriorKnow;
}

runcmd($cmdPlot);

print "Program finished succesfully\n\nFiles created:".
  "The plot of the posterior probability is ".$outputPrefix.".cont.pdf\n".
  "The contamination estimate is here ".$outputPrefix.".cont.est\n".
  "The configuration file is here ".$outputPrefix.".config\n";

