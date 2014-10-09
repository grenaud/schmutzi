#!/usr/bin/perl

#use bignum;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);

my $mock =0;

sub runcmd{
  my ($cmdtorun) = @_;

  print "running cmd ". $cmdtorun."\n";
  if($mock != 1){
    if(system($cmdtorun) != 0){
      die "system  cmd $cmdtorun failed: $?"
    }else{
    }
  }

}

sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}
my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $inserSize  = $pathdir."/insertSize";
my $approxDist = $pathdir."/approxDist.R";
my $bam2prof   = $pathdir."/bam2prof";
my $contDeam   = $pathdir."/contDeam";
my $contDeamR  = $pathdir."/posteriorDeam.R";
my $splitEndo  = $pathdir."/splitEndoVsCont/poshap2splitbam";

fileExists($insertSize);
fileExists($approxDist);
fileExists($bam2prof);
fileExists($contDeam);
fileExists($contDeamR);
fileExists($splitEndo);

#
# protocol s or d
# output prefix

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\nThis script is a wrapper that allows for the estimation of endogenous deamination and\nsubsequent contamination estimate. By default it will condition on seeing a deaminated\nbase on the 5' end to measure endogenous deamination rates on the 3' end and vice-versa.\nUsers have also the possibility of conditioning on diagnostic positions to measure\nendogenous deamination.\n\nusage:\t".$0." <options> input.bam\n\n".
"Options:\n".
"\n\t--library (single|double)\tType of library used".
"\n\t--mock\t\t\t\tDo nothing, just print the commands used\n\n".

"Output Options:\n".
"\t--out (output prefix)\t\tAll output files will share this prefix\n".
"\t--title (title)\t\t\tTitle for the graph of the posterior distribution\n".
"\t--cont (cont)\t\t\tIf you have prior knowledge about the contamination\n\t\t\t\t\trate, enter it here [0-1]\n".
"\nInput options:\n".
"\t--split (file)\t\t\tSplit endogenous/contaminant according to diagnostic positions\n".
"\t\t\t\t\tThe file must have the following format:\n".
"\t\t\t\t\t\t[coord]tab[nucleotide]tab[endo or cont]\n".
"\t\t\t\t\tWhere the coordinate is on the reference genome\n".
"\t\t\t\t\tex:\t385\tA\tendo\n".
"\nMandatory:\n".
"\t--ref (reference genome)\tThe fasta file used for alignment\n".

#"\t--help|-?".
"\n\n";
  exit;
}

my $help;
my $library        = "none";
my $outputPrefix   = "outputdeam";
my $inbam          = "none";
my $referenceFasta = "none";
my $contPriorKnow  = -1;
my $textGraph      = "Posterior probability for contamination\nusing deamination patterns";
my $splitPos          = "";

usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'split=s' => \$splitPos,'library=s' => \$library,'ref=s' => \$referenceFasta,'title=s' => \$textGraph,'cont=f' => \$contPriorKnow, 'out=s' => \$outputPrefix,'mock' => \$mock )
          or defined $help );


#die $contPriorKnow;
#print $data."\n";
#run deam

if($library ne "single" &&
   $library ne "double" ){
  die "Please enter the type of library as either double or single\n";
}

if($referenceFasta eq "none" ){
  die "Please enter the fasta reference used for mapping\n";
}

$inbam = $ARGV[ $#ARGV ];


if($splitPos ne ""){
  my $cmdBam2Prof = $bam2prof." -endo -".$library." -5p ".$outputPrefix.".5p.prof  -3p ".$outputPrefix.".3p.prof $inbam";

  runcmd($cmdBam2Prof);
} else {
  print "Examining file ".$splitPos."....\n";
  open(FILEdiag,$splitPos) or die "cannot open ".$splitPos;

  while (my $line = <FILEdiag>) {
    chomp($line);

    my @array = split("\t",$line);
    if($#array != 2){
      die "Line ".$line." does not have 3 fields\n";
    }

    if($array[2] ne "endo"  &&
       $array[2] ne "cont"  ){
      die "The third field in ".$line." is not either \"endo\" or \"cont\"\n";
    }

  }
  close(FILEdiag);
  print ".... fine\n";

  my $cmdBamSplit = $splitEndo."  ".$splitPos." ".$outputPrefix;
  runcmd($cmdBamSplit);
  #evaluate deamination
  my $cmdBam2ProfEndo = $bam2prof."  -".$library." -5p ".$outputPrefix.".endo.5p.prof  -3p ".$outputPrefix.".endo.3p.prof ".$outputPrefix."_endo.bam";
  runcmd($cmdBam2ProfEndo);

  my $cmdBam2ProfCont = $bam2prof."  -".$library." -5p ".$outputPrefix.".cont.5p.prof  -3p ".$outputPrefix.".cont.3p.prof ".$outputPrefix."_cont.bam";
  runcmd($cmdBam2ProfCont);
  
  #evaluate size
my $inserSize  = $pathdir."/insertSize";
my $approxDist = $pathdir."/approxDist.R";

}

my $cmdcontdeam = $contDeam." -deamread -deam5p ".$outputPrefix.".5p.prof  -deam3p ".$outputPrefix.".3p.prof  -log  ".$outputPrefix.".cont.deam $referenceFasta $inbam";

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


open(FILEOUT,">".$outputPrefix.".cont.est") or die "cannot write to ".$outputPrefix.".cont.est";
print FILEOUT $arrayOfValues[$maxLi]->{'cont'}."\t".$arrayOfValues[$il]->{'cont'}."\t".$arrayOfValues[$ih]->{'cont'}."\n";
close(FILEOUT);

my $cmdPlot = $contDeamR." ".$outputPrefix.".cont.deam ".$outputPrefix.".cont.pdf  \"$textGraph\" ";
if($contPriorKnow != -1){
  $cmdPlot =  $cmdPlot." ".$contPriorKnow;
}

runcmd($cmdPlot);

print "Program finished succesfully\n".
  "The plot of the posterior probability is ".$outputPrefix.".cont.pdf\n".
  "The contamination estimate is here ".$outputPrefix.".cont.est\n";

