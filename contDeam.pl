#!/usr/bin/perl

#use bignum;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);

sub runcmd{
  my ($cmdtorun) = @_;

  print "running cmd ". $cmdtorun."\n";
  if(system($cmdtorun) != 0){
    die "system  cmd $cmdtorun failed: $?"
  }else{

  }

}

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);


my $bam2prof   = $pathdir."/bam2prof";
my $contDeam   = $pathdir."/contDeam";
my $contDeamR  = $pathdir."/posteriorDeam.R";

#
# protocol s or d
# output prefix

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "This script is a wrapper that allows for the\nestimation of endogenous deamination and\nsubsequent contamination estimate.\n\n\tusage: ".$0." <options> input.bam\n\nOptions:\n\n\t--library (single|double)\n\t--out (output prefix)\n\t--cont (cont)\tIf you have prior knowledge about the contamination rate, enter it here [0-1]\n\t--ref (reference genome)\n\t--title (title)\tTitle for the graph\n\t--help|-?\n\n";
  exit;
}

my $help;
my $library        = "none";
my $outputPrefix   = "outputdeam";
my $inbam          = "none";
my $referenceFasta = "none";
my $contPriorKnow  = -1;
my $textGraph      = "Posterior probability for contamination\nusing deamination patterns";


usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'library=s' => \$library,'ref=s' => \$referenceFasta,'title=s' => \$textGraph,'cont=f' => \$contPriorKnow, 'out=s' => \$outputPrefix, )
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

my $cmdBam2Prof = $bam2prof." -endo -".$library." -5p ".$outputPrefix.".5p.prof  -3p ".$outputPrefix.".3p.prof $inbam";

runcmd($cmdBam2Prof);

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

