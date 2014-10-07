#!/usr/bin/perl


use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';

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


my $bam2prof = $pathdir."/bam2prof";
my $contDeam  = $pathdir."/contDeam";

#
#open(FILE,$ARGV[0]) or die "cannot open ".$ARGV[0];
#while(my $line = <FILE>){
#  chomp($line);
#  print $line;
#}
#close(FILE);
#
# protocol s or d
# output prefix

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "This script is a wrapper that allows for the\nestimation of endogenous deamination and\nsubsequent contamination estimate.\n\n\tusage: ".$0." [--library (single|double)] [--out (output prefix)] [--ref (reference genome)] [--help|-?] input.bam\n\n";
  exit;
}

my $help;
my $library        = "none";
my $outputPrefix   = "outputdeam";
my $inbam          = "none";
my $referenceFasta = "none";


usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'library=s' => \$library,'ref=s' => \$referenceFasta, 'out=s' => \$outputPrefix, )
          or defined $help );

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


print "Wrapper finished succesfully output available as ".$outputPrefix.".cont.deam\n";
