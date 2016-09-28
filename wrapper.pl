#!/usr/bin/perl

#use bignum;
use strict;
use warnings;
use Getopt::Long;
use Cwd 'abs_path';
use List::Util qw(min max);
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
use File::Which;                  # exports which()
use File::Which qw(which where);  # exports which() and where()

my $mock =0;

sub runcmd{
  my ($cmdtorun) = @_;

  print "running cmd ". $cmdtorun."\n";
  #if($mock != 1){
  my @argstorun = ( "bash", ,"-c", $cmdtorun );

  if(system(@argstorun) != 0){
    die "system  cmd $cmdtorun failed: $?"
  }else{
    #print "\nok";
  }
  #}
  print "\ndone";
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

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\nThis script calls shrimp and bam-rewrap for a given mitochondrial reference.\n\nusage:\t".$0." <options> [output prefix] [reference fasta] [input in bam/in fastq.gz]\n\n".
"Options:\n".
"\nYou can also define the following executables manually:\n".
  "\n\t--gmapper    [path to gmapper]\t\t\tPath to gmapper from shrimp".
"\n\t--samtools   [path to samtools]\t\t\tPath to samtools from hstlib".
"\n\t--bam-rewrap [path to bam-rewrap]\t\tPath to bam-rewrap from biohazard".
"\n\t--threads    [num threads]\t\t\tNumber of threads to use (Default: 1)".


#"\t--help|-?".
"\n\n";
  exit;
}
my $help;


my $gmapper   = "gmapper";
my $samtools  = "samtools";
my $bamrewrap = "bam-rewrap";

$gmapper    = which($gmapper);
$samtools   = which($samtools);
$bamrewrap  = which($bamrewrap);
#
#if(defined $gmapper){
#  print "gmapper used ".$gmapper."\n";
#}else{
#  die "Please install gmapper from shrimp2 http://compbio.cs.toronto.edu/shrimp/";
#}
#
#if(defined $samtools){
#  print "samtools used ".$samtools."\n";
#}else{
#  die "Please install samtools from http://www.htslib.org/";
#}
#
#if(defined $bamrewrap){
#  print "bam-rewrap used ".$bamrewrap."\n";
#}else{
#  die "Please install bam-rewrap from https://github.com/udo-stenzel/biohazard";
#}
#
#
#
#die;

my $gmapperCMDL   = "";
my $samtoolsCMDL  = "";
my $bamrewrapCMDL = "";
my $threads=1;

usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'gmapper=s' => \$gmapperCMDL, 'samtools=s' => \$samtoolsCMDL, 'bam-rewrap=s' => \$bamrewrapCMDL,'threads=i' => \$threads )
          or defined $help );

if($gmapperCMDL eq ""){
  if(defined $gmapper){
    print "gmapper used ".$gmapper."\n";
  }else{
    die "Please install gmapper from shrimp2 http://compbio.cs.toronto.edu/shrimp/";
  }

}else{
  if(-e $gmapperCMDL){
    $gmapper =  $gmapperCMDL;
  }else{
    die "Please specify an existing path for gmapper";
  }
}


if($samtoolsCMDL eq ""){
  if(defined $samtools){
    print "samtools used ".$samtools."\n";
  }else{
#  die "Please install samtools from http://www.htslib.org/";

  }

}else{
  if(-e $samtoolsCMDL){
    $samtools =  $samtoolsCMDL;
  }else{
    die "Please specify an existing path for samtools";
  }
}

if($bamrewrapCMDL eq ""){
  if(defined $bamrewrap){
    print "bamrewrap used ".$bamrewrap."\n";
  }else{
    die "Please install bam-rewrap from https://github.com/udo-stenzel/biohazard";
  }

}else{
  if(-e $bamrewrapCMDL){
    $bamrewrap =  $bamrewrapCMDL;
  }else{
    die "Please specify an existing path for bamrewrap";
  }
}

my $outPrefix      = $ARGV[ $#ARGV-2  ];
my $referenceFasta = $ARGV[ $#ARGV-1  ];
my $inBAM          = $ARGV[ $#ARGV-0 ];


my $infq=0;

if($inBAM =~ /fastq.gz$/){
    $infq=1;
}



my $id;
my $seq;
my $numIDFound=0;
open(FILE,$referenceFasta) or die "cannot open ".$referenceFasta;
while(my $line = <FILE>){
  chomp($line);
  if($line =~ /^>/){
    $id    = $line;
    $numIDFound++;
  }else{
    $seq   .= $line;
  }
}
close(FILE);


if($numIDFound != 1){
  die "We require a unique reference, found $numIDFound\n";
}
my $refLength=length($seq);
if($refLength<1000){
  die "The mitochondrial reference cannot be less than 1000bp\n";
}



#print "#".$seq."#".substr($seq,0,1000)."#\n";
my $referenceFastaWrapped = $outPrefix.".wrapped.fa";
warn "Writing to ".$referenceFastaWrapped."\n";

open(FILEOUT,">".$referenceFastaWrapped) or die "cannot open ".$referenceFastaWrapped;
print FILEOUT ">mtref\n";
print FILEOUT $seq.substr($seq,0,1000)."\n";
close(FILEOUT);

my $cmd1 = $samtools." faidx $referenceFastaWrapped";
runcmd($cmd1);

#
#my $cmd2 = $gmapper." -N $threads -o 1 --single-best-mapping --sam-unaligned --fastq --sam --no-qv-check  --qv-offset 33  <( $samtools bam2fq $inBAM   )  $referenceFastaWrapped | $samtools view -bS -   > ".$outPrefix.".1.bam";
#runcmd($cmd2);
#my $cmd3 = $samtools." fillmd -b  /dev/stdin  $referenceFastaWrapped  > ".$outPrefix.".2.bam";
#runcmd($cmd3);
#my $cmd4 = $samtools." sort -O bam  -T ".$outPrefix.".sort1 ".$outPrefix.".2.bam  ".$outPrefix.".sort1 ".$outPrefix.".3.bam";
#runcmd($cmd4);
#my $cmd5 = $bamrewrap." \"mtref:".$refLength."\"  ".$outPrefix.".3.bam > ".$outPrefix.".4.bam";
#runcmd($cmd5);
#my $cmd6 = $samtools." sort -O bam  -T ".$outPrefix.".sort2  ".$outPrefix.".4.bam   > ".$outPrefix.".bam";
#runcmd($cmd6);


my $cmd2 = $gmapper." -N $threads -o 1 --single-best-mapping --sam-unaligned --fastq --sam --no-qv-check  --qv-offset 33  ";


if($infq){
    $cmd2=$cmd2."  $inBAM     ";
}else{
    $cmd2=$cmd2." <( $samtools bam2fq $inBAM   )  ";
}

$cmd2=$cmd2."$referenceFastaWrapped | $samtools view -bS  -F4 /dev/stdin | ".$samtools." fillmd -b  /dev/stdin  $referenceFastaWrapped   | ".$samtools." sort -O bam  -T ".$outPrefix.".sort1 /dev/stdin | ".$bamrewrap." \"mtref:".$refLength."\"   | ".$samtools." sort -O bam  -T ".$outPrefix.".sort2  /dev/stdin    > ".$outPrefix.".bam";
runcmd($cmd2);


my $cmd3 = $samtools." index ".$outPrefix.".bam";
runcmd($cmd3);


# $bamrewrap "mtref:16569"
#runcmd($cmd2);
#runcmd($cmd3);
#my $cmd4 = $samtools." sort -O bam  -T ".$outPrefix.".sort1 ".$outPrefix.".2.bam  ".$outPrefix.".sort1 ".$outPrefix.".3.bam";
#runcmd($cmd4);
#my $cmd5 = $bamrewrap." \"mtref:".$refLength."\"  ".$outPrefix.".3.bam > ".$outPrefix.".4.bam";
#runcmd($cmd5);
#my $cmd6 = $samtools." sort -O bam  -T ".$outPrefix.".sort2  ".$outPrefix.".4.bam   > ".$outPrefix.".bam";
#runcmd($cmd6);

#print $cmd2."\n";


#die $contPriorKnow;
#print $data."\n";
#run deam

#if($referenceFasta eq "none" ){
#  die "Please enter the fasta reference used for mapping\n";
#}


print "Program finished succesfully\n\nSequences written to ".$outPrefix.".bam\n";
#  "The plot of the posterior probability is ".$outputPrefix.".cont.pdf\n".
#  "The contamination estimate is here ".$outputPrefix.".cont.est\n".
#  "The configuration file is here ".$outputPrefix.".config\n";

