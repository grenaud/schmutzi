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

sub writeContToLogFile{
  my ($cont,$contOutFile) = @_;

  print "Writing contamination estimate $cont to file $contOutFile\n";
  if ($mock != 1) {

    open(FILEoutcont, ">".$contOutFile) or die "cannot write to  ".$contOutFile."\n";
    print FILEoutcont $cont."\t".$cont."\t".$cont."\n";
    close(FILEoutcont);
  }
}


sub computeAmountOfPositions{
  my ($posfile) = @_;

  print "reading position file ".$posfile."\n";

  if ($mock != 1) {
    open(FILE,$posfile) or die "cannot open ".$posfile."";
    my $numberLines=0;
    while (my $line = <FILE>) {
      chomp($line);
      $numberLines++;
    }
    close(FILE);

    return $numberLines;
  }else{
    #mock, nothing to do
    return -1;
  }

}



sub computeIntervalCont{
  my ($mtcontOutputLog) = @_;

  print "reading the contamination log ".$mtcontOutputLog."\n";

  if ($mock != 1) {
    open(FILE,$mtcontOutputLog) or die "cannot open ".$mtcontOutputLog."";
    my $maxL=10;
    my $maxLi=0;

    my @arrayOfValues;
    while (my $line = <FILE>) {
      chomp($line);
      #print $line;
      my @array = split("\t",$line);
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


sub readMTCONToutlogfile{
  my ($mtcontOutputLog) = @_;

  print "reading the log file ".$mtcontOutputLog."\n";

  my $sourceFile="###";
  my $contMaxS    = -1;
  my $llikMaxS    = -1;
  my $contMaxSALL = -1;
  my $llikMaxSALL = -1;
  my $contFileInd = 0;

  if ($mock != 1) {

    open(FILEoutlogmtcont, $mtcontOutputLog) or die "cannot open ".$mtcontOutputLog."\n";

    while(my $line = <FILEoutlogmtcont>){
      my @arrayTemp = split("\t",$line);
      if($#arrayTemp != 2){
	die "Error line $line in $mtcontOutputLog does not have 3 fields";
      }

      if($arrayTemp[0] ne $sourceFile){ # new
	if($sourceFile eq "###"){ #first one
	  $sourceFile = $arrayTemp[0];
	}else{#new source
	  if( $contFileInd == 1 ){#first file
	    $contMaxSALL = $contMaxS;
	    $llikMaxSALL = $llikMaxS;
	  }else{

	    if($llikMaxSALL < $llikMaxS){
	      $contMaxSALL = $contMaxS;
	      $llikMaxSALL = $llikMaxS;
	    }

	  }
	}

	$contMaxS = $arrayTemp[1];
	$llikMaxS = $arrayTemp[2];
	$contFileInd++;

      }else{

	if($llikMaxS < $arrayTemp[2]){
	  $contMaxS = $arrayTemp[1];
	  $llikMaxS = $arrayTemp[2];
	}

      }


    }

    if($llikMaxSALL < $llikMaxS){
      $contMaxSALL = $contMaxS;
      $llikMaxSALL = $llikMaxS;
    }
    close(FILEoutlogmtcont);

    if ( $contMaxSALL == -1) {
      die "Cannot parse the contamination output file = ".$mtcontOutputLog."\n";
    }

    return $contMaxSALL;

  } else {
    return 0.5;
  }

}

sub readMTCONTRetainMostlikely{
  my ($mtcontOutputLog,$outputToWrite) = @_;

  print "reading the log file ".$mtcontOutputLog." and print most likely source to ".$outputToWrite."\n";


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


    open(FILEoutlogmtcont,      $mtcontOutputLog) or die "cannot open ".$mtcontOutputLog."\n";
    open(FILEoutlogmtcontOUT, ">".$outputToWrite) or die "cannot write to ".$outputToWrite."\n";

    while(my $line = <FILEoutlogmtcont>){
      my @arrayTemp = split("\t",$line);
      if($#arrayTemp != 2){
	die "Error line $line in $mtcontOutputLog does not have 3 fields";
      }

      if($sourceALL eq $arrayTemp[0]){
	print FILEoutlogmtcontOUT $line;
      }

    }

    close(FILEoutlogmtcont);
    close(FILEoutlogmtcontOUT);



  } else {
    #nothing to do
  }

}



sub copycmd{
  my ($source,$destination) = @_;

  print "copying  ". $source." to ".$destination."\n";

  if($mock != 1){
    copy( $source,$destination ) or die "Copy file failed: $!";
  }

}


sub makeEmptyProfFile{
  my ($filename) = @_;


  if($mock != 1){
    open(FILEprof, ">".$filename) or die "cannot write to ".$filename;

    print FILEprof "A>C\tA>G\tA>T\tC>A\tC>G\tC>T\tG>A\tG>C\tG>T\tT>A\tT>C\tT>G\n";
    for(my $i=0;$i<$lengthDeam;$i++){
      print FILEprof "0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\t0.0\n";
    }
    close(FILEprof);

  }else{
    print "writing empty deamination profile to $filename\n";
  }

}


sub readcont{
  my ($filename) = @_;
  #
  print "reading the contamination estimate file: $filename\n";
  #if ($mock != 1) {
  open(FILEcont, $filename) or die "cannot open ".$filename;

  my $line = <FILEcont>;
  my @arrayTemp = split("\t",$line);
  close(FILEcont);

  return $arrayTemp[0];
  #} else {
    #print "reading the contamination estimate file: $filename\n";
   # return 0.5;
  #}
}


sub readSizeParam{
  my ($filename) = @_;
  print "reading the insert size parameter file: $filename\n";
  if ($mock != 1) {
    open(FILEparamEndo,$filename) or die "cannot open ".$filename."\n";
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

    my $loc=$arrayEndo[0];
    my $sca=$arrayEndo[1];

    return ($loc,$sca);
  } else {
    return (1,1);
  }
}



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



sub runcmdReturnOutput {
  my $command = join ' ', @_;
  reverse ($_ = qx{$command 2>&1}, $? >> 8);
}


sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $$exeFile)) {
    if ( ( -e (($$exeFile).".exe")) ) {
      $$exeFile=(($$exeFile).".exe");
    }else{
      die "Executable ".$$exeFile." does not exist\n";
    }
  }

}

sub doesFileExists{
  my ($exeFile) = @_;

  if (!( -e $$exeFile)) {
    if ( ( -e (($$exeFile).".exe")) ) {
      return 1;
    }else{
      return 0;
    }
  }else{
    return   1;
  }

}

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $mtCont     = $pathdir."/mtCont";
my $endoCaller = $pathdir."/endoCaller";
my $insertSize = $pathdir."/insertSize";
my $countRec   = $pathdir."/countRecords";

my $approxDist = $pathdir."/approxDist.R";
my $bam2prof   = $pathdir."/bam2prof";
my $log2freq   = $pathdir."/log2freq";
my $logs2pos   = $pathdir."/logs2pos";

my $contDeam   = $pathdir."/contDeam";
my $contDeamR  = $pathdir."/posteriorDeam.R";
my $splitEndo  = $pathdir."/splitEndoVsCont/poshap2splitbam";

fileExists(\$mtCont);
fileExists(\$endoCaller);
fileExists(\$insertSize);
fileExists(\$countRec);
fileExists(\$approxDist);
fileExists(\$bam2prof);
fileExists(\$log2freq);
fileExists(\$logs2pos);
fileExists(\$contDeam);
fileExists(\$contDeamR);

if(doesFileExists(\$splitEndo)){
  $splitEndo  = $pathdir."/poshap2splitbam";
}

fileExists(\$splitEndo);


my $nameMT    = "MT";
my $nameMTc   = "MTc";
my $qualmin   = 0;
my $lengthMT=16569;
my $logindel=50;

#
# protocol s or d
# output prefix


my $multipleC=0;
#my $usepredC=0;
my $notusepredC=0;

my $numthreads=1;
my $maxIterations      = 100;
my $iterationSameCont  = 3;

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
This script is a wrapper to produce the endogenous mitochondrial
consensus genome and a contamination estimate given aligned
mitochondrial data and a set of putative contaminants. This script
takes the ouput of the contDeam.pl script and calls iteratively
the module that produces an endogenous consensus and the one that
estimates contamination. For non-hominin samples where the is no
risk of human contamination where the mitochondrial consensus is
needed, call endoCaller directly (see README.md).

\n\nusage:\t".$0." <options> [output prefix from contDeam.pl] [directory with *freq file] input.bam\n\n".

#"Mandatory parameters:\n".
#"\t[output prefix from contDeam]\tThis is the prefix of the output of contDeam.pl (option --out)\n".
#"\t[directory with *freq file]\tDirectory with all the *freq file containing the putative contaminants\n".
#"\t[reference *fasta]\tFasta file of the reference used for by aligner\n".
#"\t[input.bam]\tAligned BAM file\n\n".
#

"Options:\n".


"\n\t--iterations (# it.)\t\t\tMaximum number of iterations (Default : $maxIterations)\n".

"\t--t (threads #)\t\t\t\tUse this amount of threads (Default : $numthreads)\n".
"\t--lengthMT (bp)\t\t\t\tConsider that the real length of the MT genome is this (Default : $lengthMT)\n".
"\t--estdeam\t\t\t\tRe-estimate deamination parameters using newly computed segregating positions".
"\n\t\t\t\t\t\tThis setting is recommended if the contamination is deaminated".
"\n".

"\t--mock\t\t\t\t\tDo nothing, just print the commands used\n".
"\t--uselength\t\t\t\tUse length of the molecules as well\n".
"\t--lengthDeam (bp)\t\t\tOnly consider this about of bases to be deaminated on each end (Default : $lengthDeam)\n".
"\t--multipleC\t\t\t\tDo not assume that there is a single contaminant\n".
"\t--notusepredC\t\t\t\tIf assuming a single contaminant, do not use the predicted contaminant in the contamination estimate".

"\n\t\t\t\t\t\tThis might lead to worse results".
"\n\n".
#
"  Output Options:\n".
"\t--out (out prefix)\t\t\tOutput prefix, otherwise use the same used by contDeam.pl\n".
"\t--name (name)\t\t\t\tName of the endogenous MT genome\n".
"\t--namec (name)\t\t\t\tName of the contaminant MT genome\n".
"\t--qual (qual)\t\t\t\tMinimum quality for a base in the  MT consensus (on a PHRED scale, e.g. 50 = 1/100,000)\n".
"\t--logindel  (logindel)\t\t\tMinimum difference between logs to call an indel (on a PHRED scale default : ".$logindel.")\n".

"\t--contknown (cont)\t\t\tIf you have prior knowledge about the contamination\n\t\t\t\t\t\trate, enter it here [0-1] and it will be plotted\n".
"\t--title (title)\t\t\t\tTitle for the graph of the posterior distribution\n".
"\n".
"  Input Options:\n".
"\t--contprior (rate)\t\t\tIgnore the contamination rate prior found by contDeam.pl, use this rate instead\n".
"\t--ref (reference genome)\t\tThe fasta file used for alignment if not specified when calling contDeam.pl\n".

#"\t--out (output prefix)\t\tAll output files will share this prefix\n".

#"\t--title (title)\t\t\tTitle for the graph of the posterior distribution\n".
#"\t--cont (cont)\t\t\tIf you have prior knowledge about the contamination\n\t\t\t\t\trate, enter it here [0-1]\n".
#"\nInput options:\n".
#"\t--split (file)\t\t\tSplit endogenous/contaminant according to diagnostic positions\n".
#"\t\t\t\t\tThe file must have the following format:\n".
#"\t\t\t\t\t\t[coord]tab[nucleotide]tab[endo or cont]\n".
#"\t\t\t\t\tWhere the coordinate is on the reference genome\n".
#"\t\t\t\t\tex:\t385\tA\tendo\n".
#"\nMandatory:\n".
#"\t--ref (reference genome)\tThe fasta file used for alignment\n".

#"\t--help|-?".
"\n\n";
  exit;
}

my $starttime = time;


my $help;
my $library        = "none";
my $outputPrefix   = "outputdeam";


my $contPriorKnow         = -1;
my $contPriorKnowCMDLine  = -1;

my $textGraph      = "Posterior probability for contamination\\nusing mitochondrial positions";
my $referenceFasta = "none";
my $referenceFastaCMDL = "none";
my $splitPos       = "";
my $useLength      = 0;
my $useLengthContDEAM      = 0; #from contdeam
my $outputPrefixCMDLINE = "none";
my $estdeam = 0 ;
#my $contPriorSpecified = 0 ;
my $contPriorUser      = -1 ;

usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 'iterations=i' => \$maxIterations,'ref=s' => \$referenceFastaCMDL, 't=i' => \$numthreads, 'mock' => \$mock, 'estdeam' => \$estdeam, 'uselength' => \$useLength, 'title=s' => \$textGraph, 'out=s' => \$outputPrefixCMDLINE,'name=s' => \$nameMT, 'namec=s' => \$nameMTc, 'contknown=f' => \$contPriorKnowCMDLine, 'lengthDeam' => \$lengthDeam,'lengthMT' => \$lengthMT,'multipleC' => \$multipleC,'notusepredC' => \$notusepredC,'contprior=f' => \$contPriorUser,'qual=f' => \$qualmin,'logindel=f' => \$logindel ) #,'name=s' => \$nameMT,'namec=s' => \$nameMTc
          or defined $help );



my $prefixcontDeam = $ARGV[ $#ARGV -2 ];
my $configfiledeam = $prefixcontDeam.".deam.config";
my $freqDir        = $ARGV[ $#ARGV -1 ];
my $inbam          = $ARGV[ $#ARGV -0 ];
my $splitDeam =-1;

open(FILEcontdeam, $configfiledeam) or die "cannot open ".$configfiledeam;

while (my $line = <FILEcontdeam>) {

  if($line =~ /^library\s(\S+)$/){        $library=$1;}
  if($line =~ /^outputPrefix\s(\S+)$/){   $outputPrefix=$1;}
  #if($line =~ /^inbam\s(\S+)$/){ $inbam=$1;}
  if($line =~ /^referenceFasta\s(\S+)$/){ $referenceFasta=$1;}
  if($line =~ /^contPriorKnow\s(\S+)$/){  $contPriorKnow=$1;}
  #if($line =~ /^textGraph\s(\S+)$/){      $textGraph=$1;}
  if($line =~ /^splitPos\s(\S+)$/){       $splitPos=$1;} #file containing positions, empty if not used
  if($line =~ /^useLength\s(\S+)$/){      $useLengthContDEAM=$1;}
  if($line =~ /^splitDeam\s(\S+)$/){      $splitDeam=$1;}
}

close(FILEcontdeam);

if($outputPrefixCMDLINE  ne "none"){
  $outputPrefix=$outputPrefixCMDLINE;
}

if($referenceFastaCMDL  ne "none"){
  if($referenceFasta    ne "none" &&
    $referenceFastaCMDL ne $referenceFasta){
    print "WARNING: you are specifying the reference via the command line despite the fact that there is already one specified.\n";
  }
  $referenceFasta = $referenceFastaCMDL;
}

if($contPriorKnowCMDLine != -1){
  $contPriorKnow = $contPriorKnowCMDLine;
}

if($contPriorUser != -1){
  if($contPriorUser<0 || $contPriorUser>1){
    die "ERROR: The contamination prior specified via --contprior should be between 0 and 1.\n";
  }
}

if($numthreads <= 0){
  die "The number of threads must be greater than 1\n";
}


if($splitDeam == -1){
  die "The field splitDeam should be defined in $configfiledeam\n";
}



if($library ne "single" &&
   $library ne "double" ){
  die "Please enter the type of library as either double or single\n";
}

if($referenceFasta eq "none" ){
  die "Please enter the fasta reference used for mapping\n";
}

#read all files in $freqDir
if($freqDir !~ /\/$/){
  $freqDir .= "/";
}

my @arrayOfFreqFiles;
opendir(FREQDIR, $freqDir) || die "Cannot open directory: $freqDir\n";
while (my $freqf = readdir(FREQDIR)) {
  #print "$freqf\n";
    if($freqf =~ /\.freq$/){
      push(@arrayOfFreqFiles, $freqDir.$freqf);
    }
}
closedir(FREQDIR);
print "Found the following freq files:\n------------------\n";
foreach my $freqf (@arrayOfFreqFiles){
  print $freqf."\n";
}
print "------------------\n";
#die;
# begin loop
my $numberIteration=1;

########################
#DEAMINATION PARAMETERS#
########################
#if(! $estdeam  ){ #do not re-estimate at each iteration


copycmd( $prefixcontDeam.".endo.5p.prof" , $outputPrefix."_".$numberIteration."_endo.5p.prof" );
copycmd( $prefixcontDeam.".endo.3p.prof" , $outputPrefix."_".$numberIteration."_endo.3p.prof" );

if($splitDeam){
  copycmd( $prefixcontDeam.".cont.5p.prof" , $outputPrefix."_".$numberIteration."_cont.5p.prof" );
  copycmd( $prefixcontDeam.".cont.3p.prof" , $outputPrefix."_".$numberIteration."_cont.3p.prof" );
}else{#put dummy values for cont, the split was not done by contDeam and we will do it after the first iteration
  makeEmptyProfFile($outputPrefix."_".$numberIteration."_cont.5p.prof" );
  makeEmptyProfFile($outputPrefix."_".$numberIteration."_cont.3p.prof" );
}


#}else{
#  #we re-estimate at each iteration the deamination
  
#}

copycmd( $prefixcontDeam.".cont.est" , $outputPrefix."_".$numberIteration."_cont.est" );

########################
#  LENGTH PARAMETERS   #
########################
if($useLength){#we use the length of the molecules in the endoCaller

  if($useLengthContDEAM){ #if previously computed
    copycmd( $prefixcontDeam."_endo.size.param" , $outputPrefix."_".($numberIteration)."_split_endo.size.param" );
    copycmd( $prefixcontDeam."_cont.size.param" , $outputPrefix."_".($numberIteration)."_split_cont.size.param" );
  }else{
    #otherwise, the length will not be used for the first iteration
  }
}



print "\n\n\n############################\n\n\n";

my $previousCurrentContEst=-1;
my $previousCurrentContEstItSameVal=0;
my $notEnoughDataSize=0;


while(1){
  my $currentContEst=readcont( $outputPrefix."_".$numberIteration."_cont.est" );

  if($numberIteration == 1 ){ #first iteration
    if($contPriorUser != -1){ #user specified a prior contamination rate
      $currentContEst=$contPriorUser;
    }
  }

  my $stringToPrint1= "############################";
  my $stringToPrint2="#    ITERATION #".sprintf("%".(log10($maxIterations)+1)."s",$numberIteration)."";

  while((length($stringToPrint2)+1)<length($stringToPrint1)){
    $stringToPrint2.=" ";
  }
  $stringToPrint2.="#";
  print "\n\n\n".$stringToPrint1."\n";
  print $stringToPrint2."\n";

        #\n";
  print "#                          #\n";
  print "# Contamination rate: ".sprintf("%.2f",$currentContEst)." #\n";


  print "############################\n";


  # call endoCaller
  my $cmdEndoCaller = $endoCaller." ";

  $cmdEndoCaller .= " -seq ".$outputPrefix."_".$numberIteration."_endo.fa  ";
  $cmdEndoCaller .= " -log ".$outputPrefix."_".$numberIteration."_endo.log ";

  $cmdEndoCaller .= " -name ".$nameMT." ";
  $cmdEndoCaller .= " -qual ".$qualmin." ";
  $cmdEndoCaller .= " -logindel ".$logindel." ";
  $cmdEndoCaller .= " -deamread ";
  $cmdEndoCaller .= " -deam5p ".$outputPrefix."_".$numberIteration."_endo.5p.prof ";
  $cmdEndoCaller .= " -deam3p ".$outputPrefix."_".$numberIteration."_endo.3p.prof ";




  $cmdEndoCaller .= " -cont ".min(max($currentContEst,0.01),0.99)." ";#to allow for the possibility of having contamination

  if(!$multipleC){ #we can assume a single contaminant
  $cmdEndoCaller .= " -deam5pc ".$outputPrefix."_".$numberIteration."_cont.5p.prof ";
  $cmdEndoCaller .= " -deam3pc ".$outputPrefix."_".$numberIteration."_cont.3p.prof ";

    $cmdEndoCaller .= " -single  ";

    $cmdEndoCaller .= " -seqc ".$outputPrefix."_".$numberIteration."_cont.fa  ";
    $cmdEndoCaller .= " -logc ".$outputPrefix."_".$numberIteration."_cont.log ";
    $cmdEndoCaller .= " -namec ".$nameMTc." ";

  }

  $cmdEndoCaller .= " -l ".$lengthMT." ";

  # mol. length
  if ($useLength) {		#
    my $notDefined=0;
    my $locE=1;
    my $scaE=1;
    my $locC=1;
    my $scaC=1;

    if($numberIteration == 1){ #first iteration
      if($useLengthContDEAM){ #if previously computed
	($locE,$scaE)=   readSizeParam($outputPrefix."_".($numberIteration)."_split_endo.size.param" );
	($locC,$scaC)=   readSizeParam($outputPrefix."_".($numberIteration)."_split_cont.size.param" );

      }else{
	#do not use length for the first iteration
	$notDefined=1;
      }

    }else{

      if($notEnoughDataSize){#there was not enough data in the previous iteration
	$notDefined=1;
      }else{
	($locE,$scaE)=   readSizeParam($outputPrefix."_".($numberIteration)."_split_endo.size.param" );
	($locC,$scaC)=   readSizeParam($outputPrefix."_".($numberIteration)."_split_cont.size.param" );
      }
    }

    if(!$notDefined){
        $cmdEndoCaller .= " --loce   ".$locE." ";
        $cmdEndoCaller .= " --scalee ".$scaE." ";

        $cmdEndoCaller .= " --locc   ".$locC." ";
        $cmdEndoCaller .= " --scalec ".$scaC." ";
    }

  }



  $cmdEndoCaller .= " $referenceFasta $inbam ";
  runcmd($cmdEndoCaller);

  my @listOfFreqFiles=@arrayOfFreqFiles;

  #convert log to freq
  if(!$multipleC){ #we can assume a single contaminant

    if ( !$notusepredC  ){
      my $cmdLog2Freq =  $log2freq." ".$outputPrefix."_".$numberIteration."_cont.log >  ".$outputPrefix."_".$numberIteration."_cont.freq";
      runcmd($cmdLog2Freq);
      push(@listOfFreqFiles, $outputPrefix."_".$numberIteration."_cont.freq");
    }
  }

  # estimate cont
  my $cmdmtcont =   $mtCont. " -o ".$outputPrefix."_".$numberIteration."_mtcont.out ";
  #$cmdmtcont =   " -deam5p "
  $cmdmtcont .= " -deam5p ".$outputPrefix."_".$numberIteration."_endo.5p.prof ";
  $cmdmtcont .= " -deam3p ".$outputPrefix."_".$numberIteration."_endo.3p.prof ";
  $cmdmtcont .= " -t ".$numthreads." ";
  $cmdmtcont .= " ".$outputPrefix."_".$numberIteration."_endo.log ";
  $cmdmtcont .= " $referenceFasta $inbam ";
  $cmdmtcont .= " ".join(" ",@listOfFreqFiles). " ";
  runcmd($cmdmtcont);



  # if likelihood stable, break loop
  if($numberIteration >= $maxIterations ){
    print "Reached the maximum number of iterations $numberIteration, exiting\n";
    last;
  }



  #BEGIN read new cont esti. log
  my $mtcontOutputLog = $outputPrefix."_".$numberIteration."_mtcont.out";

  my $contFROMmtcont =  readMTCONToutlogfile($mtcontOutputLog);
  writeContToLogFile($contFROMmtcont,$outputPrefix."_".($numberIteration+1)."_cont.est");

  #END   read new cont esti. log


  if (!$multipleC) {		#we can assume a single contaminant


    if (  ($estdeam || $useLength)  ) { # we need to split if and only if we re-estimate the deamination at each iteration or use the length of the molecules

      #print seg. sites
      my $cmdLogs2Pos = $logs2pos."  ".$outputPrefix."_".$numberIteration."_endo.log ".$outputPrefix."_".$numberIteration."_cont.log  > ".$outputPrefix."_".$numberIteration."_split.pos";
      runcmd($cmdLogs2Pos);

      #count pos in .pos
      my $countPosFound=computeAmountOfPositions($outputPrefix."_".$numberIteration."_split.pos");
      if ($countPosFound<5) {
	warn "Unable to find more than 5 positions that are different in  ".$outputPrefix."_".$numberIteration."_endo.log ".$outputPrefix."_".$numberIteration."_cont.log\n".
	  "Either your sample has little or no contamination or too much, either case, we cannot split according to segrating positions.\n".
	  "WARNING: We will keep running without using --useLength nor --estdeam\n";
	$estdeam   = 0;
	$useLength = 0;
	goto SKIPESTLENGTH;
      }

      # split according to seg sites
      my $cmdBamSplit = $splitEndo."  ".$outputPrefix."_".$numberIteration."_split.pos  $inbam ".$outputPrefix."_".($numberIteration)."_split > ".$outputPrefix."_split.log 2> /dev/null ";
      #runcmd($cmdBamSplit);
      my @argsBamSplit = ( "bash", "-c", $cmdBamSplit );

      if(system(@argsBamSplit) != 0){
	warn "system  cmd $cmdBamSplit failed with code: $?";

	warn "The following command: ".$cmdBamSplit." caused a segfault\n".
	  "Please run the program again with a higher stack/heap limit.\n".
	  "WARNING: We will keep running without using --useLength nor --estdeam\n";
	$estdeam   = 0;
	$useLength = 0;
	goto SKIPESTLENGTH;

	#if ($? == -1) {
	#  die "ERROR: failed to execute: $!\n";
	#}elsif ($? & 127) {
	#  warn "The following command: ".$cmdBamSplit." caused a segfault\n".
	#    "Please run the program again with a higher stack/heap limit.\n".
	#    "WARNING: We will keep running without using --useLength nor --estdeam\n";
	#  $estdeam   = 0;
	#  $useLength = 0;
	#  goto SKIPESTLENGTH;
	#}else {
	#  die "ERROR: failed to execute: $!\n";
	#}

      }else{
      }

      #re-measure deam rateS
      if ( $estdeam  ) {   # re-estimate deamination at each iteration

	my $cmdBam2ProfEndo = $bam2prof." -length $lengthDeam -".$library." -5p ".$outputPrefix."_".($numberIteration+1)."_endo.5p.prof   -3p ".$outputPrefix."_".($numberIteration+1)."_endo.3p.prof  ".$outputPrefix."_".($numberIteration)."_split_endo.bam";
	runcmd($cmdBam2ProfEndo);

	my $cmdBam2ProfCont = $bam2prof." -length $lengthDeam -".$library." -5p ".$outputPrefix."_".($numberIteration+1)."_cont.5p.prof   -3p ".$outputPrefix."_".($numberIteration+1)."_cont.3p.prof  ".$outputPrefix."_".($numberIteration)."_split_cont.bam";
	runcmd($cmdBam2ProfCont);

	#my $cmdBam2ProfCont = $bam2prof." -length $lengthDeam -".$library." -5p ".$outputPrefix.".cont.5p.prof  -3p ".$outputPrefix.".cont.3p.prof ".$outputPrefix."_cont.bam";
	#runcmd($cmdBam2ProfCont);
      } else {

	copycmd(  $outputPrefix."_".$numberIteration."_endo.5p.prof" ,$outputPrefix."_".($numberIteration+1)."_endo.5p.prof" );
	copycmd(  $outputPrefix."_".$numberIteration."_endo.3p.prof" ,$outputPrefix."_".($numberIteration+1)."_endo.3p.prof" );
	copycmd(  $outputPrefix."_".$numberIteration."_cont.5p.prof" ,$outputPrefix."_".($numberIteration+1)."_cont.5p.prof" );
	copycmd(  $outputPrefix."_".$numberIteration."_cont.3p.prof" ,$outputPrefix."_".($numberIteration+1)."_cont.3p.prof" );

      }


      #check if enough data to measure insert size
      my ($returnCE, $readCountE) = runcmdReturnOutput($countRec,$outputPrefix."_".($numberIteration)."_split_endo.bam");
      my ($returnCC, $readCountC) = runcmdReturnOutput($countRec,$outputPrefix."_".($numberIteration)."_split_cont.bam");

      if ($returnCE !=0) {
	die "Unable to determine how many reads are in ".$outputPrefix."_".($numberIteration)."_split_endo.bam";
      }

      if ($returnCC !=0) {
	die "Unable to determine how many reads are in ".$outputPrefix."_".($numberIteration)."_split_cont.bam";
      }

      if ($readCountE < 50) {
	warn "Not enough reads are in ".$outputPrefix."_".($numberIteration)."_split_endo.bam, we will not use insert sizes for this iteration";
      }

      if ($readCountC < 50) {
	warn "Not enough reads are in ".$outputPrefix."_".($numberIteration)."_split_cont.bam, we will not use insert sizes for this iteration";
      }

      if ($useLength) { #we use the length of the molecules in the endoCaller

	if ($readCountE >= 50 &&
	    $readCountC >= 50 ) {

	  #measure insert size
	  my $cmdBam2InsertsizeEndo = $insertSize."   ".$outputPrefix."_".($numberIteration)."_split_endo.bam |gzip   > ".$outputPrefix."_".($numberIteration)."_split_endo.size.gz";
	  runcmd($cmdBam2InsertsizeEndo);
	  my $cmdBam2InsertsizeCont = $insertSize."   ".$outputPrefix."_".($numberIteration)."_split_cont.bam  |gzip  > ".$outputPrefix."_".($numberIteration)."_split_cont.size.gz";
	  runcmd($cmdBam2InsertsizeCont);

	  #Log normal fit, copying parameters for next iteration
	  my $cmdinsertsize2LognormEndo = $approxDist."  ".$outputPrefix."_".($numberIteration)."_split_endo.size.gz >  ".$outputPrefix."_".($numberIteration+1)."_split_endo.size.param";
	  runcmd($cmdinsertsize2LognormEndo);
	  my $cmdinsertsize2LognormCont = $approxDist."  ".$outputPrefix."_".($numberIteration)."_split_cont.size.gz >  ".$outputPrefix."_".($numberIteration+1)."_split_cont.size.param";
	  runcmd($cmdinsertsize2LognormCont);
	  $notEnoughDataSize=0;
	} else {
	  $notEnoughDataSize=1;
	}
      }

    } else { #if we do not use the length or re-estimate deamintion, we will re-use the deamintion profile from contDeam.pl

    SKIPESTLENGTH:
      copycmd(  $outputPrefix."_".$numberIteration."_endo.5p.prof" ,$outputPrefix."_".($numberIteration+1)."_endo.5p.prof" );
      copycmd(  $outputPrefix."_".$numberIteration."_endo.3p.prof" ,$outputPrefix."_".($numberIteration+1)."_endo.3p.prof" );
      copycmd(  $outputPrefix."_".$numberIteration."_cont.5p.prof" ,$outputPrefix."_".($numberIteration+1)."_cont.5p.prof" );
      copycmd(  $outputPrefix."_".$numberIteration."_cont.3p.prof" ,$outputPrefix."_".($numberIteration+1)."_cont.3p.prof" );

    }
  } else { #must use the previous deamination profiles from the previous iteration for the new one

    copycmd(  $outputPrefix."_".$numberIteration."_endo.5p.prof" ,$outputPrefix."_".($numberIteration+1)."_endo.5p.prof" );
    copycmd(  $outputPrefix."_".$numberIteration."_endo.3p.prof" ,$outputPrefix."_".($numberIteration+1)."_endo.3p.prof" );
    copycmd(  $outputPrefix."_".$numberIteration."_cont.5p.prof" ,$outputPrefix."_".($numberIteration+1)."_cont.5p.prof" );
    copycmd(  $outputPrefix."_".$numberIteration."_cont.3p.prof" ,$outputPrefix."_".($numberIteration+1)."_cont.3p.prof" );


  }

  # measure deam params + length


  if($mock == 1){
    die;
  }


  if($previousCurrentContEst == -1){

    $previousCurrentContEst          = $currentContEst;
    $previousCurrentContEstItSameVal = 1;

  }else{


    my $diffBetweenItContRate = sprintf("%.2f",abs($previousCurrentContEst-$currentContEst));

    if($diffBetweenItContRate <= 0.01){
      $previousCurrentContEstItSameVal++;
    }else{
      $previousCurrentContEstItSameVal=1;
    }
    $previousCurrentContEst=$currentContEst;


    if( $previousCurrentContEstItSameVal >= $iterationSameCont ){
      print "Reached the maximum number of iterations ($iterationSameCont) with stable contamination rate at iteration # $numberIteration, exiting\n";
      last;
    }

  }

  $numberIteration++;
} # go to begin loop
print "############################\n";
print "Iterations done\n\n";


#copy files
copycmd($outputPrefix."_".$numberIteration."_mtcont.out",     $outputPrefix."_final_mtcont.out");
copycmd($outputPrefix."_".$numberIteration."_endo.fa",        $outputPrefix."_final_endo.fa");
copycmd($outputPrefix."_".$numberIteration."_endo.log",       $outputPrefix."_final_endo.log");


if(!$multipleC){ #we can assume a single contaminant
  copycmd($outputPrefix."_".$numberIteration."_cont.fa",      $outputPrefix."_final_cont.fa");
  copycmd($outputPrefix."_".$numberIteration."_cont.log",     $outputPrefix."_final_cont.log");
}
readMTCONTRetainMostlikely($outputPrefix."_final_mtcont.out", $outputPrefix."_final.cont");

#Print graph
my $cmdPlot = $contDeamR." ".$outputPrefix."_final.cont ".$outputPrefix."_final.cont.pdf  \"$textGraph\" ";
if($contPriorKnow != -1){
  $cmdPlot =  $cmdPlot." ".$contPriorKnow;
}
runcmd($cmdPlot);


#report final iteration

my ($contfinal,$contfinall,$contfinalh) =  computeIntervalCont($outputPrefix."_final.cont");

open(FILEOUTCF,">".$outputPrefix."_final.cont.est") or die "cannot write to ".$outputPrefix."_final.cont.est\n";
print FILEOUTCF $contfinal."\t".$contfinall."\t".$contfinalh."\n";
close(FILEOUTCF);


print "############################\n";
print "Results:\n";
print "\tContamination estimates for all samples        : ".    $outputPrefix."_final_mtcont.out\n";
print "\tContamination estimates for most likely sample : ".    $outputPrefix."_final.cont\n";
print "\tContamination estimates with conf. intervals   : ".    $outputPrefix."_final.cont.est\n";

print "\tPosterior probability for most likely sample   : ".    $outputPrefix."_final.cont.pdf\n";

print "\tEndogenous consensus call                      : ".    $outputPrefix."_final_endo.fa\n";
print "\tEndogenous consensus log                       : ".    $outputPrefix."_final_endo.log\n";

if(!$multipleC){ #we can assume a single contaminant
print "\tContaminant consensus call                     : ".    $outputPrefix."_final_cont.fa\n";
print "\tContaminant consensus log                      : ".    $outputPrefix."_final_cont.log\n";

}



if ($mock != 1) {

  open(FILECONFIGOUT,">".$outputPrefix.".diag.config") or die "cannot write to ".$outputPrefix.".diag.config";
  print FILECONFIGOUT "lengthDeam\t$lengthDeam\n";
  print FILECONFIGOUT "nameMT\t$nameMT\n";
  print FILECONFIGOUT "nameMTc\t$nameMTc\n";
  print FILECONFIGOUT "qualmin\t$qualmin\n";
  print FILECONFIGOUT "lengthMT\t$lengthMT\n";
  print FILECONFIGOUT "multipleC\t$multipleC\n";
  print FILECONFIGOUT "numthreads\t$numthreads\n";
  print FILECONFIGOUT "maxIterations\t$maxIterations\n";
  print FILECONFIGOUT "iterationSameCont\t$iterationSameCont\n";
  print FILECONFIGOUT "library\t$library\n";
  print FILECONFIGOUT "outputPrefix\t$outputPrefix\n";
  print FILECONFIGOUT "contPriorKnow\t$contPriorKnow\n";
  print FILECONFIGOUT "contPriorKnowCMDLine\t$contPriorKnowCMDLine\n";
  print FILECONFIGOUT "textGraph\t$textGraph\n";
  print FILECONFIGOUT "referenceFasta\t$referenceFasta\n";
  print FILECONFIGOUT "referenceFastaCMDL\t$referenceFastaCMDL\n";
  print FILECONFIGOUT "splitPos\t$splitPos\n";
  print FILECONFIGOUT "useLength\t$useLength\n";
  print FILECONFIGOUT "useLengthContDEAM\t$useLengthContDEAM\n";
  print FILECONFIGOUT "estdeam\t$estdeam\n";
  print FILECONFIGOUT "contPriorUser\t$contPriorUser\n";
  close(FILECONFIGOUT);

}
print "\tConfig file                                     : ".    $outputPrefix.".diag.config\n";
print "Total runtime ".(time - $starttime)." s\n";;




















