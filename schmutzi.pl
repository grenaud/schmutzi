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

my $mock =0;
my $lengthDeam=2;


sub readMTCONToutlogfile{
  my ($mtcontOutputLog) = @_;

  print "reading the log file ".$mtcontOutputLog;

  my $sourceFile="###";
  my $contMaxS    = -1;
  my $llikMaxS    = -1;
  my $contMaxSALL = -1;
  my $llikMaxSALL = -1;
  my $contFileInd = 0;

  if ($mock != 1) {

    open(FILEoutlogmtcont, $mtcontOutputLog) or die "cannot open ".$mtcontOutputLog;

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


  }

  if( $contMaxSALL == -1){
    die "Cannot parse the contamination output file = ".$mtcontOutputLog."\n";
  }

  return $contMaxSALL;
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
  #if ($mock != 1) {
  print "reading the contamination estimate file: $filename\n";

  open(FILEcont, $filename) or die "cannot open ".$filename;

  my $line = <FILEcont>;
  my @arrayTemp = split("\t",$line);
  close(FILEcont);

  return $arrayTemp[0];
 # } else {
 #   print "reading the contamination estimate file: $filename\n";
 #   return 'X';
 # }
}


sub readSizeParam{
  my ($filename) = @_;

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




sub fileExists{
  my ($exeFile) = @_;

  if (!( -e $exeFile)) {
    die "Executable ".$exeFile." does not exist\n";
  }
}

my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);

my $mtCont     = $pathdir."/mtCont";
my $endoCaller = $pathdir."/endoCaller";
my $insertSize = $pathdir."/insertSize";
my $approxDist = $pathdir."/approxDist.R";
my $bam2prof   = $pathdir."/bam2prof";
my $log2freq   = $pathdir."/log2freq";
my $logs2pos   = $pathdir."/logs2pos";

my $contDeam   = $pathdir."/contDeam";
my $contDeamR  = $pathdir."/posteriorDeam.R";
my $splitEndo  = $pathdir."/splitEndoVsCont/poshap2splitbam";

fileExists($mtCont);
fileExists($endoCaller);
fileExists($insertSize);
fileExists($approxDist);
fileExists($bam2prof);
fileExists($log2freq);
fileExists($logs2pos);
fileExists($contDeam);
fileExists($contDeamR);
fileExists($splitEndo);


my $nameMT    = "MT";
my $nameMTc   = "MTc";
my $qualmin   = 0;
my $lengthMT=16569;

#
# protocol s or d
# output prefix


my $multipleC=0;
my $numthreads=1;

sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "\n\n
This script is a wrapper to produce the endogenous mitochondrial
consensus genome and a contamination estimate given aligned
mitochondrial data and a set of putative contaminants. This script
takes the ouput of the contDeam.pl script and calls iteratively
the module that produces an endogenous consensus and the one that
estimates contamination. For non-hominin samples without risk of
human contamination where the mitochondrial consensus is needed,
call endoCaller directly (see README.md).

\n\nusage:\t".$0." <options> [output prefix from contDeam.pl] [directory with *freq file] input.bam\n\n".

#"Mandatory parameters:\n".
#"\t[output prefix from contDeam]\tThis is the prefix of the output of contDeam.pl (option --out)\n".
#"\t[directory with *freq file]\tDirectory with all the *freq file containing the putative contaminants\n".
#"\t[reference *fasta]\tFasta file of the reference used for by aligner\n".
#"\t[input.bam]\tAligned BAM file\n\n".
#

"Options:\n".
"\n\t--t (threads #)\t\t\t\tUse this amount of threads (Default : $numthreads)\n\n".
"\n\t--lengthMT (bp)\t\t\t\tConsider that the real length of the MT genome is this (Default : $lengthMT)\n\n".

"\n\t--estdeam\tRe-estimate deamination parameters using newly computed segregating positions".
"\n\t\tThis setting is recommended if the contamination is deaminated".
"\n\t\t".

"\n\t--mock\t\t\t\tDo nothing, just print the commands used\n\n".
"\n\t--uselength\t\t\t\tUse length of the molecules as well\n\n".
"\n\t--lengthDeam (bp)\t\t\t\tOnly consider this about of bases to be deaminated on each end (Default : $lengthDeam)\n\n".
"\n\t--multipleC\tDo not assume that there is a single contaminant".
"\n\t\t\tThis might lead to worse results".

#
"Output Options:\n".
"\t--name (name)\t\tName of the endogenous MT genome\n".
"\t--namec (name)\t\tName of the contaminant MT genome\n".
"\t--qual (qual)\t\tMinimum quality for a base in the  MT consensus (PHRED scale)\n".


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

my $help;
my $library        = "none";
my $outputPrefix   = "outputdeam";

my $referenceFasta = "none";
my $contPriorKnow  = -1;
my $textGraph      = "Posterior probability for contamination\\nusing deamination patterns";
my $splitPos       = "";
my $useLength      = 0;
my $useLengthContDEAM      = 0; #from contdeam

my $estdeam = 0 ;

usage() if ( @ARGV < 1 or
	     ! GetOptions('help|?' => \$help, 't=i' => \$numthreads, 'mock' => \$mock, 'estdeam' => \$estdeam, 'uselength' => \$useLength, 'lengthDeam' => \$lengthDeam,'lengthMT' => \$lengthMT,'multipleC' => \$multipleC,'qual=f' => \$qualmin,'name=s' => \$nameMT,'namec=s' => \$nameMTc )
          or defined $help );

my $prefixcontDeam = $ARGV[ $#ARGV -2 ];
my $configfiledeam = $prefixcontDeam.".config";
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
  if($line =~ /^textGraph\s(\S+)$/){      $textGraph=$1;}
  if($line =~ /^splitPos\s(\S+)$/){       $splitPos=$1;}
  if($line =~ /^useLength\s(\S+)$/){      $useLengthContDEAM=$1;}
  if($line =~ /^splitDeam\s(\S+)$/){      $splitDeam=$1;}
}

close(FILEcontdeam);


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
my $maxIterations  =100;

########################
#DEAMINATION PARAMETERS#
########################
if(! $estdeam  ){ #do not re-estimate at each iteration


  copycmd( $prefixcontDeam.".endo.5p.prof" , $prefixcontDeam."_".$numberIteration."_endo.5p.prof" );
  copycmd( $prefixcontDeam.".endo.3p.prof" , $prefixcontDeam."_".$numberIteration."_endo.3p.prof" );

  if($splitDeam){
    copycmd( $prefixcontDeam.".cont.5p.prof" , $prefixcontDeam."_".$numberIteration."_cont.5p.prof" );
    copycmd( $prefixcontDeam.".cont.3p.prof" , $prefixcontDeam."_".$numberIteration."_cont.3p.prof" );
  }else{#put dummy values for cont, the split was not done by contDeam and we will do it after the first iteration
    makeEmptyProfFile($prefixcontDeam."_".$numberIteration."_cont.5p.prof" );
    makeEmptyProfFile($prefixcontDeam."_".$numberIteration."_cont.3p.prof" );
  }

  copycmd( $prefixcontDeam.".cont.est" , $prefixcontDeam."_".$numberIteration."_cont.est" );

}


########################
#  LENGTH PARAMETERS   #
########################
if($useLength){#
  if($useLengthContDEAM){ #if previously computed

    copycmd( $prefixcontDeam."_endo.size.param" , $prefixcontDeam."_".$numberIteration."_endo.size.param" );
    copycmd( $prefixcontDeam."_cont.size.param" , $prefixcontDeam."_".$numberIteration."_cont.size.param" );

  }else{
    #
  }
}


while(1){

  # call endoCaller
  my $cmdEndoCaller = $endoCaller." ";

  $cmdEndoCaller .= " -seq ".$prefixcontDeam."_".$numberIteration."_endo.fa  ";
  $cmdEndoCaller .= " -log ".$prefixcontDeam."_".$numberIteration."_endo.log ";

  $cmdEndoCaller .= " -name ".$nameMT." ";
  $cmdEndoCaller .= " -qual ".$qualmin." ";

  $cmdEndoCaller .= " -deamread ";
  $cmdEndoCaller .= " -deam5p ".$prefixcontDeam."_".$numberIteration."_endo.5p.prof ";
  $cmdEndoCaller .= " -deam3p ".$prefixcontDeam."_".$numberIteration."_endo.3p.prof ";
  #TODO put contamination deamination

  my $currentContEst=readcont( $prefixcontDeam."_".$numberIteration."_cont.est" );
  $cmdEndoCaller .= " -cont ".$currentContEst." ";

  if(!$multipleC){ #we can assume a single contaminant
    $cmdEndoCaller .= " -single  ";

    $cmdEndoCaller .= " -seqc ".$prefixcontDeam."_".$numberIteration."_cont.fa  ";
    $cmdEndoCaller .= " -logc ".$prefixcontDeam."_".$numberIteration."_cont.log ";
    $cmdEndoCaller .= " -namec ".$nameMTc." ";

  }

  $cmdEndoCaller .= " -l ".$lengthMT." ";

  # mol. length
  if ($useLength) {		#

    my $locE=1;
    my $scaE=1;
    my $locC=1;
    my $scaC=1;
    if($numberIteration == 1){ #first iteration
      if($useLengthContDEAM){ #if previously computed

	($locE,$scaE)=   readSizeParam($prefixcontDeam."_".$numberIteration."_endo.size.param" );
	($locC,$scaC)=   readSizeParam($prefixcontDeam."_".$numberIteration."_cont.size.param" );

      }else{
	#leave dummy values
      }

    }else{

      ($locE,$scaE)=   readSizeParam($prefixcontDeam."_".$numberIteration."_endo.size.param" );
      ($locC,$scaC)=   readSizeParam($prefixcontDeam."_".$numberIteration."_cont.size.param" );

    }
  }



  $cmdEndoCaller .= " $referenceFasta $inbam ";
  runcmd($cmdEndoCaller);

  my @listOfFreqFiles=@arrayOfFreqFiles;

  #convert log to freq
  if(!$multipleC){ #we can assume a single contaminant

    my $cmdLog2Freq =  $log2freq." ".$prefixcontDeam."_".$numberIteration."_cont.log >  ".$prefixcontDeam."_".$numberIteration."_cont.freq";
    runcmd($cmdLog2Freq);
    push(@listOfFreqFiles, $prefixcontDeam."_".$numberIteration."_cont.freq");
  }

  # estimate cont
  my $cmdmtcont =   $mtCont. " -o ".$prefixcontDeam."_".$numberIteration."_mtcont.out ";
  #$cmdmtcont =   " -deam5p "
  $cmdmtcont .= " -deam5p ".$prefixcontDeam."_".$numberIteration."_endo.5p.prof ";
  $cmdmtcont .= " -deam3p ".$prefixcontDeam."_".$numberIteration."_endo.3p.prof ";
  $cmdmtcont .= " -t ".$numthreads." ";
  $cmdmtcont .= " ".$prefixcontDeam."_".$numberIteration."_endo.log ";
  $cmdmtcont .= " $referenceFasta $inbam ";
  $cmdmtcont .= " ".join(" ",@listOfFreqFiles). " ";
  runcmd($cmdmtcont);



  # if likelihood stable, break loop
  if($numberIteration >= $maxIterations ){
    last;
  }



  #BEGIN read new cont esti. log
  my $mtcontOutputLog = $prefixcontDeam."_".$numberIteration."_mtcont.out";

  my $contFROMmtcont =  readMTCONToutlogfile($mtcontOutputLog);

  #END   read new cont esti. log


  if (!$multipleC) {		#we can assume a single contaminant

    #print seg. sites
    my $cmdLogs2Pos = $logs2pos."  ".$prefixcontDeam."_".$numberIteration."_endo.log ".$prefixcontDeam."_".$numberIteration."_cont.log  > ".$prefixcontDeam."_".$numberIteration."_split.pos";
    runcmd($cmdLogs2Pos);

    # split according to seg sites
    my $cmdBamSplit = $splitEndo."  ".$prefixcontDeam."_".$numberIteration."_split.pos  $inbam ".$prefixcontDeam."_".($numberIteration)."_split > ".$outputPrefix."_split.log 2> /dev/null ";
    runcmd($cmdBamSplit);


    #re-measure deam rateS
    if( $estdeam  ){ # re-estimate deamination at each iteration

      my $cmdBam2ProfEndo = $bam2prof." -length $lengthDeam -".$library." -5p ".$prefixcontDeam."_".($numberIteration+1)."_endo.5p.prof   -3p ".$prefixcontDeam."_".($numberIteration+1)."_endo.3p.prof  ".$prefixcontDeam."_".($numberIteration)."_split_endo.bam";
      runcmd($cmdBam2ProfEndo);

      my $cmdBam2ProfCont = $bam2prof." -length $lengthDeam -".$library." -5p ".$prefixcontDeam."_".($numberIteration+1)."_cont.5p.prof   -3p ".$prefixcontDeam."_".($numberIteration+1)."_cont.3p.prof  ".$prefixcontDeam."_".($numberIteration)."_split_cont.bam";
      runcmd($cmdBam2ProfCont);

      #my $cmdBam2ProfCont = $bam2prof." -length $lengthDeam -".$library." -5p ".$outputPrefix.".cont.5p.prof  -3p ".$outputPrefix.".cont.3p.prof ".$outputPrefix."_cont.bam";
      #runcmd($cmdBam2ProfCont);
    }else{

      copycmd(  $prefixcontDeam."_".$numberIteration."_endo.5p.prof" ,$prefixcontDeam."_".($numberIteration+1)."_endo.5p.prof" );
      copycmd(  $prefixcontDeam."_".$numberIteration."_endo.3p.prof" ,$prefixcontDeam."_".($numberIteration+1)."_endo.3p.prof" );
      #copycmd(  $prefixcontDeam."_".$numberIteration."_cont.5p.prof" ,$prefixcontDeam."_".($numberIteration+1)."_cont.5p.prof" );
      #copycmd(  $prefixcontDeam."_".$numberIteration."_cont.3p.prof" ,$prefixcontDeam."_".($numberIteration+1)."_cont.3p.prof" );

    }


    #measure insert size
    my $cmdBam2InsertsizeEndo = $insertSize."   ".$prefixcontDeam."_".($numberIteration)."_split_endo.bam |gzip  > ".$prefixcontDeam."_".($numberIteration)."_split_endo.size.gz";
    runcmd($cmdBam2InsertsizeEndo);
    my $cmdBam2InsertsizeCont = $insertSize."   ".$prefixcontDeam."_".($numberIteration)."_split_cont.bam  |gzip  > ".$outputPrefix."_cont.size.gz";
    runcmd($cmdBam2InsertsizeCont);

    #Log normal fit
    my $cmdinsertsize2LognormEndo = $approxDist."  ".$outputPrefix."_endo.size.gz >  ".$outputPrefix."_endo.size.param";
    runcmd($cmdinsertsize2LognormEndo);
    my $cmdinsertsize2LognormCont = $approxDist."  ".$outputPrefix."_cont.size.gz >  ".$outputPrefix."_cont.size.param";
    runcmd($cmdinsertsize2LognormCont);

  }else{ #must use the previous deamination profiles from the previous iteration for the new one

    copycmd(  $prefixcontDeam."_".$numberIteration."_endo.5p.prof" ,$prefixcontDeam."_".($numberIteration+1)."_endo.5p.prof" );
    copycmd(  $prefixcontDeam."_".$numberIteration."_endo.3p.prof" ,$prefixcontDeam."_".($numberIteration+1)."_endo.3p.prof" );


  }

  # measure deam params + length
  

  if($mock == 1){
    die;
  }

  $numberIteration++;
}
# go to begin loop






















