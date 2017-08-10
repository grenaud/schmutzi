#!/usr/bin/perl


use strict;
use warnings;
use Data::Dumper;
use POSIX;
use Cwd 'abs_path';

sub fileExists{
  my ($exeFile) = @_;

  return ( -e $exeFile);

}
my $threads="2";
my $library="double";




#my $dir = getcwd;
#my $abs_path = abs_path($dir);
my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);
my $skipContDeam=0;
my $skipPred=0;
my $subsample=-1;
my $nice="";

#my $installDirToFastaHaplogrep="/home/gabriel/scripts/fasta2haplogrep/";
my $installDirToFastaHaplogrep = $pathdir;
my $installDir = $pathdir;#"/home/gabriel/projects/schmutzi";
my $iterations=5;

print STDERR  "        ~~~~ Please read carefully ~~~~\n";
print STDERR  "using script path: ".$installDir."/projects/schmutzi/schmutzi.pl"."\n";
print STDERR  "using haplogred  : ".$installDirToFastaHaplogrep."/fasta2haplogrep.py"."\n";
print STDERR  "if these paths are incorrect, change them in the Perl script\n";
print STDERR  "Also, fasta2haplogrep.py uses haplogrep-cmd.jar, please hard code the path in this script\n";

my $usage= "\n\n usage:\t".$0." <options> [bam file1] [bam file2]...\n\n".
  " Options:\n".
  "\t--skippred\t\t\t\tSkip contamination estimate using prediction of the contaminant (useful for low cont. samples)\n".
  "\t--nodeam\t\t\t\tSkip contamination based on deamination, useful for UDG treated\n".
  "\t--single\t\t\t\tUse single stranded library damage (just C->T on both ends)\n".
  "\t--subsample  [XXX]\t\t\tSubsample the BAM file down to XXX, 100-300 depending on how difficult the target is are good\n".
  "\t--threads    [num]\t\t\tUse [num] of threads (default $threads)\n".
  "\t--iterations [num]\t\t\tMaximum number of iterations (default $iterations)\n".
  "\t--nice\t\t\t\t\tUse nice\n\n\n";

if($#ARGV == -1){
  die $usage;
}


my @arrayOfTargets;
my @arrayOfTargetsClean;

my $stringToPrint="";
my $i=0;

foreach my $filebam (@ARGV){
  $i++;

  if($filebam =~ /^--/){
    if($filebam eq "--nodeam"){
      $skipContDeam=1;
      next;
    }

    if($filebam eq "--skippred"){
      $skipPred=1;
      next;
    }

    if($filebam eq "--nice"){
      $nice=" nice -n 19 ";
      next;
    }

    if($filebam eq "--single"){
      $library="single";
      next;
    }

    if($filebam eq "--threads"){
      if($ARGV[$i] !~ /^\d+$/){
	die "Wrong --threads parameter ".$ARGV[$i]."\n";
      }
      $threads=$ARGV[$i];
      next;
    }


    if($filebam eq "--subsample"){
      if($ARGV[$i] !~ /^\d+$/){
	die "Wrong --subsample parameter ".$ARGV[$i]."\n";
      }
      $subsample=$ARGV[$i];
      next;
    }


    if($filebam eq "--iterations"){
      if($ARGV[$i] !~ /^\d+$/){
	die "Wrong --iterations parameter ".$ARGV[$i]."\n";
      }
      $iterations=$ARGV[$i];
      next;
    }



    die "Unknown option ".$filebam."\n";
  }

  if($filebam !~ /.bam$/){
    next;
  }

  print STDERR "Found bam file $filebam\n";
  my @filea      = split("/",$filebam);
  my $filebamname=$filea[$#filea]; #file name without path


  my $name      = substr($filebamname,0,-4);
  my $outprefix = substr($filebam,0,-4);#path without extension

  #my $prefixFilename = 
  if( $subsample != -1){
    push(@arrayOfTargets,     $outprefix."_sub.bam");
    push(@arrayOfTargetsClean,$outprefix."_sub.bam");
    $stringToPrint.="".$outprefix."_sub.bam:\n\t".$installDir."/subSampleBAM  ".$subsample." ".$outprefix.".bam  ".$outprefix."_sub.bam\n\n";
    $outprefix = $outprefix."_sub";
  }

  push(@arrayOfTargets,     $outprefix.".bai");
  push(@arrayOfTargetsClean,$outprefix.".bai");

  $stringToPrint.="".$outprefix.".bai:\n\tsamtools index ".$outprefix.".bam\n\n";

  push(@arrayOfTargets,     $outprefix.".cont.est");
  push(@arrayOfTargetsClean,$outprefix.".cont*");
  push(@arrayOfTargetsClean,$outprefix.".deam.config");
  push(@arrayOfTargetsClean,$outprefix.".endo.5p.prof");
  push(@arrayOfTargetsClean,$outprefix.".endo.3p.prof");

  if($skipContDeam == 1){
    $stringToPrint.= "".$outprefix.".cont.est: ".$outprefix.".bai\n\techo -e \"0.5\t0.49\t0.51\" > ".$outprefix.".cont.est\n\t/opt/schmutzi/bam2prof -5p ".$outprefix.".endo.5p.prof -3p ".$outprefix.".endo.3p.prof ".$outprefix.".bam\n\techo -e \"library\\t".$library."\\noutputPrefix\\t".$outprefix."\\nreferenceFasta\\t\\ncontPriorKnow\\t-1\\ntextGraph\\tPosterior probability for contamination\\\\\\nusing deamination patterns\\nsplitPos\\t\\nuseLength\\t0\\nsplitDeam\\t0\" > ".$outprefix.".deam.config\n\n";
  }else{
    $stringToPrint.= "".$outprefix.".cont.est: ".$outprefix.".bai\n\tif  ".$installDir."/contDeam.pl   --library ".$library." --lengthDeam 30 --out ".$outprefix." ".$outprefix.".bam; then echo \"command contDeam finished\"; else echo -e \"0.5\t0.49\t0.51\" > ".$outprefix.".cont.est; fi\n\n";
  }

  if (!$skipPred ) {
    push(@arrayOfTargets,     $outprefix."_wtpred_final.cont.est");
    push(@arrayOfTargetsClean,$outprefix."_wtpred*");

    $stringToPrint.= "".$outprefix."_wtpred_final.cont.est: ".$outprefix.".cont.est\n\tif $nice ".$installDir."/schmutzi.pl  --iterations $iterations                 -t $threads    --uselength   --ref ".$installDir."/refs/human_MT.fa  --out  ".$outprefix."_wtpred     ".$outprefix."      ".$installDir."/alleleFreqMT/eurasian/freqs/  ".$outprefix.".bam; then echo \"command with pred finished\"; else echo \"command with pred stopped\"; fi\n\n";
  }

  push(@arrayOfTargets,     $outprefix."_nopred_final.cont.est");
  push(@arrayOfTargetsClean,$outprefix."_nopred*");

  $stringToPrint.= "".$outprefix."_nopred_final.cont.est: ".$outprefix.".cont.est\n\tif $nice ".$installDir."/schmutzi.pl  --iterations  $iterations   --notusepredC -t $threads     --uselength   --ref ".$installDir."/refs/human_MT.fa  --out  ".$outprefix."_nopred     ".$outprefix."      ".$installDir."/alleleFreqMT/eurasian/freqs/  ".$outprefix.".bam; then echo \"command with pred finished\"; else echo \"command with pred stopped\"; fi\n\n";


  my @typeOfTargets = ("wtpred","nopred");
  if($skipPred){
    @typeOfTargets = ("nopred");
  }

  for (my $q=10;$q<=50;$q+=20) {
    foreach my $type (@typeOfTargets) {

      push(@arrayOfTargets,     $outprefix."_".$type."_final_endo.q".$q.".fa");
      push(@arrayOfTargetsClean,$outprefix."_".$type."_final_endo.q".$q.".fa");
      $stringToPrint.="".$outprefix."_".$type."_final_endo.q".$q.".fa: ".$outprefix."_".$type."_final.cont.est\n\t".
	"if [ -e ".$outprefix."_".$type."_final_endo.log ];  then  tail -n+2 ".$outprefix."_".$type."_final_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".$outprefix."_".$type."_final_endo.cov; ".$installDir."//log2fasta -name ".$name."q".$q."w -q ".$q." -indel ".$q." ".$outprefix."_".$type."_final_endo.log > ".$outprefix."_".$type."_final_endo.q".$q.".fa; else ".
        "if [ -e ".$outprefix."_".$type."_5_endo.log ];  then  tail -n+2 ".$outprefix."_".$type."_5_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".$outprefix."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q".$q."w -q ".$q." -indel ".$q." ".$outprefix."_".$type."_5_endo.log > ".$outprefix."_".$type."_final_endo.q".$q.".fa; else ".
        "if [ -e ".$outprefix."_".$type."_4_endo.log ];  then  tail -n+2 ".$outprefix."_".$type."_4_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".$outprefix."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q".$q."w -q ".$q." -indel ".$q." ".$outprefix."_".$type."_4_endo.log > ".$outprefix."_".$type."_final_endo.q".$q.".fa; else ".
	"if [ -e ".$outprefix."_".$type."_3_endo.log ];  then tail -n+2 ".$outprefix."_".$type."_3_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".$outprefix."_".$type."_final_endo.cov;  ".$installDir."/log2fasta -name ".$name."q".$q."w -q ".$q." -indel ".$q." ".$outprefix."_".$type."_3_endo.log > ".$outprefix."_".$type."_final_endo.q".$q.".fa; else ".
        "if [ -e ".$outprefix."_".$type."_2_endo.log ];  then tail -n+2 ".$outprefix."_".$type."_2_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".$outprefix."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q".$q."w -q ".$q." -indel ".$q." ".$outprefix."_".$type."_2_final_endo.log > ".$outprefix."_".$type."_final_endo.q".$q.".fa; else ".
        "if [ -e ".$outprefix."_".$type."_1_endo.log ];  then tail -n+2 ".$outprefix."_".$type."_1_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".$outprefix."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q".$q."w -q ".$q." -indel ".$q." ".$outprefix."_".$type."_1_endo.log > ".$outprefix."_".$type."_final_endo.q".$q.".fa; else echo \"error\"; fi; fi; fi; fi; fi; fi;\n\n";
    }
  }



  my @arrayOfTargetsHSD;

  for (my $q=10;$q<=50;$q+=20) {
    foreach my $type (@typeOfTargets){
      push(@arrayOfTargets,        $outprefix."_".$type."_final_endo.q".$q.".hsd");
      push(@arrayOfTargetsClean,   $outprefix."_".$type."_final_endo.q".$q.".hsd");
      push(@arrayOfTargetsHSD,     $outprefix."_".$type."_final_endo.q".$q.".hsd");

      $stringToPrint.="".$outprefix."_".$type."_final_endo.q".$q.".hsd: ".$outprefix."_".$type."_final_endo.q".$q.".fa\n\tpython ".$installDirToFastaHaplogrep."/fasta2haplogrep.py ".$outprefix."_".$type."_final_endo.q".$q.".fa > ".$outprefix."_".$type."_final_endo.q".$q.".hsd\n\n";
    }
  }

  push(@arrayOfTargets,        $name."_results.txt");
  push(@arrayOfTargetsClean,   $name."_results.txt");

  $stringToPrint.="".$name."_results.txt: ".join(" ",@arrayOfTargetsHSD)."\n\t".$installDir."/parseSchmutzi.pl ".$outprefix." ".$name."  > ".$name."_results.txt\n\n";
  #python ".$installDirToFastaHaplogrep."/fasta2haplogrep.py ".$name."_".$type."_final_endo.q".$q.".fa > ".$name."_".$type."_final_endo.q".$q.".hsd\n\n";

  
}





print "SHELL := /bin/bash\n\nall:\t".join(" ",@arrayOfTargets)."\n\nclean:\n\trm -vf ".join(" ",@arrayOfTargetsClean)."\n\n".$stringToPrint."\n";
