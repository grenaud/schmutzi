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

#my $dir = getcwd;
#my $abs_path = abs_path($dir);
my @arraycwd=split("/",abs_path($0));
pop(@arraycwd);
my $pathdir = join("/",@arraycwd);


my $installDirToFastaHaplogrep="/home/gabriel/scripts/fasta2haplogrep/";
my $installDir = $pathdir;#"/home/gabriel/projects/schmutzi";

warn "using script path: ".$installDir."/projects/schmutzi/schmutzi.pl"."\n";
warn "using haplogred  : ".$installDirToFastaHaplogrep."/fasta2haplogrep.py"."\n";

my @arrayOfTargets;

my $stringToPrint="";

foreach my $filebam (@ARGV){
  if($filebam !~ /.bam$/){
    next;
  }

  warn "Found bam file $filebam\n";
  my @filea      = split("/",$filebam);
  my $filebamname=$filea[$#filea];

  my $name=substr($filebamname,0,-4);



  push(@arrayOfTargets,substr($filebam,0,-4).".bai");
  $stringToPrint.="".substr($filebam,0,-4).".bai:\n\tsamtools index ".substr($filebamname,0,-4).".bam\n\n";

  push(@arrayOfTargets,substr($filebam,0,-4).".cont.est");
  $stringToPrint.= "".substr($filebam,0,-4).".cont.est: ".substr($filebam,0,-4).".bai\n\tif  ".$installDir."/contDeam.pl   --library double --lengthDeam 30 --out ".substr($filebamname,0,-4)." ".substr($filebam,0,-4).".bam; then echo \"command contDeam finished\"; else echo -e \"0.5\t0.49\t0.51\" > ".substr($filebam,0,-4).".cont.est; fi\n\n";

  push(@arrayOfTargets,substr($filebam,0,-4)."_wtpred_final.cont.est");
  $stringToPrint.= "".substr($filebam,0,-4)."_wtpred_final.cont.est: ".substr($filebam,0,-4).".cont.est\n\tif  ".$installDir."/schmutzi.pl  --iterations 5                 -t 2     --uselength   --ref ".$installDir."/refs/human_MT.fa  --out  ".substr($filebamname,0,-4)."_wtpred     ".substr($filebamname,0,-4)."      ".$installDir."/alleleFreqMT/eurasian/freqs/  ".substr($filebamname,0,-4).".bam; then echo \"command with pred finished\"; else echo \"command with pred stopped\"; fi\n\n";


  push(@arrayOfTargets,substr($filebam,0,-4)."_nopred_final.cont.est");
  $stringToPrint.= "".substr($filebam,0,-4)."_nopred_final.cont.est: ".substr($filebam,0,-4).".cont.est\n\tif  ".$installDir."/schmutzi.pl  --iterations 5   --notusepredC -t 2     --uselength   --ref ".$installDir."/refs/human_MT.fa  --out  ".substr($filebamname,0,-4)."_nopred     ".substr($filebamname,0,-4)."      ".$installDir."/alleleFreqMT/eurasian/freqs/  ".substr($filebamname,0,-4).".bam; then echo \"command with pred finished\"; else echo \"command with pred stopped\"; fi\n\n";



  for (my $q=10;$q<=50;$q+=20) {
    foreach my $type ("wtpred","nopred") {

      push(@arrayOfTargets,substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa");
      $stringToPrint.="".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa: ".substr($filebamname,0,-4)."_".$type."_final.cont.est\n\tif [ -e ".substr($filebamname,0,-4)."_".$type."_final_endo.log ];  then  tail -n+2 ".substr($filebamname,0,-4)."_".$type."_final_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".substr($filebamname,0,-4)."_".$type."_final_endo.cov; ".$installDir."//log2fasta -name ".$name."q10w -q ".$q." -indel ".$q." ".substr($filebamname,0,-4)."_".$type."_final_endo.log > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa; else if [ -e ".substr($filebamname,0,-4)."_".$type."_5_endo.log ];  then  tail -n+2 ".substr($filebamname,0,-4)."_".$type."_5_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".substr($filebamname,0,-4)."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q10w -q ".$q." -indel ".$q." ".substr($filebamname,0,-4)."_".$type."_5_endo.log > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa; else if [ -e ".substr($filebamname,0,-4)."_".$type."_4_endo.log ];  then  tail -n+2 ".substr($filebamname,0,-4)."_".$type."_4_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".substr($filebamname,0,-4)."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q10w -q ".$q." -indel ".$q." ".substr($filebamname,0,-4)."_".$type."_4_endo.log > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa; else if [ -e ".substr($filebamname,0,-4)."_".$type."_3_endo.log ];  then tail -n+2 ".substr($filebamname,0,-4)."_".$type."_3_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".substr($filebamname,0,-4)."_".$type."_final_endo.cov;  ".$installDir."/log2fasta -name ".$name."q10w -q ".$q." -indel ".$q." ".substr($filebamname,0,-4)."_".$type."_3_endo.log > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa; else if [ -e ".substr($filebamname,0,-4)."_".$type."_2_endo.log ];  then tail -n+2 ".substr($filebamname,0,-4)."_".$type."_2_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".substr($filebamname,0,-4)."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q10w -q ".$q." -indel ".$q." ".substr($filebamname,0,-4)."_".$type."_2_final_endo.log > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa; else if [ -e ".substr($filebamname,0,-4)."_".$type."_1_endo.log ];  then tail -n+2 ".substr($filebamname,0,-4)."_".$type."_1_endo.log |cut -f 6| awk 'BEGIN{sum=0;}{ sum+=\$\$1; }END{print sum/NR}'  > ".substr($filebamname,0,-4)."_".$type."_final_endo.cov; ".$installDir."/log2fasta -name ".$name."q10w -q ".$q." -indel ".$q." ".substr($filebamname,0,-4)."_".$type."_1_endo.log > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa; else echo \"error\"; fi; fi; fi; fi; fi; fi;\n\n";
    }
  }



  my @arrayOfTargetsHSD;

  for (my $q=10;$q<=50;$q+=20) {
    foreach my $type ("wtpred","nopred") {
      push(@arrayOfTargets,   substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".hsd");
      push(@arrayOfTargetsHSD,substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".hsd");

      $stringToPrint.="".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".hsd: ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa\n\tpython ".$installDirToFastaHaplogrep."/fasta2haplogrep.py ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".hsd\n\n";
    }
  }

  push(@arrayOfTargets,   substr($filebamname,0,-4)."_results.txt");
  $stringToPrint.="".substr($filebamname,0,-4)."_results.txt: ".join(" ",@arrayOfTargetsHSD)."\n\t".$installDir."/parseSchmutzi.pl ".substr($filebamname,0,-4)." ".substr($filebamname,0,-4)."  > ".substr($filebamname,0,-4)."_results.txt\n\n";
  #python ".$installDirToFastaHaplogrep."/fasta2haplogrep.py ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".fa > ".substr($filebamname,0,-4)."_".$type."_final_endo.q".$q.".hsd\n\n";

  
}





print "SHELL := /bin/bash\n\nall:\t".join(" ",@arrayOfTargets)."\n\nclean:\n\trm -vf ".join(" ",@arrayOfTargets)."\n\n".$stringToPrint."\n";
