=====================================================================================
  schmutzi: Bayesian maximum a posteriori contamination estimate for ancient samples
=====================================================================================

Upon sequencing ancient DNA , the DNA of the individuals involved in excavation, extraction 
of DNA and library preparation can become mixed with the actual  sample being sequenced. 
We define as endogenous the DNA pertaining to the original sample and contaminant the DNA of the experimenters. 


schmutzi is a set of programs aimed at ancient DNA data that can :

	- estimate contamination based on deamination patterns alone
	- call a mitonchondrial consensus for the endogenous genome. This consensus calls 
	  takes into account contamination and deamination.
	- estimate mitonchondrial contamination and identify the most likely contaminant 
	  from a set.

Questions :
-------------------------------------------------------------------------------------
	contact: Gabriel Renaud   
   	email:	 gabriel [dot] reno [ at sign ] gmail.com

Downloading:
-------------------------------------------------------------------------------------
Do a "git clone --recursive https://github.com/grenaud/schmutzi.git"
or download a zipped file from https://bioinf.eva.mpg.de/schmutzi/schmutzi.tar.gz

Requirements:
-------------------------------------------------------------------------------------

   - zlib.h
   - cmake to build bamtools
   - C++ compiler
   - Perl interpreter
   - R
      -  The fitdistrplus R package
      -  The MASS package

On Ubuntu, these dependencies can be resolved using :
  sudo apt-get install perl
  sudo apt-get install git
  sudo apt-get install cmake
  sudo apt-get install g++
  sudo apt-get install zlib
  sudo apt-get install zlib1g-dev
  sudo apt-get install r-base-core


Installation:
-------------------------------------------------------------------------------------


For Mac users, open a terminal please type "cmake" to check if works. If it is not installed, please install it. 
For Windows users, you would need cygwin to run schmutzi as it runs on the terminal (https://www.cygwin.com/). 

1) type "make" which should build pretty much everything. 

We tested our program on a Linux system and MacOS, please email us if you have trouble building under cygwin. 

A comment about cygwin: the high precision version of the logarithm does not seem available by default. Therefore, 
we are forced to use the low precision variant for cygwin. For higher accuracy, run our software under Linux. Also, 
some users have reported bugs with the locale used by C++ 


Running:
-------------------------------------------------------------------------------------
There are 2 main Perl scripts that drive the underlying programs written in C++.


contDeam.pl   Is a script that allows for the estimation of endogenous deamination and
              subsequent contamination estimate using those. Assuming that the contaminant
              is not deaminated, it will measure rates of deamination for the endogenous material 
	      and provide a rough contamination estimate. If say the endogenous material is 
              30% deaminated at one base but the observed deamination rate for the entire set is 15%,
              the contamination rate is therefore 0.5 . There are two ways to estimate endogenous 
              deamination:
		1) (default)
 		   Condition on one end of the read being deaminated (ex: 5') and measure 
 		   deamination rates on the other (ex: 3'). This approach is useful for 
                   all types of ancient DNA (nuclear and mitonchondrial). Please note however
                   that if the contaminant is deaminated as well, which can occur, this approach 
                   will return an underestimate whose error comensurates with rates of contaminant 
                   deamination.
                2) --split
                   Take a series of mitonchondrial diagnostic positions (positions that indicate 
                   whether the read pertains to a particular type of homonin (ex: Neandertal) or 
                   the putative contaminant homonin (ex: modern Humans) and separate the contaminant
                   from the endogenous reads. This approach is useful when prior information on the 
                   sample is available and when enough diagnostic positions are available (case of 
                   archaic homonins).


schmutzi.pl   This script performs two tasks: 
	      	   1) Call an endogenous consensus 
		   2) Using the consensus called in 1) and a set of known contaminants, compute 
                      the contamination rate for entry in this set. Return the most likely contaminant 
		      with its corresponding rate. 
	      Both 1) and 2) are called iteratively until a stable contamination rate is reached. The script 
	      produces a fasta file of the endogenous mitonchondrial genome 
                   
                       

Also, the various subparts of schmutzi can be called directly:

main modules:
contDeam  	   	Contamination estimate using deamination patterns
endoCaller	  	Consensus calling for mitochondrial data
mtCont			Contamination estimate using a set of known contaminants

sub modules:
approxDist.R		Do a maximum likelihood estimate of log-normal paramaters given molecule length data
bam2prof		Measure deamination rates and produce a matrix of deamination probabilities
insertSize		Produce all insert sizes for aligned data in BAM format
log2freq		Turn a endogenous or consensus log file into an allele frequency matrix to be used 
			as contaminant source
logs2pos		Takes two log files and produces which positions are segregating
mitoConsPDF.R		Plot various information like coverage and quality for the consensus log file
mas2freq		Turn a multiple sequence alignment file in multifasta format into an allele frequency 
			matrix to be used as contaminant source
posteriorDeam.R		Plot the posterior probability


Test Data:
-------------------------------------------------------------------------------------

To test schmutzi, we have made available emprical data with single-strand (ex: testdata/mezB9687.bam) 
and simulated data with double-stranded damage patterns (ex: testdata/simulation.bam). 


To download it, either download it manually from :
   https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam
   https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam.bai
   https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam
   https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam.bai

Or, if you have "wget" installed, just type:
    make testdata

First you need to estimate endogenous deamination rates. First create an output directory:
    mkdir outputdir/

Then run contDeam to estimation endogenous deamination rates:
    ./contDeam.pl  --library single --out outputdir/mez testdata/mezB9687.bam
or for the simulated
    ./contDeam.pl  --library double --out outputdir/sim testdata/simulation.bam    

This will produce the files:
    outputdir/[out].cont.pdf	Plot of the posterior probability for contamination based on deamination
    outputdir/[out].cont.est      Estimate for contamination based on deamination
    outputdir/[out].config	Configuration file describing the variables used


Then run the following to produce the endogenous consensus and the contamination estimate:

        ./schmutzi.pl       --uselength   --ref refs/human_MT_wrapped.fa         outputdir/mez   alleleFreqMT/197/freqs/  testdata/mezB9687.bam
        ./schmutzi.pl       --uselength   --ref refs/human_MT_wrapped.fa         outputdir/sim   alleleFreqMT/197/freqs/  testdata/simulation.bam
 	  
It will run for a few minutes and produce the following files:

For contamination:
	  [out]_final_mtcont.out        Contamination estimates for all samples in the database
	  [out]_final.cont		Contamination estimates for the sample in the database with the highest likelihood
	  [out]_final.cont.est		Contamination estimates for the most likely sample with confidence intervals 	  
	  [out]_final.cont.pdf		Posterior probability on the contamination for the most likely sample

For the respective genomes of the endogenous and contaminant:	  
	  [out]_final_endo.fa		Endogenous mitochondrial genome as fasta
	  [out]_final_endo.log		Endogenous mitochondrial genome as a log file with likelihoods on a PHRED scale
	  
	  [out]_final_cont.fa		Contaminant mitochondrial genome as fasta
	  [out]_final_cont.log		Contaminant mitochondrial genome as a log file with likelihoods on a PHRED scale



Recommended workflow for ancient samples:
-------------------------------------------------------------------------------------

- Nuclear
  The only thing schmutzi can do for nuclear data is to estimate contamination rates 
  using deamination patterns. However, use this at your own risk as you will not know
  if you have a contaminant that is deaminated thus leading to an underestimate. Simply
  use "contDeam.pl"

- MT
  1) Have your data aligned to a mitochondrial reference (see "refs/human_MT_wrapped.fa" for a wrapped reference) using 
    a sensitive mapper that produces BAM. A wrapper script is available with schmutzi (see Frequently Asked Questions)
  2) run samtools sort on your aligned bam file
  3) run samtools index on your sorted and aligned bam file
  4) If you used the wrapped reference, re-wrap your alignments exceeding the junction
     using for example https://github.com/udo-stenzel/biohazard  this step is not necessary 
     but produces equal coverage and resolution at the boundaries
  5) Estimate your initial contamination and deamination rates using "contDeam.pl"
  6) Run schmutzi.pl once with default parameters
  7) Run schmutzi.pl again  with "--usepredC"

  8) If contamination is more than a few percent re-run using the "--uselength" option
     


Frequently asked questions:
-------------------------------------------------------------------------------------

- How should I prepare the BAM file as input for endoCaller ?
   
   As we mention above the, data must be aligned and ideally, wrapped around the ends.
   We provide a little wrapper script to call the shrimp mapper and wrap the reads around
   the edge : wrapper.pl


- When should I trust the output of contDeam.pl ?

  When the deamination rates for the contaminant material are negligible. To 
  verify this, one way is to use diagnostic positions on the mitochondria. If 
  you have nuclear data, please retain the reads mapping to the mitochondria. 
  If you know in advance that the sample is either Neandertal or Denisova, you can
  use  the precomputed diagnostic positions. If you have an modern human, you 
  will have to run schmutzi.pl and use the split BAM files produced by the script. 
  Inspect, the files called *cont.5p.prof and *cont.3p.prof. If the rates of deamination
  for them is higher than 1 or 2 \%, the results provided by contDeam.pl will not 
  reliable.
  
- The contamination estimate given by contDeam.pl gives me a very different result
  than schmutzi.pl does, why is that ?

  There could be a few reasons:
  1) Your contamination is deaminated. That can be the case for museum samples.
     Check to see if the 
        - [output prefix]_[iteration]_cont.5p.prof 
        and 
        - [output prefix]_[iteration]_cont.3p.prof 
     present any sign of deamination.

  2) The algorithm diverged. This could be due to a misleading prior or an almost absence
     of endogenous material.

- How can I know if I have multiple contaminants ?
    If you have highly divergent quality values for bases in the contamination prediction,
    this could be an indication that you might have multiple contaminants. Further, if  
    the predicted contaminant does not seem to fall within a given haplogroup and has
    diagnostic positions from different haplogroups, this could be an indication that 
    your contaminant is multiple. 

- What to do if I suspect that I have multiple contaminants ?
    Do not use the --usepredC option and turn on the --multipleC option. If the contamination 
    rate is a lot lower than the one predicted by the deamination patterns, this could be an 
    indication that schmutzi did not converge and contamination is very high. Try to do a simple prediction
    using deaminated molecules (hoping that the contaminant is not deaminated) and use this as the input of 
    mtCont manually. Datasets that are heavily contaminated with multiple contaminants are very difficult 
    targets.

- Can I build my own set of putative contaminants ?

  Yes. First build a multiple sequence alignment (msa) of the various mitochondria
  including the reference used for mapping. Go to the directory where you want 
  to store the frequencies then run:
  
  ./msa2freq [msa] [name of the reference]
  
  A freqs/ directory should appear with one file per record in the msa.

- Can your software work for ancient DNA from other species ?
  
  Yes, given a reference, schmutzi can be used to call an mitochondria endogenous 
  consensus taking into account deamination. But in the case of non-hominin 
  animals, you do no need to measure contamination from humans if the 
  mitochondrion is sufficiently divergent from the human one. High divergence 
  entails that few reads from humans will align to your reference.
  
  First determine your deamination rates using bam2prof with the appropriate library type:
  
  ./bam2prof [-single|-double] -5p out5p.prof -3p out3p.prof input.bam
  
  Then call the consensus:
  ./endoCaller -cont 0 -deam5p out5p.prof -deam3p out3p.prof -seq output.fa -log outlog.log -l [length reference] /path/to/reference.fasta input.bam

- Can I use schmutzi on modern DNA for mitochondrial consensus ?
  Yes, it is the simple case where deamination does not exist. 
  Simply use endoCaller without deamination parameters and length parameters.


- The contamination estimate of contDeam.pl is very different from the one obtained from mtCont, why is that ?
  There are two possibilities:

  * schmutzi.pl did not converge 
  * The rates of deamination of the contamination are not close to zero.

  For the second option, maybe run using either --uselength or --estdeam and check the cont.prof/endo.prof files being produced.

- The program stopped and said: "Unable to find more than X positions that are different ...". Why is that ?

  This means that upon splitting into the endogenous and contaminant files, the program was not able
  to find enough positions to split. This could be due to the following:
  
  * There is either too much or too little contamination. This makes is impossible to predict the contaminant and/or the endogenous genomes. 
    If the first contamination estimate was very high, probably it was due to high contamination. Try to run endoCaller manually using higher
    contamination estimates. If there is too little contamination as indicated by an initial low contamination estimate, you could simply trust 
    the first iteration of endoCaller or rerun without the "--uselength", "--estdeam" and "--usepredC" options

  * There is very little difference between the endogenous and contaminant. You could trust the contamination rate given by the deamination 
    patterns and use the first endogenous call.

- Should I filter my BAM file for reads with high mapping quality ?
  
  * In theory, no. If the mapping quality is good proxy for the probability of mismapping, you should be fine as mapping quality is 
    incorporated


- How much can I trust the endogenous consensus call ?
 
  * That depends on two factors:
     1) Amount of contamination and coverage, you can check how good the quality of the log is. If you have low quality, not much can be done
     2) It the case of very distant mt genomes (ex: Denisova to human) the divergence might 
        be very high in certain regions and leads to misalignments. To solve this, an approach is to call the endogenous iteratively as such:
        
        ./wrapper.pl   testDenisova/endo.it1  refs/human_MT.fa testDenisova/input.bam
        ./endoCaller -seq testDenisova/endo.it1.fa -log testDenisova/endo.it1.log  refs/human_MT.fa testDenisova/endo.it1.bam
        samtools faidx testDenisova/endo.it1.fa
        ./wrapper.pl   testDenisova/endo.it2  testDenisova/endo.it1.fa testDenisova/endo.bam
        ./endoCaller -seq testDenisova/endo.it2.fa -log testDenisova/endo.it2.log testDenisova/endo.it1.fa  testDenisova/endo.it2.bam
         ...
        Until there is no difference
