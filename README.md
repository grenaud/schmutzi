=====================================================================================
  schmutzi: Bayesian maximum a posteriori contamination estimate for ancient samples
=====================================================================================

Upon sequencing ancient DNA for homonin samples, the DNA of the individuals involved in 
excavation, extraction of DNA and library preparation can become mixed with the actual 
sample being sequenced. We define as endogenous the DNA pertaining to the original sample 
and contaminant the DNA of the experimenters. schmutzi is a set of programs aimed at hominin 
ancient DNA data that can :
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


Requirements:
-------------------------------------------------------------------------------------
   - zlib.h
   - R
      -  The fitdistrplus R package
      -  The MASS package

Installation:
-------------------------------------------------------------------------------------
1) Build Bamtools first:

    cd bamtools/   
    mkdir build/   
    cd build/
    cmake ..
    make 
    cd ../..

2) type "make"


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

To test schmutzi, we have made available real data and simulated data. First you need to estimate endogenous
deamination rates. First create an output directory:
	    mkdir /output dir./

Then run contDeam to estimation endogenous deamination rates:
	    ./contDeam.pl  --library single --out /output dir./ testData/in.bam

This will produce the files:
     	       	  

Frequently asked questions:
-------------------------------------------------------------------------------------

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



