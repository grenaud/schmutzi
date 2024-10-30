
  schmutzi: Bayesian maximum a posteriori contamination estimate for ancient samples
=====================================================================================

Upon sequencing ancient DNA, the DNA of the individuals involved in excavation, extraction 
of DNA and library preparation can become mixed with the actual sample being sequenced. 
We define as endogenous the DNA pertaining to the original sample and contaminant the DNA of the experimenters. 


schmutzi is a set of programs aimed at ancient DNA data that can :

* estimate contamination based on deamination patterns alone
* call a mitochondrial consensus for the endogenous genome. This consensus calls takes into account contamination and deamination.
* estimate mitochondrial contamination and identify the most likely contaminant from a set.

Questions :
-------------------------------------------------------------------------------------
	contact: Gabriel Renaud   
	email:	 gabriel [dot] reno [ at sign ] gmail.com

Downloading:
-------------------------------------------------------------------------------------
Do a :

	git clone --recursive https://github.com/grenaud/schmutzi.git

or download a zipped file from https://bioinf.eva.mpg.de/schmutzi/schmutzi.tar.gz

Requirements:
-------------------------------------------------------------------------------------

   - zlib.h
   - cmake to build bamtools
   - git to download the submodules
   - C++ compiler
   - Perl interpreter
   - R (recent version)
      -  Rscript 
      -  The fitdistrplus R package ( install.packages("fitdistrplus") or install.packages("fitdistrplus", repos="http://R-Forge.R-project.org") on more recent versions of R)
      -  The MASS package, normally installed automatically with fitdistrplus ( install.packages("MASS") )

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


For Mac users, open a terminal please type "cmake", "git" and "R" to check if works. If it is not installed, please install it. If not, Mac Users reported that Homebrew can install those.

For Windows users, you would need cygwin to run schmutzi as it runs on the terminal (https://www.cygwin.com/). 

1) type:

	make

which should build pretty much everything. 

We tested our program on a Linux system and MacOS, please email us if you have trouble building under cygwin. 

A comment about cygwin: the high precision version of the logarithm does not seem available by default. Therefore, 
we are forced to use the low precision variant for cygwin. For higher accuracy, run our software under Linux. Also, 
some users have reported bugs with the locale used by C++ 


Running:
-------------------------------------------------------------------------------------
There are 2 main Perl scripts that drive the underlying programs written in C++.


| script       | function                                                                                    |
| ------------ | ------------------------------------------------------------------------------------------  |
|contDeam.pl   | Is a script that allows for the estimation of endogenous deamination and                    |
|              | subsequent contamination estimate using those. Assuming that the contaminant                |
|              | is not deaminated, it will measure rates of deamination for the endogenous material         |
|	       | and provide a rough contamination estimate. If say the endogenous material is               |
|              | 30% deaminated at one base but the observed deamination rate for the entire set is 15%,     |
|              | the contamination rate is therefore 0.5 . There are two ways to estimate endogenous         |
|              | deamination:                                                                                |
|	       | 1) (default)                                                                                |
|	       |     Condition on one end of the read being deaminated (ex: 5') and measure                  |
|	       |     deamination rates on the other (ex: 3'). This approach is useful for                    |
|              |     all types of ancient DNA (nuclear and mitonchondrial). Please note however              |
|              |     that if the contaminant is deaminated as well, which can occur, this approach           |
|              |     will return an underestimate whose error comensurates with rates of contaminant         |
|              |     deamination.                                                                            |
|              |  2) --split                                                                                 |
|              |     Take a series of mitonchondrial diagnostic positions (positions that indicate           |
|              |     whether the read pertains to a particular type of hominin (ex: Neanderthal) or           |
|              |     the putative contaminant hominin (ex: modern Humans) and separate the contaminant       |
|              |     from the endogenous reads. This approach is useful when prior information on the        |
|              |     sample is available and when enough diagnostic positions are available (case of         |
|              |     archaic hominins).                                                                      |
|schmutzi.pl   | This script performs two tasks:                                                             |
|	       |	   1) Call an endogenous consensus                                                   |
|	       |   2) Using the consensus called in 1) and a set of known contaminants, compute              |
|              |       the contamination rate for entry in this set. Return the most likely contaminant      |
|	       |      with its corresponding rate.                                                           |
|	       | Both 1) and 2) are called iteratively until a stable contamination rate is reached. The     |
|	       | script produces a fasta file of the endogenous mitonchondrial genome                        |
                   
                       

Also, the various subparts of schmutzi can be called directly:



| program            | function                                                                                    |
| ------------------ | ------------------------------------------------------------------------------------------  |
|contDeam  	     |	Contamination estimate using deamination patterns                                          |
|endoCaller	     |	Consensus calling for mitochondrial data                                                   |
|mtCont		     |	Contamination estimate using a set of known contaminants                                   |
|                    |                                                                                             |
|approxDist.R	     |	Do a maximum likelihood estimate of log-normal paramaters given molecule length data       |
|bam2prof	     |	Measure deamination rates and produce a matrix of deamination probabilities                |
|insertSize	     |	Produce all insert sizes for aligned data in BAM format                                    |
|log2freq	     |	Turn a endogenous or consensus log file into an allele frequency matrix to be used         |
|		     |	as contaminant source                                                                      |
|logs2pos	     |	Takes two log files and produces which positions are segregating                           |
|mitoConsPDF.R	     |	Plot various information like coverage and quality for the consensus log file              |
|mas2freq	     |	Turn a multiple sequence alignment file in multifasta format into an allele frequency      |
|		     |	matrix to be used as contaminant source                                                    |
|mas2log	     |	Transform a pairwise alignment of an endogenous sequence to the reference                  |
|		     |	into a log file for contamination estimate using mtCont                                    |
|posteriorDeam.R     |	Plot the posterior probability                                                             |
|contOut2ContEst.pl  |	Take the output of mtcont and get a point estimate                                         |


Quick start guide:
-------------------------------------------------------------------------------------

Before you start, make sure you have a BAM file that has been:

- where ALL the fragments have been mapped to the mitochondrial genome only. 
- sorted
- indexed
- fillmd/calmd has been run to fix the MD field

First, determine if you the library used was double-stranded or single-stranded. Then call:
   
    contDeam.pl --lengthDeam [length] --library [library type] --out [output prefix] [mt reference] [input bam file]

where [length] is the # of bp considered and [library type] is the library type (either double or single). If the deamination is seen for a many bases far from the 5'/3' ends, use maybe --lengthDeam 40. If you have USER treatment or single-stranded libraries and low amounts of deamination further away from the ends, you can use --lengthDeam 2 or 5.

Then run the iterative procedure without the prediction of the contaminant:

    schmutzi.pl --notusepredC --uselength --ref [mt reference] --out [out prefix]_npred [output prefix] [path to schmutzi]/share/schmutzi/alleleFreqMT/eurasian/freqs/ [input bam file]

then run it again with the prediction of the contaminant:

    schmutzi.pl               --uselength --ref [mt reference] --out [out prefix]_wpred [output prefix] [path to schmutzi]/share/schmutzi/alleleFreqMT/eurasian/freqs/ [input bam file]

Those commands will create output files with the [out prefix]_npred and [out prefix]_wpred prefixes. 

If the iterative procedure is successful, it will some files with the following suffixes:

- _final.cont.est for the final estimate present-day human contamination (format: estimate estimate_low estimate_high)
- _final.cont.pdf plot for the posterior probability
- _final_mtcont.out log for mtCont (contamination estimate) allows the most likely haplogroup for the contaminant
- _final_endo.fa for the  unfiltered fasta prediction
- _final_endo.log log file with the per-base likelihood for each position for each base of of being the endogenous base. 

We highly recommend users to use log2fasta on _final_endo.log to obtain a consensus at various quality cutoffs. 





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

        mkdir -p outputdir/

Then run contDeam to estimation endogenous deamination rates:

        ./contDeam.pl  --library single --out outputdir/mez testdata/mezB9687.bam

or for the simulated

        ./contDeam.pl  --library double --out outputdir/sim testdata/simulation.bam    

This will produce the files:

        outputdir/[out].cont.pdf	Plot of the posterior probability for contamination based on deamination
        outputdir/[out].cont.est      Estimate for contamination based on deamination
        outputdir/[out].config	Configuration file describing the variables used


Then run the following to produce the endogenous consensus and the contamination estimate:

        ./schmutzi.pl       --uselength   --ref share/schmutzi/refs/human_MT_wrapped.fa         outputdir/mez   share/schmutzi/alleleFreqMT/alleleFreqMT/197/freqs/  testdata/mezB9687.bam
        ./schmutzi.pl       --uselength   --ref share/schmutzi/refs/human_MT_wrapped.fa         outputdir/sim   share/schmutzi/alleleFreqMT/alleleFreqMT/197/freqs/  testdata/simulation.bam


        --uselength tells the program to use the length of the molecules
        --ref is for the reference

       outputdir/mez is the output from contDeam
       share/schmutzi/alleleFreqMT/197/freqs/ is the database of putative contaminants 
       testdata/mezB9687.bam is the input bam file
  
The first dataset is an empirical dataset with about 40-45% contamination and the second is a simulated dataset with 20% contamination
   
It will run for a few minutes and produce the following files:

For contamination:

| file                      | content                                                                            |
| ------------------------- | -----------------------------------------------------------------------------------|
| [out]_final_mtcont.out    | Contamination estimates for all samples in the database                            |
| [out]_final.cont          | Contamination estimates for the sample in the database with the highest likelihood |
| [out]_final.cont.est	    | Contamination estimates for the most likely sample with confidence intervals       |
| [out]_final.cont.pdf	    | Posterior probability on the contamination for the most likely sample              |

For the respective genomes of the endogenous and contaminant:

| file                      | content                                                                            |
| ------------------------- | -----------------------------------------------------------------------------------|
| [out]_final_endo.fa	    |	Endogenous mitochondrial genome as fasta                                         |
| [out]_final_endo.log	    |	Endogenous mitochondrial genome as a log file with likelihoods on a PHRED scale  | 
| [out]_final_cont.fa	    |	Contaminant mitochondrial genome as fasta                                        |
| [out]_final_cont.log	    |	Contaminant mitochondrial genome as a log file with likelihoods on a PHRED scale |



Recommended workflow for ancient samples:
-------------------------------------------------------------------------------------

- Nuclear
  I recommend DICE (http://grenaud.github.io/dice) for nuclear samples. 
  The only thing schmutzi can do for nuclear data is to estimate contamination rates 
  using deamination patterns. 
  
  However, use this at your own risk, we know three factors that lead to wrong estimates:

  1. Insufficient # of molecules (we need at least 500M)
  2. Insufficient rates of endogenous deamination (we need upwards of 5%)
  3. No or very little deamination of the contaminant fragments
  4. Since we need to condition on one end to measure deamination on the other, we need independence between 5' and 3' deamination rates 

  This method is implemented in "contDeam.pl". For mt, we can at least cross-validate the contamination estimate using "mtCont" but not for nuclear data.
  If you have a contaminant that is deaminated this will lead to an underestimate. 
  Lack of independence between the 5' and 3' deamination rates can lead to overestimated
  deamination rates for the endogenous portion and an overestimate. To test this, two programs were added as part of the package:
  jointFreqDeaminated and jointFreqDeaminatedDouble (double-stranded). For samples with low amounts of contamination, they compute the 
  joint frequency of deamination:

        ./jointFreqDeaminated  in.bam > in.freq

  Then use 

        ./jointFreqDeaminated.R in.freq 

  If you get low p-values, the method should be safe to use, you only have to worry about a deaminated contaminant.

- MT

  1. Have your data aligned to a mitochondrial reference (see "refs/human_MT_wrapped.fa" for a wrapped reference) using  a sensitive mapper that produces BAM. A wrapper script is available with schmutzi (see Frequently Asked Questions)
  2. run samtools sort on your aligned bam file
  3. run samtools calmd/fillmd on your aligned bam file
  4. run samtools index on your sorted and aligned bam file
  5. If you used the wrapped reference, re-wrap your alignments exceeding the junction using for example https://github.com/mpieva/biohazard-tools this step is not necessary but produces equal coverage and resolution at the boundaries
  6. Estimate your initial contamination and deamination rates using "contDeam.pl"
  7. Run schmutzi.pl once with default parameters
  8. Run schmutzi.pl again  with "--usepredC"
  9. If contamination is more than a few percent re-run using the "--uselength" option
     
If you have a large number of samples, I recommend using bam2makeSchmutzi.pl which creates a makefile and automates the of using contDeam.pl followed by schmutzi.pl



Interpreting the output:
-------------------------------------------------------------------------------------

Here is a list of the different meaningful output files and their meaning. 

From contDeam.pl:

| file                      | content                                                                            |
| ------------------------- | -----------------------------------------------------------------------------------|
| [out].cont.pdf            | Plot of the posterior probability for contamination based on deamination           |
| [out].cont.est            | Estimate for contamination based on deamination. format: estimate est.low est.high  |           


From schmutzi.pl for the contamination estimate:


| file                      | content                                                                            |
| ------------------------- | -----------------------------------------------------------------------------------|
| [out]_final.cont.est	    | Contamination estimates for the most likely sample with confidence intervals. format: estimate est.low est.high      |
| [out]_final_mtcont.out    | Contamination estimates for all records in the database format: record contamination log.posterior |
| [out]_final.cont          | Contamination estimates for the record in the database with the highest likelihood. format: record contamination log.posterior |
| [out]_final.cont.pdf	    | Posterior probability on the contamination for the most likely sample              |

From schmutzi.pl for the endogenous/contaminant:

| file                      | content                                                                            |
| ------------------------- | -----------------------------------------------------------------------------------|
| [out]_final_endo.fa	    |	Endogenous mitochondrial genome as fasta, this contains all the bases and has not been filtered for high-quality bases. To produce a fasta file using a quality filter, use log2fasta                                          |
| [out]_final_endo.log	    |	Endogenous mitochondrial genome as a log file with likelihoods on a PHRED scale. format: position	reference.base	predicted.base	quality(PHRED)	average.mappping.quality	coverage	support.for.predicted.base	p[base=a]	p[base=c]	p[base=g]	p[base=t] |
| [out]_final_cont.fa	    |	Contaminant mitochondrial genome as fasta, this contains all the bases and has not been filtered for high-quality bases                                         |
| [out]_final_cont.log	    |	Contaminant mitochondrial genome as a log file with likelihoods on a PHRED scale. format: position	reference.base	predicted.base	quality(PHRED)	average.mappping.quality	coverage	support.for.predicted.base	p[base=a]	p[base=c]	p[base=g]	p[base=t] |


From bam2makeSchmutzi.pl:

You will get a _results.txt, with different columns. I suggest concatenating the _results.txt files and sorting with respect to coverage and flag problematic samples. Here is the meaning of the different columns (DB=database of putative contaminants): 


| column  | shorthand ID    | meaning                                                                                                |
| ------- | --------------- | -------------------------------------------------------------------------------------------------------|
| 1       | ID              | ID of the sample                                                                                       |
| 2       | deam5p          | rate of deamination at the 5' end                                                                      |
| 3       | deam3p          | rate of deamination at the 3' end                                                                      |
| 4       | contDeam        | contamination estimate obtained via deamination patterns                                               |
| 5       | contDeamL       | contamination estimate obtained via deamination patterns (lower bound)                                 |
| 6       | contDeamH       | contamination estimate obtained via deamination patterns (upper bound)                                 |
| 7       | contWPrd        | contamination estimate obtained by predicting the contaminant+DB of putative contaminant               |
| 8       | contWPrdL       | contamination estimate obtained by predicting the contaminant+DB of putative contaminant (lower bound) |
| 9       | contWPrdH       | contamination estimate obtained by predicting the contaminant+DB of putative contaminant (upper bound) |
| 10      | contNPrd        | contamination estimate obtained by DB only                                                             | 
| 11      | contNPrdL       | contamination estimate obtained by DB only (lower bound)                                               |
| 12      | contNPrdH       | contamination estimate obtained by DB only (upper bound)                                               |
| 13      | avgCov          | average coverage                                                                                       |
| 14      | hpGrq10w        | haplogrep for endogenous consensus obtained by predicting the contaminant and QC>10 (error<1/10)       |
| 15      | hpGrq10wq       | haplogrep quality for endogenous consensus obtained by predicting the contaminant and QC>10            |
| 16      | hpGrq30w        | haplogrep for endogenous consensus obtained by predicting the contaminant and QC>30 (error<1/1000)     |
| 17      | hpGrq30wq       | haplogrep quality for endogenous consensus obtained by predicting the contaminant and QC>30            |
| 18      | hpGrq50w        | haplogrep for endogenous consensus obtained by predicting the contaminant and QC>50 (error<1/100000)   |
| 19      | hpGrq50wq       | haplogrep quality for endogenous consensus obtained by predicting the contaminant and QC>50            |
| 20      | hpGrq10n        | haplogrep for endogenous consensus obtained by DB alone and QC>10 (error<1/10)                         |
| 21      | hpGrq10nq       | haplogrep quality for endogenous consensus obtained by DB alone and QC>10                              |
| 22      | hpGrq30n        | haplogrep for endogenous consensus obtained by DB alone and QC>30 (error<1/1000)                       |
| 23      | hpGrq30nq       | haplogrep quality for endogenous consensus obtained by DB alone and QC>30                              | 
| 24      | hpGrq50n        | haplogrep for endogenous consensus obtained by DB alone and QC>50 (error<1/100000)                     | 
| 25      | hpGrq50nq       | haplogrep quality for endogenous consensus obtained by DB alone and QC>50                              |



Frequently asked questions:
-------------------------------------------------------------------------------------

- How should I prepare the BAM file as input for endoCaller?
   
   As we mentioned above the, data must be aligned and ideally, wrapped around the ends.
   We provide a little wrapper script to call the shrimp mapper and wrap the reads around
   the edge : wrapper.pl


- When should I trust the output of contDeam.pl?

  When the deamination rates for the contaminant material are negligible. To 
  verify this, one way is to use diagnostic positions on the mitochondria. If 
  you have nuclear data, please retain the reads mapping to the mitochondria. 
  If you know in advance that the sample is either Neanderthal or Denisova, you can
  use  the precomputed diagnostic positions. If you have a modern human, you 
  will have to run schmutzi.pl and use the split BAM files produced by the script. 
  Inspect, the files called *cont.5p.prof and *cont.3p.prof. If the rates of deamination
  for them is higher than 1 or 2 \%, the results provided by contDeam.pl will not 
  reliable.
  
- The contamination estimate given by contDeam.pl gives me a very different result
  than schmutzi.pl does, why is that?

  There could be a few reasons:
  1. Your contamination is deaminated. That can be the case for museum samples.
     Check to see if the 

        - [output prefix]_[iteration]_cont.5p.prof 
     and 

        - [output prefix]_[iteration]_cont.3p.prof 

     present any sign of deamination.

  2. The algorithm diverged. This could be due to a misleading prior or an almost absence
     of endogenous material.

- How can I know if I have multiple contaminants?

    If you have highly divergent quality values for bases in the contamination prediction,
    this could be an indication that you might have multiple contaminants. Further, if  
    the predicted contaminant does not seem to fall within a given haplogroup and has
    diagnostic positions from different haplogroups, this could be an indication that 
    your contaminant is multiple. 

- What to do if I suspect that I have multiple contaminants?

    Do not use the --usepredC option and turn on the --multipleC option. If the contamination 
    rate is a lot lower than the one predicted by the deamination patterns, this could be an 
    indication that schmutzi did not converge and contamination is very high. Try to do a simple prediction
    using deaminated molecules (hoping that the contaminant is not deaminated) and use this as the input of 
    mtCont manually. Datasets that are heavily contaminated with multiple contaminants are very difficult 
    targets.

- How can I just use mtCont using a fasta file or a fasta file from a proxy mitochondria?
  
  Given you have a fasta file with the endogenous mitochondria called endo.fa and using the reference ref/human_MT.fa. You can run mafft and pipe into msa2log:

          cat endo.fa refs/human_MT.fa |mafft --auto /dev/stdin  |./msa2log /dev/stdin  mtref  > endo.log

  Then you have a log file for mtCont.

- Can I build my own set of putative contaminants?

  Yes. First, build a multiple sequence alignment (msa) of the various mitochondria
  including the reference used for mapping. Go to the directory where you want 
  to store the frequencies then run:
  
          ./msa2freq [msa] [name of the reference]
  
  A freqs/ directory should appear with one file per record in the msa.

- Can your software work for ancient DNA from other species?
  
  Yes, given a reference, schmutzi can be used to call a mitochondria endogenous 
  consensus taking into account deamination. But in the case of non-hominin 
  animals, you do no need to measure contamination from humans if the 
  mitochondrion is sufficiently divergent from the human one. High divergence 
  entails that few reads from humans will align to your reference.
  
  First determine your deamination rates using bam2prof with the appropriate library type:
  
           ./bam2prof [-single|-double] -5p out5p.prof -3p out3p.prof input.bam
  
  Then call the consensus:

            ./endoCaller -cont 0 -deam5p out5p.prof -deam3p out3p.prof -seq output.fa -log outlog.log -l [length reference] /path/to/reference.fasta input.bam

- Can I use schmutzi on modern DNA for mitochondrial consensus?

  Yes, it is a simple case where deamination does not exist. 
  Simply use endoCaller without deamination parameters and length parameters.


- The contamination estimate of contDeam.pl is very different from the one obtained from mtCont, why is that?
  There are two possibilities:

  1. schmutzi.pl did not converge 
  2. The rates of deamination of the contamination are not close to zero.

  For the second option, maybe run using either --uselength or --estdeam and check the cont.prof/endo.prof files being produced.

- The program stopped and said: "Unable to find more than X positions that are different ...". Why is that ?

  This means that the program could not find enough positions to split upon splitting into the endogenous and contaminant files. This could be due to the following:
  
  * There is either too much or too little contamination. This makes is impossible to predict the contaminant and/or the endogenous genomes. 
    If the first contamination estimate was very high, it was probably due to high contamination. Try to run endoCaller manually using higher
    contamination estimates. If there is too little contamination as indicated by an initial low contamination estimate, you could simply trust 
    the first iteration of endoCaller or rerun without the "--uselength", "--estdeam" and "--usepredC" options

  * There is very little difference between the endogenous and contaminant. You could trust the contamination rate given by the deamination 
    patterns and use the first endogenous call.

- I am getting a very high contamination rate (98%-99%), why?
  
  * This is likely a fluke, first make sure that you're not running with a prediction of the contaminant (using option --notusepredC). Not using this option only works if you have very high contamination rates. The second thing is to make sure that you have sufficient coverage (at least 10-15X on the mitochondrial genome). If not, that means that you cannot infer the endogenous consensus properly, however, we have a near-perfect resolution of the contaminant because we have the database. This means that the most likely explanation is that everything is contaminated. We highly suggest to disregard the contamination estimates for samples with less than 10x and doubt the estimate of those between 10-20X. 

- Should I filter my BAM file for reads with high mapping quality?
  
  * In theory, no. If the mapping quality is a good proxy for the probability of mismapping, you should be fine as mapping quality is 
    incorporated. However, most aligners do not compute the mismapping probability properly. We have found that sometimes applying a cutoff
    of 30 or 35 bp reduces the chance of mismapping.


- How much can I trust the endogenous consensus call?
 
  * That depends on two factors:
     1. Amount of contamination and coverage, you can check how good the quality of the log is. If you have low quality, not much can be done
     2. It the case of very distant mt genomes (ex: Denisova to human) the divergence might 
        be very high in certain regions and lead to misalignments. To solve this, an approach is to call the endogenous iteratively as such:
        
              ./wrapper.pl   testDenisova/endo.it1  refs/human_MT.fa testDenisova/input.bam
              ./endoCaller -seq testDenisova/endo.it1.fa -log testDenisova/endo.it1.log  refs/human_MT.fa testDenisova/endo.it1.bam
              samtools faidx testDenisova/endo.it1.fa
              ./wrapper.pl   testDenisova/endo.it2  testDenisova/endo.it1.fa testDenisova/endo.bam
              ./endoCaller -seq testDenisova/endo.it2.fa -log testDenisova/endo.it2.log testDenisova/endo.it1.fa  testDenisova/endo.it2.bam
              ...

        Until there is no difference


- I run contDeam.pl on my BAM file and get nan errors:

        ./contDeam.pl  --library double --out outputdir/test testdata/test.bam 
        running cmd bam2prof -length  2 -endo -double -5p outputdir/test.endo.5p.prof  -3p outputdir/test.endo.3p.prof testdata/test.bam
        running cmd contDeam  -deamread -deam5p outputdir/test.endo.5p.prof  -deam3p outputdir/test.endo.3p.prof  -log  outputdir/test.cont.deam   testdata/test.bam
        Utils.cpp: destringify() Unable to convert string="-nan"
        system  cmd contDeam  -deamread -deam5p outputdir/test.endo.5p.prof  -deam3p outputdir/test.endo.3p.prof  -log  outputdir/test.cont.deam   testdata/test.bam failed: 256 at ./contDeam.pl line 22.

  * This could be due to two factors:
     1. The number of aDNA fragments is too low to get an estimate of misincorporation patterns due to deamination. Less than 1X for example.
     2. There is no deamination present hence no measurable misincorporation patterns.

  Problem #1 cannot be solved per se. It is possible that the iterative procedure will not work either. Problem #2 requires a bit of thinking. I suggest running mapDamage and checking if you have fragmentation and misincorporation patterns. If so, was USER/UDG treatment used hence explaining the lack of deamination patterns? If you are confident that you have no deamination pattern, please note that is it very difficult to infer the endogenous base and therefore to estimate present-day contamination in the absence of deamination patterns. I recommend using endoCaller and calling mtCont manually just once instead of iteratively.   
     
- When endoCaller is run, get: Query reference base is not the same for read XXXX pos 3107
 
  Make sure you have used "samtools calmd" (or samtools fillmd for earlier versions).


- What is share/schmutzi/alleleFreqMT/197/freqs/

  This is a very small putative contaminant database of 197 human mitogenomes using a worldwide sampling. For samples collected in Eurasia, we recommend the eurasian/ database. However, this database might be better if you suspect that your contaminant is African in origin. As ancient DNA only expanded into Africa relatively recently, we are currently making a more comprehensive curated mitogenome database with a worldwide distribution.


- When I run bam2makeSchmutzi.pl, I get NA columns, why?

  If you're getting NA, it is likely due to a failure in one of the programs. I recommend running the programs manually (contDeam+schmutzi.pl) using the commands in the Makefile and check what you're getting. Also, check that you have all the R packages installed.


