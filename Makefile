

all: src/endoCaller src/splitEndoVsCont/poshap2splitbam 

install: all
	cp -fv src/posteriorDeam.R src/approxDist.R src/schmutzi.pl src/endoCaller src/mtCont src/damage2profile src/log2freq src/log2fasta src/contDeam src/contDeam.pl src/msa2freq src/bam2prof src/insertSize src/logs2pos src/countRecords  src/jointFreqDeaminated src/jointFreqDeaminatedDouble src/msa2log src/log2ConsensusLog src/logdiff src/logmask src/filterlog src/msa2singlefreq src/subSampleBAM src/addXACircular src/bam2makeSchmutzi.pl /usr/bin/
	cp -fv src/splitEndoVsCont/poshap2splitbam /usr/bin/
	cp -rfv share/schmutzi /usr/share/



splitEndoVsCont/poshap2splitbam:
	make -C src/splitEndoVsCont/

src/endoCaller:
	make -C src/

testdata: testdata/mezB9687.bam testdata/simulation.bam


testdata/mezB9687.bam:
	wget -O testdata/mezB9687.bam       https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam
	wget -O testdata/mezB9687.bam.bai   https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam.bai

testdata/simulation.bam:
	wget -O testdata/simulation.bam     https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam
	wget -O testdata/simulation.bam.bai https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam.bai


clean :
	make -C src/ clean
	make -C src/splitEndoVsCont/ clean
	make -C lib/libgab/ clean
	rm -fv /usr/bin/posteriorDeam.R /usr/bin/approxDist.R /usr/bin/schmutzi.pl  /usr/bin/endoCaller /usr/bin/mtCont /usr/bin/damage2profile /usr/bin/log2freq /usr/bin/log2fasta /usr/bin/contDeam /usr/bin/contDeam.pl /usr/bin/msa2freq /usr/bin/bam2prof /usr/bin/insertSize /usr/bin/logs2pos /usr/bin/countRecords /usr/bin/jointFreqDeaminated /usr/bin/jointFreqDeaminatedDouble /usr/bin/msa2log /usr/bin/log2ConsensusLog /usr/bin/logdiff /usr/bin/logmask /usr/bin/filterlog /usr/bin/msa2singlefreq /usr/bin/subSampleBAM /usr/bin/addXACircular /usr/bin/bam2makeSchmutzi.pl
	rm -fv /usr/bin/poshap2splitbam
	rm -rfv /usr/share/schmutzi/


