CXX      = /usr/bin/g++ #-g -pg
BAMTOOLS= bamtools/
LIBGAB   = libgab/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a -lpthread -lm -lz



all: endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq bam2prof insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords libgab/utils.o bamtools/lib/libbamtools.a jointFreqDeaminated jointFreqDeaminatedDouble msa2log log2ConsensusLog logdiff logmask filterlog msa2singlefreq subSampleBAM

splitEndoVsCont/poshap2splitbam:
	make -C splitEndoVsCont/

libgab/utils.h:
	rm -rf libgab/
	git clone --recursive https://github.com/grenaud/libgab.git


libgab/utils.o: bamtools/lib/libbamtools.a  libgab/utils.h
	make -C libgab

bamtools/src/api/BamAlignment.h:
	rm -rf bamtools/
	git clone --recursive https://github.com/pezmaster31/bamtools.git


bamtools/lib/libbamtools.a: bamtools/src/api/BamAlignment.h
	cd bamtools/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../..


miscfunc.o:	libgab/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp

insertSize.o:	libgab/utils.o insertSize.cpp 
	${CXX} ${CXXFLAGS} insertSize.cpp

insertSize: insertSize.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

countRecords.o:	libgab/utils.o countRecords.cpp 
	${CXX} ${CXXFLAGS} countRecords.cpp

countRecords: countRecords.o  ${LIBGAB}utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

jointFreqDeaminated.o:	libgab/utils.o jointFreqDeaminated.cpp 
	${CXX} ${CXXFLAGS} jointFreqDeaminated.cpp

jointFreqDeaminatedDouble.o:	libgab/utils.o jointFreqDeaminatedDouble.cpp 
	${CXX} ${CXXFLAGS} jointFreqDeaminatedDouble.cpp

jointFreqDeaminated: jointFreqDeaminated.o  ${LIBGAB}utils.o ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

jointFreqDeaminatedDouble: jointFreqDeaminatedDouble.o  ${LIBGAB}utils.o ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

mtCont.o:	libgab/utils.o mtCont.cpp 
	${CXX} ${CXXFLAGS} mtCont.cpp

mtCont:	libgab/utils.o miscfunc.o mtCont.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

bam2prof.o:	libgab/utils.o bam2prof.cpp 
	${CXX} ${CXXFLAGS} bam2prof.cpp

bam2prof:	libgab/utils.o miscfunc.o bam2prof.o ${LIBGAB}utils.o  ${LIBGAB}/ReconsReferenceBAM.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

subSampleBAM.o:	libgab/utils.o subSampleBAM.cpp 
	${CXX} ${CXXFLAGS} subSampleBAM.cpp

subSampleBAM:	libgab/utils.o subSampleBAM.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

contDeam.o:	libgab/utils.o contDeam.cpp 
	${CXX} ${CXXFLAGS} contDeam.cpp

contDeam:	libgab/utils.o ${LIBGAB}/ReconsReferenceBAM.o miscfunc.o contDeam.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


endoCaller.o:	libgab/utils.o  endoCaller.cpp
	${CXX} ${CXXFLAGS} endoCaller.cpp

endoCaller:	libgab/utils.o ${LIBGAB}/ReconsReferenceBAM.o endoCaller.o miscfunc.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

log2freq.o:	libgab/utils.o log2freq.cpp
	${CXX} ${CXXFLAGS} log2freq.cpp

log2freq:	libgab/utils.o log2freq.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

log2fasta.o:	libgab/utils.o log2fasta.cpp
	${CXX} ${CXXFLAGS} log2fasta.cpp

log2fasta:	libgab/utils.o log2fasta.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterlog.o:	libgab/utils.o filterlog.cpp
	${CXX} ${CXXFLAGS} filterlog.cpp

filterlog:	libgab/utils.o filterlog.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logs2pos.o:	libgab/utils.o logs2pos.cpp
	${CXX} ${CXXFLAGS} logs2pos.cpp

logs2pos:	libgab/utils.o logs2pos.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

log2ConsensusLog.o:	libgab/utils.o log2ConsensusLog.cpp
	${CXX} ${CXXFLAGS} log2ConsensusLog.cpp

log2ConsensusLog:	libgab/utils.o log2ConsensusLog.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logdiff.o:	libgab/utils.o logdiff.cpp
	${CXX} ${CXXFLAGS} logdiff.cpp

logdiff:	libgab/utils.o logdiff.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logmask.o:	libgab/utils.o logmask.cpp
	${CXX} ${CXXFLAGS} logmask.cpp

logmask:	libgab/utils.o logmask.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


msa2freq.o:	libgab/utils.o msa2freq.cpp
	${CXX} ${CXXFLAGS} msa2freq.cpp

msa2freq:	libgab/utils.o msa2freq.o ${LIBGAB}utils.o  libgab/FastQParser.o libgab/FastQObj.o libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

msa2singlefreq.o:	libgab/utils.o msa2singlefreq.cpp
	${CXX} ${CXXFLAGS} msa2singlefreq.cpp

msa2singlefreq:	libgab/utils.o msa2singlefreq.o ${LIBGAB}utils.o  libgab/FastQParser.o libgab/FastQObj.o libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


msa2log.o:	libgab/utils.o msa2log.cpp
	${CXX} ${CXXFLAGS} msa2log.cpp

msa2log:	libgab/utils.o msa2log.o ${LIBGAB}utils.o  libgab/FastQParser.o libgab/FastQObj.o libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



damage2profile.o:	libgab/utils.o damage2profile.cpp
	${CXX} ${CXXFLAGS} damage2profile.cpp

damage2profile:	libgab/utils.o damage2profile.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testdata: testdata/mezB9687.bam testdata/simulation.bam


testdata/mezB9687.bam:
	wget -O testdata/mezB9687.bam       https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam
	wget -O testdata/mezB9687.bam.bai   https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam.bai

testdata/simulation.bam:
	wget -O testdata/simulation.bam     https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam
	wget -O testdata/simulation.bam.bai https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam.bai


clean :
	rm -f *.o endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq msa2log bam2prof subSampleBAM insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords libgab/utils.o bamtools/lib/libbamtools.a jointFreqDeaminated jointFreqDeaminatedDouble log2ConsensusLog logdiff logmask filterlog msa2singlefreq
	make -C splitEndoVsCont/ clean
	make -C libgab/ clean


