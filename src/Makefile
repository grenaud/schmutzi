CXX      = /usr/bin/g++ #-g -pg
BAMTOOLS = ../lib/bamtools/
LIBGAB   = ../lib/libgab/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/build/src/api/libbamtools.a -lpthread -lm -lz


all: endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq bam2prof insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords ${LIBGAB}/utils.o ${BAMTOOLS}/build/src/api/libbamtools.a jointFreqDeaminated jointFreqDeaminatedDouble msa2log log2ConsensusLog logdiff logmask filterlog msa2singlefreq subSampleBAM addXACircular

splitEndoVsCont/poshap2splitbam:
	make -C splitEndoVsCont/

${LIBGAB}/utils.h:
	rm -rf ${LIBGAB}
	git clone --recursive https://github.com/grenaud/libgab.git ${LIBGAB}

${LIBGAB}/utils.o: ${BAMTOOLS}/build/src/api/libbamtools.a   ${LIBGAB}/utils.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}
	git clone --recursive https://github.com/pezmaster31/bamtools.git ${BAMTOOLS}/ && cd ${BAMTOOLS}/ &&  git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

${BAMTOOLS}/build/src/api/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && cmake .. && make && cd ../../../src/

miscfunc.o:	${LIBGAB}/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp

insertSize.o:	${LIBGAB}/utils.o insertSize.cpp 
	${CXX} ${CXXFLAGS} insertSize.cpp

insertSize: insertSize.o  ${LIBGAB}/utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

countRecords.o:	${LIBGAB}/utils.o countRecords.cpp 
	${CXX} ${CXXFLAGS} countRecords.cpp

countRecords: countRecords.o  ${LIBGAB}/utils.o
	${CXX} -o $@ $^ $(LDLIBS) 

jointFreqDeaminated.o:	${LIBGAB}/utils.o jointFreqDeaminated.cpp 
	${CXX} ${CXXFLAGS} jointFreqDeaminated.cpp

jointFreqDeaminatedDouble.o:	${LIBGAB}/utils.o jointFreqDeaminatedDouble.cpp 
	${CXX} ${CXXFLAGS} jointFreqDeaminatedDouble.cpp

jointFreqDeaminated: jointFreqDeaminated.o  ${LIBGAB}/utils.o ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

jointFreqDeaminatedDouble: jointFreqDeaminatedDouble.o  ${LIBGAB}/utils.o ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

mtCont.o:	${LIBGAB}/utils.o mtCont.cpp 
	${CXX} ${CXXFLAGS} mtCont.cpp

mtCont:	${LIBGAB}/utils.o miscfunc.o mtCont.o ${LIBGAB}/utils.o    ${LIBGAB}/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

bam2prof.o:	${LIBGAB}/utils.o bam2prof.cpp 
	${CXX} ${CXXFLAGS} bam2prof.cpp

bam2prof:	${LIBGAB}/utils.o miscfunc.o bam2prof.o ${LIBGAB}/utils.o  ${LIBGAB}/ReconsReferenceBAM.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

subSampleBAM.o:	${LIBGAB}/utils.o subSampleBAM.cpp 
	${CXX} ${CXXFLAGS} subSampleBAM.cpp

subSampleBAM:	${LIBGAB}/utils.o subSampleBAM.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

addXACircular.o:	${LIBGAB}/utils.o addXACircular.cpp 
	${CXX} ${CXXFLAGS} addXACircular.cpp

addXACircular:	${LIBGAB}/utils.o addXACircular.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

contDeam.o:	${LIBGAB}/utils.o contDeam.cpp 
	${CXX} ${CXXFLAGS} contDeam.cpp

contDeam:	${LIBGAB}/utils.o ${LIBGAB}/ReconsReferenceBAM.o miscfunc.o contDeam.o ${LIBGAB}/utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


endoCaller.o:	${LIBGAB}/utils.o  endoCaller.cpp
	${CXX} ${CXXFLAGS} endoCaller.cpp

endoCaller:	${LIBGAB}/utils.o ${LIBGAB}/ReconsReferenceBAM.o endoCaller.o miscfunc.o ${LIBGAB}/utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

log2freq.o:	${LIBGAB}/utils.o log2freq.cpp
	${CXX} ${CXXFLAGS} log2freq.cpp

log2freq:	${LIBGAB}/utils.o log2freq.o ${LIBGAB}/utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

log2fasta.o:	${LIBGAB}/utils.o log2fasta.cpp
	${CXX} ${CXXFLAGS} log2fasta.cpp

log2fasta:	${LIBGAB}/utils.o log2fasta.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterlog.o:	${LIBGAB}/utils.o filterlog.cpp
	${CXX} ${CXXFLAGS} filterlog.cpp

filterlog:	${LIBGAB}/utils.o filterlog.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logs2pos.o:	${LIBGAB}/utils.o logs2pos.cpp
	${CXX} ${CXXFLAGS} logs2pos.cpp

logs2pos:	${LIBGAB}/utils.o logs2pos.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

log2ConsensusLog.o:	${LIBGAB}/utils.o log2ConsensusLog.cpp
	${CXX} ${CXXFLAGS} log2ConsensusLog.cpp

log2ConsensusLog:	${LIBGAB}/utils.o log2ConsensusLog.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logdiff.o:	${LIBGAB}/utils.o logdiff.cpp
	${CXX} ${CXXFLAGS} logdiff.cpp

logdiff:	${LIBGAB}/utils.o logdiff.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logmask.o:	${LIBGAB}/utils.o logmask.cpp
	${CXX} ${CXXFLAGS} logmask.cpp

logmask:	${LIBGAB}/utils.o logmask.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


msa2freq.o:	${LIBGAB}/utils.o msa2freq.cpp
	${CXX} ${CXXFLAGS} msa2freq.cpp

msa2freq:	${LIBGAB}/utils.o msa2freq.o ${LIBGAB}/utils.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

msa2singlefreq.o:	${LIBGAB}/utils.o msa2singlefreq.cpp
	${CXX} ${CXXFLAGS} msa2singlefreq.cpp

msa2singlefreq:	${LIBGAB}/utils.o msa2singlefreq.o ${LIBGAB}/utils.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


msa2log.o:	${LIBGAB}/utils.o msa2log.cpp
	${CXX} ${CXXFLAGS} msa2log.cpp

msa2log:	${LIBGAB}/utils.o msa2log.o ${LIBGAB}/utils.o  ${LIBGAB}/FastQParser.o ${LIBGAB}/FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



damage2profile.o:	${LIBGAB}/utils.o damage2profile.cpp
	${CXX} ${CXXFLAGS} damage2profile.cpp

damage2profile:	${LIBGAB}/utils.o damage2profile.o ${LIBGAB}/utils.o  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



clean :
	rm -f *.o endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq msa2log bam2prof subSampleBAM addXACircular insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords ${LIBGAB}/utils.o bamtools/lib/libbamtools.a jointFreqDeaminated jointFreqDeaminatedDouble log2ConsensusLog logdiff logmask filterlog msa2singlefreq
	make -C splitEndoVsCont/ clean
	make -C ${LIBGAB}/ clean

