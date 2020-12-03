CXX      = /usr/bin/g++ #-g -pg
BAMTOOLS = ../lib/bamtools/
LIBGAB   = ../lib/libgab/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/build/src/api/libbamtools.a -lpthread -lm -lz


all: endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq bam2prof insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords ${LIBGAB}/libgab.a ${BAMTOOLS}/build/src/api/libbamtools.a jointFreqDeaminated jointFreqDeaminatedDouble msa2log log2ConsensusLog logdiff logmask filterlog msa2singlefreq subSampleBAM addXACircular

splitEndoVsCont/poshap2splitbam:
	make -C splitEndoVsCont/

${LIBGAB}/libgab.h:
	rm -rf ${LIBGAB}
	git clone --recursive https://github.com/grenaud/libgab.git ${LIBGAB}

${LIBGAB}/libgab.a: ${BAMTOOLS}/build/src/api/libbamtools.a   ${LIBGAB}/libgab.h
	make -C ${LIBGAB}

${BAMTOOLS}/src/api/BamAlignment.h:
	rm -rf ${BAMTOOLS}
	git clone --recursive https://github.com/pezmaster31/bamtools.git ${BAMTOOLS}/ && cd ${BAMTOOLS}/ &&  git reset --hard d24d850de17134fe4e7984b26493c5c0a1844b35

${BAMTOOLS}/build/src/api/libbamtools.a: ${BAMTOOLS}/src/api/BamAlignment.h
	cd ${BAMTOOLS}/ && mkdir -p build/  && cd build/ && if cmake ..; then echo ""; else if cmake3 ..; then echo ""; else echo "cmake failed, please install cmake v3"; fi  fi && make && cd ../../../src/

miscfunc.o:	${LIBGAB}/libgab.a miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp

insertSize.o:	${LIBGAB}/libgab.a insertSize.cpp 
	${CXX} ${CXXFLAGS} insertSize.cpp

insertSize: insertSize.o  ${LIBGAB}/libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

countRecords.o:	${LIBGAB}/libgab.a countRecords.cpp 
	${CXX} ${CXXFLAGS} countRecords.cpp

countRecords: countRecords.o  ${LIBGAB}/libgab.a
	${CXX} -o $@ $^ $(LDLIBS) 

jointFreqDeaminated.o:	${LIBGAB}/libgab.a jointFreqDeaminated.cpp 
	${CXX} ${CXXFLAGS} jointFreqDeaminated.cpp

jointFreqDeaminatedDouble.o:	${LIBGAB}/libgab.a jointFreqDeaminatedDouble.cpp 
	${CXX} ${CXXFLAGS} jointFreqDeaminatedDouble.cpp

jointFreqDeaminated: jointFreqDeaminated.o  ${LIBGAB}/libgab.a ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

jointFreqDeaminatedDouble: jointFreqDeaminatedDouble.o  ${LIBGAB}/libgab.a ${LIBGAB}/ReconsReferenceBAM.o
	${CXX} -o $@ $^ $(LDLIBS) 

mtCont.o:	${LIBGAB}/libgab.a mtCont.cpp 
	${CXX} ${CXXFLAGS} mtCont.cpp

mtCont:	${LIBGAB}/libgab.a miscfunc.o mtCont.o ${LIBGAB}/libgab.a    ${LIBGAB}/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

bam2prof.o:	${LIBGAB}/libgab.a bam2prof.cpp 
	${CXX} ${CXXFLAGS} bam2prof.cpp

bam2prof:	${LIBGAB}/libgab.a miscfunc.o bam2prof.o ${LIBGAB}/libgab.a  ${LIBGAB}/ReconsReferenceBAM.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

subSampleBAM.o:	${LIBGAB}/libgab.a subSampleBAM.cpp 
	${CXX} ${CXXFLAGS} subSampleBAM.cpp

subSampleBAM:	${LIBGAB}/libgab.a subSampleBAM.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

addXACircular.o:	${LIBGAB}/libgab.a addXACircular.cpp 
	${CXX} ${CXXFLAGS} addXACircular.cpp

addXACircular:	${LIBGAB}/libgab.a addXACircular.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

contDeam.o:	${LIBGAB}/libgab.a contDeam.cpp 
	${CXX} ${CXXFLAGS} contDeam.cpp

contDeam:	${LIBGAB}/libgab.a ${LIBGAB}/ReconsReferenceBAM.o miscfunc.o contDeam.o ${LIBGAB}/libgab.a    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


endoCaller.o:	${LIBGAB}/libgab.a  endoCaller.cpp
	${CXX} ${CXXFLAGS} endoCaller.cpp

endoCaller:	${LIBGAB}/libgab.a ${LIBGAB}/ReconsReferenceBAM.o endoCaller.o miscfunc.o ${LIBGAB}/libgab.a    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

log2freq.o:	${LIBGAB}/libgab.a log2freq.cpp
	${CXX} ${CXXFLAGS} log2freq.cpp

log2freq:	${LIBGAB}/libgab.a log2freq.o ${LIBGAB}/libgab.a  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

log2fasta.o:	${LIBGAB}/libgab.a log2fasta.cpp
	${CXX} ${CXXFLAGS} log2fasta.cpp

log2fasta:	${LIBGAB}/libgab.a log2fasta.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

filterlog.o:	${LIBGAB}/libgab.a filterlog.cpp
	${CXX} ${CXXFLAGS} filterlog.cpp

filterlog:	${LIBGAB}/libgab.a filterlog.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logs2pos.o:	${LIBGAB}/libgab.a logs2pos.cpp
	${CXX} ${CXXFLAGS} logs2pos.cpp

logs2pos:	${LIBGAB}/libgab.a logs2pos.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

log2ConsensusLog.o:	${LIBGAB}/libgab.a log2ConsensusLog.cpp
	${CXX} ${CXXFLAGS} log2ConsensusLog.cpp

log2ConsensusLog:	${LIBGAB}/libgab.a log2ConsensusLog.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logdiff.o:	${LIBGAB}/libgab.a logdiff.cpp
	${CXX} ${CXXFLAGS} logdiff.cpp

logdiff:	${LIBGAB}/libgab.a logdiff.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

logmask.o:	${LIBGAB}/libgab.a logmask.cpp
	${CXX} ${CXXFLAGS} logmask.cpp

logmask:	${LIBGAB}/libgab.a logmask.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


msa2freq.o:	${LIBGAB}/libgab.a msa2freq.cpp
	${CXX} ${CXXFLAGS} msa2freq.cpp

msa2freq:	${LIBGAB}/libgab.a msa2freq.o ${LIBGAB}/libgab.a  ${LIBGAB}/FastQParser.o ${LIBGAB}/FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

msa2singlefreq.o:	${LIBGAB}/libgab.a msa2singlefreq.cpp
	${CXX} ${CXXFLAGS} msa2singlefreq.cpp

msa2singlefreq:	${LIBGAB}/libgab.a msa2singlefreq.o ${LIBGAB}/libgab.a  ${LIBGAB}/FastQParser.o ${LIBGAB}/FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)


msa2log.o:	${LIBGAB}/libgab.a msa2log.cpp
	${CXX} ${CXXFLAGS} msa2log.cpp

msa2log:	${LIBGAB}/libgab.a msa2log.o ${LIBGAB}/libgab.a  ${LIBGAB}/FastQParser.o ${LIBGAB}/FastQObj.o ${LIBGAB}/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



damage2profile.o:	${LIBGAB}/libgab.a damage2profile.cpp
	${CXX} ${CXXFLAGS} damage2profile.cpp

damage2profile:	${LIBGAB}/libgab.a damage2profile.o ${LIBGAB}/libgab.a  ${LIBGAB}/gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



clean :
	rm -f *.o endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq msa2log bam2prof subSampleBAM addXACircular insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords ${LIBGAB}/libgab.a bamtools/lib/libbamtools.a jointFreqDeaminated jointFreqDeaminatedDouble log2ConsensusLog logdiff logmask filterlog msa2singlefreq
	make -C splitEndoVsCont/ clean
	make -C ${LIBGAB}/ clean


