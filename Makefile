CXX      = g++ #-g -pg
BAMTOOLS= bamtools/
LIBGAB   = libgab/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a -lpthread -lm -lz



all: endoCaller  mtCont damage2profile log2freq log2fasta contDeam msa2freq bam2prof insertSize splitEndoVsCont/poshap2splitbam logs2pos countRecords libgab/utils.o bamtools/lib/libbamtools.a

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
	cd bamtools/ && mkdir build/  && cd build/ && cmake .. && make && cd ../..


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

mtCont.o:	libgab/utils.o mtCont.cpp 
	${CXX} ${CXXFLAGS} mtCont.cpp

mtCont:	libgab/utils.o miscfunc.o mtCont.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

bam2prof.o:	libgab/utils.o bam2prof.cpp 
	${CXX} ${CXXFLAGS} bam2prof.cpp

bam2prof:	libgab/utils.o miscfunc.o bam2prof.o ${LIBGAB}utils.o  ${LIBGAB}/ReconsReferenceBAM.o  ${LIBGAB}gzstream/libgzstream.a
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

logs2pos.o:	libgab/utils.o logs2pos.cpp
	${CXX} ${CXXFLAGS} logs2pos.cpp

logs2pos:	libgab/utils.o logs2pos.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



msa2freq.o:	libgab/utils.o msa2freq.cpp
	${CXX} ${CXXFLAGS} msa2freq.cpp

msa2freq:	libgab/utils.o msa2freq.o ${LIBGAB}utils.o  libgab/FastQParser.o libgab/FastQObj.o libgab/gzstream/libgzstream.a 
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



damage2profile.o:	libgab/utils.o damage2profile.cpp
	${CXX} ${CXXFLAGS} damage2profile.cpp

damage2profile:	libgab/utils.o damage2profile.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)

testdata: testdata/mezB9687.bam testdata/simulation.bam


testdata/mezB9687.bam:
	wget -O mezB9687.bam       https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam
	wget -O mezB9687.bam.bai   https://bioinf.eva.mpg.de/schmutzi/testData/mezB9687.bam.bai

testdata/simulation.bam:
	wget -O simulation.bam     https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam
	wget -O simulation.bam.bai https://bioinf.eva.mpg.de/schmutzi/testData/simulation.bam.bai


clean :
	rm -f mtCont.o mtCont contDeam.o contDeam endoCaller endoCaller.o damage2profile.o damage2profile miscfunc.o log2freq log2fasta msa2freq contDeam insertSize countRecords logs2pos
	make -C splitEndoVsCont/ clean
	make -C libgab/ clean


