CXX      = g++ #-g -pg
BAMTOOLS= bamtools/
LIBGAB   = libgab/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a -lpthread -lm -lz

#all: libgab/utils.o mitochondrialDeam  mitochondrialDeamCorrection mitochondrialDeamCorrectionSingle  mtCont damage2profile log2freq contDeam msa2freq

all: libgab/utils.o endoCaller  mtCont damage2profile log2freq contDeam msa2freq bam2prof insertSize splitEndoVsCont/poshap2splitbam logs2pos

splitEndoVsCont/poshap2splitbam:
	make -C splitEndoVsCont/

libgab/utils.o:
	make -C libgab/

miscfunc.o:	libgab/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp

insertSize.o:	libgab/utils.o insertSize.cpp 
	${CXX} ${CXXFLAGS} insertSize.cpp

insertSize: insertSize.o  ${LIBGAB}utils.o
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

#
#mitochondrialDeam.o:	libgab/utils.o mitochondrialDeam.cpp
#	${CXX} ${CXXFLAGS} mitochondrialDeam.cpp
#
#mitochondrialDeam:	libgab/utils.o mitochondrialDeam.o miscfunc.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
#	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)
#
#mitochondrialDeamCorrection.o:	libgab/utils.o  mitochondrialDeamCorrection.cpp
#	${CXX} ${CXXFLAGS} mitochondrialDeamCorrection.cpp
#
#mitochondrialDeamCorrection:	libgab/utils.o ${LIBGAB}/ReconsReferenceBAM.o mitochondrialDeamCorrection.o miscfunc.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
#	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)
#
#
#mitochondrialDeamCorrectionSingle.o:	libgab/utils.o  mitochondrialDeamCorrectionSingle.cpp
#	${CXX} ${CXXFLAGS} mitochondrialDeamCorrectionSingle.cpp
#
#mitochondrialDeamCorrectionSingle:	libgab/utils.o ${LIBGAB}/ReconsReferenceBAM.o mitochondrialDeamCorrectionSingle.o miscfunc.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
#	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)
#

endoCaller.o:	libgab/utils.o  endoCaller.cpp
	${CXX} ${CXXFLAGS} endoCaller.cpp

endoCaller:	libgab/utils.o ${LIBGAB}/ReconsReferenceBAM.o endoCaller.o miscfunc.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

log2freq.o:	libgab/utils.o log2freq.cpp
	${CXX} ${CXXFLAGS} log2freq.cpp

log2freq:	libgab/utils.o log2freq.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
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



clean :
	rm -f mtCont.o mtCont contDeam.o contDeam endoCaller endoCaller.o damage2profile.o damage2profile miscfunc.o log2freq msa2freq contDeam insertSize logs2pos
	make -C splitEndoVsCont/ clean
#	rm -f mtCont.o mtCont contDeam.o contDeam mitochondrialDeamCorrection mitochondrialDeamCorrectionSingle.o mitochondrialDeamCorrectionSingle mitochondrialDeamCorrection.o mitochondrialDeam mitochondrialDeam.o damage2profile.o damage2profile miscfunc.o log2freq msa2freq contDeam

