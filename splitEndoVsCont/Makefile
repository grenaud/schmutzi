
CXX      = g++
BAMTOOLS= /mnt/solexa/bin/bamtools-2.2.2/
LIBGAB   = /home/gabriel_renaud/lib/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I/home/gabriel_renaud/projects/divergenceCompute/ -I${LIBGAB}/VCFparser/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a -lm -lz

all: poshap2splitbam 

poshap2splitbam.o:	poshap2splitbam.cpp
	${CXX} ${CXXFLAGS} poshap2splitbam.cpp

poshap2splitbam:	poshap2splitbam.o ${LIBGAB}utils.o   /home/gabriel_renaud/projects/divergenceCompute/GenomicRange.o ${LIBGAB}VCFparser/gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)



clean :
	rm -f poshap2splitbam.o poshap2splitbam 

