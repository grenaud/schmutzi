
CXX      = g++
BAMTOOLS= bamtools/
LIBGAB   = libgab/

CXXFLAGS  = -lm -O3 -Wall -I${LIBGAB} -I${LIBGAB}/gzstream/ -I${BAMTOOLS}/include/  -I${BAMTOOLS}/src/ -c
LDLIBS   +=    ${BAMTOOLS}/build/src/utils/CMakeFiles/BamTools-utils.dir/*cpp.o  ${BAMTOOLS}/lib/libbamtools.a -lm -lz

all: libgab/utils.o mitochondrialDeam  mtCont damage2profile


libgab/utils.o:
	make -C libgab/

miscfunc.o:	libgab/utils.o miscfunc.cpp
	${CXX} ${CXXFLAGS} miscfunc.cpp


mtCont.o:	libgab/utils.o mtCont.cpp 
	${CXX} ${CXXFLAGS} mtCont.cpp

mtCont:	libgab/utils.o miscfunc.o mtCont.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)


mitochondrialDeam.o:	libgab/utils.o mitochondrialDeam.cpp
	${CXX} ${CXXFLAGS} mitochondrialDeam.cpp

mitochondrialDeam:	libgab/utils.o mitochondrialDeam.o miscfunc.o ${LIBGAB}utils.o    ${LIBGAB}gzstream/libgzstream.a
	${CXX} $(LDFLAGS) -o $@ $^ $(LDLIBS)

damage2profile.o:	libgab/utils.o damage2profile.cpp
	${CXX} ${CXXFLAGS} damage2profile.cpp

damage2profile:	libgab/utils.o damage2profile.o ${LIBGAB}utils.o  ${LIBGAB}gzstream/libgzstream.a
	${CXX} -o $@ $^ $(LDLIBS) $(LDFLAGS)



clean :
	rm -f mtCont.o mtCont mitochondrialDeam mitochondrialDeam.o damage2profile.o damage2profile miscfunc.o

