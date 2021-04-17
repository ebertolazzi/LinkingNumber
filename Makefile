LIBS="-lm"

GCCF="-g0 -O -Wall" 
CXXF="-w0 -g0 -O -std strict_ansi"
KCCF="-g -O --strict"

cc:
	make CXX="c++ -Isrc -std=c++11" CXXFLAGS=${GCCF} all

gcc:
	make CXX="g++ -Isrc -std=c++11" CXXFLAGS=${GCCF} all

all:
	$(CXX) -Isrc tests/test1.cc src/linking_number.cc -o t1 $(LIBS)
	$(CXX) -Isrc tests/test2.cc src/linking_number.cc -o t2 $(LIBS)

clean:
	rm -f t1 t2
