CXX=g++
CPPFLAGS=-O3

all: pairfq
	

pairfq: pairfq.cpp
	$(CXX) $(CPPFLAGS) $? -o $@

coverage: pairfq.cpp
	$(CXX) -coverage -O0 $? -o $@
	./coverage data/r1.fq data/r2.fq out/out.fq out/unpaired.fq
	gcov $?

check: pairfq
	mkdir -p out
	./pairfq data/r1.fq data/r2.fq out/out.fq out/unpaired.fq
	diff data/out.fq out/out.fq
	diff data/unpaired.fq out/unpaired.fq

test: check
	
